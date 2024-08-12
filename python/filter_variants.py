#!/usr/bin/env python3

import argparse
import logging
import math
import pickle

import pysam
import tqdm


def load_model(path):
    with open(path, "rb") as f:
        return pickle.load(f)


def pct_strand_imbalance(combined_cnt, fwd_strand_cnt):
    try:
        deviation = 0.5 - (fwd_strand_cnt / combined_cnt)
        return round(deviation * 100, 2)
    except ZeroDivisionError:
        return 0


def parse_sample_fields(variant, sidx, out_key_prefix):
    data = {k: variant.samples[sidx][k] for k in variant.format.keys()}
    data = {k: round(v, 2) if isinstance(v, float) else v for k, v in data.items()}

    result = {
        "DEPTH": data["DP"],
        "REF_COUNT": data["AD"][0],
        "ALT_COUNT": data["AD"][1],
        "REF_RATIO": 0 if data["DP"] == 0 else round(data["AD"][0] / data["DP"], 3),
        "ALT_RATIO": 0 if data["DP"] == 0 else round(data["AD"][1] / data["DP"], 3),
        "REF_PCT_STRAND_BIAS": pct_strand_imbalance(data["AD"][0], data["ADF"][0]),
        "ALT_PCT_STRAND_BIAS": pct_strand_imbalance(data["AD"][1], data["ADF"][1]),
        "PCT_HQ_READS_IN_WINDOW": 0 if math.isinf(data["PRF"]) else round(data["PRF"] * 100, 2),
    }

    # Round 0.0 to 0, since VAF is rounded to 3 decimal points
    result["REF_RATIO"] = 0 if result["REF_RATIO"] == 0 else result["REF_RATIO"]
    result["ALT_RATIO"] = 0 if result["ALT_RATIO"] == 0 else result["ALT_RATIO"]

    for k in ("RAQS", "AAQS", "RMQS", "AMQS", "RAPDS", "AAPDS"):
        minimum, median, maximum, absdev = data[k]
        out_prefix = k[:-1]
        result.update(
            {
                f"{out_prefix}_MINIMUM": minimum,
                f"{out_prefix}_MEDIAN": median,
                f"{out_prefix}_MAXIMUM": maximum,
                f"{out_prefix}_ABSDEV": absdev,
            }
        )

    return {out_key_prefix + k: v for k, v in result.items()}


def get_variant_state(v):
    return "SHARED" if v.info.get("SHARED") else "NORMAL" if v.info.get("NORMAL") else "TUMOR"


def build_variant_info(variant):
    result = {
        "VARIANT_STATE": get_variant_state(variant),
        "VARIANT_TYPE": variant.info.get("TYPE"),
        "KMER_LENGTH": variant.info.get("KMERLEN"),
        "VARIANT_LENGTH": variant.info.get("LENGTH"),
        "SOMATIC_FET_SCORE": round(variant.qual, 2),
        "VARIANT_NEAR_STR": variant.info.get("STR"),
    }

    result.update(parse_sample_fields(variant, 0, "NML_"))
    result.update(parse_sample_fields(variant, 1, "TMR_"))

    vaf_diff = round(result["TMR_ALT_RATIO"] - result["NML_ALT_RATIO"], 2)
    result.update(ABS_TMR_NML_VAF_DIFF=abs(vaf_diff))

    return result


def phred_score(error_probability):
    return round(-10 * math.log10(error_probability), 2)


def main(raw_vcf_path, model_path):
    # Setup basic logging configuration
    msg_fmt = "%(asctime)s | %(levelname)s | %(message)s"
    dt_fmt = "%Y-%m-%d %H:%M:%S"
    logging.basicConfig(format=msg_fmt, level=logging.INFO, datefmt=dt_fmt)
    logging.info(f"Starting to filter and score input VCF {raw_vcf_path}")

    logging.info(f"Building variants dataframe for ML model evaluation")
    variants_itr = tqdm.tqdm(pysam.VariantFile(raw_vcf_path, threads=2))
    df = [tuple(build_variant_info(v).values()) for v in variants_itr]
    logging.info(f"Done building dataframe with {len(df)} variants")

    logging.info("Loading somatic ML model into memory")
    ml_model = load_model(model_path)
    logging.info("Done loading somatic ML model into memory")

    logging.info("Applying somatic ML model to the variants dataframe")
    predictions = ml_model.predict(df)
    probabilities = ml_model.predict_proba(df)
    logging.info("Done applying somatic ML model to the variants dataframe")

    logging.info("Writing filtered and scored variants in VCF format to stdout")
    invcf = pysam.VariantFile(raw_vcf_path, threads=2)
    outvcf = pysam.VariantFile("-", "w", header=invcf.header)

    for idx, v in tqdm.tqdm(enumerate(invcf)):
        class_scores = probabilities[idx]
        error_probability = sorted(class_scores)[0]
        score = phred_score(error_probability)
        is_pass_variant = predictions[idx]
        v.qual = score

        is_snv = len(v.ref) == 1 and len(v.alts[0]) == 1
        is_indel = len(v.ref) != 1 or len(v.alts[0]) != 1

        if is_pass_variant and is_snv and score >= 10:
            v.filter.add("PASS") # SNV needs > 90% probability
        elif is_pass_variant and is_indel and score >= 13:
            v.filter.add("PASS") # InDel needs > 95% probability

        outvcf.write(v)

    invcf.close()
    outvcf.close()
    logging.info(f"Done writing filtered and scored variants to output VCF")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="filter_variants.py", description="Filter and score Lancet2 variants")
    parser.add_argument("lancet_raw_vcf", help="Path to raw lancet VCF to be filtered")
    parser.add_argument("ml_model", help="Path to ML model to be used for filtering")
    args = parser.parse_args()
    main(args.lancet_raw_vcf, args.ml_model)
