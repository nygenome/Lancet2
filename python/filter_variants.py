#!/usr/bin/env python3

import argparse
import logging
import math
import pickle
import urllib.request

import pysam
import tqdm

SOMATIC_URL = "https://storage.googleapis.com/lancet-ml-models/somatic_ebm_clf.v1.pkl"


def download_and_load_model(link):
    with urllib.request.urlopen(link) as response:
        return pickle.loads(response.read())


def phred_score(error_probability):
    return round(-10 * math.log10(error_probability), 2)


def build_variant_data(variant, normal_sample_idx=0, tumor_sample_idx=1):
    smpls = [i for i in variant.samples.itervalues()]

    NML = normal_sample_idx
    TMR = tumor_sample_idx
    REF = 0
    ALT = 1

    nml_ref_cnt = smpls[NML]["AD"][REF]
    nml_fwd_ref_cnt = smpls[NML]["ADF"][REF]

    nml_alt_cnt = smpls[NML]["AD"][ALT]
    nml_fwd_alt_cnt = smpls[NML]["ADF"][ALT]

    tmr_ref_cnt = smpls[TMR]["AD"][REF]
    tmr_fwd_ref_cnt = smpls[TMR]["ADF"][REF]

    tmr_alt_cnt = smpls[TMR]["AD"][ALT]
    tmr_fwd_alt_cnt = smpls[TMR]["ADF"][ALT]

    data = [
        # REF_ALLELE_LENGTH
        len(variant.ref),
        # ALT_ALLELE_LENGTH
        len(variant.alts[0]),
        # IS_SHARED_VARIANT
        variant.info.get("SHARED", False),
        # IS_TUMOR_ONLY_VARIANT
        variant.info.get("TUMOR", False),
        # VARIANT_LENGTH
        variant.info.get("LENGTH"),
        # KMER_LENGTH
        variant.info.get("KMERLEN"),
        # IS_STR_VARIANT
        variant.info.get("STR", False),
        # STR_REGION_LENGTH
        variant.info.get("STR_LEN", 0),
        # STR_MOTIF_LENGTH
        len(variant.info.get("STR_MOTIF", "")),
        # NORMAL_DEPTH
        smpls[NML]["DP"],
        # NORMAL_REF_COUNT
        nml_ref_cnt,
        # NORMAL_ALT_COUNT
        nml_alt_cnt,
        # NORMAL_REF_PCT_STRAND_IMBALANCE
        round(abs((0.5 - (0.5 if nml_ref_cnt == 0 else nml_fwd_ref_cnt / nml_ref_cnt)) * 100), 2),
        # NORMAL_ALT_PCT_STRAND_IMBALANCE
        round(abs((0.5 - (0.5 if nml_alt_cnt == 0 else nml_fwd_alt_cnt / nml_alt_cnt)) * 100), 2),
        # NORMAL_PCT_FAIL_READS_IN_WINDOW
        round(100 - (100 * smpls[NML]["PRF"]), 2),
        # NORMAL_VAF
        smpls[NML]["VAF"],
        # NORMAL_REF_ALLELE_QUALITY_MINIMUM
        smpls[NML]["RAQS"][0],
        # NORMAL_ALT_ALLELE_QUALITY_MINIMUM
        smpls[NML]["AAQS"][0],
        # NORMAL_REF_ALLELE_QUALITY_MEDIAN
        smpls[NML]["RAQS"][1],
        # NORMAL_ALT_ALLELE_QUALITY_MEDIAN
        smpls[NML]["AAQS"][1],
        # NORMAL_REF_ALLELE_QUALITY_MAXIMUM
        smpls[NML]["RAQS"][2],
        # NORMAL_ALT_ALLELE_QUALITY_MAXIMUM
        smpls[NML]["AAQS"][2],
        # NORMAL_REF_ALLELE_QUALITY_ABSDEV
        smpls[NML]["RAQS"][3],
        # NORMAL_ALT_ALLELE_QUALITY_ABSDEV
        smpls[NML]["AAQS"][3],
        # NORMAL_REF_MAPPING_QUALITY_MINIMUM
        smpls[NML]["RMQS"][0],
        # NORMAL_ALT_MAPPING_QUALITY_MINIMUM
        smpls[NML]["AMQS"][0],
        # NORMAL_REF_MAPPING_QUALITY_MEDIAN
        smpls[NML]["RMQS"][1],
        # NORMAL_ALT_MAPPING_QUALITY_MEDIAN
        smpls[NML]["AMQS"][1],
        # NORMAL_REF_MAPPING_QUALITY_MAXIMUM
        smpls[NML]["RMQS"][2],
        # NORMAL_ALT_MAPPING_QUALITY_MAXIMUM
        smpls[NML]["AMQS"][2],
        # NORMAL_REF_MAPPING_QUALITY_ABSDEV
        smpls[NML]["RMQS"][3],
        # NORMAL_ALT_MAPPING_QUALITY_ABSDEV
        smpls[NML]["AMQS"][3],
        # NORMAL_REF_ALN_PCT_DIFF_MINIMUM
        smpls[NML]["RAPDS"][0],
        # NORMAL_ALT_ALN_PCT_DIFF_MINIMUM
        smpls[NML]["AAPDS"][0],
        # NORMAL_REF_ALN_PCT_DIFF_MEDIAN
        smpls[NML]["RAPDS"][1],
        # NORMAL_ALT_ALN_PCT_DIFF_MEDIAN
        smpls[NML]["AAPDS"][1],
        # NORMAL_REF_ALN_PCT_DIFF_MAXIMUM
        smpls[NML]["RAPDS"][2],
        # NORMAL_ALT_ALN_PCT_DIFF_MAXIMUM
        smpls[NML]["AAPDS"][2],
        # NORMAL_REF_ALN_PCT_DIFF_ABSDEV
        smpls[NML]["RAPDS"][3],
        # NORMAL_ALT_ALN_PCT_DIFF_ABSDEV
        smpls[NML]["AAPDS"][3],
        # TUMOR_DEPTH
        smpls[TMR]["DP"],
        # TUMOR_REF_COUNT
        tmr_ref_cnt,
        # TUMOR_ALT_COUNT
        tmr_alt_cnt,
        # TUMOR_REF_PCT_STRAND_IMBALANCE
        round(abs((0.5 - (0.5 if tmr_ref_cnt == 0 else tmr_fwd_ref_cnt / tmr_ref_cnt)) * 100), 2),
        # TUMOR_ALT_PCT_STRAND_IMBALANCE
        round(abs((0.5 - (0.5 if tmr_alt_cnt == 0 else tmr_fwd_alt_cnt / tmr_alt_cnt)) * 100), 2),
        # TUMOR_PCT_FAIL_READS_IN_WINDOW
        round(100 - (100 * smpls[TMR]["PRF"]), 2),
        # TUMOR_VAF
        smpls[TMR]["VAF"],
        # TUMOR_REF_ALLELE_QUALITY_MINIMUM
        smpls[TMR]["RAQS"][0],
        # TUMOR_ALT_ALLELE_QUALITY_MINIMUM
        smpls[TMR]["AAQS"][0],
        # TUMOR_REF_ALLELE_QUALITY_MEDIAN
        smpls[TMR]["RAQS"][1],
        # TUMOR_ALT_ALLELE_QUALITY_MEDIAN
        smpls[TMR]["AAQS"][1],
        # TUMOR_REF_ALLELE_QUALITY_MAXIMUM
        smpls[TMR]["RAQS"][2],
        # TUMOR_ALT_ALLELE_QUALITY_MAXIMUM
        smpls[TMR]["AAQS"][2],
        # TUMOR_REF_ALLELE_QUALITY_ABSDEV
        smpls[TMR]["RAQS"][3],
        # TUMOR_ALT_ALLELE_QUALITY_ABSDEV
        smpls[TMR]["AAQS"][3],
        # TUMOR_REF_MAPPING_QUALITY_MINIMUM
        smpls[TMR]["RMQS"][0],
        # TUMOR_ALT_MAPPING_QUALITY_MINIMUM
        smpls[TMR]["AMQS"][0],
        # TUMOR_REF_MAPPING_QUALITY_MEDIAN
        smpls[TMR]["RMQS"][1],
        # TUMOR_ALT_MAPPING_QUALITY_MEDIAN
        smpls[TMR]["AMQS"][1],
        # TUMOR_REF_MAPPING_QUALITY_MAXIMUM
        smpls[TMR]["RMQS"][2],
        # TUMOR_ALT_MAPPING_QUALITY_MAXIMUM
        smpls[TMR]["AMQS"][2],
        # TUMOR_REF_MAPPING_QUALITY_ABSDEV
        smpls[TMR]["RMQS"][3],
        # TUMOR_ALT_MAPPING_QUALITY_ABSDEV
        smpls[TMR]["AMQS"][3],
        # TUMOR_REF_ALN_PCT_DIFF_MINIMUM
        smpls[TMR]["RAPDS"][0],
        # TUMOR_ALT_ALN_PCT_DIFF_MINIMUM
        smpls[TMR]["AAPDS"][0],
        # TUMOR_REF_ALN_PCT_DIFF_MEDIAN
        smpls[TMR]["RAPDS"][1],
        # TUMOR_ALT_ALN_PCT_DIFF_MEDIAN
        smpls[TMR]["AAPDS"][1],
        # TUMOR_REF_ALN_PCT_DIFF_MAXIMUM
        smpls[TMR]["RAPDS"][2],
        # TUMOR_ALT_ALN_PCT_DIFF_MAXIMUM
        smpls[TMR]["AAPDS"][2],
        # TUMOR_REF_ALN_PCT_DIFF_ABSDEV
        smpls[TMR]["RAPDS"][3],
        # TUMOR_ALT_ALN_PCT_DIFF_ABSDEV
        smpls[TMR]["AAPDS"][3],
    ]

    return data


def make_somatic_variant_list(raw_vcf_path):
    data = list()
    for v in tqdm.tqdm(pysam.VariantFile(raw_vcf_path)):
        # Skip adding NORMAL only variants,
        # since they don't need to be eval'ed
        # by the somatic machine learning model
        if not v.info.get("NORMAL", False):
            data.append(build_variant_data(v))
    return data


def main(raw_vcf_path):
    # Setup basic logging configuration
    msg_fmt = "%(asctime)s | %(levelname)s | %(message)s"
    dt_fmt = "%Y-%m-%d %H:%M:%S"
    logging.basicConfig(format=msg_fmt, level=logging.INFO, datefmt=dt_fmt)

    logging.info("Building SHARED/TUMOR variants dataframe for further evaluation")
    variants = make_somatic_variant_list(raw_vcf_path)
    logging.info(f"Done building dataframe with {len(variants)} SHARED/TUMOR variants")

    logging.info("Loading somatic machine learning model into memory")
    ml_model = download_and_load_model(SOMATIC_URL)

    logging.info("Applying machine learning model to the variants dataframe")
    preds = ml_model.predict(variants)
    probs = ml_model.predict_proba(variants)
    logging.info("Done applying machine learning model to the variants dataframe")

    logging.info("Writing final filtered and scored output VCF")
    invcf = pysam.VariantFile(raw_vcf_path)
    outvcf = pysam.VariantFile("-", "w", header=invcf.header)

    idx_ml_model = -1
    for v in tqdm.tqdm(invcf):
        if v.info.get("NORMAL", False):
            continue

        idx_ml_model += 1
        # Use lowest probablity class of machine learning binary
        # classification as error probability for phred score
        score = phred_score(sorted(probs[idx_ml_model])[0])
        is_pass_variant = preds[idx_ml_model]
        if is_pass_variant:
            v.qual = score
            v.filter.add("PASS")
            outvcf.write(v)

    invcf.close()
    outvcf.close()
    logging.info("Done writing final filtered output VCF")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(prog="filter_variants.py", description="Filter and score Lancet2 variants")
    parser.add_argument("lancet_raw_vcf", help="Path to raw lancet VCF to be filtered")
    args = parser.parse_args()
    main(args.lancet_raw_vcf)
