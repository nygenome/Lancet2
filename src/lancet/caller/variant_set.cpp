#include "lancet/caller/variant_set.h"

#include <algorithm>
#include <array>
#include <cstdint>
#include <limits>
#include <string>
#include <string_view>
#include <vector>

#include "absl/container/flat_hash_map.h"
#include "absl/types/span.h"
#include "spoa/graph.hpp"
#include "lancet/base/assert.h"
#include "lancet/base/types.h"
#include "lancet/caller/raw_variant.h"
#include "lancet/core/window.h"

namespace {

using lancet::caller::RawVariant;

// =========================================================================================
// STRICT SEQUENCE CORE LENGTH CALCULATOR
// -----------------------------------------------------------------------------------------
// Calculates the biological length of a mutation independent of VCF padding requirements.
// Without extracting the core, MNP lengths are vastly artificially inflated by the shared 
// anchored padding bases bounding multi-allelic clusters!
// =========================================================================================
inline auto CalculateVariantLength(std::string_view ref, std::string_view alt, RawVariant::Type vtype) -> i64 {
  if (vtype == RawVariant::Type::SNV) return 1;

  const auto ref_len = static_cast<i64>(ref.length());
  const auto alt_len = static_cast<i64>(alt.length());
  const auto diff = alt_len - ref_len;

  // Net structural variance mathematically aligns trivially with string length discrepancies
  if (vtype == RawVariant::Type::INS || vtype == RawVariant::Type::DEL || vtype == RawVariant::Type::CPX) {
    return diff;
  }

  // For strict MNPs (diff == 0), the biological length IS exactly the Sequence Core!
  // E.g., REF="ATGC", ALT="ACCC". Squeezing `start=1` ('A') and `end=1` ('C') bounds 
  // explicitly extracts the pure mutation core `TG`->`CC` (length 4 - 1 - 1 = 2).
  usize start_match = 0;
  while (start_match < ref_len && start_match < alt_len && ref[start_match] == alt[start_match]) {
    start_match++;
  }

  usize end_match = 0;
  while (end_match < (ref_len - start_match) && end_match < (alt_len - start_match) &&
         ref[ref_len - 1 - end_match] == alt[alt_len - 1 - end_match]) {
    end_match++;
  }

  return alt_len - static_cast<i64>(start_match) - static_cast<i64>(end_match);
}

}  // namespace

static constexpr usize REF_HAP_IDX = 0;

namespace lancet::caller {

// =========================================================================================
// THE INTUITIVE PRIMER: NATIVE TOPOLOGICAL BUBBLE SINKING
// =========================================================================================
// Rather than relying on rigid 2D String Matrices (MSA) which are mathematically bloated, or 
// pairwise mapping which erroneously artificially fractures overlapping Multiallelic loci, 
// we view biological sequencing natively as tracking Pointers traversing a structured DAG graph.
//
// -> 1. BIALLELIC SNVs (Paths traverse independent topological nodes before natively reuniting):
//
//                                [REF]
//                          .--> (T)[3] --.
//                         /               \
//   Anchor: (A)[2] ------+                 +-----> Target: (G)[5] (CONVERGED!)
//                         \               /
//                          `--> (C)[4] --'
//                                [ALT]
//
// -> 2. DELETIONS (Independent sequences shortcut massive blocks of graph architecture directly):
//
//                                            [REF]
//                          .--> (T)[3] --> (C)[4] --> (G)[5] --.
//                         /                                     \
//   Anchor: (A)[2] ------+                                       +-----> Target: (T)[6] (CONVERGED!)
//                         \                                     /
//                          `-----------------------------------'
//                                            [ALT]
//
// -> 3. MULTIALLELIC COMPLEXES (Massive concurrent branching geometries across the locus):
//       This requires sweeping ALL pointer coordinates entirely in unison natively. 
//
//                               [ALT 1]
//                          .--> (C)[3] --> (A)[4] ---------.
//                         /                                 \
//   Anchor: (T)[2] ------+----> (T)[5] --> (G)[6] --> (C)[7] +-----> Target: (A)[10] (CONVERGED!)
//                         \     [REF]                       /
//                          `--> (G)[8] --> (A)[9] ---------'
//                               [ALT 2]
//
// HOW THE SWEEP EXTRACTOR ALGORITHM WORKS: ("Greedy Sink")
// 1. We track an array of active pointers representing exactly where every path is locally.
// 2. If the pointers lose uniform consensus (meaning they point to different topological DAG nodes), 
//    a mathematical "Bubble" has triggered!
// 3. To perfectly map the bubble without traversing paths unevenly, we evaluate every active pointer's 
//    `Rank` (its globally unique, topologically sorted 5'-to-3' graph depth indicator).
// 4. We continually advance EXCLUSIVELY the pointer that possesses the MATHEMATICALLY LOWEST rank, 
//    strictly extracting its characters into a string buffer natively. 
// 5. This causes computationally "lagging" pointers to march forward until mathematically EVERY path  
//    inherently syncs identically onto the exact same `Rank` (the unified convergence target!).
// 6. We universally collect the buffered sequence strings across all paths, route them through the 
//    `VariantBubble` for robust VCF multi-allelic trimming parsimony, and register the payload!
// =========================================================================================

// Isolates algorithmic string Math and left-align computations fully away from graph logic headers
class VariantBubble {
public:
    // DATA STRUCTURE CHOICE: We map the variant string -> array of integers (haplotype ids)
    // Abseil maintains robust O(1) retrieval dynamically mapping haplotypes concurrently. 
    absl::flat_hash_map<std::string, std::vector<usize>> alt_alleles_to_haps; // 40B+ natively
    std::string ref_allele;                                                   // 24B
    std::vector<usize> hap_starts;                                            // 24B
    usize genome_start_pos;                                                   // 8B

    // Normalizing a Multiallelic block demands slicing bounding boxes universally based on ALL 
    // participating alleles simultaneously. You cannot over-trim the REF if one strict ALT block
    // needs bounding anchor integrity (e.g. an Indel becoming empty).
    void NormalizeVcfParsimony() {
        if (alt_alleles_to_haps.empty() || ref_allele.empty()) return;

        // Right Trim (eg. REF: "ATCG", ALTS: ["ACCG", "AGGG"] => "ATC", ["ACC", "AGG"])
        // Erases rightward bloat first allowing indels to left-align against leftward structural boundaries. 
        ApplyUnifiedTrim(
            [](std::string_view r, std::string_view a) { return a.length() > 1 && a.back() == r.back(); },
            [](std::string& t) { t.pop_back(); }
        );

        // Left Trim (e.g. REF: "TTC", ALTS: ["TGC"] => "TC", ["GC"])
        const usize initial_ref_len = ref_allele.length();
        ApplyUnifiedTrim(
            [](std::string_view r, std::string_view a) { return a.length() > 1 && a.front() == r.front(); },
            [](std::string& t) { t.erase(0, 1); }
        );
        
        // Correcting Start Pos! Every left character dropped dynamically sweeps Genomic Start `+1`.
        genome_start_pos += (initial_ref_len - ref_allele.length()); 
    }

private:
    // Lambda generically isolates unified string validation.
    // Takes `can_trim` (checks boundary char identically across all alleles utilizing 
    // `string_view` for speed) and `do_trim` (executes the heavy string mutation).
    template <typename CanTrimFunc, typename DoTrimFunc>
    void ApplyUnifiedTrim(CanTrimFunc can_trim, DoTrimFunc do_trim) {
        while (ref_allele.length() > 1) { // Guard ensuring bounding box doesn't evaporate completely!
            bool unanimously_shared = true;
            std::string_view ref_view(ref_allele);
            for (const auto& [alt_seq, _] : alt_alleles_to_haps) {
                if (!can_trim(ref_view, std::string_view(alt_seq))) {
                    unanimously_shared = false;
                    break; // Any deviation stops the trimming completely across ALL alleles globally
                }
            }
            if (!unanimously_shared) break;
            
            // Once verified unconditionally, physically mutate all tracking memories uniformly 
            do_trim(ref_allele);
            
            // Flat Hash Maps utilize const `Keys` protecting structural hash-index legitimacy. 
            // We cannot arbitrarily mutate strings inside keys in-place. We allocate a blank
            // hashmap natively and populate it exploiting `std::move` to skip memory reallocations
            // of the heavy `std::vector<usize>` properties! 
            absl::flat_hash_map<std::string, std::vector<usize>> rehashed_map;
            for (auto& [alt_seq, haps] : alt_alleles_to_haps) {
                std::string new_alt = alt_seq;
                do_trim(new_alt);
                rehashed_map[std::move(new_alt)] = std::move(haps);
            }
            alt_alleles_to_haps = std::move(rehashed_map);
        }
    }
};

// =========================================================================================
// VariantExtractor Class Object (Finite State Machine Array Tracking)
// -----------------------------------------------------------------------------------------
// WHY A CLASS? The looping logic relies profoundly heavily on 6+ sweeping memory tracking
// variables (`active_ptrs_`, `node_to_rank_`, `current_ref_pos_`). Emulating these dynamically
// locally inside an ultra massive routine produces extreme repetitive clutter.
// Encapsulating state cleanly minimizes internal functions into robust English pseudo-code.
// =========================================================================================
class VariantExtractor {
public:
    VariantExtractor(const spoa::Graph& graph, const core::Window& win, usize anchor_start)
        : graph_(graph), win_(win), ref_anchor_start_(anchor_start), current_ref_pos_(anchor_start) {
        
        num_seqs_ = graph_.sequences().size();
        if (num_seqs_ < 2) return; // Natively nothing to discover against isolated references
        
        current_hap_pos_.assign(num_seqs_, 0); // Initializes array accurately bounding paths

        // ===========================================================================
        // RANK LOOKUP INITIALIZATION:  O(N) Inverse Topological Indexing
        // ---------------------------------------------------------------------------
        // WHY? In `spoa`, a node's physical `id` is dynamically assigned based exclusively
        // on chronological graph insertion. A downstream reference base might be Node #5,
        // while an upstream variant inserted significantly later could miraculously be Node #500.
        // Utilizing native `.id` to evaluate chronological progression generates catastrophic bugs.
        // 
        // We utilize the `spoa` strictly validated 5'-to-3' topological sorting vector array.
        // We inverse it here, meaning if you pass ANY native `.id`, `node_to_rank_[id]`
        // retrieves its true mathematical biological left-to-right progressive rank value `O(1)`.
        // ===========================================================================
        const auto& topological_order = graph_.rank_to_node();
        // Instantiate graph memory maps matching natively biological pathways.
        active_ptrs_.assign(num_seqs_, nullptr);
        node_to_rank_.assign(graph_.nodes().size(), std::numeric_limits<u32>::max());
        for (u32 rank = 0; rank < topological_order.size(); ++rank) {
            node_to_rank_[topological_order[rank]->id] = rank;
        }

        // Cache initial coordinate origins representing all continuous biological pathways 
        for (usize i = 0; i < num_seqs_; ++i) active_ptrs_[i] = graph_.sequences()[i];
    }

    // ===================================================================================================
    // VARIANT EXTRACTOR FLOWCHART: Biological String Unraveling Topology
    // ===================================================================================================
    // The VariantExtractor continuously pushes active pointers 5'-to-3' topologically across the POA graph.
    // When pathways diverge (a variant occurs), the Extractor halts uniform sweeping and triggers a Bubble.
    //
    //                              (C)  (Prior Match Node! All Pointers Converged Here)
    //                            /     \
    //                  (ALT)    /       \   (REF)
    //                 (C)->(T) .         . (T)->(G)
    //                           \       /
    //                            \     /
    //                              (G)      (Target Convergence Node! Paths reunite here)
    // 
    // HOW ARE PRIVATE HELPERS UTILIZED NATIVELY?
    // 1. `InitializeBubbleAnchor`    : Before diving into the divergence, it retroactively graphs the prior
    //                                  Match Node base `(C)` uniformly appending it to structural sequences.
    // 2. `SinkPointers`              : The Hot Loop Engine! It identifies the computationally lowest active 
    //    |-- `FindLowestActiveRank`  : topological Rank across all diverging paths.
    //    |-- `ConsumePathsAtRank`    : It advances strictly those bottlenecking paths sequentially, looping 
    //                                  the matrix until ALL pointers jump perfectly simultaneously to `(G)`.
    // 3. `CreateNormalizedBubble`    : Resolves the aggregated string bundles and routes them structurally 
    //                                  into `VariantBubble::NormalizeVcfParsimony` permanently removing bloat!
    // 4. `AssembleMultiallelicVariant`: Evaluates each isolated pure sequence core structurally by calling 
    //                                  `RawVariant::ClassifyVariant` and emitting `RawVariant` sets seamlessly!
    // ===================================================================================================
    // Single public interface endpoint parsing variants continuously directly into the tracking payload container
    void SearchAndExtractTo(absl::btree_set<RawVariant>& out_variants) {
        if (num_seqs_ < 2) return;

        while (true) {
            if (AreAllPathsConverged()) {
                // If everyone identically evaluates native `nullptr`, biological tracking universally concludes!
                if (active_ptrs_[REF_HAP_IDX] == nullptr) break; 
                AdvanceConvergedPaths();
            } else {
                EatTopologicalBubble(out_variants);
            }
        }
    }

private:
    // Memory Alignment: Large containers (24B vector objects) top-loaded, followed by
    // uniformly 8B sized references, integers, and raw pointers perfectly eliminating waste.
    std::vector<u32> node_to_rank_;
    std::vector<const spoa::Graph::Node*> active_ptrs_;
    std::vector<usize> current_hap_pos_; // Dynamically tracks each active path's exact physical array offset
    const spoa::Graph& graph_;
    const core::Window& win_;
    usize ref_anchor_start_;
    usize current_ref_pos_;
    usize num_seqs_;
    // Mandatory dynamically updated VCF Anchor boundary
    const spoa::Graph::Node* prev_match_node_ = nullptr;

    // O(M_paths) inline evaluating if NO multiallelic divergence exists currently locally
    inline bool AreAllPathsConverged() const {
        const auto* target = active_ptrs_[REF_HAP_IDX];
        for (usize i = 1; i < num_seqs_; ++i) {
            if (active_ptrs_[i] != target) return false;
        }
        return true;
    }

    // Step continuous matched blocks homogeneously concurrently 
    inline void AdvanceConvergedPaths() {
        prev_match_node_ = active_ptrs_[REF_HAP_IDX];
        for (usize i = 0; i < num_seqs_; ++i) {
            if (active_ptrs_[i]) {
                active_ptrs_[i] = active_ptrs_[i]->Successor(i);
                current_hap_pos_[i]++; // Dynamically advance sequence coordinate
            }
        }
        current_ref_pos_++;
    }


    // Controller natively consuming open bubbles and emitting clean Multiallelic Set payloads!
    inline void EatTopologicalBubble(absl::btree_set<RawVariant>& out_variants) {
        std::vector<std::string> raw_alleles(num_seqs_, ""); // Initialize empty base sequences uniformly
        std::vector<usize> bubble_hap_starts(num_seqs_, 0);
        
        // 1. Anchor 
        usize exact_start_pos = InitializeBubbleAnchor(absl::MakeSpan(raw_alleles), bubble_hap_starts);

        // 2. Engage Hot Loop Vector Engine to sweep topologies
        SinkPointers(absl::MakeSpan(raw_alleles));

        // 3. Normalize into a clean VCF mathematical block
        VariantBubble bubble = CreateNormalizedBubble(exact_start_pos, std::move(raw_alleles));
        bubble.hap_starts = std::move(bubble_hap_starts);

        if (!bubble.alt_alleles_to_haps.empty()) {
            out_variants.insert(AssembleMultiallelicVariant(std::move(bubble)));
        }
    }

    // BIOLOGICAL VCF ANCHORING REQUIREMENT:
    // Complex Indels universally mandate a shared prefixed Match anchor natively appended uniformly.
    inline auto InitializeBubbleAnchor(absl::Span<std::string> raw_alleles, std::vector<usize>& out_hap_starts) -> usize {
        const bool has_prev = (prev_match_node_ != nullptr);
        const usize anchor_offset = static_cast<usize>(has_prev);

        // Branchlessly shifts the universal genomic coordinate logically backward by exactly 1 if an anchor exists
        const usize bubble_start_pos = current_ref_pos_ - anchor_offset;

        // We retroactively extract the LAST confirmed unified anchor base to prefix the sequences!
        if (has_prev) {
            const char decoder_val = static_cast<char>(graph_.decoder(prev_match_node_->code));
            for (auto& al : raw_alleles) al += decoder_val;
        }

        // Branchlessly record matrix boundaries uniformly
        for (usize i = 0; i < num_seqs_; ++i) {
            out_hap_starts[i] = current_hap_pos_[i] - anchor_offset;
        }

        return bubble_start_pos;
    }

    // ===========================================================================
    // TOPOLOGICAL SINK: The Biological Sweeping Engine 
    // Using `absl::Span` ensures `raw_alleles` memory modifications are written precisely 
    // natively back to the source string vectors without executing vector copying internally.
    // ===========================================================================
    // Pushes active pointers topologically forward until universal convergence is re-established natively
    inline void SinkPointers(absl::Span<std::string> raw_alleles) {
        while (!AreAllPathsConverged()) {
            const u32 min_rank = FindLowestActiveRank();
            if (min_rank == std::numeric_limits<u32>::max()) break;
            ConsumePathsAtRank(min_rank, raw_alleles);
        }
    }

    // Phase 1: Determine the lowest active topological boundary amongst traversing sets
    inline auto FindLowestActiveRank() const -> u32 {
        u32 min_rank = std::numeric_limits<u32>::max();
        for (const auto* p : active_ptrs_) {
            if (p != nullptr) min_rank = std::min(min_rank, node_to_rank_.at(p->id));
        }
        return min_rank;
    }

    // Phase 2: Selectively consume and roll exclusively paths pinned precisely against that lowest rank
    inline void ConsumePathsAtRank(u32 target_rank, absl::Span<std::string> raw_alleles) {
        for (usize i = 0; i < num_seqs_; ++i) {
            if (active_ptrs_[i] != nullptr && node_to_rank_.at(active_ptrs_[i]->id) == target_rank) {
                raw_alleles[i] += static_cast<char>(graph_.decoder(active_ptrs_[i]->code));
                active_ptrs_[i] = active_ptrs_[i]->Successor(i);
                
                current_hap_pos_[i]++; // Tracking local path lengths natively mathematically 
                
                // Critical Update! Only the native biological Reference Index modifies universal REF coordinates.
                if (i == REF_HAP_IDX) current_ref_pos_++;
            }
        }
    }

    // Phase 3: Cluster disjoint sequence strings upon convergence organically, and precisely purge parsimony padding.
    inline auto CreateNormalizedBubble(usize genome_start_pos, std::vector<std::string> raw_alleles) -> VariantBubble {
        VariantBubble bubble;
        bubble.genome_start_pos = genome_start_pos;
        bubble.ref_allele = std::move(raw_alleles[REF_HAP_IDX]);
        
        // Exclude the pure identical reference tracks. Everything else goes into the variant parser
        for (usize alt_idx = 1; alt_idx < num_seqs_; ++alt_idx) {
            if (raw_alleles[alt_idx] != bubble.ref_allele) { 
                bubble.alt_alleles_to_haps[raw_alleles[alt_idx]].push_back(alt_idx);
            }
        }
        
        bubble.NormalizeVcfParsimony();
        return bubble;
    }

    // Phase 4: Mathematically squeeze sequence cores against variant classes and emit normalized VCF multiallelic payloads.
    inline auto AssembleMultiallelicVariant(VariantBubble bubble) -> RawVariant {
        // Instantiate unified multiallelic VCF bucket correctly capturing complex fields.
        RawVariant multiallelic_var;
        multiallelic_var.mChromIndex = win_.ChromIndex();
        multiallelic_var.mChromName = win_.ChromName();
        multiallelic_var.mGenomeChromPos1 = bubble.genome_start_pos;
        multiallelic_var.mLocalRefStart0Idx = bubble.hap_starts[REF_HAP_IDX]; // Statically locks to REF coordinate flawlessly 
        multiallelic_var.mRefAllele = std::move(bubble.ref_allele);

        // Populate disjoint variant structures
        for (auto& [normalized_alt, haps] : bubble.alt_alleles_to_haps) {
            RawVariant::AltAllele sub_alt;
            sub_alt.mSequence = std::move(normalized_alt);
            sub_alt.mType = RawVariant::ClassifyVariant(multiallelic_var.mRefAllele, sub_alt.mSequence);
            sub_alt.mLength = CalculateVariantLength(multiallelic_var.mRefAllele, sub_alt.mSequence, sub_alt.mType);

            for (const usize hap_id : haps) {
                sub_alt.mLocalHapStart0Idxs.emplace(hap_id, bubble.hap_starts[hap_id]);
            }
            
            multiallelic_var.mAlts.push_back(std::move(sub_alt));
        }

        // `RawVariant` utilizes std::vector equality natively, so sequence elements MUST
        // strictly align uniformly to generate completely identical hashes reliably inside the set arrays later on.
        std::sort(multiallelic_var.mAlts.begin(), multiallelic_var.mAlts.end());
        return multiallelic_var;
    }
};

void VariantSet::ExtractVariantsFromGraph(const spoa::Graph& graph, const core::Window& win, usize ref_anchor_start) {
    if (graph.sequences().size() < 2) return;
    VariantExtractor extractor(graph, win, ref_anchor_start);
    extractor.SearchAndExtractTo(this->mResultVariants);
}

} // namespace lancet::caller
