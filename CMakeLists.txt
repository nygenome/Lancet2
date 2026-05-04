# ═══════════════════════════════════════════════════════════════════════════════
# Lancet2 — Micro-assembly based variant caller
#
# Build dependency chain (each layer depends on everything above it):
#   lancet_base  → lancet_hts → lancet_cbdg → lancet_caller → lancet_core → lancet_cli → Lancet2
#
# Within each library, sources are ordered:
#   1. Header-only files (types, constants, data structures)
#   2. Implementation pairs (.cpp + .h), from foundational → orchestrator
# ═══════════════════════════════════════════════════════════════════════════════
message(STATUS "Current CMake version: ${CMAKE_VERSION}")
cmake_minimum_required(VERSION 3.25 FATAL_ERROR)

file(STRINGS "${CMAKE_CURRENT_SOURCE_DIR}/VERSION.txt" LANCET_VERSION_FROM_FILE LIMIT_COUNT 1)
string(STRIP "${LANCET_VERSION_FROM_FILE}" LANCET_VERSION_FROM_FILE)
project(lancet VERSION ${LANCET_VERSION_FROM_FILE} LANGUAGES C CXX
		HOMEPAGE_URL https://github.com/nygenome/Lancet2
		DESCRIPTION "Micro-assembly based variant caller")

# ── Version & Options ─────────────────────────────────────────────────────────
include(cmake/version.cmake)
include(CMakeDependentOption)
message(STATUS "Lancet v${PROJECT_VERSION}")
option(LANCET_TESTS "Builds the Lancet tests executable" "OFF")
option(LANCET_BENCHMARKS "Builds the Lancet benchmarks executable" "OFF")
option(LANCET_PROFILE_MODE "Enables CPU profiling for the Lancet executable" "OFF")
option(LANCET_BUILD_STATIC "Builds a statically linked Lancet executable" "ON")
option(LANCET_NATIVE_BUILD "Optimizes for host CPU. Faster but may not run on different machine" "OFF")
option(LANCET_VERBOSE_LOGS "Enables TRACE level logging with --verbose flag" "OFF")
option(LANCET_SANITIZE_BUILD "Builds with an LLVM sanitizer (ASan/MSan/TSan/UBSan)" "OFF")
option(LANCET_QUIET_DEPS "Suppresses cmake configure output for third-party libraries" "ON")
cmake_dependent_option(LANCET_ENABLE_CLOUD_IO "Enables native S3 and GCS streaming via libcurl" OFF "NOT LANCET_BUILD_STATIC" OFF)

set_property(CACHE LANCET_TESTS PROPERTY STRINGS "OFF" "ON")
set_property(CACHE LANCET_BENCHMARKS PROPERTY STRINGS "OFF" "ON")
set_property(CACHE LANCET_PROFILE_MODE PROPERTY STRINGS "OFF" "ON")
set_property(CACHE LANCET_BUILD_STATIC PROPERTY STRINGS "OFF" "ON")
set_property(CACHE LANCET_NATIVE_BUILD PROPERTY STRINGS "OFF" "ON")
set_property(CACHE LANCET_VERBOSE_LOGS PROPERTY STRINGS "OFF" "ON")
set_property(CACHE LANCET_SANITIZE_BUILD PROPERTY STRINGS "OFF" "ON")
set_property(CACHE LANCET_QUIET_DEPS PROPERTY STRINGS "OFF" "ON")
if (NOT LANCET_BUILD_STATIC)
	set_property(CACHE LANCET_ENABLE_CLOUD_IO PROPERTY STRINGS "OFF" "ON")
endif ()

# ── Language Standard ─────────────────────────────────────────────────────────
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ── Build Configuration ───────────────────────────────────────────────────────
message(STATUS "CMAKE_HOST_SYSTEM_NAME: ${CMAKE_HOST_SYSTEM_NAME}")
message(STATUS "CMAKE_HOST_SYSTEM_PROCESSOR: ${CMAKE_HOST_SYSTEM_PROCESSOR}")
MESSAGE(STATUS "CMAKE_C_COMPILER: ${CMAKE_C_COMPILER}")
MESSAGE(STATUS "CMAKE_CXX_COMPILER: ${CMAKE_CXX_COMPILER}")

include(cmake/platform_checks.cmake)
include(cmake/compiler_flags.cmake)

# ── System Dependencies ───────────────────────────────────────────────────────
# BZip2 and LibLZMA are required by HTSlib for compressed I/O.
# CURL and OpenSSL are optional, enabled via -DLANCET_ENABLE_CLOUD_IO=ON
# for native S3/GCS streaming in dynamic-linked builds.
find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)

if (LANCET_ENABLE_CLOUD_IO)
	find_package(CURL REQUIRED)
	find_package(OpenSSL REQUIRED)
endif ()

# ── Vendored Dependencies ────────────────────────────────────────────────────
include(cmake/dependencies.cmake)

# ═══════════════════════════════════════════════════════════════════════════════
# Layer 1: lancet_base — types, math utilities, sequence analysis primitives
# ═══════════════════════════════════════════════════════════════════════════════
add_library(lancet_base STATIC
		# ── Header-only: type aliases, macros, utilities ──────────────────
		src/lancet/base/types.h
		src/lancet/base/assert.h
		src/lancet/base/version.h
		src/lancet/base/logging.h
		src/lancet/base/timer.h
		src/lancet/base/memory.h
		src/lancet/base/rev_comp.h
		src/lancet/base/compute_stats.h
		src/lancet/base/sliding.h
		src/lancet/base/polar_coords.h
		src/lancet/base/longdust_scorer.h
		# ── Implementation pairs: foundational → derived ──────────────────
		src/lancet/base/hash.cpp src/lancet/base/hash.h
		src/lancet/base/repeat.cpp src/lancet/base/repeat.h
		src/lancet/base/eta_timer.cpp src/lancet/base/eta_timer.h
		src/lancet/base/crash_handler.cpp src/lancet/base/crash_handler.h
		src/lancet/base/sequence_complexity.cpp src/lancet/base/sequence_complexity.h
		src/lancet/base/gzip_ostream.cpp src/lancet/base/gzip_ostream.h
		src/lancet/base/tar_gz_writer.cpp src/lancet/base/tar_gz_writer.h)
target_link_libraries(lancet_base PRIVATE absl::flat_hash_set zlibstatic
		PUBLIC spdlog::spdlog absl::span absl::fixed_array absl::strings absl::time)
target_include_directories(lancet_base PUBLIC "${CMAKE_SOURCE_DIR}/src" "${CMAKE_BINARY_DIR}/generated")
set_target_properties(lancet_base PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

# ═══════════════════════════════════════════════════════════════════════════════
# Layer 2: lancet_hts — HTSlib/BAM I/O, read extraction, CIGAR parsing
#
# Data flow:
#
#   ┌─────────────────────────────────┐
#   │ cigar_unit, cigar_utils, aux_tag│  BAM field types
#   └────────────────┬────────────────┘
#                    ▼
#   ┌────────────────────────────────┐
#   │            reference           │  FASTA reader
#   └────────────────┬───────────────┘
#                    ▼
#   ┌────────────────────────────────┐
#   │      sam_flag → alignment      │  read representation
#   └────────────────┬───────────────┘
#                    ▼
#   ┌────────────────────────────────┐
#   │      iterator → extractor      │  BAM region traversal
#   └────────────────────────────────┘
# ═══════════════════════════════════════════════════════════════════════════════
add_library(lancet_hts STATIC
		# ── Header-only: BAM field types ──────────────────────────────────
		src/lancet/hts/cigar_unit.h
		src/lancet/hts/cigar_utils.h
		src/lancet/hts/aux_tag.h
		src/lancet/hts/mate_info.h
		# ── Implementation pairs: reference → alignment → extraction ─────
		src/lancet/hts/bgzf_ostream.cpp src/lancet/hts/bgzf_ostream.h
		src/lancet/hts/phred_quality.cpp src/lancet/hts/phred_quality.h
		src/lancet/hts/reference.cpp src/lancet/hts/reference.h
		src/lancet/hts/sam_flag.cpp src/lancet/hts/sam_flag.h
		src/lancet/hts/alignment.cpp src/lancet/hts/alignment.h
		src/lancet/hts/iterator.cpp src/lancet/hts/iterator.h
		src/lancet/hts/extractor.cpp src/lancet/hts/extractor.h
		src/lancet/hts/uri_utils.cpp src/lancet/hts/uri_utils.h)
add_dependencies(lancet_hts htslib)
set_target_properties(lancet_hts PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
target_include_directories(lancet_hts SYSTEM PUBLIC ${HTSLIB_ROOT_DIR})
target_include_directories(lancet_hts PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(lancet_hts PUBLIC lancet_base absl::statusor absl::fixed_array absl::flat_hash_set ${LIB_HTS}
		PRIVATE absl::flat_hash_map INTERFACE BZip2::BZip2 LibLZMA::LibLZMA zlibstatic libdeflate_static)

if (LANCET_ENABLE_CLOUD_IO)
	target_link_libraries(lancet_hts INTERFACE CURL::libcurl OpenSSL::Crypto)
endif ()

# ═══════════════════════════════════════════════════════════════════════════════
# Layer 3: lancet_cbdg — Colored de Bruijn Graph assembly engine
#
# Data flow:
#
#   ┌──────────────────────────────────────────────────────┐
#   │                 label, graph_params                  │  configuration
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                 sample_mask → kmer                   │  k-mer encoding + sample tagging
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                 node → edge → path                   │  graph topology (path + walk)
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │              cycle_finder, traversal_index           │  graph algorithms
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │               max_flow, graph_complexity             │  flow + complexity metrics
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                  component_result                    │  per-component assembly output
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │   probe_index → probe_tracker → probe_results_writer │  truth variant k-mer tracking
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │ dot_layers / dot_walk_layers / dot_overlay_factories │  layered DOT renderer
#   │    dot_plan / dot_renderer / dot_snapshot_buffer     │  (deferred output via buffer)
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                        graph                         │  assembly orchestrator
#   └──────────────────────────────────────────────────────┘
# ═══════════════════════════════════════════════════════════════════════════════
add_library(lancet_cbdg STATIC
		# ── Header-only: graph configuration ──────────────────────────────
		src/lancet/cbdg/label.h
		src/lancet/cbdg/graph_params.h
		src/lancet/cbdg/edge.h
		src/lancet/cbdg/read.h
		# ── Implementation pairs: kmer → node → path → algorithms → graph ─
		src/lancet/cbdg/sample_mask.cpp src/lancet/cbdg/sample_mask.h
		src/lancet/cbdg/kmer.cpp src/lancet/cbdg/kmer.h
		src/lancet/cbdg/node.cpp src/lancet/cbdg/node.h
		src/lancet/cbdg/path.cpp src/lancet/cbdg/path.h
		src/lancet/cbdg/cycle_finder.cpp src/lancet/cbdg/cycle_finder.h
		src/lancet/cbdg/traversal_index.cpp src/lancet/cbdg/traversal_index.h
		src/lancet/cbdg/max_flow.cpp src/lancet/cbdg/max_flow.h
		src/lancet/cbdg/graph_complexity.cpp src/lancet/cbdg/graph_complexity.h
		# ── Per-component assembly output: paths + walks + complexity ─────
		src/lancet/cbdg/component_result.cpp src/lancet/cbdg/component_result.h
		# ── Probe diagnostics: index → tracker → results writer ───────────
		src/lancet/cbdg/probe_index.cpp src/lancet/cbdg/probe_index.h
		src/lancet/cbdg/probe_results_writer.cpp src/lancet/cbdg/probe_results_writer.h
		src/lancet/cbdg/probe_tracker.cpp src/lancet/cbdg/probe_tracker.h
		# ── DOT visualization: layered overlays + deferred per-k buffer ───
		src/lancet/cbdg/dot_layers.h
		src/lancet/cbdg/dot_plan.h
		src/lancet/cbdg/dot_walk_layers.h
		src/lancet/cbdg/dot_overlay_factories.h
		src/lancet/cbdg/dot_renderer.cpp src/lancet/cbdg/dot_renderer.h
		src/lancet/cbdg/dot_snapshot_buffer.cpp src/lancet/cbdg/dot_snapshot_buffer.h
		# ── Assembly orchestrator ─────────────────────────────────────────
		src/lancet/cbdg/graph.cpp src/lancet/cbdg/graph.h)
target_include_directories(lancet_cbdg PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(lancet_cbdg PUBLIC lancet_hts absl::cord absl::inlined_vector absl::flat_hash_map absl::synchronization)
set_target_properties(lancet_cbdg PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

# ═══════════════════════════════════════════════════════════════════════════════
# Layer 4: lancet_caller — variant calling: MSA, genotyping, VCF construction
#
# Data flow:
#
#   ┌──────────────────────────────────────────────────────┐
#   │               alt_allele, raw_variant                │  data structures
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │          variant_bubble, variant_extractor           │  SPOA DAG bubble detection
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                     variant_set                      │  extraction orchestrator
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │            local_scorer → combined_scorer            │  alignment scoring
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │         genotype_likelihood, posterior_bq            │  statistical models
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │           variant_support, support_array             │  per-allele evidence
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │          sample_format_data, vcf_formatter           │  per-sample FORMAT fields
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                     variant_call                     │  VCF record assembly
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                      genotyper                       │  read-to-haplotype orchestrator
#   └──────────────────────────────────────────────────────┘
# ═══════════════════════════════════════════════════════════════════════════════
add_library(lancet_caller STATIC
		# ── Header-only: constants, enums, data structures ────────────────
		src/lancet/caller/allele_scoring_types.h
		src/lancet/caller/alt_allele.h
		src/lancet/caller/per_allele_data.h
		src/lancet/caller/scoring_constants.h
		src/lancet/caller/vcf_formatter.h
		# ── MSA: SPOA graph construction ──────────────────────────────────
		src/lancet/caller/msa_builder.cpp src/lancet/caller/msa_builder.h
		# ── Variant discovery: data → bubble → extractor → set ────────────
		src/lancet/caller/raw_variant.cpp src/lancet/caller/raw_variant.h
		src/lancet/caller/variant_bubble.cpp src/lancet/caller/variant_bubble.h
		src/lancet/caller/variant_extractor.cpp src/lancet/caller/variant_extractor.h
		src/lancet/caller/variant_set.cpp src/lancet/caller/variant_set.h
		# ── Alignment scoring: local → combined ───────────────────────────
		src/lancet/caller/local_scorer.cpp src/lancet/caller/local_scorer.h
		src/lancet/caller/combined_scorer.cpp src/lancet/caller/combined_scorer.h
		# ── Statistical models: genotype likelihoods + base quality ───────
		src/lancet/caller/genotype_likelihood.cpp src/lancet/caller/genotype_likelihood.h
		src/lancet/caller/posterior_base_qual.cpp src/lancet/caller/posterior_base_qual.h
		# ── Evidence aggregation: per-allele support → sample FORMAT ──────
		src/lancet/caller/variant_support.cpp src/lancet/caller/variant_support.h
		src/lancet/caller/support_array.cpp src/lancet/caller/support_array.h
		src/lancet/caller/sample_format_data.cpp src/lancet/caller/sample_format_data.h
		# ── VCF output: variant call assembly + genotyper orchestrator ────
		src/lancet/caller/variant_call.cpp src/lancet/caller/variant_call.h
		src/lancet/caller/genotyper.cpp src/lancet/caller/genotyper.h)
add_dependencies(lancet_caller minimap2)
target_include_directories(lancet_caller PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_include_directories(lancet_caller SYSTEM PUBLIC "${MM2_ROOT_DIR}")
target_link_libraries(lancet_caller PUBLIC lancet_cbdg spoa absl::btree ${LIB_MM2} PRIVATE zlibstatic)
set_target_properties(lancet_caller PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

# ═══════════════════════════════════════════════════════════════════════════════
# Layer 5: lancet_core — windowing, read collection, variant building, pipeline execution
#
# Data flow:
#
#   ┌──────────────────────────────────────────────────────┐
#   │                  window, sample_info                 │  genomic region + sample metadata
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                  input_spec_parser                   │  --normal/--tumor/--sample parsing
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                sample_header_reader                  │  BAM header SM extraction
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                      bed_parser                      │  BED file streaming (HTSlib)
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                    window_builder                    │  region partitioning
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                active_region_detector                │  mutation evidence pre-scan
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                    read_collector                    │  BAM read extraction per window
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                   probe_diagnostics                  │  truth variant MSA analysis
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                   variant_builder                    │  assembly + calling per window
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │         variant_annotator → variant_store            │  annotation + dedup + VCF output
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                     async_worker                     │  thread pool executor
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                  pipeline_executor                   │  batch-fed thread pool orchestration
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                  tar_gz_shard_merger                 │  cat-merge per-worker tar.gz shards
#   └──────────────────────────────────────────────────────┘
# ═══════════════════════════════════════════════════════════════════════════════
add_library(lancet_core STATIC
		# ── Header-only: window + sample metadata ─────────────────────────
		src/lancet/core/window.h
		src/lancet/core/sample_info.h
		# ── Implementation pairs: parsing → windowing → reads → build → diagnostics → store → execute → merge ─
		src/lancet/core/input_spec_parser.cpp src/lancet/core/input_spec_parser.h
		src/lancet/core/sample_header_reader.cpp src/lancet/core/sample_header_reader.h
		src/lancet/core/bed_parser.cpp src/lancet/core/bed_parser.h
		src/lancet/core/active_region_detector.cpp src/lancet/core/active_region_detector.h
		src/lancet/core/window_builder.cpp src/lancet/core/window_builder.h
		src/lancet/core/read_collector.cpp src/lancet/core/read_collector.h
		src/lancet/core/probe_diagnostics.cpp src/lancet/core/probe_diagnostics.h
		src/lancet/core/variant_builder.cpp src/lancet/core/variant_builder.h
		src/lancet/core/variant_annotator.cpp src/lancet/core/variant_annotator.h
		src/lancet/core/variant_store.cpp src/lancet/core/variant_store.h
		src/lancet/core/async_worker.cpp src/lancet/core/async_worker.h
		src/lancet/core/pipeline_executor.cpp src/lancet/core/pipeline_executor.h
		src/lancet/core/tar_gz_shard_merger.cpp src/lancet/core/tar_gz_shard_merger.h)
target_include_directories(lancet_core PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(lancet_core PUBLIC lancet_caller absl::synchronization concurrentqueue PRIVATE absl::hash)
set_target_properties(lancet_core PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

# ═══════════════════════════════════════════════════════════════════════════════
# Layer 6: lancet_cli — command-line interface, VCF header construction, pipeline entry
#
# Data flow:
#
#   ┌──────────────────────────────────────────────────────┐
#   │                     cli_params                       │  user-facing option definitions
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                  vcf_header_builder                  │  VCF v4.5 header construction
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                    pipeline_runner                   │  validate → I/O → delegate
#   └──────────────────────────┬───────────────────────────┘
#                              ▼
#   ┌──────────────────────────────────────────────────────┐
#   │                    cli_interface                     │  CLI11 argument parser + validators
#   └──────────────────────────────────────────────────────┘
# ═══════════════════════════════════════════════════════════════════════════════
add_library(lancet_cli STATIC
		# ── Header-only: CLI parameter definitions ────────────────────────
		src/lancet/cli/cli_params.h
		# ── Implementation pairs: VCF header → pipeline runner → CLI parser ─
		src/lancet/cli/vcf_header_builder.cpp src/lancet/cli/vcf_header_builder.h
		src/lancet/cli/pipeline_runner.cpp src/lancet/cli/pipeline_runner.h
		src/lancet/cli/cli_interface.cpp src/lancet/cli/cli_interface.h)
target_include_directories(lancet_cli PUBLIC "${CMAKE_SOURCE_DIR}/src")
target_link_libraries(lancet_cli PUBLIC lancet_core CLI11)
set_target_properties(lancet_cli PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)

# ═══════════════════════════════════════════════════════════════════════════════
# Lancet2 — main executable
# ═══════════════════════════════════════════════════════════════════════════════
add_executable(Lancet2 src/lancet/main.cpp)
target_include_directories(Lancet2 PRIVATE "${CMAKE_SOURCE_DIR}")

# ── Sanitizer support ─────────────────────────────────────────────────────────
if (LANCET_SANITIZE_BUILD)
	target_link_libraries(Lancet2 PRIVATE lancet_cli absl::cleanup)
	target_compile_definitions(Lancet2 PRIVATE LANCET_SANITIZE_BUILD=1)
else ()
	target_link_libraries(Lancet2 PRIVATE lancet_cli mimalloc-static absl::cleanup)
endif ()
# ──────────────────────────────────────────────────────────────────────────────

set_target_properties(Lancet2 PROPERTIES CXX_STANDARD 20 CXX_STANDARD_REQUIRED ON)
set_property(TARGET Lancet2 PROPERTY $<$<ENABLE_LTO>:INTERPROCEDURAL_OPTIMIZATION TRUE>)
install(TARGETS Lancet2 DESTINATION bin)

if (LANCET_VERBOSE_LOGS)
	add_compile_definitions(LANCET_VERBOSE_LOGS)
	message(STATUS "Lancet2: Verbose logs ON")
endif ()

if (CMAKE_BUILD_TYPE MATCHES Debug)
	add_compile_definitions(LANCET_DEBUG_MODE)
	message(STATUS "Lancet2: Debug mode ON")
endif ()

# ── Optional: CPU profiling (RelWithDebInfo only) ─────────────────────────────
# platform_checks.cmake enforces RelWithDebInfo when LANCET_PROFILE_MODE=ON.
# https://github.com/google/pprof
# https://gperftools.github.io/gperftools/cpuprofile.html
if (LANCET_PROFILE_MODE AND ${CMAKE_BUILD_TYPE} MATCHES RelWithDebInfo)
	add_compile_definitions(LANCET_PROFILE_MODE)
	message(STATUS "Lancet2: Profile mode ON")
	add_dependencies(lancet_core gperftools)
	add_dependencies(lancet_cli gperftools)
	target_include_directories(lancet_core SYSTEM PRIVATE ${GPERFTOOLS_INC_DIR})
	target_include_directories(lancet_cli SYSTEM PRIVATE ${GPERFTOOLS_INC_DIR})
	target_link_libraries(lancet_core PRIVATE ${LIB_PROFILER})
	target_link_libraries(lancet_cli PRIVATE ${LIB_PROFILER})
endif ()

# ── Optional: Tests & Benchmarks ──────────────────────────────────────────────
# https://github.com/catchorg/Catch2/blob/v3.14.0/docs/Readme.md
if (LANCET_TESTS)
	add_subdirectory(tests)
endif ()

# https://github.com/google/benchmark/blob/main/docs/user_guide.md
if (LANCET_BENCHMARKS)
	add_subdirectory(benchmarks)
endif ()
