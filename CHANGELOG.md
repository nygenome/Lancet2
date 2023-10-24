
<a name="v2.6.0"></a>
## [v2.6.0](https://github.com/nygenome/Lancet2/compare/v2.5.0...v2.6.0) (2023-10-24)

### Bug Fixes

* ensure min chunk length for fuzzy matching
* use 90% of chunk length for fuzzy match cutoff
* count longer indels better with using min of var and read len

### Maintenance/Refactoring

* remove extra commands in version extract


<a name="v2.5.0"></a>
## [v2.5.0](https://github.com/nygenome/Lancet2/compare/v2.4.0...v2.5.0) (2023-10-23)

### Bug Fixes

* reserve vector length before use
* return untouched ref and qry
* resolve typo in logic

### Maintenance/Refactoring

* Bump version to v2.5.0
* ensure we don't go over max kmer len
* update deps

### New Features

* add fuzzy matching for genotyping long indels


<a name="v2.4.0"></a>
## [v2.4.0](https://github.com/nygenome/Lancet2/compare/v2.3.0...v2.4.0) (2023-10-19)

### Bug Fixes

* update vcf header text
* simplify and remove max normal filter
* resolve typo in if statement
* reduce flank to 3 exact match
* remove str specific filter
* check for non gap alignment of flanks instead of match
* reduce AS/XS diff to 1 percent
* reduce AS/XS diff to 5 percent
* reduce AS/XS diff to 10 percent
* increase min difference between AS and XS
* use same value for ends and flank match
* revert back remove short link
* reduce min uniq seq len before trim links
* resolve bug which fails to compile with latest btree
* resolve bug in pre loop value
* increment kmer len by cli step len value
* re-run docs action on change in config

### Maintenance/Refactoring

* Bump version to v2.4.0
* rename function
* set ends non gap len to 11
* revert back run config
* update deps
* rename actions config
* remove and re-add yml to debug
* Update github-action-deploy.yml
* add styles.module.css
* Add back index.module.css
* Add docker info to README
* auto-use git tag version

### New Features

* genotype insertions where read is not fully contained
* change default kmer min, max and step
* add indel specific fisher score threshold
* add strict normal vaf filter for STRs
* double min odds ratio for str variants
* use 1% nml vaf threshold in STR
* add AS/XS filter with 10% difference threshold
* add flank match check to improve precision
* add back short link removal with very lenient params
* add back AS-XS, XT/XA tag read filters
* add kmer step size cli flag

### Reverts

* feat: add back AS-XS, XT/XA tag read filters

### Pull Requests

* Merge pull request [#8](https://github.com/nygenome/Lancet2/issues/8) from nygenome/web_docs


<a name="v2.3.0"></a>
## [v2.3.0](https://github.com/nygenome/Lancet2/compare/v2.2.0...v2.3.0) (2023-08-30)

### Maintenance/Refactoring

* Bump version to v2.3.0

### Performance Improvements

* optimize graph build by removing cord, turn off secure mimalloc


<a name="v2.2.0"></a>
## [v2.2.0](https://github.com/nygenome/Lancet2/compare/v2.1.0...v2.2.0) (2023-08-30)

### Bug Fixes

* reduce default min tumor cov filter

### Maintenance/Refactoring

* Bump version to v2.2.0
* remove redundant flag from run config
* remove runtime stats option
* comment out profiling code from release build
* update deps
* use 10mb non gap window for testing
* update dependencies
* update dependencies and docker file
* update deps and dockerfile
* add explicit preset option in genotyper, defaulting to sr
* use max per sample coverage of 500
* rename image name
* update cmake project version
* change cloud build project
* remove redundant cmake configure
* always use commit group type

### Performance Improvements

* optimize queue operations
* throttle counter and reduce default window & overlap


<a name="v2.1.0"></a>
## [v2.1.0](https://github.com/nygenome/Lancet2/compare/v2.0.0.rc2-main-a055ecdecf...v2.1.0) (2023-07-06)

### Bug Fixes

* use VCFv4.3 and skip adding MQ0 reads

### Maintenance/Refactoring

* Bump version to v2.1.0


<a name="v2.0.0.rc2-main-a055ecdecf"></a>
## [v2.0.0.rc2-main-a055ecdecf](https://github.com/nygenome/Lancet2/compare/v2.0.0.rc1...v2.0.0.rc2-main-a055ecdecf) (2023-06-25)

### Bug Fixes

* resolve undercounting ref for long alleles
* resolve score overflow issue
* resolve typo in while check
* count end gaps and handle it correctly in var ranges
* use more conservative assembly to ref scoring
* create parent path of output vcf if needed
* revert back to v1.14 libdeflate because latest version breaks build
* resolve typo in spdlog usage
* skip regions which have 0 or more than 3 tokens
* include cmake install line
* use debian as default base image
* use relative path to inject right require on build Former-commit-id: ec10bd90fe3df428a08f7076fdd2c3635650f00e Former-commit-id: 40c9b35f7aa984b7b79c6614064cd03ce1aaec5d Former-commit-id: 74cf99ef548050f677de5cce6da1d9ebb9bea6a6
* change config to include static assets Former-commit-id: d818cc677ca0ce1529041ba804ec6e2d4ba67a29 Former-commit-id: c65a78edf57ae926af2bee2ee65f13913b8a32a0 Former-commit-id: e8e292b5af811cd57aa7e2f4536244f10113ab1d
* use public api to refer to transcipt member variable Former-commit-id: f8e1ebc84f3450960bbe7da1a1175204b3ca4dac Former-commit-id: d51fb61953e975c2ed71a5f9c3d031346e6bed41 Former-commit-id: 9cc4265b74796c1f7f3eae08e397890ed611b518
* use node fill colors consistent with v1 Former-commit-id: f698fe76eccb7858781a6941b0fd9e7640dc8913 Former-commit-id: 1691140cbaf41b10c3649af78fa0b458945c0d3b Former-commit-id: 073c09532824d6984acea58824197e7ecc65216d
* alternate index file missing file extension dot Former-commit-id: 9fbc4ef335efbbbbf6cd8599aa1bd492c292d964 Former-commit-id: 248c5500dbd3141600d610f207cf8099f30f112d Former-commit-id: 240652b1a022c224731387819a7ff0d58a95ee48

### Maintenance/Refactoring

* fix typo in dockerfile
* add dockerfile fi
* update dependencies and build tasks
* update local build settings & gitignore
* update CMP0135 policy
* update CMP0135 policy
* update deps
* keep only gcc build Former-commit-id: a4f8b4aeb1353e97ddac1eb08995ba40e5a05944 Former-commit-id: 6d01d6a27eb1c7fdf65dfe33bdbae974bf4d810b Former-commit-id: 4dcb3a240bde0daf58d98d139bd0fb3c9263fe7e
* update graphviz image
* remove ghcr from build to prevent push 403 errors Former-commit-id: 09ffe64accc7479a16c46a4dcdbfe23c0adc3028 Former-commit-id: a8fd94a3e496a1bc24fd319466b4501e054a38fe Former-commit-id: fae57a5d155e73bc566face970b44d53e6f0cefa
* update readme Former-commit-id: 0f920b34997abc4f2bc023d2485dd2db5e05d9ab Former-commit-id: f4f01d634fd6c48fc75cfae8be8dbd5b74ff54b0 Former-commit-id: cc929dc6fc954ca834fa14f32a1ff97cad9b0edf
* Update docs Former-commit-id: 06b10918bcfbf4ec9e7a6a9cf14f68cc61db1ef5 Former-commit-id: 80f708b507a345305a5c006b4d6f6514e1c8f65a Former-commit-id: a1302cd2c95236f53c7a5162c0a9caea9d52ea37
* Update docs Former-commit-id: 155eca1f9665713c1f18c2516862f057cf496260 Former-commit-id: ce5d5a2c9756fe8909d914c61d5ec24f2f0b71a3 Former-commit-id: 93ea17c7d3a5d9ba273c582b2752239a88a05dbb
* Update docs Former-commit-id: 7cd131e9340febfdf68f3275a24c55a344a7fb6c Former-commit-id: e565f6d9889b716cdd122e3435f7ca633541d134 Former-commit-id: 9d0c32e191be8fef0af955acd649b40e5c2348b3
* Update docs Former-commit-id: 4f4f66aba54f151ecb21c8ee124136b5d60283dd Former-commit-id: 5d4a4b1c7645a3b89a758686ba2af59601dfdcd5 Former-commit-id: abd3d761420e5e6d0c004f5ed366c5969e14f901
* Update docs Former-commit-id: 0d716afec1717f698bb1c733b5ae5417f7d52a65 Former-commit-id: 33dbeccd68ef8e375f7fa2ae59c2a20c6f964490 Former-commit-id: 9735acbb78ea328068b708be28cd27cc49caa9da
* Update docs Former-commit-id: e4250b703098de7422cd81c51c4160874fd5eae7 Former-commit-id: 3a1e40f66fd9a3145632940761ae3c43dc866dfb Former-commit-id: 6a4a699221bc73671d7df2a8607cda6e9fd32cdc
* Update docs Former-commit-id: 29971e42f81f179b458ac99a28d06edd719e45a8 Former-commit-id: d47e8ce4a3cbf063556104d481c6438b5e529d29 Former-commit-id: ce29c3f75a97617c61c41edbd33902d1606b7d7c
* Update docs Former-commit-id: 432bfeefe6f13dda0db719a17825ad55a7b89d2c Former-commit-id: f7dbf5412c23e2d7f45740bb8300a5a1acdf5c38 Former-commit-id: bfa255ddaf79995a1d126cd011697cdf12cf003e
* Update docs Former-commit-id: 363da3fe3a87798fd7a25a379fd8ddaa0d9f3e26 Former-commit-id: 20fcca4ef3717bfd7c4fb8852e95fe57e91eeeae Former-commit-id: 5446c5fec6669b4d92ed9213a07ba36dd4daf848
* Update docs Former-commit-id: 3e7ca93dcde23520f35ff7031bfe802944ac9404 Former-commit-id: 25255c194fe0cd248490ddff72db792a54f779e5 Former-commit-id: 1b66b168849637c3ddd7037773a69581379024dc
* Update docs Former-commit-id: 3278677702c49493024d0e682b46ca924f7b78f2 Former-commit-id: 54e99f29258b47b3cef70fa7d9aeef39ddd54038 Former-commit-id: 5ede0eda990754840c590afcac290ba85a26ccab
* Update docs Former-commit-id: 61b6b2e8163eb283929c740955ba69c3175016cd Former-commit-id: 30d5b3e2a4c763040dd3323137f3186286f790c1 Former-commit-id: 26c3254e281e5e1ff71d857420bf3ed0af15db46
* Update docs Former-commit-id: 3b0940c7cf98adb162969448f61f184f8d994484 Former-commit-id: a02f16c46d0020d6b1ed5ac653c8deedba3d8016 Former-commit-id: f0cf954c65c0ece715a28c54e54b283acb2f02c0
* Update docs Former-commit-id: 26435f116131c7d88897d57a5365c1beab8fd28f Former-commit-id: c68d6163558c27a60906615c6a5ff33afcc3b781 Former-commit-id: 3d14e32fcf3e83950d1e5534df9ee3a45d828537
* Update docs Former-commit-id: 2fee87afe551b62c140b2bc4a596d9b47684cbbf Former-commit-id: a1c6cdb20ad154faeaac9cc07ec5bf583dd9ed75 Former-commit-id: 50b0f7c062d0f7306833d460426474e9904ad70b
* update docs Former-commit-id: daaba8c6dcce07bd7685e0c90f2211477c1c1a25 Former-commit-id: 7190b770e35badc7730974caf37c5d764147ff2d Former-commit-id: b23a08bb6d9ace976ac42d4c1d89ce600d0a3eb7
* update docs config Former-commit-id: d13e0716a13205f46cee7c687952f612f537d8ea Former-commit-id: c205583da1f63b476e552eff0bb74747fb54b9c9 Former-commit-id: 9614b69eac2c6b1f9a1694cd651c308601ee7116
* use permalinks Former-commit-id: e77f6c3d2512e3cbe3e03bab5679f058fadbc2d4 Former-commit-id: fc3ffc309b50c02fa62c01ab61825da09eb2841f Former-commit-id: ae0b1c6e0a42011f10f92deaa8750142472f2fcd
* update visualization section in docs Former-commit-id: db225350de4da2a125694e03e01dbc006a9865fc Former-commit-id: d97e23e1f3e4e8ef7d493577cdebf959d1dff559 Former-commit-id: 866a20e845164ff7cbcb0a6564cf1ba134769e81
* use 0 component for raw graph Former-commit-id: e817b66b39179a7307011141582d39c4f68ba14f Former-commit-id: c8c60bb1bdc3fd95b610cdfac6dc08680935a6a8 Former-commit-id: 5c8977e2ba9dc06547219b7469213f9084f53172
* also build dot file before low cov removal Former-commit-id: db2029fd2c4cab5b0563ff8ade81126e43e519e8 Former-commit-id: ab3449a1a319a9f2db05452f14abff3ddf833079 Former-commit-id: 490be6ee4bafe791ad62233b2fd348ad0cf9d953
* ensure non ref alt allele Former-commit-id: 89a8b155e949ba1add275ee95416a0457aa1a6aa Former-commit-id: bdcba126ca1593b7fdcd9673a9ceffe83cac3558 Former-commit-id: 62098deb02f19bae5a4488a181d4b3a7b9df6b26

### New Features

* add rname hash to prevent double counting
* add median
* Fix scoring to call variants longer than read length

### Reverts

* Updating website


<a name="v2.0.0.rc1"></a>
## v2.0.0.rc1 (2022-08-10)

### Bug Fixes

* use separate kmer sizes for ref and alt
* resolve typo which removes all normal reads
* use uniq kmer for path seq
* use bq pass only for low quality singleton normal alt snvs
* check if all low qual bases in normal, skip left/right flank for indels
* resolve logic for when no flank is needede
* resolve logic error in if block
* use bq filter for normal only when one low qual SNV alt base found
* use bq filter for snvs and tumor reads
* use bq filter only for snvs
* use bq filter for all single allele haplotypes
* 0 based ref window excludes end, so remove +1
* use right ref and qry indexes after trimming alignment ([#3](https://github.com/nygenome/Lancet2/issues/3))
* only check high qual snv bases in tumor
* skip adding to kmer coverage for snvs in low qual bases
* add header lines for base qual filters
* remove redundant param
* remove overlap check for read pairs
* always use overlapping reads to increase recall, skip adding empty readInfos
* use right binary path
* go back to full left flank to undo regression in perf
* use small fixed left flank to improve recall
* remove isSomatic ? 0 blocks to fix logic
* use index (0-based) instead of pos (1-based) to pick the spanning node
* break when best path is confirmed
* dont skip nodes in component in remove tips
* swap data with temp to force release of hash table memory
* use vector and .at() to remove undefined behavior
* use padded region end to build windows
* handle unmapped alignments where read and/or mate is unmapped
* use flush ref window instead of index
* handle all window exceptions to prevent deadlock in main thread
* skip on invalid iterator before checking for node ptr
* adjust anchor start on trimming, ensure non bases are not added
* skip early to prevent operating on invalid nodes
* skip returning erased/rehashed node ptrs
* increment index to flush even when no variants are flushed
* allow for merged data source/sink in checks
* trim end gaps in alignments if they exist
* use right positons to point to mismatches (0-based aln start)
* remove -march=native to allow cross compatability
* remove alpine and use ubuntu/debian based builds to fix SIGILLs
* add mutex to log calls to prevent interleaved messages
* move log message before removing
* return from microassembler thread on error
* use shared_ptr without going out of scope
* use right syntax to get value
* use right syntax to get value
* impose strict weak ordering in sort
* add option to skip truncated reference windows

### Maintenance/Refactoring

* use unique path kmer
* use bioinformatics link Former-commit-id: 7033a21a28235d471868e1827d47a01dbb505094 Former-commit-id: 2043ca4c05c59bbb57b13ad1dc82b1b53eaf7acc Former-commit-id: 434ccc77009e54405fb7ec14c04d0053dd91f317
* skip low qual normal snvs
* use raw cov for normal always
* change normal read filters
* minor refactoring and edits
* use old merging logic for node labels
* always use bq pass filter for both tmr and nml read snvs
* use overlapping reads
* use contained reads
* make method static
* add new base coverage class
* use non inlined import
* minor edits and refactoring
* add helper canonical hash helper
* replace variants with higher fet score
* remove use contained flag
* remove alt qual from vcf
* reorder loop for perf
* only filter tmr alt quals
* minor edits to fix clang tidy warnings, log trace on alignment erorr
* remove float specifier
* bump to release candidate 1
* use caught exception directly
* add cmake ide profiles to vcs
* update ide files
* add status badge and change workflow name
* fix case for docker image tag
* add default dockerfile and change scope
* add packages write permission
* fix tag names in actions
* push only to ghcr
* run only on release tags and nightly
* build only amd64 images
* remove extra copies from typo
* remove old badges
* update actions yaml to build and push docker images
* remove typo
* update readme
* update dockerfiles to newer compiler versions
* fix clang-12 compilation error
* fix clang-12 compilation error
* fix clang-12 compilation error
* fix clang-12 compilation error
* fix clang-12 compilation error
* add clang tidy fixes
* add license badge and update docs url Former-commit-id: 40d6386f06c45e158aee48d775b8fd39e7a6e37a Former-commit-id: 3fe3a7abd1aded64e727ebbc7155de6cb7a14edf Former-commit-id: a624022b7bd50cfc367251d7d1cf5db9d933b0ec
* update to BSD license Former-commit-id: 1002421774521f079eb6642a8de94ce93f24dbe2 Former-commit-id: d810ed6a2f7766c88f0be7455b7688b89dfc7b71 Former-commit-id: 92e0c257074d24c61b5739e25443b942f948f8ff
* add gh pages config
* update dependencies and fix related issues for upgrade
* update ide config files
* ignore variable length check
* always return true when fractionToKeep is 1
* use cflags env var directly
* add local cpp env dockerfile
* add local editor files
* add local editor files
* run configure scripts in bash
* add clion cfg files
* add vscode prefs
* fix diagnostics branch case
* revert concurrentqueue version to fix gcc build error
* temporarily ignore gcc warning
* use newly added absl::StatusOr directly
* use htslib's fisher exact impl
* update deps
* cleanup mark ends and print path flow only if present
* use std::fs and allow existing graphs dir
* minor edits
* update dependencies
* add link to post comparing perf of mutex vs spinlocks
* minor edits and revisions
* add version info on startup
* add --verbose option instead of compile time level selection
* add some context for failed assertion
* add log macros
* minor renaming and edits
* use string vector for regions
* use semicolon separated string as regions input
* revert back to uncompressed outvcf for immediate flushing
* write bgzip compressed output vcf
* flush after roughly every 10000bp
* increase required buffer windows
* skip stripping build files in sanitizer images
* remove build files from docker image, run GH action everyday
* edit file watcher outdir
* fix gcc build errors
* remove LTO and remove deps message on build
* rename queue tokens to clarify usage
* reduce IO operations and reuse file handles per thread
* edit log level for messages
* catch align exceptions
* add LTO for gcc builds
* clarify log message
* add terminate handler to prevent deadlock and add more context in logs
* add log message end of microassembler
* use thread id and single letter level fields in logs
* use env vars to get mimalloc stats
* add another failing complex test case
* fix format string
* remove milliseconds from time format
* output uncompressed vcf and check done bits before breaking
* always skip windows with truncated reference sequence
* tag gcc images as latest
* rehash variant store to ensure consistent memory load
* get vcf prefix and always bgzip output vcf
* allow sanitizer build regardless of build type
* remove unused file from build
* write uncompressed vcf by default
* handle non-ref ends in ProcessPath, skip missing nodes
* check revcomp explicitly and ignore missing buddy neighbours
* continue when neighbour is nullptr
* check in inner remaining loop before merging
* remove alignment gap assert
* skip compressing data source and sink nodes
* add more debug asserts
* continue on errors to debug
* add sanitizer blacklist file
* fix typo and reorder deps
* use alpine linux for release builds after SIGILL fix
* add sha and ref tags only for clang release builds
* include integration test data set
* resolve typo in config
* add thread annotations for tsan
* fix tags for docker images
* run CI/CD only manually
* re-add fixed actions config
* remove config to fix actions trigger issue
* point to right repo name
* add ghub actions config
* cleanup lists and remove unused packages
* change sanitizer build compiler flags and add dockerfiles
* remove lto, remove mimalloc in sanitizer builds
* fix clang builds
* add extra check in debug build
* fix typo in condition
* throw exception to debug state
* flush variants more frequently
* move to debug messages
* use reference order in header or fasta
* add informative log and error message
* log stderr to show feedback without buffer
* initial commit with c++17 codebase

### New Features

* also check for non repeat kmer in refseq
* do some basic read filtering on normal
* add left and right flank haplotypes for coverage count
* require atleast 25% overlap to include read in window
* use avg base qual for both snvs and indels
* reduce min tumor vaf callable
* add low avg base qual filter and info fields
* add contig header lines to vcf output
* write bgzf compressed vcf and auto index at end
* use x86 and arm specific pause/yield instructions
* add flag to use only fully contained reads
* skip adding reads not overlapping region
* use weighted harmonic mean to merge node counts
* write per window contigs fasta
* write one path flow dot file with different hues
* add option for cpu profiling
* use spinlock instead of mutex & use relaxed lock when blocking is not needed
* add option to use overlapping reads to build graph
* catch and rethrow thread exceptions to show failures
* use spdlog instead of custom logger
* use concurrent queue for better window result and progress
* add spython converted singularity recipes
* use concurrent queue for pull-based worker model

### Performance Improvements

* batch variants/results to reduce lock contention
* use bool array to reduce allocs
* use inlined vector to reduce allocs for small path builders
* add consumer token when dequeueing windows
* use shared_ptrs to reduce window copies

