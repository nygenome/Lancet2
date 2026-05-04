# Native Cloud Streaming

Lancet2 supports streaming reads, reference genomes, BED windows, and writing output VCFs directly from/to cloud buckets on AWS (Amazon S3) and GCP (Google Cloud Storage), as well as standard HTTP(S) and FTP(S) endpoints. This eliminates the need for local file copies or data staging.

> [!TIP]
> Pre-built packages from [prefix.dev](https://prefix.dev/channels/lancet2) include cloud I/O support out of the box. See [Installation](installation.md#pre-built-packages-recommended) for install commands — no source build required.

## Enabling Cloud I/O Support

When building from source, cloud streaming is an opt-in feature because it requires dynamic linking against the system `CURL` and `OpenSSL` libraries.

Configure and build Lancet2 with the `-DLANCET_ENABLE_CLOUD_IO=ON` flag and dynamic linkage (`LANCET_BUILD_STATIC=OFF`) enabled:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DLANCET_BUILD_STATIC=OFF -DLANCET_ENABLE_CLOUD_IO=ON ..
make -j$(nproc)
```

## Authentication

Authentication is handled by `libcurl` and `htslib`'s AWS/GCP remote protocol APIs. You only need to export valid environmental credentials.

### AWS (S3 Pipelines)

To analyze BAM/CRAM files located directly at `s3://bucket/data.bam`:

```bash
export AWS_ACCESS_KEY_ID="your_access_key"
export AWS_SECRET_ACCESS_KEY="your_secret_key"
export AWS_DEFAULT_REGION="us-east-1"

Lancet2 pipeline \
    --normal s3://bucket/normal.cram \
    --tumor s3://bucket/tumor.cram \
    --reference s3://bucket/GRCh38.fa \
    --bed-file s3://bucket/regions.bed \
    --out-vcfgz s3://bucket/results.vcf.gz
```

### Google Cloud (GCS Pipelines)

For `gs://` pipelines, authenticate via standard Google Application Credentials. Due to a limitation in `libcurl`/`htslib` where static tokens expire after 1 hour, **you should not use `GCS_OAUTH_TOKEN` for workloads > 1 hour.**

Instead, the standard workaround (validated in [samtools/htslib#803](https://github.com/samtools/htslib/issues/803)) uses a background script and a named UNIX pipe (`mkfifo`) routed through `HTS_AUTH_LOCATION` to automatically refresh expired tokens.

To enable persistent GCS streaming, run this script alongside the pipeline:

```bash
# 1. Create a transient Named Pipe
mkfifo /tmp/gcp_token_fifo

# 2. Launch background token regeneration daemon
( while true ; do gcloud auth application-default print-access-token > /tmp/gcp_token_fifo ; done ) &
export TOKEN_DAEMON_PID=$!

# 3. Ensure cleanup on exit (signal or crash)
trap 'kill $TOKEN_DAEMON_PID 2>/dev/null; rm -f /tmp/gcp_token_fifo' EXIT INT TERM

# 4. Point htslib to the pipe for automatic token refresh
export HTS_AUTH_LOCATION="/tmp/gcp_token_fifo"

Lancet2 pipeline \
    --normal gs://bucket/normal.cram \
    --tumor gs://bucket/tumor.cram \
    --reference gs://bucket/GRCh38.fa \
    --bed-file gs://bucket/regions.bed \
    --out-vcfgz gs://bucket/results.vcf.gz
```

> [!TIP]
> **Cloud Compute Engines**
> If you are running Lancet2 inside a Google Compute Engine (GCE), Vertex AI, or GKE pod with an attached Service Account, you can replace the `gcloud` command with a faster metadata endpoint request:
> `curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/gcp_token_fifo`

### HTTP & FTP Authentication

For standard `http://`, `https://`, `ftp://`, or `ftps://` web endpoints, `libcurl` handles Basic Authentication via URL-embedded credentials (`protocol://user:password@host/path/to/file`):

```bash
Lancet2 pipeline \
    --tumor https://username:password@data.example.com/datasets/tumor.bam \
    --reference ftps://user:pass@ftp.example.org/GRCh38.fa \
    ...
```

## Cloud Authentication Pre-Validation

When the `--out-vcfgz` path points to a cloud bucket, Lancet2 performs a **zero-byte HTTP PUT** immediately at startup before processing any windows. This upfront authentication check validates that your credentials are valid and the output bucket is writable.

!!! tip "Why this matters"
    Cloud storage APIs (AWS S3, GCS) use **5 MB+ multipart chunk caching** over `libcurl`. Because the VCF header is smaller than this threshold, `libcurl` defers all HTTP network handshakes until the final flush at the end of pipeline execution. Without the pre-validation check, a misconfigured credential or bucket permission error would only surface **after the entire pipeline completes** — potentially 40+ hours into a WGS run.

    The zero-byte PUT bypasses the multipart cache and forces an immediate API round-trip, catching authentication failures within the first second of execution.

## Index Streaming UX

When the argument `--out-vcfgz` points to a cloud bucket, Lancet2 supports streaming uploads for the compressed VCF payload.

Immediately following the upload, Lancet2 compiles a `tabix` index (`.tbi`) by re-reading the uploaded VCF from the cloud to compute the index blocks, then streams the `.tbi` file back to the same bucket.

