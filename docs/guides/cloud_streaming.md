# Native Cloud Streaming (S3 & GCS)

Lancet2 supports streaming reads, reference genomes, BED windows, and writing output VCFs natively from/to private cloud buckets on AWS (Amazon S3) and GCP (Google Cloud Storage), as well as supporting regular HTTP(S) URLs. This enables highly efficient workflows natively within cloud environments without requiring massive local storage provisioning or cumbersome data localization orchestration.

## Enabling Cloud I/O Support

Since Lancet2 prioritizes maximal portability for standalone source builds by employing rigid static linkage, native cloud streaming is an opt-in feature. Cloud networking intrinsically requires dynamically linking the underlying system `CURL` and `OpenSSL` networking stacks.

Configure and build Lancet2 with the `-DLANCET_ENABLE_CLOUD_IO=ON` flag and dynamic linkage (`LANCET_BUILD_STATIC=OFF`) enabled:

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DLANCET_BUILD_STATIC=OFF -DLANCET_ENABLE_CLOUD_IO=ON ..
make -j$(nproc)
```

> **Note:** If you pull our pre-built official Docker containers (Conda packaging is actively in-development), these dependencies are inherently resolved dynamically and pre-configured out-of-the-box!

## Authentication

Authentication protocols are transparently managed by `libcurl` and `htslib`'s AWS/GCP remote protocol APIs deep within Lancet2. You only need to export valid environmental credentials.

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

For `gs://` pipelines, authenticate via standard Google Application Credentials. Due to a native limitation in `libcurl`/`htslib` on GCP where static tokens naturally rigidly expire after 1 hour, **you should not use `GCS_OAUTH_TOKEN` for workloads > 1 hour!** 

Instead, the universally accepted native streaming workaround (validated in [samtools/htslib#803](https://github.com/samtools/htslib/issues/803)) uses a lightweight background script and a named UNIX pipe (`mkfifo`) specifically routed through `HTS_AUTH_LOCATION` to seamlessly negotiate dynamic token regeneration asynchronously behind the scenes!

To enable infinite resilient GCS streaming, run this quick script alongside the pipeline:

```bash
# 1. Create a transient Named Pipe
mkfifo /tmp/gcp_token_fifo

# 2. Launch background token regeneration daemon
( while true ; do gcloud auth application-default print-access-token > /tmp/gcp_token_fifo ; done ) &
export TOKEN_DAEMON_PID=$!

# 3. Inform Lancet2 to poll the pipe whenever connections natively drop
export HTS_AUTH_LOCATION="/tmp/gcp_token_fifo"

Lancet2 pipeline \
    --normal gs://bucket/normal.cram \
    --tumor gs://bucket/tumor.cram \
    --reference gs://bucket/GRCh38.fa \
    --bed-file gs://bucket/regions.bed \
    --out-vcfgz gs://bucket/results.vcf.gz

# 4. Clean up the daemon securely afterwards
kill $TOKEN_DAEMON_PID && rm /tmp/gcp_token_fifo
```

> [!TIP]
> **Cloud Compute Engines**  
> If you are running Lancet2 inside a Google Compute Engine (GCE), Vertex AI, or GKE pod with an attached Service Account, you can optimize the daemon massively by replacing the `gcloud` subshell command entirely with a wildly fast localized metadata payload request:
> `curl -s -H "Metadata-Flavor: Google" http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token > /tmp/gcp_token_fifo`

### HTTP & FTP Authentication

For standard `http://`, `https://`, `ftp://`, or `ftps://` web endpoints, Lancet2 relies on `libcurl`'s native URL parsing to handle Basic Authentication seamlessly. 

If your web server requires basic HTTP user authentication, simply embed your credentials directly into the URL using standard web schema formatting (`protocol://user:password@host/path/to/file`):

```bash
Lancet2 pipeline \
    --tumor https://username:password@data.example.com/datasets/tumor.bam \
    --reference ftps://user:pass@ftp.example.org/GRCh38.fa \
    ...
```

## Index Streaming UX

When the argument `--out-vcfgz` points to a cloud bucket, Lancet2 natively supports multipart/streaming uploads for the compressed VCF payload in real-time.

Immediately following the upload execution, Lancet2 guarantees data-integrity by compiling a `tabix` index (`.tbi`). This triggers a brief, automated secondary read sequence of the fresh VCF from the cloud to compute the background index blocks, synchronously streaming the `.tbi` file straight back to your private bucket as well.
