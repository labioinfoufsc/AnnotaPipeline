# Configuration

Copy `annotapipeline.yaml` and edit it for your run. Software paths are
pre-filled for the conda environment; replace them with full paths if you
use different installations.

Configure the reference database paths before execution. You can remove the
Kallisto and Comet sections entirely when transcriptomics or proteomics are
not used.

Accepted values are documented inline in `annotapipeline.yaml`.
`config_example.yaml` shows a complete legacy configuration.

The Nextflow workflow uses a JSON `-params-file` instead; see
[`docs/nextflow.md`](../docs/nextflow.md).
