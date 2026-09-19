process fastqc {

    // container ""
    conda "bioconda::fastqc=0.12.1"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip
    tuple val(sample_id), path("*_fastqc.html"), emit: html

    // fastqc will process a list of files sequentially
    script:
    """
    fastqc ${reads}
    """
}