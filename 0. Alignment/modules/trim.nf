process trim{
    
    // container ""
    conda "bioconda::trim-galore=2.3.0"

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), path("*{_trimmed_,_val_}{1,2}.fq.gz") , emit: trimmed_reads
    tuple val(sample_id), path("*_trimming_report.txt"), emit: trimming_reports
    tuple val(sample_id), path("*_fastqc.html"), emit: html
    tuple val(sample_id), path("*_fastqc.zip"), emit: zip

    script:
    """
    trim_galore --fastqc --paired ${reads}
    """
}