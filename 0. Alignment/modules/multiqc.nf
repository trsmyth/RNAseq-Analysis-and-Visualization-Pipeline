process multiqc{

    // container ""
    conda "bioconda::multiqc=1.35"

    input:
    path("inputs/*")
    val output_name

    output:
    path("multiqc_report.html"), emit: report
    path("multiqc_report_data"), emit: data

    script:
    """
    multiqc inputs/ --filename "multiqc_report.html"
    """
}