process align {
    // container ""
    conda "bioconda::kallisto=0.52.0"

    input:
    tuple val(sample_id), path(reads)
    path transcriptome

    output:
    tuple val(sample_id), path("${sample_id}_kallisto/abundance.tsv"), emit: tsv
    tuple val(sample_id), path("${sample_id}_kallisto/abundance.h5") , emit: h5
    tuple val(sample_id), path("${sample_id}_kallisto/${sample_id}_run_info.json"), emit: json

    // kallisto quant : Calls Kallisto
    // -i ${transcriptome} : Loads the transcriptome reference file defined in workflow/params
    // -o ${prefix}_kallisto : Specifies output
    // ${args} : Injects specified args, or ""
    // ${read_string} : Insert the path(s)

    script:
    def args = task.ext.args ?: "" // If args is not defined in configs, set to ""
    def read_string = "${reads[0]} ${reads[1]}" // Since workflow uses paired, provide string 1 and 2

    """
    kallisto quant \
    -i ${transcriptome} \
    -o ${sample_id}_kallisto \
    ${args} \
    ${read_string}

    mv ${sample_id}_kallisto/run_info.json ${sample_id}_kallisto/${sample_id}_run_info.json
    """
}