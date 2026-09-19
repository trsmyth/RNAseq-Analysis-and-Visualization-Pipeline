process count_estimate{
    
    // container ""
    conda "bioconda::pytximport=0.13.0"

    input:
    tuple val(sample_id), path(reads)
    path t2g

    output:
    path("./counts/${sample_id}_counts.csv"), emit: count_estimate

    script:
    """
    pytximport -i ${reads} -t kallisto -m ${t2g} -o ./counts/${sample_id}_counts.csv
    sed -i '1s/.*/Gene,${sample_id}/' ./counts/${sample_id}_counts.csv
    """
}