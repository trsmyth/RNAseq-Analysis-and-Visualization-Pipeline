process make_index {

    conda "bioconda::kallisto=0.52.0"

    input:
    path reference
    val index_str

    output:
    path "./processed_index.idx", emit: processed_file

    script:
    def target_file = file(index_str)
    """
    if [ ! -f ${target_file} ]; then
        echo "Index not found. Generating index."
        kallisto index -i ${target_file} ${reference}
        cp "${target_file}" ./processed_index.idx
    else
        echo "Index found. Skipping index generation."
        cp "${target_file}" ./processed_index.idx
    fi
    """
}