process merge_counts{

    conda "conda-forge::csvkit=2.2.0"

    input:
    path reads

    output:
    path "./merged_output.csv", emit: merged_counts

    script:
    """
    csvjoin -c "Gene" *.csv > ./merged_output.csv
    """
}