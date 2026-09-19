process map_t2g {

    input:
    path map_file_path
    val mapped_file_path

    output:
    path "./t2g.tsv", emit: processed_file

    script:
    def target_file = file(mapped_file_path) 
    """
    if [ ! -f "${target_file}" ]; then

        echo "Transcript to gene map not found. Generating map."

        # Write the Python code to a temporary file using a Heredoc
        cat << 'EOF' > generate_map.py

    #!/usr/bin python3
    import gzip
    import csv
    import re

    t2g = {"transcript_id": [], "gene_id": []} # Create storage dict

    with gzip.open("${map_file_path}", 'rt') as f_in:
        for line in f_in:
            if line.startswith('#'):
                continue # Skip lines that start with hash

            columns = line.strip().split('\t') # Clear new line and split at tab

            # Homo_sapiens.GRCh38.116.gtf.gz has third column of 'feature' (transcript, gene, exon, etc)
            # and 'attributes' at ninth column which contains transcript and gene information
            if len(columns) > 8 and columns[2] == 'transcript': 
                t_id_match = re.search(r'transcript_id "([^"]+)"', columns[8])
                g_id_match = re.search(r'gene_id "([^"]+)"', columns[8])
                
                # If both are found, add them to the storage dict
                if t_id_match and g_id_match:
                    t2g["transcript_id"].append(t_id_match.group(1))
                    t2g["gene_id"].append(g_id_match.group(1))

    pairs = list(zip(t2g["transcript_id"], t2g["gene_id"])) # Pair the matching transcript_id and gene_id
    unique_pairs = list(dict.fromkeys(pairs)) # Find unique pairs
    t2g = {"transcript_id": [p[0] for p in unique_pairs], "gene_id": [p[1] for p in unique_pairs]} # Remake the storage dict

    # Export the results as a .tsv file
    with open("${target_file}", "w", newline = "", encoding = "utf-8") as f:
        writer = csv.writer(f, delimiter = "\t")
        writer.writerow(["transcript_id", "gene_id"])
        for t, g in zip(t2g["transcript_id"], t2g["gene_id"]):
            writer.writerow([t, g])

    EOF

        python3 generate_map.py

        cp "${target_file}" ./t2g.tsv

    else
        echo "Transcript to gene map found. Skipping map generation."
        cp "${target_file}" ./t2g.tsv
    fi
    """
}