include {make_index} from './modules/make_index.nf'
include {map_t2g} from './modules/map_t2g.nf'
include {fastqc} from './modules/fastqc.nf'
include {trim} from './modules/trim.nf'
include {align} from './modules/align.nf'
include {count_estimate} from './modules/count_estimate.nf'
include {merge_counts} from './modules/merge_counts.nf'
include {multiqc} from './modules/multiqc.nf'

workflow{

    main:

    read_ch = channel.fromFilePairs("${params.input}/*R{1,2}_001.fastq.gz") // Current data format is *R{1,2}_001.fastq.gz. Replace as needed.

    ch_ref = channel.fromPath(params.reference, checkIfExists: true) // Create single channel for reference file
    target_path_ch = channel.value(params.index) // Create single value channel for index path
    make_index(ch_ref, target_path_ch) // Check if kallisto reference index already exists and create it if not. Copy the file and pass the path as output
    ch_index = make_index.out.processed_file.collect() // Collect the output for downstream alignment

    ///////////////////////

    ch_t_g_map = channel.fromPath(params.t2g_mapping, checkIfExists: true) // Create single channel for transcript to gene mapping file
    target_path_ch = channel.value(params.mapped_t2g) // Create single value channel for path of mapped file
    map_t2g(ch_t_g_map, target_path_ch) // Check if transcript to gene map already exists and create it if not. Copy the file and pass the path as output
    ch_map = map_t2g.out.processed_file.collect() // Collect the output for downstream alignment

    ///////////////////////
    
    fastqc(read_ch)
    trim(read_ch)

    align(trim.out.trimmed_reads, ch_index)

    count_estimate(align.out.tsv, ch_map)
    all_counts = count_estimate.out.count_estimate.collect()

    merge_counts(all_counts)

    multiqc_files_ch = channel.empty() // Create an empty channel and mix previous outputs
        .mix(
            fastqc.out.zip, 
            fastqc.out.html, 
            trim.out.trimming_reports,
            trim.out.html,
            trim.out.zip,
            align.out.json
        )

    multiqc_files_list = multiqc_files_ch.groupTuple() // Group files as tup(sample_id, [[files_output_x], [files_output_y]])
        .map{nested_files ->  return nested_files[1].flatten()} // Map tup(sample_id, [[files_output_x], [files_output_y]]) to [output_files]
        .collect() // Collect (Sample_X to Sample_Y) [output_files] to [output_files]

    multiqc(multiqc_files_list, params.report_id)

    ///////////////////////

    publish:

    align_tsv = align.out.tsv
    align_log = align.out.json

    estimated_counts = count_estimate.out.count_estimate
    counts = merge_counts.out.merged_counts

    multiqc_report = multiqc.out.report
    multiqc_data = multiqc.out.data

}

output{

    align_tsv{
        path 'align'
        mode 'copy'
    }
    align_log{
        path 'align'
        mode 'copy'
    }

    ////////////////////////

    estimated_counts{
        path 'counts'
        mode 'copy'
    }

    ////////////////////////

    counts{
        path 'merged_counts'
        mode 'copy'
    }

    ////////////////////////

    multiqc_report{
        path 'multiqc'
        mode 'copy'
    }
    multiqc_data{
        path 'multiqc'
        mode 'copy'
    }
}
