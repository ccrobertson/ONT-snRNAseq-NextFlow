#!/usr/bin/env nextflow

import groovy.json.JsonBuilder
import nextflow.util.BlankSeparatedList;
nextflow.enable.dsl = 2

include { fastq_ingress; xam_ingress } from './lib/ingress'
include { preprocess } from './subworkflows/preprocess'
include { process_bams } from './subworkflows/process_bams'


process getVersions {
    label "singlecell"
    cpus 1
    memory "1 GB"
    output:
        path "versions.txt"
    script:
    """
    python -c "import pysam; print(f'pysam,{pysam.__version__}')" >> versions.txt
    python -c "import parasail; print(f'parasail,{parasail.__version__}')" >> versions.txt
    python -c "import pandas; print(f'pandas,{pandas.__version__}')" >> versions.txt
    python -c "import rapidfuzz; print(f'rapidfuzz,{rapidfuzz.__version__}')" >> versions.txt
    python -c "import sklearn; print(f'scikit-learn,{sklearn.__version__}')" >> versions.txt
    fastcat --version | sed 's/^/fastcat,/' >> versions.txt
    minimap2 --version | sed 's/^/minimap2,/' >> versions.txt
    samtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    bedtools --version | head -n 1 | sed 's/ /,/' >> versions.txt
    gffread --version | sed 's/^/gffread,/' >> versions.txt
    seqkit version | head -n 1 | sed 's/ /,/' >> versions.txt
    stringtie --version | sed 's/^/stringtie,/' >> versions.txt
    gffcompare --version | head -n 1 | sed 's/ /,/' >> versions.txt
    """
}


process getParams {
    label "singlecell"
    cpus 1
    memory "1 GB"
    output:
        path "params.json"
    script:
        def paramsJSON = new JsonBuilder(params).toPrettyString()
    """
    # Output nextflow params object to JSON
    echo '$paramsJSON' > params.json
    """
}


process parse_kit_metadata {
    label "singlecell"
    cpus 1
    memory "1 GB"
    input:
        path 'sample_ids'
        path sc_sample_sheet
        path kit_config, stageAs: 'kit_config.csv'
    output:
        path "merged.csv"
    script:
    if (sc_sample_sheet.name != "OPTIONAL_FILE"){
        """
        workflow-glue parse_kit_metadata from_sheet \
            --user_config ${sc_sample_sheet} \
            --kit_config kit_config.csv \
            --sample_ids sample_ids \
            --output merged.csv
        """
    }else{
        """
        workflow-glue parse_kit_metadata from_cli \
            --kit_config kit_config.csv \
            --kit "$params.kit" \
            --expected_cells $params.expected_cells \
            --sample_ids $sample_ids \
            --output merged.csv
        """
    }
}



// workflow module
workflow pipeline {
    take:
        chunks
        ref_genome_dir
    main:
        ref_genome_fasta = file("${params.ref_genome_dir}/fasta/genome.fa", checkIfExists: true)
        ref_genome_idx = file("${params.ref_genome_dir}/fasta/genome.fa.fai", checkIfExists: true)
        ref_genes_gtf = file("${params.ref_genome_dir}/genes/genes.gtf", checkIfExists: true)
        ref_gtfdb = file(params.gtfdb, checkIfExists: true)
        software_versions = getVersions()
        workflow_params = getParams()

        bc_longlist_dir = file("${projectDir}/data", checkIfExists: true)

        preprocess(
            chunks.map{meta, fastq, stats -> [meta, fastq]},
            bc_longlist_dir,
            ref_genome_fasta,
            ref_genome_idx,
            ref_genes_gtf)
        // preprocess emits:
        // bam_sort = call_adapter_scan.out.bam_sort
        // bam_stats = call_adapter_scan.out.bam_stats
        // read_tags = call_adapter_scan.out.read_tags
        // high_qual_bc_counts = call_adapter_scan.out.barcode_counts
        // adapter_summary = call_adapter_scan.out.adapter_summary

        process_bams(
            preprocess.out.bam_sort,
            preprocess.out.read_tags,
            preprocess.out.high_qual_bc_counts.groupTuple(),
            ref_genes_gtf,
            ref_gtfdb,
            ref_genome_fasta,
            ref_genome_idx)

}


// entrypoint workflow
WorkflowMain.initialise(workflow, params, log)
workflow {

    ref_genome_dir = file(params.ref_genome_dir, checkIfExists: true)


    if (params.kit_config){
        kit_configs_file = file(params.kit_config, checkIfExists: true)
    }else{
        kit_configs_file = file("${projectDir}/kit_configs.csv", checkIfExists: true)
    }

    if (params.fastq) {
        samples = fastq_ingress([
                "input":params.fastq,
                "sample":params.sample,
                "sample_sheet":params.sample_sheet,
                "fastq_chunk": params.fastq_chunk,
                "stats": true,
                "per_read_stats": false])
    } else {
        samples = xam_ingress([
                "input":params.bam,
                "sample":params.sample,
                "sample_sheet":params.sample_sheet,
                "fastq_chunk": params.fastq_chunk,
                "keep_unaligned": true,
                "return_fastq": true,
                "stats": true,
                "per_read_stats": false])

    }
    if (!params.single_cell_sample_sheet) {
        sc_sample_sheet = file("$projectDir/data/OPTIONAL_FILE")
    } else {
        // Read single_cell_sample_sheet
        sc_sample_sheet = file(params.single_cell_sample_sheet, checkIfExists: true)
    }

    fastqingress_ids = samples.map {meta, file, stats -> meta.alias }.unique().collectFile(newLine: true)
    // Get [sample_id, kit_meta]
    kit_meta = parse_kit_metadata(fastqingress_ids, sc_sample_sheet, kit_configs_file)
        .splitCsv(header:true)
        .map {it -> [it['sample_id'], it]}
    // Merge the kit metadata onto the sample metadata
    sample_and_kit_meta = kit_meta
        .cross(samples
            // Put sample_id as first element for join
            .map {meta, chunk, stats -> [meta.alias, meta, chunk, stats]})
        // Extract the joined sample and kit info from the cross results
        .map {kit, sample -> [ sample[1] + kit[1], sample[2], sample[3]]}
        // we never need the chunk index for merging items so discard it
        .map {meta, chunk, stats ->
            def new_meta = meta.clone()
            new_meta.remove('group_index')
            [new_meta, chunk, stats]}

    pipeline(
        sample_and_kit_meta,
        ref_genome_dir)
}

