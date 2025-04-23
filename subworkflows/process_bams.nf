process assign_barcodes{
    label "singlecell"
    cpus 1
    memory "15 GB"
    input:
         tuple val(meta),
               path("high_qual_bc_counts/?_barcode.tsv"),
               path("extract_barcodes.tsv")
    output:
        tuple val(meta),
              path("extract_barcodes_with_bc.tsv"),
              emit: tags
    """
    custom-barcode-correction.py --tag-file extract_barcodes.tsv --umi-counts high_qual_bc_counts/* --method phred > extract_barcodes_with_bc.tsv
    """
}


process concat_barcodes {
    label "singlecell"
    cpus 1
    memory "2 GB"
    publishDir "${params.out_dir}/${meta.alias}/barcode-corrections", mode: 'copy'
    input:
        tuple val(meta),
              path("corrected_barcodes/?_extract_barcodes_with_bc.tsv")
    output:
        tuple val(meta),
              path("extract_barcodes_with_bc.tsv")
    
    """
    cat corrected_barcodes/*_extract_barcodes_with_bc.tsv | awk 'NR==1' > header.txt
    cp header.txt extract_barcodes_with_bc.tsv
    cat corrected_barcodes/*_extract_barcodes_with_bc.tsv | grep -v -f header.txt >> extract_barcodes_with_bc.tsv
    rm header.txt
    """
}


process merge_bams {
    // Combine all BAMs derived from the initial chunking into per sample files
    cpus params.threads
    memory "8 GB"
    container 'library://porchard/default/general:20220107'

    input:
        tuple val(meta),
            path('bams/*aln.bam'),
            path('bams/*aln.bam.bai')
    output:
        tuple val(meta),
              path("merged.sorted.bam"),
              path("merged.sorted.bam.bai"),
              emit: merged_bam

    """
    samtools merge -@ ${task.cpus -1} --write-index -o "merged.sorted.bam##idx##merged.sorted.bam.bai" bams/*.bam
    """
}


process isoquant {

    cpus params.threads
    memory '40 GB'
    publishDir "${params.out_dir}/${meta.alias}/isoquant", mode: 'copy'
    container 'docker://porchard/isoquant:20240911'
    time '24h'

    input:
    tuple val(meta), path(bam), path(bam_index), path(fasta), path(gtfdb)

    output:
    tuple val(meta), path("${library}/${library}.read_assignments.tsv.gz"), emit: read_assignments

    script:
    library = meta.alias

    """
    # HOME is reset to current working directory because isoquant tries to create a directory in and write in $HOME, which results in an error when run in a read-only singularity container
    mkdir -p genedb
    export HOME=. && isoquant.py --output . --bam $bam --data_type ont --reference $fasta --no_model_construction --complete_genedb --genedb $gtfdb --stranded none --no_secondary --labels $library --bam_tags CR,UR --threads ${task.cpus - 1} --prefix $library --genedb_output genedb/
    """

}


process assign_reads_to_transcripts {

    memory '100 GB'
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    container 'library://porchard/default/general:20220107'
    time '5h'
    label 'largemem'

    input:
    tuple val(meta), path(isoquant)

    output:
    tuple val(meta), path("${meta.alias}.transcript-assignments.txt")

    """
    isoquant-read-assignments-to-transcript-assignment.py $isoquant > ${meta.alias}.transcript-assignments.txt
    """
    
}


process assign_reads_to_genes {

    memory '7 GB'
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    container 'library://porchard/default/general:20220107'
    time '10h'

    input:
    tuple val(meta), path(bam), path(bam_index), path(gtf)

    output:
    tuple val(meta), path("${meta.alias}.gene-assignments.txt")

    """
    assign-reads-to-genes.py --strandedness forward --gtf $gtf --bam $bam > ${meta.alias}.gene-assignments.txt
    """

}


process correct_umis {

    memory '150 GB'
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    time '24h'
    label 'largemem'
    container 'docker://ontresearch/wf-single-cell:sha0fcdf10929fbef2d426bb985e16b81153a88c6f4'

    input:
    tuple val(meta), path(corrected_barcodes), path(gene_assignments)

    output:
    tuple val(meta), path("${meta.alias}.corrected-umis.txt")

    """
    correct-umis.py --corrected-barcodes $corrected_barcodes --read-to-gene-assignments $gene_assignments > ${meta.alias}.corrected-umis.txt
    """

}


process tag_bam_and_make_count_matrices_and_calculate_qc_metrics {

    memory '100 GB'
    publishDir "${params.out_dir}/${meta.alias}", mode: 'copy'
    time '24h'
    label 'largemem'
    container 'library://porchard/default/general:20220107'

    input:
    tuple val(meta), path(barcodes), path(gene_assignments), path(transcript_assignments), path(bam), path(bam_index), path(gtf)

    output:
    tuple val(meta), path("${meta.alias}.*")
    tuple val(meta), path("${meta.alias}.genes.matrix.mtx"), path("${meta.alias}.genes.features.tsv"), path("${meta.alias}.genes.barcodes.tsv"), emit: gene_matrices
    tuple val(meta), path("${meta.alias}.transcripts.matrix.mtx"), path("${meta.alias}.transcripts.features.tsv"), path("${meta.alias}.transcripts.barcodes.tsv"), emit: transcript_matrices
    tuple val(meta), path("${meta.alias}.tagged.bam"), emit: bam

    """
    tag-bam-and-make-count-matrices-and-calculate-qc-metrics.py --gene-assignments $gene_assignments --transcript-assignments $transcript_assignments --barcodes $barcodes --gtf $gtf --bam-in $bam --bam-out ${meta.alias}.tagged.bam --count-matrix-prefix ${meta.alias}. --qc-metrics-out ${meta.alias}.qc.txt
    """

}

// changed to use trim-for-scafe-fixed.py
// old command:
// trim-for-scafe.py --input-bam $bam --output-bam ${meta.alias}.trimmed.unsorted.bam --trim-to 100 --ggg-mismatches-allowed 1
process trim_for_scafe {

    publishDir "${params.out_dir}/${meta.alias}/trim-for-scafe", mode: 'copy'
    container 'docker://porchard/general:20220406125608'
    time '10h'
    memory '5 GB'

    input:
    tuple val(meta), path(bam)

    output:
    tuple val(meta), path("${meta.alias}.trimmed.bam"), path("${meta.alias}.trimmed.bam.bai")

    """
    trim-for-scafe-fixed.py --max-softclipping 4 --input-bam $bam --output-bam ${meta.alias}.trimmed.unsorted.bam
    samtools sort -m 3G -o ${meta.alias}.trimmed.bam ${meta.alias}.trimmed.unsorted.bam
    samtools index ${meta.alias}.trimmed.bam
    """

}


process cellbender {

    cpus 1
    memory '40 GB'
    publishDir "${params.out_dir}/${meta.alias}/cellbender"
    container 'docker://porchard/cellbender:0.3.0'
    time '24h'
    errorStrategy 'ignore'
    label 'gpu'

    input:
    tuple val(meta), path('matrix.mtx'), path('genes.tsv'), path('barcodes.tsv'), val(feature_type)

    output:
    path("${meta.alias}*")
    path("${meta.alias}*.h5"), emit: h5_files

    """
    cellbender remove-background --cuda --epochs 150 --fpr 0.01 0.05 0.1 --input . --output ./${meta.alias}.${feature_type}.cellbender.h5
    cp .command.log ${meta.alias}.${feature_type}.log
    """

}


workflow process_bams {
    take:
        bam
        extracted_barcodes
        high_qual_bc_counts
        gtf
        gtfdb
        ref_genome_fasta
        ref_genome_idx
    main:


        // run IsoQuant
        // get gene counts (custom python script)
        // correct barcodes (custom python script)
        // correct UMIs
        // now have gene, transcript, CB, UMI for all reads
        // create the count matrices
        // tag the bam file

        assign_barcodes(
            high_qual_bc_counts
            .cross(extracted_barcodes)
            .map {it ->
                meta = it[0][0]
                counts = it[0][1]
                barcodes = it[1][1]
                [meta, counts, barcodes]})

        concat_barcodes(assign_barcodes.out.tags.groupTuple())

        // Combine the BAM chunks per-sample
        merge_bams(bam.groupTuple())

        // get read --> transcript assignments
        isoquant(merge_bams.out.merged_bam.combine(Channel.fromPath(ref_genome_fasta)).combine(Channel.fromPath(gtfdb)))
        assign_reads_to_transcripts(isoquant.out.read_assignments)
        assign_reads_to_genes(merge_bams.out.merged_bam.combine(Channel.fromPath(gtf)))
        
        correct_umis(concat_barcodes.out.combine(assign_reads_to_genes.out, by: 0))
        tagged_bam = tag_bam_and_make_count_matrices_and_calculate_qc_metrics(correct_umis.out.combine(assign_reads_to_genes.out, by: 0).combine(assign_reads_to_transcripts.out, by: 0).combine(merge_bams.out.merged_bam, by: 0).combine(Channel.fromPath(gtf)))

        trim_for_scafe(tagged_bam.bam)
        tagged_bam.gene_matrices.map({it -> it + ['genes']}).mix(tagged_bam.transcript_matrices.map({it -> it + ['transcripts']})) | cellbender


    emit:
        merged_bam = merge_bams.out.merged_bam
}
