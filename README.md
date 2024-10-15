# NextFlow pipeline for 10X snRNA-seq ONT data

## Dependencies
Singularity (v. 3) and NextFlow (>= v. 20.10.0).

## Configuration

### Input files

Input files should be organized as follows:

Fasta/GTF reference files should be placed in a `ref_genome_dir`, with structure:

```
--ref_genome_dir
  |
  |--fasta
  |  |
  |  |--genome.fa
  |  |--genome.fa.fai
  |
  |--genes
     |
     |-- genes.gtf
```

Fastq files should be organized within a `fastq` directory. This directory should have one subdirectory per library (subdirectory name is the library name), and the subdirectory for a library should contain all the fastq files (could be a single fastq file or many fastq files), with names ending in `*.fastq.gz`, `*.fq.gz`, `*.fastq`, or `*.fq`. For example:

```
─── fastq_directory
    ├── sample_1
    │   ├── reads0.fastq
    │   └── reads1.fastq
    ├── sample_2
    │   ├── reads0.fastq
    │   ├── reads1.fastq
    │   └── reads2.fastq
    └── sample_3
        └── reads0.fastq
```

You'll also need a gffutils DB file. If you don't have one of these, you can generate it (will takes several hours) with:

```
#!/usr/bin/env python

import gffutils
db = gffutils.create_db('/path/to/gtf', 'genome.db', disable_infer_genes=True, disable_infer_transcripts=True)
```

When launching the pipeline, as shown in the `nextflow` command below, you'll need to provide the following parameters:

1. 10X kit used (3prime:v2, 3prime:v3, 3prime:v4, or multiome:v1; e.g. `--kit 'multiome:v1`)
2. Threads available (if running on greatlakes, `--threads 24` should be fine)
3. The path to the fastq directory (`--fastq /path/to/fastq_dir/fastq`)
4. The path to the ref_genome_dir directory (`--ref_genome_dir /path/to/ref-genome-dir`)
5. The path to the gffutils DB file (`--gtfdb /path/to/genome.db`)


## Running
Once you have all of the above information, you can run the pipeline as follows (the results will be generated in the directory in which the pipeline is run):

```bash
mkdir -p my_results && cd my_results && nextflow run -resume --gtfdb /scratch/scjp_root/scjp0/porchard/2023-HSM-ONT/work/gffutils-db/results/db/genes-no-ercc.db --fastq /path/to/fastq_dir/fastq --kit 'multiome:v1' --threads 24 --ref_genome_dir /path/to/ref-genome-dir -profile singularity /path/to/main.nf
```


## Pipeline steps

1. Adapter identification, fused read splitting and stranding (this is taken from the [original ONT pipeline](https://github.com/epi2me-labs/wf-single-cell)).
2. Mapping of reads to genomic reference using minimap2 (also from the original ONT pipeline)
3. Barcode correction (custom python script)
4. Quantification of transcript expression (using a reference GTF and IsoQuant)
5. Quantification of gene expression (using a reference GTF and custom python script)
6. UMI correction  (also from the original ONT pipeline)
7. Calculation of QC statistics (custom python script)
8. Trimming of reads for compatibility w/ SCAFE (custom python script)
9. Removal of ambient RNA (cellbender)


## Output
The important output files are:

* `output/{sample_name}/{sample_name}.genes.{matrix.mtx,features.tsv,barcodes.tsv}`: Gene-level count matrices
* `output/{sample_name}/{sample_name}.transcripts.{matrix.mtx,features.tsv,barcodes.tsv}`: Transcript-level count matrices
* `output/{sample_name}/{sample_name}.qc.txt`: QC metrics
* `output/{sample_name}/{sample_name}.tagged.bam`: BAM file
* `output/{sample_name}/cellbender`: Cellbender results (gene level and transcript level)
* `output/{sample_name}/trim-for-scafe`: BAM file prepped for SCAFE 5' end processing


Other intermediate files will be present but can generally be ignored:
* `output/{sample_name}/{sample_name}.corrected-umis.txt`
* `output/{sample_name}/{sample_name}.gene-assignments.txt`
* `output/{sample_name}/{sample_name}.transcript-assignments.txt`
* `output/{sample_name}/isoquant`