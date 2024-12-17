#!/usr/bin/env nextflow

nextflow.enable.dsl=2

def input_fasta = file(params.input)
def blast_db = file(params.db)
def output_dir = file(params.outdir)

workflow {
    RUN_BLAST()

    blast_results.view {
        println "BLAST results written to: ${it.name}"
    }
}

process RUN_BLAST {
    container "ncbi/blast:latest"

    input:
    path query_file from input_fasta
    path db_files from blast_db

    output:
    path "*.out" into blast_results

    script:
    """
    blastn -query $query_file \
           -db $db_files \
           -out results.out \
           -outfmt 6
    """
}
