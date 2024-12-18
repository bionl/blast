#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    def input_fasta = file(params.input)
    def blast_db = file(params.db)

    RUN_BLAST(input_fasta, blast_db)
        | view {
            println("BLAST results written to: ${it.name}")
        }
}

process RUN_BLAST {
    container "ncbi/blast:latest"

    input:
    path query_file
    path db_files

    output:
    path "*.out"

    script:
    """
    blastn -query ${query_file} -db ${db_files} -out results.out -outfmt 6
    """
}
