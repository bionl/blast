#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    def input_fasta = file(params.input)
    def blast_db = file(params.db_files)

    RUN_BLAST(input_fasta, blast_db, params.db_name)
        | view {
            println("BLAST results written to: ${it.name}")
        }
}

process RUN_BLAST {
    container "ncbi/blast:latest"
    publishDir params.outdir, mode: 'copy'

    input:
    path query_file
    path db_files
    val db_name

    output:
    path "*.txt"

    script:
    """
    ulimit -n 65536
    blastn -query ${query_file} -db ${db_files}/${db_name} -out results.txt -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore" -num_threads 4
    """
}
