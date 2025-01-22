#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

workflow {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true)
        .set { all_samples }

    all_samples
        .filter { it.fasta.endsWith('.fastq') }
        .set { fastq_samples }

    all_samples
        .filter { it.fasta.endsWith('.fasta') }
        .set { fasta_samples }

    CONVERT_FASTQ_TO_FASTA(fastq_samples)

    fasta_samples.mix(CONVERT_FASTQ_TO_FASTA.out) | (COUNT_READS & FILTER_QUALITY)

    RUN_BLAST(FILTER_QUALITY.out, params.db_files)

    MAP_TAXONOMIES(RUN_BLAST.out, params.taxonomy)

    DECISION_MAKING(MAP_TAXONOMIES.out, params.script).collect() | MERGE_FEATURES
}


// 1. (Python)
process COUNT_READS {
    container "pegi3s/biopython"
    publishDir "${params.outdir}/counts", mode: 'copy'

    input:
    tuple val(sample_name), path(fasta_path)

    output:
    tuple val(sample_name), path("*.txt")

    script:
    """
    python3 -c "
    from Bio import SeqIO
    format_type = 'fastq' if '${fasta_path}'.endswith('.fastq') else 'fasta'
    count = sum(1 for _ in SeqIO.parse('${fasta_path}', format_type))
    print(count)
    " \
    > ${sample_name}.txt
    """
}

// 2. (seqtk)
process CONVERT_FASTQ_TO_FASTA {
    container "staphb/seqtk"

    input:
    tuple val(sample_name), path(inputFile)

    output:
    tuple val(sample_name), path("converted_*.fasta")

    script:
    """
    seqtk seq -A ${inputFile} > converted_${sample_name}.fasta
    """
}

// 3. (Python)
process FILTER_QUALITY {
    container "pegi3s/biopython"

    input:
    tuple val(sample_name), path(fastaFile)

    output:
    tuple val(sample_name), path("filtered_*.fasta")

    script:
    """
    python3 -c "
    from Bio import SeqIO
    import numpy as np

    input_file = '${fastaFile}'
    output_file = 'filtered_${sample_name}.fasta'

    # Compute sequence lengths
    lengths = [len(record.seq) for record in SeqIO.parse(input_file, 'fasta')]

    # Define threshold dynamically
    median_length = np.median(lengths)
    threshold = 1200 if median_length >= 1200 else 150
    mode = 'FULL' if threshold == 1200 else 'PARTIAL'

    print(f'Detected mode: {mode} (Threshold: {threshold} bp)')

    # Filter sequences based on detected mode
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        for record in SeqIO.parse(infile, 'fasta'):
            if len(record.seq) >= threshold:
                SeqIO.write(record, outfile, 'fasta')

    print(f'Filtering complete: Keeping reads >= {threshold} bp')
    "
    """
}

// 4. (Blast)
process RUN_BLAST {
    container "ncbi/blast:latest"
    publishDir "${params.outdir}/${params.blast_results_dir}", mode: 'copy'

    input:
    tuple val(sample_name), path(filteredFasta)
    path db_files

    output:
    tuple val(sample_name), path("*.txt")

    script:
    """
    ulimit -n 65536
    blastn \
        -task megablast \
        -evalue 0.01 \
        -query ${filteredFasta} \
        -db ${db_files}/${params.db_name} \
        -out blast_results_${sample_name}.txt \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -outfmt "6 qseqid sseqid qlen evalue pident length qstart qend sstart send" \
        -num_threads 4
    """
}

// 5. (Python)
process MAP_TAXONOMIES {
    container "biocontainers/pandas:1.5.1_cv1"

    input:
    tuple val(sample_name), path(blastResults)
    path taxonomyFile

    output:
    tuple val(sample_name), path("results_with_taxonomy_${sample_name}.txt")

    script:
    """
    python -c "
    import pandas as pd
    print('Loading BLAST results for ${sample_name}...')
    blast_df = pd.read_csv('${blastResults}', sep='\t', header=None, names=['qseqid', 'sseqid', 'qlen', 'evalue', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send'])
    
    print('Calculating query coverage...')
    blast_df['query_coverage (%)'] = (blast_df['length'] / blast_df['qlen']) * 100
    
    print('Loading taxonomy file...')
    taxonomy_df = pd.read_csv('${taxonomyFile}', sep='\t', names=['sseqid', 'taxonomy'])
    
    print('Merging BLAST results with taxonomy...')
    merged_df = pd.merge(blast_df, taxonomy_df, on='sseqid', how='left')
    
    merged_df.to_csv('results_with_taxonomy_${sample_name}.txt', sep='\t', index=False)
    print(f'Taxonomy-mapped results saved for ${sample_name}')
    "
    """
}

// 6. (Python) -- custom container?
process DECISION_MAKING {
    container "biocontainers/pandas:1.5.1_cv1"
    publishDir "${params.outdir}/${params.decision_output_dir}", mode: 'copy'

    input:
    tuple val(sample_name), path(mappedResults)
    path scriptFile

    output:
    path "${sample_name}/*.txt"

    script:
    """
    mkdir -p ${sample_name}
    cd ${sample_name}
    python ../${scriptFile} ../${mappedResults}
    """
}

// 7. (Python)
process MERGE_FEATURES {
    container "biocontainers/pandas:1.5.1_cv1"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path "**/features.txt"

    output:
    path "merged_features.txt"

    script:
    """
    python -c "
    import pandas as pd
    import glob
    import os

    # Collect all features.txt files
    feature_files = glob.glob('**/features.txt')

    # Initialize empty list to store dataframes
    dfs = []

    for file in feature_files:
        # Extract sample name correctly by taking the last directory name before 'features.txt'
        sample_name = os.path.basename(os.path.dirname(file))

        # Read file with first row as header
        df = pd.read_csv(file, sep='\\t', header=0)
        
        df = df.set_index('Taxon ID')  # Set taxonomy as index
        dfs.append(df)

    # Merge all dataframes by Taxon ID (outer join to keep all taxonomies)
    merged_df = pd.concat(dfs, axis=1, sort=False).fillna(0).astype(int)

    # Save as TSV without unintended extra columns
    merged_df.to_csv('merged_features.txt', sep='\\t')

    print(f'Merged features table saved as merged_features.txt')
    "
    """
}
