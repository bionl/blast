#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

params.input = "samples.tsv"
params.taxonomy = "taxonomy.tsv"
params.script = "blastdm.py"
params.db_name_prefix = "reference_db"
params.db_title = "Demo BLAST DB"
params.blast_results_dir = "blast_results"
params.decision_output_dir = "blastdm_outputs"
params.final_merged_features = "merged_features.txt"
params.read_type = "PARTIAL"

workflow {
    Channel
        .fromPath(params.input)
        .splitCsv(header: true, sep: '\t')
        .set { sample_info }

    sample_info | COUNT_READS

    sample_info
        | convertFastqToFasta
        | filterQuality
        | MAKE_BLAST_DB
        | RUN_BLAST
        | mapTaxonomies
        | decisionMaking
        | mergeFeatures
}
// Set to "FULL" for full-length sequences


// 1. (Python)
process COUNT_READS {
    input:
    tuple val(sample_name), path(inputFile)

    output:
    tuple val(sample_name), val("read_counts_*.txt")

    script:
    """
    python -c "
    from Bio import SeqIO
    format_type = 'fastq' if '${inputFile}'.endswith('.fastq') else 'fasta'
    count = sum(1 for _ in SeqIO.parse('${inputFile}', format_type))
    print(f'{sample_name}\\t{count}')
    "
    > read_counts_${sample_name}.txt
    """
}

// 2. (seqtk)
process convertFastqToFasta {
    input:
    tuple val(sample_name), path(inputFile)

    output:
    tuple val(sample_name), path("converted_${sample_name}.fasta")

    when:
    inputFile.endsWith('.fastq')

    script:
    """
    seqtk seq -A ${inputFile} > converted_${sample_name}.fasta
    """
}

// 3. (Python)
process filterQuality {
    input:
    tuple val(sample_name), path(fastaFile)

    output:
    tuple val(sample_name), path("filtered_${sample_name}.fasta")

    script:
    """
    python -c "
    from Bio import SeqIO
    import numpy as np

    input_file = '${fastaFile}'
    output_file = 'filtered_${sample_name}.fasta'

    # Compute sequence lengths
    lengths = [len(record.seq) for record in SeqIO.parse(input_file, 'fasta')]

    # Define threshold dynamically
    median_length = np.median(lengths)
    threshold = 1200 if median_length >= 1200 else 150  # FULL if median length >= 1200, else PARTIAL

    print(f'Detected mode: {"FULL" if threshold == 1200 else "PARTIAL"} (Threshold: {threshold} bp)')

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
process MAKE_BLAST_DB {
    container "ncbi/blast:latest"

    input:
    tuple val(sample_name), path(filteredFasta)

    output:
    tuple val(sample_name), path("${params.db_name_prefix}_${sample_name}.nhr")

    script:
    """
    makeblastdb \
        -in ${filteredFasta} \
        -dbtype nucl \
        -out ${params.db_name_prefix}_${sample_name} \
        -parse_seqids \
        -title "${params.db_title}_${sample_name}"
    """
}

// 5. (Blast)
process RUN_BLAST {
    container "ncbi/blast:latest"

    input:
    tuple val(sample_name), path(filteredFasta)

    output:
    tuple val(sample_name), path("${params.blast_results_dir}/blast_results_${sample_name}.txt")

    script:
    """
    blastn \
        -task megablast \
        -evalue 0.01 \
        -query ${filteredFasta} \
        -db ${params.db_name_prefix}_${sample_name} \
        -out ${params.blast_results_dir}/blast_results_${sample_name}.txt \
        -max_target_seqs 10 \
        -max_hsps 1 \
        -outfmt "6 qseqid sseqid qlen evalue pident length qstart qend sstart send"
        -num_threads 4
    """
}

// 6. (Python)
process mapTaxonomies {
    input:
    tuple val(sample_name), path(blastResults)
    path taxonomyFile

    output:
    tuple val(sample_name), path("results_with_taxonomy_${sample_name}.txt")

    script:
    """
    python -c "
    import pandas as pd;
    print('Loading BLAST results for ${sample_name}...');
    blast_df = pd.read_csv('${blastResults}', sep='\t', header=None, names=['qseqid', 'sseqid', 'qlen', 'evalue', 'pident', 'length', 'qstart', 'qend', 'sstart', 'send']);
    
    print('Calculating query coverage...');
    blast_df['query_coverage (%)'] = (blast_df['length'] / blast_df['qlen']) * 100;
    
    print('Loading taxonomy file...');
    taxonomy_df = pd.read_csv('${taxonomyFile}', sep='\t', names=['sseqid', 'taxonomy']);
    
    print('Merging BLAST results with taxonomy...');
    merged_df = pd.merge(blast_df, taxonomy_df, on='sseqid', how='left');
    
    merged_df.to_csv('results_with_taxonomy_${sample_name}.txt', sep='\t', index=False);
    print(f'Taxonomy-mapped results saved for ${sample_name}');
    "
    """
}

// 7. (Python) -- custom container?
process decisionMaking {
    input:
    tuple val(sample_name), path(mappedResults)
    path scriptFile

    output:
    tuple val(sample_name), path("${params.decision_output_dir}/${sample_name}/features.txt")

    script:
    """
    mkdir -p ${params.decision_output_dir}/${sample_name}
    python ${scriptFile} ${mappedResults} --output ${params.decision_output_dir}/${sample_name}
    """
}

// 8. (Python)
process mergeFeatures {
    input:
    // path(featureFiles) from "${params.decision_output_dir}/*/features.txt"
    path featureFiles

    output:
    path params.final_merged_features

    script:
    """
    python -c "
    import pandas as pd
    import glob
    import os

    # Collect all features.txt files
    feature_files = glob.glob('${params.decision_output_dir}/*/features.txt')

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
    merged_df.to_csv('${params.final_merged_features}', sep='\\t')

    print(f'Merged features table saved as ${params.final_merged_features}')
    "
    """
}
