process MergeFastq {

    label 'merge_label'
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    //publishDir params.merged_dir

    input:
        tuple val(sample_id), path(fastq_files), val(pair)
    
    output:
        tuple(val(sample_id), path("*fq.gz")) 

    script:
    
    def is_multiple = (fastq_files instanceof List && fastq_files.size() > 1)

    // 2. Generate the Bash script based on the result
    if ( is_multiple ) {
        """
        echo "Merging files..."
        cat ${fastq_files} > ${sample_id}_${pair}.fq.gz
        """
    } else {
        """
        echo "Linking single file..."
        ln -s ${fastq_files} ${sample_id}_${pair}.fq.gz
        """
    }
}


process FastQC {
    label 'FastQC_label'
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    errorStrategy 'ignore'

    input:
        val trimmed_file

    script:
    """
    mkdir -p ${params.fastqc_path}
    fastqc -t ${task.cpus} -o ${params.fastqc_path} $trimmed_file

    """
}

process CUTadapt {
    label 'cutadapt'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/pipeline'

    input:
        tuple(val(sample_id), path(fastq))

    output:
        tuple(val(sample_id), path("trimmed_*fq.gz"))

    script:
        """
        cutadapt -a "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA" -A "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT" -o trimmed_${sample_id}_1.fq.gz -p trimmed_${sample_id}_2.fq.gz --max-n 15 --max-ee 24 --cores=${task.cpus} ${fastq}
        """
}

