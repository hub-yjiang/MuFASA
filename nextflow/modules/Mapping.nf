process BWAMapping {
    
    label 'bwa_label'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    errorStrategy 'finish'

    input:
        tuple(val(sample_id), path(fastq))

    output:
        tuple(val(sample_id), path("${sample_id}.bam"))

    script:
        def bwa_readgroup = "\"@RG\\tID:${sample_id}\\tSM:${sample_id}\\tPL:ILLUMINA\\tLB:${sample_id}\\tPU:UNKNOWN\""
        """
        bwa mem -R $bwa_readgroup -K 10000000 -M -t ${task.cpus} ${params.ref} ${fastq} | \
samtools view -b --threads ${task.cpus} -o ${sample_id}.bam

        """
}
   
process Picard_cleansam {
    label 'picard_cleansam'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/samtools_picard'
    errorStrategy 'finish'

    input:
        tuple(val(sample_id), path(bam))

    output:
        tuple(val(sample_id), path("${sample_id}_cleaned.bam"))

    script:
        """
        picard CleanSam INPUT=${bam} OUTPUT=${sample_id}_cleaned.bam USE_JDK_DEFLATER=True USE_JDK_INFLATER=True
        """
}


process Samtools_processing {
    label 'Samtools_processing'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/samtools_picard'
    publishDir params.mapping_dir, mode: 'copy'
    errorStrategy 'finish'

    input:
        tuple(val(sample_id), path(cleaned_bam))

    output:
        //tuple(val(sample_id), path("${sample_id}_samtools_fixed.bam"))
        tuple(val(sample_id), path("${sample_id}_samtools_markdup.bam"), path("${sample_id}_samtools_markdup.bai"))
    script:
        """
        samtools fixmate -m -@ ${task.cpus} -O BAM ${cleaned_bam} ${sample_id}_samtools_fixed.bam
        
        samtools sort -O BAM --threads ${task.cpus} -o ${sample_id}_sorted.bam ${sample_id}_samtools_fixed.bam
        rm -f ${sample_id}_samtools_fixed.bam

        samtools markdup -@ ${task.cpus} -O BAM ${sample_id}_sorted.bam ${sample_id}_samtools_markdup.bam
        rm -f ${sample_id}_sorted.bam

        samtools index -@ ${task.cpus} -o ${sample_id}_samtools_markdup.bai ${sample_id}_samtools_markdup.bam
        
        """
}



process Samtools_fixmate {
    label 'samtools_fixmate'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/samtools_picard'
    errorStrategy 'finish'

    input:
        tuple(val(sample_id), path(cleaned_bam))

    output:
        tuple(val(sample_id), path("${sample_id}_samtools_fixed.bam"))

    script:
        """
        samtools fixmate -m -@ ${task.cpus} -O BAM ${cleaned_bam} ${sample_id}_samtools_fixed.bam
        """
}

process Samtools_sort {

    label 'samtools_sort'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/samtools_picard'
    errorStrategy 'finish'

    input:
        tuple(val(sample_id), path(fixed_bam))

    output:
        tuple(val(sample_id), path("${sample_id}_sorted.bam"))

    script:
        """
        samtools sort -O BAM --threads ${task.cpus} -o ${sample_id}_sorted.bam ${fixed_bam}

        """
}



process Samtools_markdup {
    label 'samtools_markdup'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/samtools_picard'
    publishDir params.mapping_dir, mode: 'copy'
    errorStrategy 'finish'
    
    input:
        tuple(val(sample_id), path(sorted_bam))

    output:
        tuple(val(sample_id), path("${sample_id}_samtools_markdup.bam"))

    script:
        """
        samtools markdup -@ ${task.cpus} -O BAM ${sorted_bam} ${sample_id}_samtools_markdup.bam
        """
}

process Samtools_index {
    label 'samtools_index'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/samtools_picard'
    publishDir params.mapping_dir, mode: 'copy'
    errorStrategy 'finish'
    
    input:
        tuple(val(sample_id), path(unindexed_bam))

    output:
        tuple(val(sample_id), path("${sample_id}_samtools_markdup.bam"), path("${sample_id}_samtools_markdup.bai"))

    script:
        """
        samtools index -@ ${task.cpus} -o ${sample_id}_samtools_markdup.bai ${unindexed_bam}
        """
}

process Samtools_bamtocram {
    label 'samtools_sort'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/Yang/miniforge3_YJ/envs/samtools_1.19.1'
    publishDir params.mapping_dir, mode: 'copy'
    errorStrategy 'finish'

    input:
        tuple(val(sample_id), path(dupmarked_bam))
   
    output:
        tuple(val(sample_id), path("${sample_id}_markdup.cram"))

    script:
    """
    samtools view -@ ${task.cpus} -C -T ${params.ref} -o ${sample_id}_markdup.cram ${dupmarked_bam}
    samtools index ${sample_id}_markdup.cram
    """
}