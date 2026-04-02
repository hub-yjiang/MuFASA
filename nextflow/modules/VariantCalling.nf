process Strelka {

    label 'strelka'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/strelka'
    publishDir "${params.strelka}", mode: 'copy'
    executor 'slurm'

    input:
        tuple val(tumorID), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai)

    output:
        val(tumorID) 

    script:
      """
      rm -rf ${params.strelka}/${tumorID}_runworkflow/
      configureStrelkaSomaticWorkflow.py --normalBam ${normal_bam} --tumorBam ${tumor_bam} --referenceFasta ${params.ref} --outputCallableRegions --runDir ${params.strelka}/${tumorID}_runworkflow/
      python2 ${params.strelka}/${tumorID}_runworkflow/runWorkflow.py -j ${task.cpus} -m local -g 10
      """
}

process handleStrelka {
    label 'mutect_flag'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    publishDir "${params.strelka}", mode: 'copy'
    executor 'slurm'

  input:
    val(sampleid)
    
  output:
    tuple val(sampleid), path("*indels.vcf.gz"), path("*snvs.vcf.gz")

  script:
      """
      cp ${params.strelka}/${sampleid}_runworkflow/results/variants/somatic.indels.vcf.gz ${sampleid}.somatic.indels.vcf.gz
      cp ${params.strelka}/${sampleid}_runworkflow/results/variants/somatic.snvs.vcf.gz ${sampleid}.somatic.snvs.vcf.gz
      """
}

process Mutect2 {

  label 'mutect'
  shell = ['/bin/bash', '-euo', 'pipefail']
  conda '/groups/group-garaycoechea/linda/envs/gatk4'
//  publishDir "${params.mutect}"  
  executor 'slurm'

  input:
    //[[a1d2, a1, a1d2_markdups.CRAM, a1_markdups.CRAM, 12]  
    tuple val(tumorid), val(normalid), val(tumor_bam), val(normal_bam), val(chrID)

  output:
    tuple val(tumorid),path("${tumorid}_${chrID}.vcf.gz"), val(chrID)

  script:
      """
      gatk Mutect2 --native-pair-hmm-threads 4 -R ${params.ref} --callable-depth ${params.callable_depth} --max-mnp-distance 0 -I ${tumor_bam} -normal ${normalid} -I ${normal_bam} -L chr${chrID} -O ${tumorid}_${chrID}.vcf.gz
      """
}

process Mutect2_flag {

  label 'mutect_flag'
  shell = ['/bin/bash', '-euo', 'pipefail']
  conda '/groups/group-garaycoechea/linda/envs/gatk4'
  publishDir "${params.mutect}", mode: 'copy'
  executor 'slurm'

  input:
    tuple val(tumorid),val(variants),val(chrID)

  output:
    tuple val(tumorid),path("${tumorid}_${chrID}.flagged.vcf.gz"), val(chrID)

  script:
      """
      gatk FilterMutectCalls -R ${params.ref} -V ${variants} -O ${tumorid}_${chrID}.flagged.vcf.gz
      """
}

process Mutect2_concat {
  label 'mutect_flag'
  shell = ['/bin/bash', '-euo', 'pipefail']
  conda '/groups/group-garaycoechea/linda/envs/pipeline'
  publishDir "${params.mutect}", mode: 'copy'
  executor 'slurm'

  input:
    tuple val(tumorid), path(flagged_files), val(chrIDs)

  output:
    tuple val(tumorid), path("*indels.vcf.gz"), path("*snvs.vcf.gz")

  script:
      """
      bcftools concat -O v -o ${tumorid}_final.vcf ${flagged_files}
      bcftools view -v snps -O z -o ${tumorid}_Mutect.snvs.vcf.gz ${tumorid}_final.vcf
      bcftools view -v indels -O z -o ${tumorid}_Mutect.indels.vcf.gz ${tumorid}_final.vcf
      """
}


process Gridss {
  label 'gridss_run'
  container = 'docker://gridss/gridss:2.13.2'
  shell = ['/bin/bash', '-euo', 'pipefail']
  publishDir "${params.gridss}", mode: 'copy'
  executor 'slurm'

  input:
    tuple( val(sample_id), val(normal_id), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai) )

  output:
    tuple( val(sample_id), path("${sample_id}.gridss.vcf.gz"))//, path("${sample_id}.gridss.vcf.gz.tbi"))

  script:
    """
    gridss --jvmheap 50g -s preprocess,assemble,call -w . -o ${sample_id}.gridss.vcf.gz -r ${params.ref} -t ${task.cpus} --labels ${normal_id},${sample_id} --maxcoverage 200 ${normal_bam} ${tumor_bam}
    """
  

}

//process gridss_diff {
//  label 'gridss'
//  container = 'docker://gridss/gridss:2.13.2'
//  shell = ['/bin/bash', '-euo', 'pipefail']
//  publishDir "${params.gridss_output}", mode: 'copy'
//  executor 'slurm'

//  input:
//    tuple( val(sample_id), path(bam), path(bai), val(type) )

//  output:
//    tuple( val(sample_id), path("${sample_id}.gridss.vcf.gz"))

// script:
//    """
//    gridss --jvmheap 50g -s preprocess,assemble,call -o ${sample_id}.gridss.vcf.gz -r ${params.ref} -t ${task.cpus} --labels ${type} --maxcoverage 200 ${bam}
//    """
//}

process Gridss_extract {
    label 'gridss'
    container = 'docker://gridss/gridss:2.13.2'
    shell = ['/bin/bash', '-euo', 'pipefail']
    publishDir "${params.gridss}", mode: 'copy'
    executor 'slurm'

    input:
      tuple( val(sample_id), path(bam), path(bai), path(bed) )

    output:
      tuple( val(sample_id), path("${sample_id}.regions.bam") )

    script:
      """

        if [ ! -s "${bed}" ]; then
        echo "WARNING: BED file ${bed} is empty for ${sample_id}. gridss_extract skipped"
        
        # Option A: Create a dummy BAM 
        samtools view -H ${bam} | samtools view -b - > ${sample_id}.regions.bam
        
        exit 0
    fi

        gridss_extract_overlapping_fragments -t ${task.cpus} -w ${params.gridss} -o ${sample_id}.regions.bam --targetbed ${bed} ${bam}
      
      """
}

process Gridss_index {
    label 'samtools_index'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/samtools_picard'
    publishDir params.gridss, mode: 'copy'
    errorStrategy 'finish'
    
    input:
        tuple(val(sample_id), path(unindexed_bam))

    output:
        tuple(val(sample_id), path("${sample_id}.regions.bam"), path("${sample_id}.regions.bai"))

    script:
        """
        samtools index -@ ${task.cpus} -o ${sample_id}.regions.bai ${unindexed_bam}
        """


}

process Manta {
    label 'manta'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/miniforge3/envs/manta'
    publishDir "${params.manta}", mode: 'copy'
    executor 'slurm'
    errorStrategy 'ignore'

    input:
        tuple( val(tumor_id), val(normal_sample_id), path(tumor_bam), path(normal_bam), path(tumor_bai), path(normal_bai) )

    output:
        val(tumor_id)

    script:
      """
      rm -rf ${params.manta}/${tumor_id}_runworkflow/
      python /groups/group-garaycoechea/miniforge3/envs/manta/bin/configManta.py --normalBam ${normal_bam} --tumorBam ${tumor_bam} --referenceFasta ${params.ref} --runDir ${params.manta}/${tumor_id}_runworkflow/
      python ${params.manta}/${tumor_id}_runworkflow/runWorkflow.py -j ${task.cpus} -g 30
      """
}



