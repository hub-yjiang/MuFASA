process PASSfilter {
  label 'mutect_flag'
  shell = ['/bin/bash', '-euo', 'pipefail']
  conda '/groups/group-garaycoechea/linda/envs/pipeline'
  publishDir "${params.filtered_dir}", mode: 'copy'  
  executor 'slurm'

  input:
    tuple(val(sample_id), path(raw_vcf), val(tool), val(type)) //-> change to path..

  output:
    tuple(val(sample_id),path("*_PASS.vcf.g*"),val(tool), val(type))

  script:
      """
      bcftools view --threads ${task.cpus} -O z -o ${sample_id}.${tool}.${type}_PASS.vcf.gz -f PASS ${raw_vcf}
      bcftools index ${sample_id}.${tool}.${type}_PASS.vcf.gz 
      """

}

process ISEC_overlap {
    label 'isec_overlap'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    publishDir "${params.filtered_dir}", mode: 'copy'
    executor 'slurm'

    input:
      tuple(val(sample_id),path(strelka),path(strelka_index),path(mutect), path(mutect_index), val(type))

    output:
      tuple(val(sample_id),path("${sample_id}_shared_${type}.vcf"), val(type))

    script:
      """
      bcftools isec --threads ${task.cpus} -p isec_${sample_id} ${mutect} ${strelka}
      cp isec_${sample_id}/0002.vcf ${sample_id}_shared_${type}.vcf
      """
}

process PoN_filter {
    label 'isec_overlap'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    publishDir "${params.filtered_dir}", mode: 'copy'
    executor 'slurm'

    input:
      tuple(val(sample_id),path(shared),val(type))

    output:
      tuple(val(sample_id),path("${sample_id}_pon_${type}.vcf"), val(type))

    script:
      """
      if [ ! -f ${shared}.gz.csi ]; then
        bgzip ${shared}
        bcftools index ${shared}.gz
      fi

      bcftools isec --threads ${task.cpus} -O v -p isec_${sample_id} -C ${shared}.gz ${params.pon}
      cp isec_${sample_id}/0000.vcf ${sample_id}_pon_${type}.vcf

      """
}

process FiNGS {
    label 'fings'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/fings'
    publishDir "${params.snvs_filtered}", mode: 'copy'
    executor 'slurm'

    input:
         tuple(val(sample_id),path(pon),path(tumorbam),path(tumorbai),path(normalbam), path(normalbai))

    output:
        val(sample_id)

    script:
      """
        fings -j ${task.cpus} --overwrite -r ${params.ref} -p ${params.fings} -n ${normalbam} -t ${tumorbam} -v ${pon} -d ${params.snvs_filtered}/${sample_id}_fings/ --PASSonlyin --PASSonlyout

      """
}

process handleFings {
    label 'copy'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    publishDir "${params.snvs_filtered}", mode: 'copy'
    executor 'local'

    input:
        val(sample_id)
      
    output:
        val(sample_id)

    script:
      """
        cp ${params.snvs_filtered}/${sample_id}_fings/*.filtered.vcf ${params.snvs_filtered}/${sample_id}_fings.vcf
        cp ${params.snvs_filtered}/${sample_id}_fings/plots.pdf ${params.pdf_path}/${sample_id}_fings.pdf
      """
}

process snvsVAF {
    label 'vaf'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/miniforge3/envs/scripts'
    publishDir "${params.snvs_filtered}", mode: 'copy'
    executor 'local'

    input:
         val(sample_id)

    output:
        tuple(val(sample_id),path("${sample_id}.vaf.vcf")) 

    script:
      """
      python ${params.script_dir}/filter_mutect2_allelic_frequency.py ${params.snvs_filtered}/${sample_id}_fings.vcf ${sample_id}.ALLvaf.tab ${sample_id}.vaf.vcf ${sample_id}.PASSvaf.tab ${params.minvaf} ${params.maxvaf}

      """
}

process indelFiltering {
    label 'indel'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/miniforge3/envs/scripts'
    publishDir "${params.indel_filtered}", mode: 'copy'
    executor 'local'

    input:
        tuple val(sample_id),path(pon),val(type) 

    output:
        tuple(val(sample_id),path("${sample_id}.NoAlt.vcf"),path("${sample_id}.SBtolerant.vcf"), path("${sample_id}.SBzero.vcf"), path("${sample_id}.RF.vcf"),path("${sample_id}.AF.vcf"),path("${sample_id}.AF.csv"))

    script:
      """
      python ${params.script_dir}/filter_shared_indels.py ${pon} ${sample_id} ${params.minvaf} ${params.maxvaf}

      """

}
//  output gridss tuple( val(tumor_sample_id), path("${tumor_sample_id}.gridss.vcf.gz"), path("${tumor_sample_id}.gridss.vcf.gz.tbi"))
// conda activate Renv
// Rscript --vanilla gridss_somatic_filter.R --ref BSgenome.Hsapiens.UCSC.hg38 --input M1R7_FU3_C3.gridss.driver.vcf --output out.vcf --plotdir . --scriptdir .  
// zcat out.vcf.bgz > out.vcf
//[JG0163m, /groups/group-garaycoechea/linda/MuFASA/SBS17_Human_diploid/nextflow/work/3f/9d6c6ee70a4dbd02128ea819581c7e/JG0163m.gridss.vcf.gz, /groups/group-garaycoechea/linda/MuFASA/SBS17_Human_diploid/nextflow/work/3f/9d6c6ee70a4dbd02128ea819581c7e/JG0163m.gridss.vcf.gz.tbi]
process Gridss_filter {
    label 'gridss_filter'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/miniforge3/envs/Renv'
    publishDir "${params.gridss}", mode: 'copy'  
    executor 'slurm'
    //errorStrategy 'ignore'

    input:
        tuple val(sample_id), path(input) 

    output:
        tuple val(sample_id),path("${sample_id}.somatic.filter.vcf.bgz")

    script:
      """
      Rscript --vanilla ${params.script_dir}/gridss_filter/gridss_somatic_filter_LB.R --ref ${params.gridss_ref} --input ${input} --output ${sample_id}.somatic.filter.vcf --plotdir ${params.gridss} --scriptdir ${params.script_dir}/gridss_filter/  

      """
}

process GridssManta_validate {
    label 'indel'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/miniforge3/envs/scripts'
    publishDir "${params.structural}", mode: 'copy'
    executor 'slurm'
    errorStrategy 'ignore'

    input:
        tuple val(sample_id), path(gridss), path(manta) 

    output:
        tuple val(sample_id), 
        path("${sample_id}.gridss_manta_unmatched.bed"),
        path("${sample_id}.all_structural_overlap.tsv"),
        path("${sample_id}.final.bedpe"),
        path("${sample_id}.StructuralVariants.csv"),
        path("${sample_id}.structural.final.vcf")

    script:
      """
      zcat ${gridss} > ${sample_id}.gridss.filtered.vcf
      python ${params.script_dir}/Validate_Manta_Gridss_intersect.py ${sample_id} ${manta} ${sample_id}.gridss.filtered.vcf

      """


}


process Manta_pass {
    label 'pass'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    publishDir "${params.manta}"
    executor 'slurm'

    input:
        val(tumor_id)
    output:
        tuple val(tumor_id), path("${tumor_id}.Manta_somaticSV.PASS.vcf")

    script:
      """
      bcftools view --threads ${task.cpus} -O v -o ${tumor_id}.Manta_somaticSV.PASS.vcf -f PASS ${params.manta}/${tumor_id}_runworkflow/results/variants/somaticSV.vcf.gz

      """
}

process Manta_encode {
    label 'pass'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    publishDir "${params.manta}"
    executor 'slurm'

    input:
        tuple val(tumor_id), path(input)

    output:
        tuple val(tumor_id), path("${tumor_id}.Manta.bl.filtered.vcf")

    script:
      """
      bcftools view --threads ${task.cpus} -T ^${params.blacklist} -O v -o ${tumor_id}.Manta.bl.filtered.vcf ${input}
      """
}

process Manta_filtering {
    label 'indel'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/miniforge3/envs/scripts'
    publishDir "${params.structural}", mode: 'copy'
    //  publishDir "${params.manta}"
    executor 'slurm'

    input:
        tuple val(sample_id), path(bl_input) 

    output:
        tuple val(sample_id), path("${sample_id}.Manta.filtered.vcf") , path("${sample_id}.Manta.regions.bed")

    script:
      def chrs = params.chrs.join(',')
      """
      python ${params.script_dir}/Manta_filtering.py ${bl_input} ${sample_id} ${chrs}

      """

}


process Concat_variants {
    label 'indel'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/pipeline'
    errorStrategy 'ignore'
    publishDir "${params.snvs_filtered}", mode: 'copy'
    //  publishDir "${params.manta}"
    executor 'slurm'

    input:
        tuple val(sample_id), path(id_final), path(snv_final) 

    output:
        tuple val(sample_id), path("${sample_id}.concatenated.vcf")

    script:
      
      """
      bgzip ${id_final} 
      bgzip ${snv_final}

      bcftools index ${id_final}.gz
      bcftools index ${snv_final}.gz
      
      bcftools concat -a -Ov -o ${sample_id}.concatenated.vcf ${id_final}.gz ${snv_final}.gz

      """

}


