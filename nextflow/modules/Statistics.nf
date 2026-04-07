
process Depth {
    label 'depth_label'
    shell = ['/bin/bash', '-euo', 'pipefail']
    //publishDir '/hpc/hub_garayco/projects/hepg2_OC2_redo/Stats/'
    conda '/groups/group-garaycoechea/linda/envs/pipeline'

    input:
        tuple(val(sample_id), path(input))

    output:
        tuple(val(sample_id), path("${sample_id}_stats.txt"))

    shell:
        '''
        samtools depth -a !{input} | awk '{sum+=$3; sumsq+=$3*$3} END { print "Average = ",sum/NR; print "Stdev = ",sqrt(sumsq/NR - (sum/NR)**2)}' > !{sample_id}_stats.txt
        '''
}

process WgsMetrics {
    label 'wgs_metrics'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/linda/envs/gatk4'
    publishDir params.fastqc_path, mode: 'copy'
    errorStrategy 'ignore'

    input:
        tuple(val(sample_id), path(input))

    output:
        tuple val(sample_id), path("${sample_id}_collect_wgs_metrics.txt")


    script:
        """
        picard CollectWgsMetrics R=${params.ref} I=${input} O=${sample_id}_collect_wgs_metrics.txt VALIDATION_STRINGENCY=SILENT
    
        """
}

process AlignmentSummary {
    label 'Align_sum'
    shell = ['/bin/bash', '-euo', 'pipefail']    
    conda '/groups/group-garaycoechea/linda/envs/gatk4'
    publishDir params.fastqc_path, mode: 'copy'
    errorStrategy 'ignore'

    input:
        tuple(val(sample_id), path(input))

    output:
         tuple(val(sample_id), path("${sample_id}_collect_align_metrics.txt"))

    script:
        """
        picard CollectAlignmentSummaryMetrics R=${params.ref} I=${input} O=${sample_id}_collect_align_metrics.txt VALIDATION_STRINGENCY=SILENT
        
        """

}

process CoverageOverview {
    label 'plot_vaf'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/miniforge3/envs/scripts'
    errorStrategy 'ignore'
    publishDir "${params.pdf_path}", mode: 'copy'
    //executor 'local'

    input:
        val(samples)

    output:
        path("${project_name}_picard_wgs_plots.pdf")

    script:
        """
        python ${params.script_dir}/coverage_plots.py ${params.project_name} "${samples}" "${params.fastqc_path}" "${params.pdf_path}"
        """
}



process PlotVAF {
    label 'plot_vaf'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/miniforge3/envs/scripts'
    errorStrategy 'ignore'
    publishDir "${params.pdf_path}", mode: 'copy'

    input:
        tuple(val(sample_id), path(input), val(type), val(color))

    output:
        tuple(val(sample_id), path("${sample_id}_${type}_VAFplot.pdf"))

    script:
        """
        python ${params.script_dir}/plot_VAF_mutect.vcf.py ${input} ${color} ${type} ${sample_id}
        """


}

process VariantCounts {
    label 'plot_vaf'
    shell = ['/bin/bash', '-euo', 'pipefail']
    conda '/groups/group-garaycoechea/miniforge3/envs/scripts'
    errorStrategy 'ignore'
    publishDir "${params.project_dir}", mode: 'copy'

    input:
        val(sample_id)

    output:
        path("counts.tsv")

    script:
        """
        cd ${params.project_dir}
        ${params.script_dir}/count_one_sample.sh ${sample_id}
        """
}


//process MatrixGenerator {
  //  label 'plot_vaf'
  //  shell = ['/bin/bash', '-euo', 'pipefail']
  //  conda '/groups/group-garaycoechea/miniforge3/envs/alexandrov'
  //  errorStrategy 'ignore'
  //  publishDir "/groups/group-garaycoechea/linda/MuFASA/SBS17_Human_diploid/Matrix/", mode: 'copy'

  //  output:
    //    path("${project_dir}/Matrix/")
//}
