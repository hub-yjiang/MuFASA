
include { MergeFastq; FastQC; CUTadapt} from '/groups/group-garaycoechea/Yang/MuFASA_YJ/nextflow/modules/FASTQprocessing.nf'
include { BWAMapping; Samtools_index; Samtools_sort; Samtools_fixmate; Picard_cleansam; Samtools_markdup; Samtools_bamtocram; Samtools_processing} from '/groups/group-garaycoechea/Yang/MuFASA_YJ/nextflow/modules/Mapping.nf'
include { WgsMetrics; AlignmentSummary; PlotVAF; VariantCounts } from '/groups/group-garaycoechea/Yang/MuFASA_YJ/nextflow/modules/Statistics.nf'
include { preVariantCalling; add_indel_snvs_tags; defineOrderISEC } from '/groups/group-garaycoechea/Yang/MuFASA_YJ/nextflow/modules/processing.nf'
include { Strelka; handleStrelka; Mutect2; Mutect2_flag; Mutect2_concat; Manta; Gridss_extract; Gridss_index; Gridss } from '/groups/group-garaycoechea/Yang/MuFASA_YJ/nextflow/modules/VariantCalling.nf'
include { PASSfilter; ISEC_overlap; PoN_filter; FiNGS; snvsVAF; handleFings; indelFiltering; Manta_pass; Manta_encode; Manta_filtering; Gridss_filter;  GridssManta_validate; Concat_variants} from '/groups/group-garaycoechea/Yang/MuFASA_YJ/nextflow/modules/VariantFiltering.nf'

params.cram = false
params.bam_only = false
params.start_from = "fastq"
params.ext_bam = null

main: workflow {
	
	if (params.start_from == 'fastq' && params.bam_only) {
    	
		names_paths = Channel.fromPath(params.samples)
				   		 .splitCsv()//.view()
	
		pair1 = names_paths.map{ it -> [it[0],params.fastq_path + it[1],"1"]}.groupTuple(by: 0)
		pair2 = names_paths.map{ it -> [it[0],params.fastq_path + it[2],"2"]}.groupTuple(by: 0)
	
		merge_input = pair1.concat(pair2).map{ it -> [it[0],it[1],it[2][0]]}
	
		// ------------FASTQ processing
		MergeFastq(merge_input) 
		merge_paired = MergeFastq.out.groupTuple(by:0)
		CUTadapt(merge_paired)
		fastqc_input = CUTadapt.out.flatMap{ sample_id, file_list -> file_list}
		FastQC(fastqc_input)
	
		// --------------MAPPING-------------------
		BWAMapping(CUTadapt.out) 
		Picard_cleansam(BWAMapping.out)
		Samtools_processing(Picard_cleansam.out)
		//Samtools_fixmate(Picard_cleansam.out)
		//Samtools_sort(Samtools_fixmate.out)
		//Samtools_markdup(Samtools_sort.out)
		bam_indexed = Samtools_processing.out
	

	
		// -------------MAPPING METRICS -----------
		
		//input_Metrics = Samtools_processing.out.map { sample_id, bam, bai -> tuple(sample_id, bam)}

		//WgsMetrics(input_Metrics)
		//AlignmentSummary(input_Metrics)

    } else {


		if (params.start_from == 'fastq') {
	
	
			names_paths = Channel.fromPath(params.samples)
					   		 .splitCsv()//.view()
		
			pair1 = names_paths.map{ it -> [it[0],params.fastq_path + it[1],"1"]}.groupTuple(by: 0)
			pair2 = names_paths.map{ it -> [it[0],params.fastq_path + it[2],"2"]}.groupTuple(by: 0)
		
			merge_input = pair1.concat(pair2).map{ it -> [it[0],it[1],it[2][0]]}
		
			// ------------FASTQ processing
			MergeFastq(merge_input) 
			merge_paired = MergeFastq.out.groupTuple(by:0)
			CUTadapt(merge_paired)
			fastqc_input = CUTadapt.out.flatMap{ sample_id, file_list -> file_list}
			FastQC(fastqc_input)
		
			// --------------MAPPING-------------------
			BWAMapping(CUTadapt.out) 
			Picard_cleansam(BWAMapping.out)
			Samtools_processing(Picard_cleansam.out)
			//Samtools_fixmate(Picard_cleansam.out)
			//Samtools_sort(Samtools_fixmate.out)
			//Samtools_markdup(Samtools_sort.out)
			bam_indexed = Samtools_processing.out
	

	
			// -------------MAPPING METRICS -----------
		
			//input_Metrics = Samtools_processing.out.map { sample_id, bam, bai -> tuple(sample_id, bam)}

			//WgsMetrics(input_Metrics)
			//AlignmentSummary(input_Metrics)
	
	
			if (params.cram) {
    			// Only run CRAM conversion and stop
    		    Samtools_bamtocram(Samtools_markdup.out)
    			//error "Pipeline stopped intentionally as cram = true"
    			}
		
		} else if (params.start_from == 'bam') {
			
			bam_indexed = Channel.fromPath(params.samples)
    				.splitCsv()
    				//.map { row -> [ row.sampleID, file(row.bam), file(row.bai) ]}
		}	
			
		

		// -------------VARIANT CALLING -----------
		// prepare channels for variant calling:
		chrs = Channel.from(params.chrs)//.view()
		sample_pairs = Channel
				.fromPath(params.strelka_conf)
				.splitCsv()//.view()
	
		if (params.ext_bam) {
			extNormal = Channel
				.fromPath(params.ext_bam)
				.splitCsv()


			bam_indexed = bam_indexed.concat(extNormal)//.view()  
		
		}

	
	
		


		paired = bam_indexed.combine(sample_pairs,by:0).groupTuple(by: 4, size: 2)//.view()
		//ensure correct order: tumor normal
		ordered_paired = preVariantCalling(paired)
		
		// STRELKA variant calling//
		strelka_input = ordered_paired.map{ it -> [it[0][0],it[1][0],it[1][1],it[2][0],it[2][1]]}
		Strelka(strelka_input)
		strelka_out = handleStrelka(Strelka.out)
	
		// MUTECT //
		mutect_input = ordered_paired.map{ it -> [it[0][0],it[0][1],it[1][0],it[1][1]]}
	
		mut_paired = mutect_input.combine(chrs)
		Mutect2(mut_paired) //--> TO DO: seperate this, so Chr1, Chr2 are done with different time limit
		// Mutect2.out.view()
	
		Mutect2_flag(Mutect2.out) // --> out: tumorID, *flagged.vcf.gz, chrID
		concat_input = Mutect2_flag.out.groupTuple(by: 0, size: params.nrChrs) // only start concatenate when all Chromosomes are done.
		Mutect2_concat(concat_input)//.view() // --> out: tumorID, indel.vcf.gz, snvs.vcf.gz
		
		// -------------VARIANT FILTERING -----------
	
		// WITHOUT MUTECT //
		//transposed = strelka_out.map{ it -> [it[0],it.drop(1)]}.transpose()
		
		// WITH MUTECT //
		// All output from Mutect and Strelka are combined to be able to go through filtering.	
		all_raw_variants = strelka_out.concat(Mutect2_concat.out)//.view()
		transposed = all_raw_variants.map{ it -> [it[0],it.drop(1)]}.transpose()// sampleid, path(.vcf.gz) both strelka and Mutect
	
		pass_input = add_indel_snvs_tags(transposed)
		PASSfilter(pass_input)
		isec_grouped = PASSfilter.out.groupTuple(by: [0,3], size: 2)//.view()
		  
		// rearrange for isec: strelka_mutect combinations and correct order: [sampleid, strelka_*, strelka_*.bai, mutect_*, mutect_*.bai, indel/snvs]
		isec_input = defineOrderISEC(isec_grouped)
		ISEC_overlap(isec_input)
		PoN_filter(ISEC_overlap.out)
	
		// ------------ SNV filtering ---------
	
		pon_snvs = PoN_filter.out.filter{ it.contains("snvs")}
		fings_input = pon_snvs.combine(strelka_input,by:0).map{ it -> [it[0],it[1],it[3],it[5],it[4],it[6]]}// -> id, vcffile, tumorbam, tumorbai, normalbam, normalbai
		FiNGS(fings_input)//.view()
		handleFings(FiNGS.out) // --> sample id, _fings.vcf
		snvsVAF(handleFings.out)  // --> out: sample_id, ${sample_id}.vaf.vcf
		
		toplot = pon_snvs.map{ listItem -> [listItem[0],listItem[1] , "PoN", "red"]}
		  .concat(handleFings.out.map{ listItem -> [listItem,params.snvs_filtered + "/" + listItem + "_fings.vcf", "FiNGS", "green"]})
		  .concat(snvsVAF.out.map{ listItem -> listItem + ["VAF", "blue"]})
		
		PlotVAF(toplot)
	
		// ------------ INDEL filtering ------
	
		pon_indels = PoN_filter.out.filter{ it.contains("indel")}
		indelFiltering(pon_indels)
	
		// ------------ MANTA Structural Variant Calling ------
	
		manta_input = ordered_paired.map{ it -> [it[0][0],it[0][1],it[1][0],it[1][1],it[2][0],it[2][1]]}
		Manta(manta_input)
		Manta_pass(Manta.out)
		Manta_encode(Manta_pass.out)	
		manta_out = Manta_filtering(Manta_encode.out)
	
		// ----------- Gridss Structural Variant Calling ------
	
		gridss_input = manta_input.combine(manta_out,by:0).map{ it -> [it[0],it[2],it[4],it[7]]}
		// gridss_input.view()
	
	
    	Gridss_extract(gridss_input)
    	region_index = Gridss_index(Gridss_extract.out)
    	gridss_in = manta_input.combine(region_index,by:0).map{ it -> [it[0],it[1],it[6],it[3],it[7],it[5]]}//.view()
    	Gridss(gridss_in)  
    	Gridss_filter(Gridss.out)
    	input_intersect = Gridss_filter.out.combine(manta_out,by:0).map{ it -> [it[0],it[1],it[2]]} // -> id. gridss, manta
	
		// ---------- Intersect Manta & Gridss output ------
	
    	GridssManta_validate(input_intersect)

    	id_final = indelFiltering.out.map { sample_id, noalt, sbtol, sbzero, rf, af_vcf, af_csv -> tuple(sample_id, af_vcf)}
    	to_concat = id_final.join(snvsVAF.out)

    	Concat_variants(to_concat)
	
		VariantCounts(FiNGS.out)
	}

}