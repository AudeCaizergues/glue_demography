# Uses ANGSD to estimate SFS and genotype likelihoods across the genome using all samples
# as input (i.e., a global dataset). Only really used to perform PCA of global samples. 


###############
#### SETUP ####
###############


rule create_bam_list_finalSamples_withoutRelated:
	"""
	Create text file with paths to BAMs, excluding samples with high alignment error rates
	"""
	input:
		bams = ancient(get_all_bams_withoutRelated(BAM_DIR)),
		highErr = rules.create_samples_to_remove.output.error_df
	output:
		'{0}/bam_lists/withoutRelated/finalSamples_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
	log: LOG_DIR + '/create_bam_list/withoutRelated/finalSamples_{site}_bam_list.log'
	run:
		import os
		import re
		import pandas as pd
		from contextlib import redirect_stderr
		with open(log[0], 'w') as stderr, redirect_stderr(stderr):
			bad_samples = pd.read_table(input.highErr, header=None).iloc[:,0].tolist()
			with open(output[0], 'w') as f:
				for bam in input.bams:
					search = re.search('^(.+)(?=_\w)', os.path.basename(bam))
					sample = search.group(1)
					if sample not in bad_samples:
						f.write('{0}\n'.format(bam))

		
##############################
#### GENOTYPE LIKELIHOODS ####
##############################
		
rule angsd_gl_allSamples_alldegenerates_withoutRelated:
	"""
	Estimate genotype likelihoods for all samples separately for each of 16 chromosomes using ANGSD.
	"""
	input:
		bams = rules.create_bam_list_finalSamples_withoutRelated.output,
		ref = rules.unzip_reference.output,
		sites = rules.split_angsd_sites_byChrom.output,
		idx = rules.angsd_index_sites_byChrom.output
	output:
		gls = '{0}/gls/allSamples_alldegenerates/{{site}}/withoutRelated/{{chrom}}/{{chrom}}_{{site}}_withoutRelated.beagle.gz'.format(ANGSD_DIR),
		mafs = '{0}/gls/allSamples_alldegenerates/{{site}}/withoutRelated/{{chrom}}/{{chrom}}_{{site}}_withoutRelated.mafs.gz'.format(ANGSD_DIR),
	log: LOG_DIR + '/angsd_gl_allSamples_alldegenerates/withoutRelated/{chrom}_{site}_angsd_gl_withoutRelated.log'
	container: 'library://james-s-santangelo/angsd/angsd:0.933' 
	params:
		out = '{0}/gls/allSamples_alldegenerates/{{site}}/withoutRelated/{{chrom}}/{{chrom}}_{{site}}_withoutRelated'.format(ANGSD_DIR),
		max_dp = ANGSD_MAX_DP
	threads: 12
	resources:
		mem_mb = lambda wildcards, attempt: attempt * 30000,
		time = '12:00:00'
	shell:
		"""
		NUM_IND=$( wc -l < {input.bams} );
		MIN_IND=$(( NUM_IND*50/100 ))
		angsd -GL 1 \
			-out {params.out} \
			-nThreads {threads} \
			-doGlf 2 \
			-doMajorMinor 1 \
			-SNP_pval 1e-6 \
			-doMaf 1 \
			-doCounts 1 \
			-setMaxDepth {params.max_dp} \
			-baq 2 \
			-ref {input.ref} \
			-minInd $MIN_IND \
			-minQ 20 \
			-minMapQ 30 \
			-sites {input.sites} \
			-r {wildcards.chrom} \
			-bam {input.bams} 2> {log}
		"""
 
rule concat_angsd_alldegenerates_gl_withoutRelated:
	"""
	Concatenated GLs from all 16 chromosomes into single file. Done separately for each site type.
	"""
	input:
		get_angsd_alldegenerates_gl_toConcat_withoutRelated
	output:
		'{0}/gls/allSamples/{{site}}/withoutRelated/allChroms_{{site}}.beagle.gz'.format(ANGSD_DIR)
	log: LOG_DIR + '/concat_angsd_alldegenerates_gl/withoutRelated/allSamples_{site}_concat.log'
	container: 'library://james-s-santangelo/angsd/angsd:0.933' 
	shell:
		"""
		first=1
		for f in {input}; do
			if [ "$first"  ]; then
				zcat "$f"
				first=
			else
				zcat "$f"| tail -n +2
			fi
		done | bgzip -c > {output} 2> {log}
		"""

rule concat_angsd_alldegenerates_mafs_withoutRelated:
	"""
	Concatenate MAF files for each of 16 chromosomes into single file. Done separately for each site type.
	"""
	input:
		get_angsd_alldegenerates_maf_toConcat_withoutRelated
	output:
		'{0}/gls/allSamples/{{site}}/withoutRelated/allChroms_{{site}}.mafs.gz'.format(ANGSD_DIR)
	log: LOG_DIR + '/concat_angsd_alldegenerates_mafs/withoutRelated/allSamples_{site}_concat.log'
	container: 'library://james-s-santangelo/angsd/angsd:0.933' 
	shell:
		"""
		first=1
		for f in {input}; do
			if [ "$first"  ]; then
				zcat "$f"
				first=
			else
				zcat "$f"| tail -n +2
			fi
		done | bgzip -c > {output} 2> {log}
		"""

rule angsd_allSamples_alldegenerates_withoutRelated_done:
    """
    Generate empty flag file signalling successful completion of GL estimation across all samples. 
    """
    input:
        expand(rules.concat_angsd_alldegenerates_gl_withoutRelated.output, site=['4fold']),
        expand(rules.concat_angsd_alldegenerates_mafs_withoutRelated.output, site=['4fold'])
    output:
        '{0}/angsd_allSamples_alldegenerates_withoutRelated.done'.format(ANGSD_DIR)
    shell:
        "touch {output}"

