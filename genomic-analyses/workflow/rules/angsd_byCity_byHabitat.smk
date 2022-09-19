# Rules for estimating SFS (1D and 2D), summary stats, and GLs for urban and rural habitats within cities

###############################
#### SFS AND SUMMARY STATS ####
###############################

rule create_samples_to_remove:
    """
    Writes file with sample names for those with high alignment error rates.
    Thresholds were assessed through exploratory analysis of QC data. 
    """
    input:
        qc_data = ancient(config['qualimap_bamqc'])
    output: 
        error_df = '{0}/highErrorRate_toRemove.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        import pandas as pd
        qc_data = pd.read_table(input.qc_data, sep = '\t')
        qc_data['sample'] = qc_data['Sample'].str.extract('(\w+_\d+_\d+)')
        cols = ['sample', 'mean_coverage', 'general_error_rate']
        qc_data = qc_data[cols]
        # Samples with high ealignment errors have error rates > 0.04
        highError_samples = qc_data[qc_data['general_error_rate'] >= 0.03]
        highError_samples['sample'].to_csv(output.error_df, header = None, index = None)
	
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

rule create_bam_list_byCity_byHabitat_withoutRelated:
    """
    Create text file with paths to BAM files in each habitat by city. 
    """
    input:
        rules.create_bam_list_finalSamples_withoutRelated.output
    output:
        '{0}/bam_lists/by_city/withoutRelated/{{city}}/{{city}}_{{habitat}}_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/create_bam_list_byCity_byHabitat_withoutRelated/withoutRelated/{city}_{habitat}_{site}_concat.log' 
    run:
        import os
        import pandas as pd
        df = pd.read_table(config['samples_noRelated'], sep = '\t')
        df_sub = df[(df['city'] == wildcards.city) & (df['site'] == wildcards.habitat)]
        samples_city_habitat = df_sub['sample'].tolist()
        bams = open(input[0], 'r').readlines()
        with open(output[0], 'w') as f:
            for bam in bams:
                search = re.search('^(.+)(?=_\w)', os.path.basename(bam))
                sample = search.group(1)
                if sample in samples_city_habitat:
                    f.write('{0}'.format(bam))

                    
rule angsd_saf_likelihood_byCity_byHabitat:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each habitat in each city using ANGSD. 
    Uses only 4fold sites.
    """
    input:
        unpack(get_files_for_saf_estimation_byCity_byHabitat)
    output:
        saf = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.saf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_saf_likelihood_byCity_byHabitat/{city}_{habitat}_{site}_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}'.format(ANGSD_DIR)
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '03:00:00'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*60/100 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doMajorMinor 4 \
            -baq 2 \
            -ref {input.ref} \
            -sites {input.sites} \
            -minInd $MIN_IND \
            -minMapQ 30 \
            -minQ 20 \
            -doSaf 1 \
            -anc {input.ref} \
            -bam {input.bams} 2> {log}
        """

rule angsd_estimate_joint_sfs_byCity:
    """
    Estimated folded, two-dimensional urban-rural SFS for each city using realSFS. Uses 4fold sites.
    """
    input:
        saf = get_habitat_saf_files_byCity,
        sites = "/scratch/projects/trifolium/glue/demography/glue_demography/results/program_resources/angsd_sites/Trepens_4fold.sites",
        idx = "/scratch/projects/trifolium/glue/demography/glue_demography/results/program_resources/angsd_sites/Trepens_4fold.sites.idx",
    output:
        '{0}/sfs/by_city/{{city}}/{{city}}_{{site}}_r_u.2dsfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_2dsfs_byCity/{city}_{site}.2dsfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 10000,
        time = '03:00:00'
    shell:
        """
        realSFS {input.saf} \
            -sites {input.sites} \
            -maxIter 2000 \
            -seed 42 \
            -fold 1 \
            -P {threads} > {output} 2> {log}
        """

rule angsd_estimate_sfs_byCity_byHabitat:
    """
    Estimate folded SFS separately for each habitat in each city (i.e., 1D SFS) using realSFS. 
    """
    input:
        saf = rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx,
        sites = "/scratch/projects/trifolium/glue/demography/glue_demography/results/program_resources/angsd_sites/Trepens_4fold.sites",
        idx = "/scratch/projects/trifolium/glue/demography/glue_demography/results/program_resources/angsd_sites/Trepens_4fold.sites.idx",
    output:
        '{0}/sfs/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.sfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_sfs_byCity_byHabitat/{city}_{habitat}_{site}_sfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5000,
        time = '03:00:00'
    shell:
        """
        realSFS {input.saf} \
            -sites {input.sites} \
            -P {threads} \
            -fold 1 \
            -maxIter 2000 \
            -seed 42 > {output} 2> {log}
        """

#######################
#### FST AND THETA ####
#######################

rule angsd_fst_index:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation.
    """
    input: 
        saf_idx = get_habitat_saf_files_byCity,
        joint_sfs = rules.angsd_estimate_joint_sfs_byCity.output,
        sites = "/scratch/projects/trifolium/glue/demography/glue_demography/results/program_resources/angsd_sites/Trepens_4fold.sites"
    output:
        fst = '{0}/summary_stats/hudson_fst/{{city}}/{{city}}_{{site}}_r_u.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/hudson_fst/{{city}}/{{city}}_{{site}}_r_u.fst.idx'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_fst_index/{city}_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/summary_stats/hudson_fst/{{city}}/{{city}}_{{site}}_r_u'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf_idx} \
            -sites {input.sites} \
            -sfs {input.joint_sfs} \
            -fold 1 \
            -P {threads} \
            -whichFst 1 \
            -fstout {params.fstout} 2> {log}
        """

rule angsd_fst_readable:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_fst_index.output.idx
    output:
        '{0}/summary_stats/hudson_fst/{{city}}/{{city}}_{{site}}_r_u_readable.fst'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_fst_readable/{city}_{site}_readable_fst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """

rule angsd_estimate_thetas_byCity_byHabitat:
    """
    Generate per-site thetas in each habitat for each city from 1DSFS
    """
    input:
        saf_idx = rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx,
        sfs = rules.angsd_estimate_sfs_byCity_byHabitat.output,
        sites = "/scratch/projects/trifolium/glue/demography/glue_demography/results/program_resources/angsd_sites/Trepens_4fold.sites"
    output:
        idx = '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.thetas.idx'.format(ANGSD_DIR),
        thet = '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.thetas.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_thetas_byCity_byHabitat/{city}_{habitat}_{site}_thetas.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    params:
        out = '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}'.format(ANGSD_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        realSFS saf2theta {input.saf_idx} \
            -sites {input.sites} \
            -P {threads} \
            -fold 1 \
            -sfs {input.sfs} \
            -outname {params.out} 2> {log}
        """

rule angsd_diversity_neutrality_stats_byCity_byHabitat:
    """
    Estimate pi, Waterson's theta, Tajima's D, etc. in each habitat in each city.
    """
    input:
        rules.angsd_estimate_thetas_byCity_byHabitat.output.idx
    output:
       '{0}/summary_stats/thetas/by_city/{{city}}/{{city}}_{{habitat}}_{{site}}.thetas.idx.pestPG'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_diversity_neutrality_stats_byCity_byHabitat/{city}_{habitat}_{site}_diversity_neutrality.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    shell:
        """
        thetaStat do_stat {input} 2> {log}
        """

##############
#### POST ####
##############

rule angsd_byCity_byHabitat_done:
    """
    Generate empty flag file signalling successful completion of SFS, summary stat and GL estimation 
    for habitats within cities
    """
    input:
        expand(rules.angsd_saf_likelihood_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold']),
        expand(rules.angsd_estimate_joint_sfs_byCity.output, city=CITIES, site=['4fold']),
        expand(rules.angsd_estimate_sfs_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold']),
        expand(rules.angsd_fst_readable.output, city=CITIES, site=['4fold']),
        expand(rules.angsd_diversity_neutrality_stats_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold'])
    output:
        '{0}/angsd_byCity_byHabitat.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """

