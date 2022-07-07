
#####################
#### RELATEDNESS ####
#####################

# Estimate pariwise relatedness among all individuals within a city

rule pruned_degenerate_angsd_format:
    """
    Create ANGSD sites-formatted file with position of LD-pruned 4fold sites
    """
    input:
        rules.prune_degenerateSNPs_forPopStructure.output
    output:
        '{0}/angsd_sites/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.sites'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/pruned_degenerate_angsd_format/{chrom}_{site}_maf{maf}.log'
    shell:
        """
        sed 's/:/\t/g' {input} > {output} 2> {log}
        """

rule concat_pruned_degenerate_angsd_format:
	"""
	Concatenated ANGSD sites-formatted file with position of LD-pruned 4fold sites from all 16 chromosomes 
    into single file. Done separately for each site type.
	"""
	input:
		get_pruned_angsd_files_toConcat
	output:
		'{0}/angsd_sites/{{site}}/{{site}}_maf{{maf}}_all_chrom_pruned.sites'.format(PROGRAM_RESOURCE_DIR)
	log: LOG_DIR + '/concat_pruned_degenerate_angsd_format/allchrom_pruned_{site}_maf{maf}_concat.log'
	shell:
		"""
		first=1
		for f in {input}; do
				cat "$f"
			fi
		done | bgzip -c > {output} 2> {log}
		"""

rule angsd_index_prunedSNPs:
    input:
        rules.concat_pruned_degenerate_angsd_format.output
    output:
        binary = '{0}/angsd_sites/{{site}}/{{site}}_maf{{maf}}_all_chrom_pruned.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{site}}/{{site}}_maf{{maf}}_all_chrom_pruned.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    log: LOG_DIR + '/angsd_index_prunedSNPs/{site}_maf{maf}_all_chom_prunedIndex.log'
    shell:
        """
        angsd sites index {input} 2> {log}
        """
        
rule angsd_gl_forNGSrelate:
    input:
        bams = rules.remove_lowCovSamples_forPCA_byCity.output,
        ref = rules.unzip_reference.output,
        sites = rules.concat_pruned_degenerate_angsd_format.output,
        idx = rules.angsd_index_prunedSNPs.output
    output:
        gls = '{0}/gls/ngsrelate/{{city}}/{{site}}_maf{{maf}}_forNGSrelate.glf.gz'.format(ANGSD_DIR),
        mafs = '{0}/gls/ngsrelate/{{city}}/{{site}}_maf{{maf}}_forNGSrelate.mafs.gz'.format(ANGSD_DIR),
        pos = temp('{0}/gls/ngsrelate/{{city}}/{{site}}_maf{{maf}}_forNGSrelate.glf.pos.gz'.format(ANGSD_DIR))
    log: LOG_DIR + '/angsd_gl_forNGSrelate/{city}/{city}_{site}_maf{maf}.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/gls/ngsrelate/{{city}}/{{site}}_maf{{maf}}_forNGSrelate'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '3:00:00'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*50/100 ))
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 3 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -setMaxDepth {params.max_dp} \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -sites {input.sites} \
            -minQ 20 \
            -minMapQ 30 \
            -minMaf {wildcards.maf} \
            -bam {input.bams} 2> {log}
        """

rule convert_freq_forNGSrelate:
    input:
        rules.angsd_gl_forNGSrelate.output.mafs
    output:
        '{0}/gls/ngsrelate/{{city}}/{{site}}_maf{{maf}}_forNGSrelate.freqs'.format(ANGSD_DIR)
    log: LOG_DIR + '/convert_freq_forNGSrelate/{city}/{city}_{site}_maf{maf}_convert_freqs.log'
    shell:
        """
        zcat {input} | cut -f6 | sed 1d > {output} 2> {log}
        """

rule ngsrelate:
    input:
        bams = rules.remove_lowCovSamples_forPCA_byCity.output,
        gls = rules.angsd_gl_forNGSrelate.output.gls,
        freq = rules.convert_freq_forNGSrelate.output
    output:
        '{0}/ngsrelate/{{city}}/{{city}}_{{site}}_maf{{maf}}_NGSrelate.out'.format(POP_STRUC_DIR)
    log: LOG_DIR + '/ngsrelate/{city}/{city}_{site}_maf{maf}.log'
    container: 'library://james-s-santangelo/ngsrelate/ngsrelate:2.0' 
    threads: 10
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '02:00:00'
    shell:
        """
        N=$( wc -l < {input.bams} );
        ngsRelate -f {input.freq} \
            -O {output} \
            -g {input.gls} \
            -p {threads} \
            -z {input.bams}\
            -n $N 2> {log}
        """
        
        
##############
#### POST ####
##############

rule Relatedness_done:
    """
    Generate empty flag file signaling successful completion of LD and relatedness analyses
    """
    input:
        expand(rules.angsd_gl_forNGSrelate.output, site = '4fold', maf = ['0.05'], city = CITIES),
        expand(rules.ngsrelate.output, site = '4fold', maf = ['0.05'], city = CITIES)
    output:
        '{0}/Relatedness.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """