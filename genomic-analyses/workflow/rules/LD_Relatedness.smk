###############
#### SETUP ####
###############

# Estimate LD among 4fold SNPs with MAF > 0.05. 
# Prune these SNPs within 20 Kb using r-squared cutoff of 0.2

rule create_pos_file_for_ngsLD:
    input:
        rules.angsd_gl_allSamples_alldegenerates.output.mafs
    output:
        '{0}/ngsld_pos/{{chrom}}_{{site}}_maf{{maf}}.pos'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/create_pos_file_for_ngsLD/{chrom}_{site}_maf{maf}_pos.log'
    shell:
        """
        zcat {input} | cut -f 1,2 | tail -n +2 > {output} 2> {log}
        """

rule ngsLD_degenerateSites:
    input:
        pos = rules.create_pos_file_for_ngsLD.output,
        gls = rules.angsd_gl_allSamples_alldegenerates.output.gls
    output:
        '{0}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}.ld.gz'.format(NGSLD_DIR)
    log: LOG_DIR + '/ngsld/{chrom}_{site}_maf{maf}_calc_ld.log'
    container: 'library://james-s-santangelo/ngsld/ngsld:1.1.1'
    threads: 8
    params:
        n_ind = len(FINAL_SAMPLES)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '1:00:00'
    shell:
        """
        ( NUM_SITES=$(cat {input.pos} | wc -l) &&
          ngsLD --geno {input.gls} \
            --pos {input.pos} \
            --n_ind {params.n_ind} \
            --n_sites $NUM_SITES \
            --probs \
            --n_threads {threads} \
            --max_kb_dist 20 | gzip --best > {output} ) 2> {log}
        """

rule prune_degenerateSNPs_forPopStructure:
    input:
        rules.ngsLD_degenerateSites.output
    output:
        '{0}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.id'.format(NGSLD_DIR)
    log: LOG_DIR + '/prune_degenerateSNP_forPopStructure/{chrom}_{site}_maf{maf}_prune_ld.log'
    container: 'library://james-s-santangelo/ngsld/ngsld:1.1.1'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '3:00:00'
    shell:
        """
        ( zcat {input} | perl /opt/bin/prune_graph.pl \
            --max_kb_dist 20 \
            --min_weight 0.2 | sort -V > {output} ) 2> {log}
        """

rule pruneGLs_degenerateSNPs:
    input:
        gls = rules.angsd_gl_allSamples_alldegenerates.output.gls,
        pos = rules.prune_degenerateSNPs_forPopStructure.output
    output:
        '{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.beagle.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/pruneGLs_degenerateSNPs/{chrom}_{site}_maf{maf}_pruneGLs.log'
    params:
        out = '{0}/gls/allSamples/{{site}}/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.beagle'.format(ANGSD_DIR)
    shell:
        """
        ( zgrep 'marker' {input.gls} > {params.out} &&
                sed 's/:/_/g' {input.pos} | zgrep -w -f - {input.gls} >> {params.out} &&
                gzip {params.out} ) 2> {log}
        """

rule concat_angsd_gl_pruned:
    """
    Concatenated GLs from all 16 chromosomes into single file. Done separately for each site type.
    """
    input:
    	lambda wildcards: expand(rules.pruneGLs_degenerateSNPs.output, chrom=CHROMOSOMES, site=wildcards.site, maf=wildcards.maf)
    output:
        '{0}/gls/allSamples/{{site}}/allChroms_{{site}}_maf{{maf}}_pruned.beagle.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/concat_angsd_gl_pruned/allSamples_{site}_{maf}_concat.log'
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

rule angsd_index_prunedSNPs:
    input:
        rules.pruned_degenerate_angsd_format.output
    output:
        binary = '{0}/angsd_sites/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.sites.bin'.format(PROGRAM_RESOURCE_DIR),
        idx = '{0}/angsd_sites/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.sites.idx'.format(PROGRAM_RESOURCE_DIR)
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    log: LOG_DIR + '/angsd_index_prunedSNPs/{chrom}_{site}_maf{maf}_prunedIndex.log'
    shell:
        """
        angsd sites index {input} 2> {log}
        """
     
rule angsd_gl_forNGSrelate:
    input:
        bams = rules.angsd_gl_byCity.output,
        ref = rules.unzip_reference.output,
        sites = rules.pruned_degenerate_angsd_format.output,
        idx = rules.angsd_index_prunedSNPs.output
    output:
        gls = temp('{0}/gls/ngsrelate/{{city}}/{{city}}_{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate.glf.gz'.format(ANGSD_DIR)),
        mafs = temp('{0}/gls/ngsrelate/{{city}}/{{city}}_{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate.mafs.gz'.format(ANGSD_DIR)),
        pos = temp('{0}/gls/ngsrelate/{{city}}/{{city}}_{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate.glf.pos.gz'.format(ANGSD_DIR))
    log: LOG_DIR + '/angsd_gl_forNGSrelate/{{city}}/{{city}}_{chrom}_{site}_maf{maf}.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/gls/ngsrelate/{{city}}/{{city}}_{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate'.format(ANGSD_DIR),
        max_dp = ANGSD_MAX_DP,
        min_dp_ind = ANGSD_MIN_DP_IND_GL
    threads: 8
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '3:00:00'
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND*80/100 ))
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 3 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -setMinDepthInd {params.min_dp_ind} \
            -setMaxDepth {params.max_dp} \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -sites {input.sites} \
            -minQ 20 \
            -minMapQ 30 \
            -minMaf {wildcards.maf} \
            -r {wildcards.chrom} \
            -bam {input.bams} 2> {log}
        """

rule convert_freq_forNGSrelate:
    input:
        rules.angsd_gl_forNGSrelate.output.mafs
    output:
        '{0}/gls/ngsrelate/{{city}}/{{city}}_{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_forNGSrelate.freqs'.format(ANGSD_DIR)
    log: LOG_DIR + '/convert_freq_forNGSrelate/{{city}}/{{city}}_{chrom}_{site}_maf{maf}_convert_freqs.log'
    shell:
        """
        zcat {input} | cut -f6 | sed 1d > {output} 2> {log}
        """

rule ngsrelate:
    input:
        bams = rules.angsd_gl_byCity.output,
        gls = rules.angsd_gl_forNGSrelate.output.gls,
        freq = rules.convert_freq_forNGSrelate.output
    output:
        '{0}/ngsrelate/{{city}}/{{city}}_{{chrom}}_{{site}}_maf{{maf}}_NGSrelate.out'.format(POP_STRUC_DIR)
    log: LOG_DIR + '/ngsrelate/{{city}}/{{city}}_{chrom}_{site}_maf{maf}.log'
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
            -n $N 2> {log}
        """
        
        
##############
#### POST ####
##############

rule LD_Relatedness_done:
    """
    Generate empty flag file signaling successful completion of PCAngsd
    """
    input:
        expand(rules.ngsrelate.output, site = '4fold', maf = ['0.05'], chrom=CHROMOSOMES, city = CITIES)
    output:
        '{0}/LD_Relatedness.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """
