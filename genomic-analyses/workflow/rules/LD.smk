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
        n_ind = 2058
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
        ( sed 's/:/_/g' {input.pos} > {input.pos}.temp
          zgrep 'marker' {input.gls} > {params.out} &&
                zgrep -w -f {input.pos}.temp {input.gls} >> {params.out} &&
                gzip {params.out}
                rm {input.pos}.temp ) 2> {log}
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


##############
#### POST ####
##############

rule LD_done:
    """
    Generate empty flag file signaling successful completion of LD and relatedness analyses
    """
    input:
        expand(rules.ngsLD_degenerateSites.output, site = '4fold', maf = ['0.05'], chrom=CHROMOSOMES),
        expand(rules.concat_angsd_gl_pruned.output, site = '4fold', maf = ['0.05'], chrom=CHROMOSOMES)
    output:
        '{0}/LD.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
