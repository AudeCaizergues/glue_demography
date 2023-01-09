 # Population structure analyses. 

###############
#### SETUP ####
###############

rule pruneGLs_degenerateSNPs_withoutRelated:
    input:
        gls = rules.angsd_gl_allSamples_alldegenerates_withoutRelated.output.gls,
        pos = rules.prune_degenerateSNPs_forPopStructure.output
    output:
        '{0}/gls/allSamples_alldegenerates/{{site}}/withoutRelated/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.beagle.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/pruneGLs_degenerateSNPs/withoutRelated/{chrom}_{site}_maf{maf}_pruneGLs.log'
    params:
        out = '{0}/gls/allSamples_alldegenerates/{{site}}/withoutRelated/{{chrom}}/{{chrom}}_{{site}}_maf{{maf}}_pruned.beagle'.format(ANGSD_DIR)
    shell:
        """
        ( sed 's/:/_/g' {input.pos} > {input.pos}.temp
          zgrep 'marker' {input.gls} > {params.out} &&
                zgrep -w -f {input.pos}.temp {input.gls} >> {params.out} &&
                gzip {params.out}
                rm {input.pos}.temp ) 2> {log}
        """
        
rule concat_angsd_gl_pruned_withoutRelated:
    """
    Concatenated GLs from all 16 chromosomes into single file. Done separately for each site type.
    """
    input:
    	lambda wildcards: expand(rules.pruneGLs_degenerateSNPs_withoutRelated.output, chrom=CHROMOSOMES, site=wildcards.site, maf=wildcards.maf)
    output:
        '{0}/gls/allSamples/{{site}}/withoutRelated/allChroms_{{site}}_maf{{maf}}_pruned_withoutRelated.beagle.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/concat_angsd_gl_pruned/withoutRelated/allSamples_{site}_{maf}_concat_withoutRelated.log'
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
      

#by city
rule pruneGLs_byCity_degenerateSNPs_withoutRelated:
    input:
        gls = rules.angsd_gl_byCity_withoutRelated.output.gls,
        pos = rules.prune_degenerateSNPs_forPopStructure.output
    output:
        '{0}/gls/by_city/{{site}}/withoutRelated/{{city}}/{{chrom}}/{{city}}_{{chrom}}_{{site}}_maf{{maf}}_pruned_withoutRelated.beagle.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/pruneGLs_degenerateSNPs_by_city/withoutRelated/{city}/{chrom}/{city}_{chrom}_{site}_maf{maf}_pruneGLs_withoutRelated.log'
    params:
        out = '{0}/gls/by_city/{{site}}/withoutRelated/{{city}}/{{chrom}}/{{city}}_{{chrom}}_{{site}}_maf{{maf}}_pruned_withoutRelated.beagle'.format(ANGSD_DIR)
    priority: 50
    shell:
        """
        ( sed 's/:/_/g' {input.pos} > {input.pos}.temp
          zgrep 'marker' {input.gls} > {params.out} &&
                zgrep -w -f {input.pos}.temp {input.gls} >> {params.out} &&
                gzip {params.out} ) 2> {log}
        """
       
rule concat_angsd_byCity_gl_pruned_withoutRelated:
    """
    Concatenated GLs from all 16 chromosomes into single file. Done separately for each site type.
    """
    input:
    	lambda w: expand(rules.pruneGLs_byCity_degenerateSNPs_withoutRelated.output, chrom=CHROMOSOMES, city=w.city, site=w.site, maf=w.maf)
    output:
        '{0}/gls/byCity/{{site}}/withoutRelated/allChrom_{{city}}_{{site}}_maf{{maf}}_pruned_withoutRelated.beagle.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/concat_angsd_gl_pruned/withoutRelated/by_city/{city}_{site}_{maf}_concat_withoutRelated.log'
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
## gls devrait avoir les noms d'indiv en premiere ligne, mais c'est du binaire, regarder sur les fichiers chrom si ça vient de la à l'origine
 
#######################
### PCA & ADMIXTURE ###
#######################

# PCA & admixture analysis using LD-pruned 4fold SNPs from above (MAF > 0.05)

####################
#### GLOBAL PCA ####
####################

rule pcangsd_allSamples_withoutRelated:
    """
    Perform PCA using genome-wide 4fold dengenerate sites using all samples from all cities.
    """
    input:
        rules.concat_angsd_gl_pruned_withoutRelated.output
    output:
        '{0}/pcangsd/allSamples/pcangsd_withoutRelated/allSamples_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
    log: LOG_DIR + '/pcangsd_allSamples/pcangsd_withoutRelated/allSamples_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 10
    params:
        out = '{0}/pcangsd/allSamples/pcangsd_withoutRelated/allSamples_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 12000,
        time = '03:00:00'
    shell:
        """
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            -iter 5000 \
            > {log}
        """

#####################
#### PCA BY CITY ####
#####################


rule pcangsd_byCity_withoutRelated:
    """
    Perform PCA by city using genome-wide 4fold dengenerate sites and estimate admixture proportions
    """
    input:
        rules.concat_angsd_byCity_gl_pruned_withoutRelated.output
    output:
        cov = '{0}/pcangsd/by_city/pcangsd_withoutRelated/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.cov'.format(POP_STRUC_DIR),
        Q = '{0}/pcangsd/by_city/pcangsd_withoutRelated/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd.admix.Q.npy'.format(POP_STRUC_DIR)
    log: LOG_DIR + '/pcangsd_byCity/pcangsd_withoutRelated/{city}_{site}_maf{maf}_pcangsd.log'
    container: 'library://james-s-santangelo/pcangsd/pcangsd:0.99'
    threads: 6
    params:
        out = '{0}/pcangsd/by_city/pcangsd_withoutRelated/{{city}}/{{city}}_{{site}}_maf{{maf}}_pcangsd'.format(POP_STRUC_DIR)
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 8000,
        time = '02:00:00'
    shell:
        """
        python3 /opt/pcangsd-v.0.99/pcangsd.py \
            -beagle {input} \
            -o {params.out} \
            -minMaf {wildcards.maf} \
            -threads {threads} \
            -admix \
            -admix_seed 42 \
            -iter 5000 \
            > {log}
        """
###########################
#### ADMIXTURE BY CITY ####
###########################

rule ngsadmix:
    input:
        rules.concat_angsd_byCity_gl_pruned_withoutRelated.output
    output:
        fopt = '{0}/ngsadmix/{{city}}/K{{k}}/ngsadmix_{{city}}_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}.fopt.gz'.format(POP_STRUC_DIR),
        qopt = '{0}/ngsadmix/{{city}}/K{{k}}/ngsadmix_{{city}}_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}.qopt'.format(POP_STRUC_DIR),
        lf = '{0}/ngsadmix/{{city}}/K{{k}}/ngsadmix_{{city}}_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}.log'.format(POP_STRUC_DIR)
    log: LOG_DIR + '/ngsadmix/{city}_{site}_maf{maf}_K{k}_seed{seed}_ngsadmix.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    threads: 10
    params:
        out = '{0}/ngsadmix/{{city}}/K{{k}}/ngsadmix_{{city}}_{{site}}_maf{{maf}}_K{{k}}_seed{{seed}}'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '02:00:00'
    shell:
        """
        NGSadmix -likes {input} \
            -K {wildcards.k} \
            -seed {wildcards.seed} \
            -P {threads} \
            -outfiles {params.out} 2> {log}
        """

rule logfile_for_clumpak:
    """
    Create Inputfile for CLUMPAK containing Log likelihood values of NGSadmix runs for each K
    """
    input:
        lambda w: expand(rules.ngsadmix.output.lf, city=w.city, site='4fold', maf='0.05', k=NGSADMIX_K, seed=NGSADMIX_SEEDS)
    output:
        '{0}/clumpak/ngsadmix_logfile_for_clumpak_{{city}}.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        import re
        with open(output[0], 'w') as fout:
            for lf in input:
                # Get K
                m1 = re.search('(?<=_K)(\d+)', lf)
                k = m1.group(1)
                # Get likelihood
                line = open(lf, 'r').readlines()[-1]  # Likelihood always on last line
                m2 = re.search('(?<=like=)(-?\d+.\d+)', line)
                like = m2.group(1)
                fout.write('{0}\t{1}\n'.format(k, like))

rule clumpak_best_k_by_evanno:
    """
    Find optimal K value by city using Evanno method, as implemented in CLUMPAK
    """
    input:
        ancient(rules.logfile_for_clumpak.output)
    output:
        directory('{0}/bestKbyEvanno/{{city}}/'.format(POP_STRUC_DIR))
    log: LOG_DIR + '/clumpak_best_k_by_evanno/evanno_{city}.log'
    container: 'library://james-s-santangelo/clumpak/clumpak:1.1'
    params:
        outdir = '{0}/bestKbyEvanno/{{city}}/'.format(POP_STRUC_DIR)
    resources:
        mem_mb = 1000,
        time = '01:00:00'
    shell:
        """
        perl /opt/bin/BestKByEvanno.pl --id clumpak_best_k_out \
            --d {params.outdir} \
            --f {input} \
            --inputtype lnprobbyk 2>&1 > {log}
        """
        

##########################
#### ADMIXTURE GLOBAL ####
##########################

rule ngsadmix_global:
    input:
        rules.concat_angsd_gl_pruned_withoutRelated.output
    output:
        fopt = '{0}/ngsadmix/global/K{{k2}}/ngsadmix_{{site}}_maf{{maf}}_K{{k2}}_seed{{seed}}.fopt.gz'.format(POP_STRUC_DIR),
        qopt = '{0}/ngsadmix/global/K{{k2}}/ngsadmix_{{site}}_maf{{maf}}_K{{k2}}_seed{{seed}}.qopt'.format(POP_STRUC_DIR),
        lf = '{0}/ngsadmix/global/K{{k2}}/ngsadmix_{{site}}_maf{{maf}}_K{{k2}}_seed{{seed}}.log'.format(POP_STRUC_DIR)
    log: LOG_DIR + '/ngsadmix/global_{site}_maf{maf}_K{k2}_seed{seed}_ngsadmix.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    threads: 10
    params:
        out = '{0}/ngsadmix/global/K{{k2}}/ngsadmix_{{site}}_maf{{maf}}_K{{k2}}_seed{{seed}}'.format(POP_STRUC_DIR)
    wildcard_constraints:
        site = '4fold'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '02:00:00'
    shell:
        """
        NGSadmix -likes {input} \
            -K {wildcards.k2} \
            -seed {wildcards.seed} \
            -P {threads} \
            -outfiles {params.out} 2> {log}
        """

rule logfile_for_clumpak_global:
    """
    Create Inputfile for CLUMPAK containing Log likelihood values of NGSadmix runs for each K
    """
    input:
        lambda w: expand(rules.ngsadmix_global.output.lf, site='4fold', maf='0.05', k2=NGSADMIX_K2, seed=NGSADMIX_SEEDS)
    output:
        '{0}/clumpak/ngsadmix_logfile_for_clumpak_global.txt'.format(PROGRAM_RESOURCE_DIR)
    run:
        import re
        with open(output[0], 'w') as fout:
            for lf in input:
                # Get K
                m1 = re.search('(?<=_K)(\d+)', lf)
                k = m1.group(1)
                # Get likelihood
                line = open(lf, 'r').readlines()[-1]  # Likelihood always on last line
                m2 = re.search('(?<=like=)(-?\d+.\d+)', line)
                like = m2.group(1)
                fout.write('{0}\t{1}\n'.format(k, like))

rule clumpak_best_k_by_evanno_global:
    """
    Find optimal K value by city using Evanno method, as implemented in CLUMPAK
    """
    input:
        rules.logfile_for_clumpak_global.output
    output:
        directory('{0}/bestKbyEvanno/global/'.format(POP_STRUC_DIR))
    log: LOG_DIR + '/clumpak_best_k_by_evanno/evanno_global.log'
    container: 'library://james-s-santangelo/clumpak/clumpak:1.1'
    params:
        outdir = '{0}/bestKbyEvanno/global/'.format(POP_STRUC_DIR)
    resources:
        mem_mb = 1000,
        time = '01:00:00'
    shell:
        """
        perl /opt/bin/BestKByEvanno.pl --id clumpak_best_k_out \
            --d {params.outdir} \
            --f {input} \
            --inputtype lnprobbyk 2>&1 > {log}
        """


##############
#### POST ####
##############

rule population_structure_withoutRelated:
    """
    Generate empty flag file signaling successful completion of PCAngsd
    """
    input:
        expand(rules.pruneGLs_byCity_degenerateSNPs_withoutRelated.output, site = '4fold', maf = '0.05', city = CITIES, chrom=CHROMOSOMES),
        expand(rules.pcangsd_allSamples_withoutRelated.output, site = '4fold', maf = ['0.05']),
        expand(rules.pcangsd_byCity_withoutRelated.output, site = '4fold', maf = '0.05', city = CITIES),
        expand(rules.ngsadmix.output, site = '4fold', maf = '0.05', city = CITIES, seed=NGSADMIX_SEEDS, k=NGSADMIX_K),
        expand(rules.clumpak_best_k_by_evanno.output, city = CITIES),
        expand(rules.ngsadmix_global.output, site = '4fold', maf = '0.05', seed=NGSADMIX_SEEDS, k2=NGSADMIX_K2),
        expand(rules.clumpak_best_k_by_evanno_global.output)
    output:
        '{0}/population_structure_withoutRelated.done'.format(POP_STRUC_DIR)
    shell:
        """
        touch {output}
        """
