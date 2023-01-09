# Rules for estimating SFS (1D and 2D), summary stats, and GLs for urban and rural habitats within cities

###############################
#### SFS AND SUMMARY STATS ####
###############################

rule create_bam_list_byCity_byPop:
    """
    Create text file with paths to BAM files in each habitat by city. 
    """
    input:
        rules.create_bam_list_finalSamples.output
    output:
        '{0}/bam_lists/by_city_byPop/{{city_pop}}._{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    run:
        import os
        import pandas as pd
        df = pd.read_table(config['samples'], sep = '\t')
        city_name = wildcards.city_pop.split('_.')[0]
        pop_number = wildcards.city_pop.split('_.')[1]
        df_sub = df[(df['city'] == city_name)  & (df['pop'] == int(pop_number))]
        samples_city_habitat = df_sub['sample'].tolist()
        bams = open(input[0], 'r').readlines()
        with open(output[0], 'w') as f:
            for bam in bams:
                search = re.search('^(.+)(?=_\w)', os.path.basename(bam))
                sample = search.group(1)
                if sample in samples_city_habitat:
                    f.write('{0}'.format(bam))

rule angsd_saf_likelihood_byCity_byPop:
    """
    Generate Site Allele Frequency (SAF) likelihood file for each habitat in each city using ANGSD. 
    Uses only 4fold sites.
    """
    input:
        unpack(get_files_for_saf_estimation_byCity_byPop)
    output:
        saf = '{0}/sfs/by_city/byPop/{{city_pop}}._{{site}}.saf.gz'.format(ANGSD_DIR),
        saf_idx = '{0}/sfs/by_city/byPop/{{city_pop}}._{{site}}.saf.idx'.format(ANGSD_DIR),
        saf_pos = '{0}/sfs/by_city/byPop/{{city_pop}}._{{site}}.saf.pos.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_saf_likelihood_byCity_byHabitat_byPop/{city_pop}_{site}_saf.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    params:
        out = '{0}/sfs/by_city/byPop/{{city_pop}}._{{site}}'.format(ANGSD_DIR)
    threads: 4
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
        
rule angsd_estimate_joint_sfs_byCity_byPop:
    """
    Estimated folded, two-dimensional pop1-pop2 SFS for each city using realSFS. Uses 4fold sites.
    """
    input:
        unpack(get_city_population_saf_files)
    output:
        '{0}/sfs/by_city/byPop/{{pop_combo}}_{{site}}.2dsfs'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_estimate_2dsfs_byCity_byPop/{pop_combo}_{site}.2dsfs.log'
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


##############
#### FST  ####
##############

rule angsd_fst_index_byCity_byPop:
    """
    Estimate per-site alphas (numerator) and betas (denominator) for Fst estimation.
    """
    input: 
        unpack(get_city_population_saf_sfs_files)
    output:
        fst = '{0}/summary_stats/byPop/hudson_fst/{{pop_combo}}_{{site}}.fst.gz'.format(ANGSD_DIR),
        idx = '{0}/summary_stats/byPop/hudson_fst/{{pop_combo}}_{{site}}.fst.idx'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_fst_index/{pop_combo}_{site}_index.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    threads: 4
    resources:
        mem_mb = 4000,
        time = '01:00:00'
    params:
        fstout = '{0}/summary_stats/byPop/hudson_fst/{{pop_combo}}_{{site}}'.format(ANGSD_DIR)
    shell:
        """
        realSFS fst index {input.saf} \
            -sites {input.sites} \
            -sfs {input.sfs} \
            -fold 1 \
            -P {threads} \
            -whichFst 1 \
            -fstout {params.fstout} 2> {log}
        echo {input.saf_gz1}
        echo {input.saf_gz2}
        """

rule angsd_fst_readable_byCity_byPop:
    """
    Create readable Fst files. Required due to format of realSFS fst index output files. 
    """
    input:
        rules.angsd_fst_index_byCity_byPop.output.idx
    output:
        '{0}/summary_stats/byPop/hudson_fst/{{pop_combo}}_{{site}}_readable.fst'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_fst_readable/byPop/{{pop_combo}}_{{site}}_readable_fst.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    shell:
        """
        realSFS fst print {input} > {output} 2> {log}
        """


##############
#### POST ####
##############

rule angsd_byCity_byPop_done:
    """
    Generate empty flag file signalling successful completion of SFS, summary stat and GL estimation 
    for populations within cities
    """
    input:
        expand(rules.angsd_estimate_joint_sfs_byCity_byPop.output, city_pop=CITY_POP, site=['4fold'], pop_combo=POP_COMBO),
        expand(rules.angsd_fst_index_byCity_byPop.output, site=['4fold'], pop_combo=POP_COMBO)
    output:
        '{0}/angsd_byCity_byPop.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """