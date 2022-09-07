# Rules to get genotype likelihoods across all individuals per city. 
# Used for population structure analyses. 

###############
#### SETUP ####
###############

# REMOVE IT WHEN RUNNIN angsd_byCity_byHabitat to avoid doublons
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
#

rule concat_habitat_bamLists_withinCities_noRelated:
    """
    Concatenate urban and rural sample BAM lists within cities. Generates a single file with
    the paths to all of the BAM files for samples within a city
    """
    input:
        get_bamLists_toConcat_withoutRelated
    output:
        '{0}/bam_lists/by_city/withoutRelated/{{city}}/{{city}}_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/concat_habitat_bamLists_withinCities/withoutRelated/{city}_{site}_concat.log'
    shell:
        """
        cat {input} > {output} 2> {log}
        """

rule remove_lowCovSamples_forPCA_byCity_withoutRelated:
    """
    Removes samples with mean coverage lower than `params.cov` from within-city BAM lists.
    """
    input:
        bams = rules.concat_habitat_bamLists_withinCities_noRelated.output,
        qc_data = config['qualimap_bamqc']
    output:
        '{0}/bam_lists/by_city/withoutRelated/{{city}}/{{city}}_{{site}}_lowCovRemoved_bams.list'.format(PROGRAM_RESOURCE_DIR)
    log: LOG_DIR + '/remove_lowCovSamples_forPCA_byCity/withoutRelated/{city}_{site}_remove_lowCovSamples.log'
    params:
        cov = 0.2
    run:
        import logging
        import re
        import pandas as pd
        logging.basicConfig(filename=log[0], level=logging.DEBUG)
        try:
            qc_data = pd.read_table(input.qc_data, sep = '\t')
            qc_data['sample'] = qc_data['Sample'].str.extract('(\w+_\d+_\d+)')
            cols = ['sample', 'mean_coverage']
            qc_data = qc_data[cols]
            samples_lowCovRemoved = qc_data[qc_data['mean_coverage'] >= float(params.cov)]['sample'].tolist()
            bams = open(input.bams[0], 'r').readlines()
            with open(output[0], 'w') as fout:
                for bam in bams:
                    match = re.search('^(.+)(?=_\w)', os.path.basename(bam))
                    sample = match.group(1)
                    if sample in samples_lowCovRemoved:
                        fout.write(bam)
        except:
            logging.exception("An error occured!") 
            raise

##############################
#### GENOTYPE LIKELIHOODS ####
##############################

rule angsd_gl_byCity_withoutRelated:
    """
    Estimate Beagle genotype likelihoods jointly for all samples within city.
    """
    input:
        bams = rules.remove_lowCovSamples_forPCA_byCity_withoutRelated.output,
        sites = rules.convert_sites_for_angsd.output,
        sites_idx = rules.angsd_index_degenerate_sites.output, 
        ref = rules.unzip_reference.output,
        chroms = config['chromosomes']
    output:
        gls = '{0}/gls/by_city/withoutRelated/{{city}}/{{city}}_{{site}}_maf{{maf}}.beagle.gz'.format(ANGSD_DIR),
        mafs = '{0}/gls/by_city/withoutRelated/{{city}}/{{city}}_{{site}}_maf{{maf}}.mafs.gz'.format(ANGSD_DIR)
    log: LOG_DIR + '/angsd_gl_byCity_beagle/withoutRelated/{city}_{site}_maf{maf}_beagleGL.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933' 
    params:
        out = '{0}/gls/by_city/withoutRelated/{{city}}/{{city}}_{{site}}_maf{{maf}}'.format(ANGSD_DIR)
    threads: 6
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 5000,
        time = '08:00:00' 
    shell:
        """
        NUM_IND=$( wc -l < {input.bams} );
        MIN_IND=$(( NUM_IND * 50/100 ));
        angsd -GL 1 \
            -out {params.out} \
            -nThreads {threads} \
            -doGlf 2 \
            -doMajorMinor 1 \
            -SNP_pval 1e-6 \
            -doMaf 1 \
            -doCounts 1 \
            -baq 2 \
            -ref {input.ref} \
            -minInd $MIN_IND \
            -minQ 20 \
            -minMapQ 30 \
            -minMaf {wildcards.maf} \
            -sites {input.sites} \
            -rf {input.chroms} \
            -bam {input.bams} 2> {log}
        """

##############
#### POST ####
##############

rule angsd_byCity_withoutRelated_done:
    input:
        expand(rules.angsd_gl_byCity_withoutRelated.output, city=CITIES, site='4fold', maf='0.05')
    output:
        '{0}/angsd_byCity_withoutRelated.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """