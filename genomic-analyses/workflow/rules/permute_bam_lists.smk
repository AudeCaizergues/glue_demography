# Rules to permute urban and rural samples by city to estimate null distribution of pi and Fst

###############
#### SETUP ####
###############

rule create_random_bam_list_byCity_byHabitat:
    """
    Create text files with paths to BAM files for each city. 'urban' and 'rural' samples are randomly
    selected from all possible bams for a city. Sample sizes are the same as those in the observed samples.
    """
    input:
        unpack(get_urban_rural_bam_lists)
    output:
        urban = '{0}/bam_lists/by_city/{{city}}/randomized/{{city}}_randU_seed{{seed}}_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR),
        rural = '{0}/bam_lists/by_city/{{city}}/randomized/{{city}}_randR_seed{{seed}}_{{site}}_bams.list'.format(PROGRAM_RESOURCE_DIR)
    priority: 1
    run:
        import random
        urban_bams = list(open(input.urban_bams, 'r'))
        rural_bams = list(open(input.rural_bams, 'r'))
        urban_n = len(urban_bams)
        rural_n = len(rural_bams)
        all_bams = urban_bams + rural_bams
        random.seed(int(wildcards.seed))
        randU = random.sample(all_bams, urban_n)
        randR = [bam for bam in all_bams if not bam in randU] 
        with open(output.urban, 'w') as uout:
            for bam in randU:
                uout.write(bam)
        with open(output.rural, 'w') as rout:
            for bam in randR:
                rout.write(bam)
                
                
                
##############
#### POST ####
##############

rule permuted_bam_list_done:
    """
    Generate empty flag file signalling successful completion of pairwise pi and Fst analysis
    """
    input:
        expand(rules.create_random_bam_list_byCity_byHabitat.output, city=CITIES, habitat=HABITATS, site=['4fold'], seed=BOOT_SEEDS)
    output:
        '{0}/permuted_bam_list.done'.format(ANGSD_DIR)
    shell:
        """
        touch {output}
        """
