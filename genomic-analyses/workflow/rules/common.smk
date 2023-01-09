# Python functions used throughout snakemake workflow

def get_all_bams(BAM_DIR):
    """
    Returns list with paths to GLUE bams 
    """
    bams = expand(BAM_DIR + '/{sample}_{site}.bam', sample=SAMPLES, site='4fold')
    return bams

def get_all_bams_withoutRelated(BAM_DIR):
    """
    Returns list with paths to GLUE bams 
    """
    bams = expand(BAM_DIR + '/{sample}_{site}.bam', sample=FINAL_SAMPLES, site='4fold')
    return bams

def get_bed(wildcards):
    """
    Get correct BED file for conversion to ANGSD sites format
    """
    bed = expand(rules.get_fourfold_zerofold.output, site=['4fold', '0fold'])
    bed_out = [x for x in bed if wildcards.site in os.path.basename(x)]
    return bed_out

def get_angsd_alldegenerates_gl_toConcat(wildcards):
    """
    Returns list with correct genotype likelihood files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    out = expand(rules.angsd_gl_allSamples_alldegenerates.output.gls, chrom=CHROMOSOMES, site=wildcards.site)
    return out

def get_angsd_alldegenerates_maf_toConcat(wildcards):
    """
    Returns list with correct minor allele frequency files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    out = expand(rules.angsd_gl_allSamples_alldegenerates.output.mafs, chrom=CHROMOSOMES, site=wildcards.site)
    return out
    
def get_angsd_alldegenerates_gl_toConcat_withoutRelated(wildcards):
    """
    Returns list with correct genotype likelihood files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    out = expand(rules.angsd_gl_allSamples_alldegenerates_withoutRelated.output.gls, chrom=CHROMOSOMES, site=wildcards.site)
    return out

def get_angsd_alldegenerates_maf_toConcat_withoutRelated(wildcards):
    """
    Returns list with correct minor allele frequency files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    out = expand(rules.angsd_gl_allSamples_alldegenerates_withoutRelated.output.mafs, chrom=CHROMOSOMES, site=wildcards.site)
    return out

def get_angsd_gl_toConcat(wildcards):
    """
    Returns list with correct genotype likelihood files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    out = expand(rules.angsd_gl_allSamples.output.gls, chrom=CHROMOSOMES, site=wildcards.site)
    return out

def get_pruned_angsd_files_toConcat(wildcards):
    """
    Returns list with correct ANGSD sites-formatted file with position of LD-pruned 4fold sites for concatenation, 
    depending on combination of "chrom", "site", and "maf" wildcard values
    """
    out = expand(rules.pruned_degenerate_angsd_format.output, chrom=CHROMOSOMES, site=wildcards.site, maf=wildcards.maf)
    return out

def get_angsd_alldegenerates_maf_toConcat(wildcards):
    """
    Returns list with correct minor allele frequency files for concatenation, depending on combination
    of "sample_set", "site", and "maf" wildcard values
    """
    out = expand(rules.angsd_gl_allSamples_alldegenerates.output.mafs, chrom=CHROMOSOMES, site=wildcards.site)
    return out

def get_files_for_saf_estimation_byCity_byHabitat(wildcards):
    """
    Get files to estimate SAF likelihhods for urban and rural habitats by city.
    """
    sites_idx = expand(rules.angsd_index_degenerate_sites.output.idx, site=wildcards.site)
    sites = expand(rules.convert_sites_for_angsd.output, site=wildcards.site)
    ref = rules.unzip_reference.output
    bams = expand(rules.create_bam_list_byCity_byHabitat_withoutRelated.output, city=wildcards.city, habitat=wildcards.habitat, site = wildcards.site, sample=FINAL_SAMPLES)
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_files_for_saf_estimation_byCity_byPop(wildcards):
    """
    Get files to estimate SAF likelihhods for urban and rural habitats by city.
    """
    sites_idx = expand(rules.angsd_index_degenerate_sites.output.idx, site=wildcards.site)
    sites = expand(rules.convert_sites_for_angsd.output, site=wildcards.site)
    ref = rules.unzip_reference.output
    all_bams = expand(rules.create_bam_list_byCity_byPop.output, city_pop=wildcards.city_pop, site = wildcards.site, sample=FINAL_SAMPLES)
    city_pop_name = wildcards.city_pop
    bams = [x for x in all_bams if '{0}'.format(city_pop_name) in os.path.basename(x)]
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }


def get_files_for_alleleFreq_estimation_byCity_byHabitat(wildcards):
    """
    Get files to estimate Allele Frequencies for urban and rural habitats by city.
    """
    sites_idx = expand(rules.angsd_index_city_snps.output.idx, site=wildcards.site, city=wildcards.city)
    sites = expand(rules.snps_forAlleleFreqs_byCity_byHabitat.output, site=wildcards.site, city=wildcards.city)
    ref = rules.unzip_reference.output
    bams = expand(rules.create_bam_list_byCity_byHabitat.output, city=wildcards.city, habitat=wildcards.habitat, site = wildcards.site, sample=FINAL_SAMPLES)
    chroms = config['chromosomes']
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_habitat_saf_files_byCity(wildcards):
    """
    Returns list with 4fold urban and rural SAF files by city
    """
    city_saf_files = expand(rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx, city=wildcards.city, habitat=HABITATS, site=wildcards.site)
    return city_saf_files
    
def get_habitat_saf_files_byCity_byPop(wildcards):
    """
    Returns list with 4fold SAF files by city by population
    """
    city_saf_files = expand(rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx, city=wildcards.city, habitat=HABITATS, site=wildcards.site)
    return city_saf_files
    
def get_city_population_saf_files(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byCity_byPop.output.saf_idx, city_pop=CITY_POP, site='4fold')
    pop1 = wildcards.pop_combo.split('.')[0] + "_." + str(wildcards.pop_combo.split('.')[1])
    pop2 = wildcards.pop_combo.split('.')[0] + "_." + str(wildcards.pop_combo.split('.')[2])
    saf1 = [x for x in all_saf_files if '{0}'.format(pop1) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '{0}'.format(pop2) in os.path.basename(x)]
    sites = expand(rules.convert_sites_for_angsd.output, site=wildcards.site)
    return {'saf' : saf1 + saf2, 'sites' : sites}


def get_city_population_saf_sfs_files(wildcards):
    all_saf_files = expand(rules.angsd_saf_likelihood_byCity_byPop.output.saf_idx, city_pop=CITY_POP, site='4fold')
    all_saf_gz_files = expand(rules.angsd_saf_likelihood_byCity_byPop.output.saf, city_pop=CITY_POP, site='4fold')
    pop1 = wildcards.pop_combo.split('.')[0] + "_." + str(wildcards.pop_combo.split('.')[1]) + "."
    pop2 = wildcards.pop_combo.split('.')[0] + "_." + str(wildcards.pop_combo.split('.')[2]) + "."
    saf1 = [x for x in all_saf_files if '{0}'.format(pop1) in os.path.basename(x)]
    saf2 = [x for x in all_saf_files if '{0}'.format(pop2) in os.path.basename(x)]
    saf_gz1 = [x for x in all_saf_gz_files if '{0}'.format(pop1) in os.path.basename(x)] 
    saf_gz2 = [x for x in all_saf_gz_files if '{0}'.format(pop2) in os.path.basename(x)] 
    all_sfs_files = expand(rules.angsd_estimate_joint_sfs_byCity_byPop.output, pop_combo=wildcards.pop_combo, site='4fold')
    combo = wildcards.pop_combo
    sfs = [x for x in all_sfs_files if '{0}'.format(combo) in os.path.basename(x)]
    sites = expand(rules.convert_sites_for_angsd.output, site=wildcards.site)
    return {'saf' : saf1 + saf2, 'sfs' : sfs, 'sites' : sites, 'saf_gz1' : saf_gz1 , 'saf_gz2' : saf_gz2 }



def get_urban_rural_bam_lists(wildcards):
    """
    Collect files with paths to urban and rural bams by City. Return as dictionary. 
    """
    urban = expand(rules.create_bam_list_byCity_byHabitat_withoutRelated.output, city=wildcards.city, habitat='u', site=wildcards.site)[0]
    rural = expand(rules.create_bam_list_byCity_byHabitat_withoutRelated.output, city=wildcards.city, habitat='r', site=wildcards.site)[0]
    return { 'urban_bams' : urban, 'rural_bams' : rural }

def get_files_for_permuted_saf_estimation(wildcards):
    """
    Get files to estimate SAF likelihoods for permuted versions of "urban" and "rural" populations
    """
    sites_idx = expand(rules.angsd_index_degenerate_sites.output.idx, site=wildcards.site)
    sites = expand(rules.convert_sites_for_angsd.output, site=wildcards.site)
    ref = rules.unzip_reference.output
    if wildcards.habitat == 'u':
        bams = expand(rules.create_random_bam_list_byCity_byHabitat.output.urban, city=wildcards.city, seed=wildcards.seed, site=wildcards.site)
    elif wildcards.habitat == 'r':
        bams = expand(rules.create_random_bam_list_byCity_byHabitat.output.rural, city=wildcards.city, seed=wildcards.seed, site=wildcards.site)
    return { 'bams' : bams, 'sites_idx' : sites_idx , 'sites' : sites, 'ref' : ref }

def get_habitat_saf_files_byCity_permuted(wildcards):
    """
    Returns list with 4fold urban and rural SAF files by city
    """
    city_saf_files = expand(rules.angsd_permuted_saf_likelihood_byCity_byHabitat.output.saf_idx, city=wildcards.city, habitat=HABITATS, site=wildcards.site, seed=wildcards.seed)
    return city_saf_files

def get_bamLists_toConcat(wildcards):
    """
    Collect text files with paths to urban and rural bams by city
    """
    all_bam_lists = expand(rules.create_bam_list_byCity_byHabitat.output, city = wildcards.city, habitat = HABITATS, site = wildcards.site)
    return all_bam_lists
    
def get_bamLists_toConcat_withoutRelated(wildcards):
    """
    Collect text files with paths to urban and rural bams by city
    """
    all_bam_lists = expand(rules.create_bam_list_byCity_byHabitat_withoutRelated.output, city = wildcards.city, habitat = HABITATS, site = wildcards.site)
    return all_bam_lists    

def get_bams_for_read_counts(wildcards):
    """
    Returns the correct GLUE or Toronto BAM file
    """
    tor_bams = expand(rules.glue_dnaSeqQC_downsample_toronto_bam.output, sample=TOR_SAMPLES)
    glue_bams = expand(rules.glue_dnaSeqQC_samtools_markdup.output.bam, sample=SAMPLES)
    glue_bams = [bam for bam in glue_bams if not os.path.basename(bam).startswith('s_')]
    return tor_bams + glue_bams

def get_dadi_sfs_input_files(wildcards):
    hab1 = '_u_'
    hab2 = '_r_'
    saf_files = expand(rules.angsd_saf_likelihood_byCity_byHabitat.output.saf_idx, city=wildcards.city, habitat=HABITATS, site='4fold')  
    sfs_files = expand(rules.angsd_estimate_sfs_byCity_byHabitat.output, city=wildcards.city, habitat=HABITATS, site='4fold') 
    saf_urban = [x for x in saf_files if '{0}'.format(hab1) in os.path.basename(x)]
    saf_rural = [x for x in saf_files if '{0}'.format(hab2) in os.path.basename(x)]
    sfs_urban = [x for x in sfs_files if '{0}'.format(hab1) in os.path.basename(x)]
    sfs_rural = [x for x in sfs_files if '{0}'.format(hab2) in os.path.basename(x)]
    ref = rules.unzip_reference.output
    return { 'saf_urban' : saf_urban , 'saf_rural' : saf_rural, 'sfs_urban' : sfs_urban, 'sfs_rural' : sfs_rural, 'ref' : ref }

