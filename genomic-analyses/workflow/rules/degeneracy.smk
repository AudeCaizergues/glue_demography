rule download_reference:
    """
    Download white clover reference genome from NCBI's FTP directory
    """
    output:
        '{0}/GCA_005869975.1_AgR_To_v5_genomic.fna.gz'.format(REF_DIR)
    log: 'logs/download_reference/download_reference.log'
    params:
        outpath = '{0}/'.format(REF_DIR)
    shell:
        """
        wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/005/869/975/GCA_005869975.1_AgR_To_v5/GCA_005869975.1_AgR_To_v5_genomic.fna.gz -P {params.outpath} 2> {log}
        """

rule unzip_reference:
    """
    Unzip reference genome
    """
    input:
        rules.download_reference.output
    output:
        '{0}/GCA_005869975.1_AgR_To_v5_genomic.fna'.format(REF_DIR)
    log: 'logs/unzip_reference/unzip_reference.log'
    shell:
        """
        gunzip {input} 2> {log}
        """

rule clone_degeneracy:
    """
    Clone Degeneracy GitHub repo for getting 4fold and 0fold sites.
    """
    output:
        temp(directory('Degeneracy'))
    log: 'logs/clone_degeneracy/clone_degeneracy.log'
    shell:
        """
        git clone https://github.com/James-S-Santangelo/Degeneracy.git 2> {log}
        """

rule get_fourfold_zerofold:
    """
    Uses get_4fold_sites.sh from Degeneracy to get 4fold and 0fold sites across white clover
    genome from reference sequence (FASTA) and annotation file (GFF).
    """
    input:
        degen = rules.clone_degeneracy.output,
        ref = rules.unzip_reference.output,
        gff = GFF_FILE
    output:
        expand('{0}/4fold_0fold/Trepens_{{site}}.bed'.format(PROGRAM_RESOURCE_DIR), site=['0fold','4fold'])
    log: 'logs/4fold_0fold/get_fourfold_zerofold.log'
    conda: '../envs/degeneracy.yaml'
    params:
        outpath = '{0}/4fold_0fold/'.format(PROGRAM_RESOURCE_DIR)
    resources:
        mem_mb = 4000,
        time = '06:00:00'
    shell:
        """
        OUT_ABS=$( readlink -f {params.outpath} );
        REF_ABS=$( readlink -f {input.ref} );
        GFF_ABS=$( readlink -f {input.gff} )
        ( cd {input.degen} &&
        bash get_4fold_sites.sh $GFF_ABS $REF_ABS $OUT_ABS ) 2> {log}
        """

rule degeneracy_done:
    """
    Write empty flag file signalling successful completion of Degenerate site extraction
    """
    input:
        expand(rules.get_fourfold_zerofold.output, site=['4fold', '0fold'])
    output:
        '{0}/degeneracy.done'.format(PROGRAM_RESOURCE_DIR)
    shell:
        """
        touch {output}
        """
