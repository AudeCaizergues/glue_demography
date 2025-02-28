# Rules to perform demographic modelling using dadi

rule dadi_sfs:
    input:
        unpack(get_dadi_sfs_input_files)
    output:
        '{0}/{{city}}_dadi_fromAngsd.sfs'.format(DADI_DIR),
    log: LOG_DIR + '/dadi_sfs/{city}_dadi_sfs.log'
    container: 'library://james-s-santangelo/angsd/angsd:0.933'
    resources:
        mem_mb = lambda wildcards, attempt: attempt * 4000,
        time = '01:00:00'
    threads: 6
    shell:
        """
        realSFS dadi \
            {input.saf_urban} {input.saf_rural} \
            -sfs {input.sfs_urban} \
            -sfs {input.sfs_rural} \
            -ref {input.ref} \
            -anc {input.ref} \
            -seed 42 -P {threads} -maxIter 2000 > {output} 2> {log}
        """

#rule format_dadi_sfs:
#    input:
#        rules.dadi_sfs.output
#    output:
#        '{0}/{{city}}_dadi_formatted.sfs'.format(DADI_DIR)
#    params:
#        pop1_n = 41,
#        pop2_n = 41
#    shell:
#        """
#        
#        perl scripts/perl/realsfs2dadi.pl {input} {params.pop1_n} {params.pop2_n} > {output}
#        """

rule dadi_sfs_done:
    input:
        expand(rules.dadi_sfs.output, city=CITIES)
    output:
        '{0}/sfs_for_dadi.done'.format(DADI_DIR)
    shell:
        """
        touch {output}
        """
