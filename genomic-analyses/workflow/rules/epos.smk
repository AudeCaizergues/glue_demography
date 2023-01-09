# Rules to estimate populations size for each city and habitat using EPOS


###############
#### SETUP ####
###############


rule file_epos_format:
    """
    Convert sfs from angsd to EPOS format.
    """
    input:
        rules.angsd_estimate_sfs_byCity_byHabitat.output
    output:
        '{0}/epos/sfs/{{city}}_{{habitat}}_{{site}}.epos.sfs'.format(EPOS_DIR)
    shell:
        """
        touch {output}
        python3 scripts/python/angsd_sfs_to_epos.py {input} {output}

rule epos:
    """
    Run EPOS by city by habitat.
    """
    input:
        rules.file_epos_format.output
    output:
        '{0}/epos/run/{{city}}_{{habitat}}_{{site}}.epos.dat'.format(EPOS_DIR)
    shell:
        """
        bootSfs -i 1000 {input} | epos -u 1.2e-8| epos2plot > {output}
        """

rule recap_current_ne:
    """
    Built a file with all current Ne.
    """
    input:
        rules.epos_byCity_byHabitat.output
    output:
        '{0}/epos/run/ne_current_data.txt'.format(EPOS_DIR)
    shell:
        """
        ls *r_4fold.epos.dat > list_r
        ls *u_4fold.epos.dat > list_u

       for line in $(cat list_r)
        do
            echo -e "$(sed -n 2p $line)""\t""$(echo $line | cut -d'_' -f 1)""\t""rural" >> ne_current_data.txt
        done

        for line in $(cat list_u)
        do
           echo -e "$(sed -n 2p $line)""\t""$(echo $line | cut -d'_' -f 1)""\t""urban" >> ne_current_data.txt
        done
        """

rule recap_100ya_ne:
    """
    Build a file with Ne 100 years ago.
    """
    input:
        rules.epos_byCity_byHabitat.output
    output:
        '{0}/epos/run/ne_100ya_data.txt'.format(EPOS_DIR)
    shell:
        """   
        IFS=$'\n'
        ls *r_4fold.epos.dat > list_r
        for file in $(cat list_r) ; do
            N=0
            for line in $(cat $file) ; do
                N=$(($N +1))
                if [ "$N" > 2 ] ; then
                    A="$(echo $line | awk '{print $1}')"
                    if [ $A -gt 99 ]; then
                        echo -e "100""\t""$line""\t""$(echo $file | cut -d'_' -f 1)""\t""rural" >> {output}
                        break
                    fi
                fi	
            done
        done

        ls *u_4fold.epos.dat > list_u
        for file in $(cat list_u) ; do
            N=0
            for line in $(cat $file) ; do
                N=$(($N +1))
                if [ "$N" > 2 ] ; then
                    A="$(echo $line | awk '{print $1}')"
                    if [ $A -gt 99 ]; then
                        echo -e "100""\t""$line""\t""$(echo $file | cut -d'_' -f 1)""\t""urban" >> {output}
                        break
                    fi
                fi	
            done
        done
        """

rule recap_250ya_ne:
    """
    Build a file with Ne 250 years ago.
    """
    input:
        rules.epos_byCity_byHabitat.output
    output:
        '{0}/epos/run/ne_250ya_data.txt'.format(EPOS_DIR)
    shell:
        """   
    IFS=$'\n'
    ls *r_4fold.epos.dat > list_r
    for file in $(cat list_r) ; do
        N=0
        for line in $(cat $file) ; do
            N=$(($N +1))
            if [ "$N" > 2 ] ; then
                A="$(echo $line | awk '{print $1}')"
                if [ $A -gt 249 ]; then
                    echo -e "250""\t""$line""\t""$(echo $file | cut -d'_' -f 1)""\t""rural" >> {output}
                    break
                fi
            fi	
        done
    done

    ls *u_4fold.epos.dat > list_u
    for file in $(cat list_u) ; do
        N=0
        for line in $(cat $file) ; do
            N=$(($N +1))
            if [ "$N" > 2 ] ; then
                A="$(echo $line | awk '{print $1}')"
                if [ $A -gt 249 ]; then
                    echo -e "250""\t""$line""\t""$(echo $file | cut -d'_' -f 1)""\t""urban" >> {output}
                    break
                fi
            fi	
        done
    done
        """
rule recap_500ya_ne:
    """
    Build a file with Ne 500 years ago.
    """
    input:
        rules.epos_byCity_byHabitat.output
    output:
        '{0}/epos/run/ne_500ya_data.txt'.format(EPOS_DIR)
    shell:
        """   
        IFS=$'\n'
        ls *r_4fold.epos.dat > list_r
        for file in $(cat list_r) ; do
            N=0
            for line in $(cat $file) ; do
                N=$(($N +1))
                if [ "$N" > 2 ] ; then
                    A="$(echo $line | awk '{print $1}')"
                    if [ $A -gt 499 ]; then
                        echo -e "500""\t""$line""\t""$(echo $file | cut -d'_' -f 1)""\t""rural" >> {output}
                        break
                    fi
                fi	
            done
        done

        ls *u_4fold.epos.dat > list_u
        for file in $(cat list_u) ; do
            N=0
            for line in $(cat $file) ; do
                N=$(($N +1))
                if [ "$N" > 2 ] ; then
                    A="$(echo $line | awk '{print $1}')"
                    if [ $A -gt 499 ]; then
                        echo -e "500""\t""$line""\t""$(echo $file | cut -d'_' -f 1)""\t""urban" >> {output}
                        break
                    fi
                fi	
            done
        done
        """

##############
#### POST ####
##############

rule epos:
    """
    Generate empty flag file signaling successful completion of PCAngsd
    """
    input:
        expand(rules.file_epos_format.output, site = '4fold', city = CITIES, habitat=HABITATS),
        expand(rules.recap_current_ne.output, site = '4fold', city = CITIES, habitat=HABITATS)
    output:
        '{0}/epos.done'.format(EPOS_DIR)
    shell:
        """
        touch {output}
        """
