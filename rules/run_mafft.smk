rule run_mafft:
    input:
        qry = "{query}",
        ref = REF_FASTA
    output:
        "{query}.mafft"
    conda:
        "../envs/mafft.yml"
    threads:
        10
    shell:
        r"""
            cat
            mafft --maxiterate 1000 --localpair --quiet --nuc {input.ref} {input.qry} > {output}
        """
