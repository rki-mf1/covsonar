rule pairwise_align:
	input:
		qry = os.path.join(config['seqdir'], "{query}.fna"),
	output:
		fna = temp(os.path.join(config['algndir'] + "{query}")),
		algn = os.path.join(config['algndir'] + "{query}.fna")
	conda:
		"../envs/mafft.yml"
	threads: 10
	params:
		ref = ">sonar reference \n" + snr.refseq
	shell:
		r"""
			cp {input.qry} {output.fna}
			echo "{params.ref}" >> {output.fna}
			mafft --auto --quiet --nuc {output.fna} > {output.algn}
		"""
