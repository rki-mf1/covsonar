rule align:
	input:
		aseq = os.path.join(config['seq_dir'], '{subdir}', '{fname}' + '.aseq'),
		bseq = os.path.join(config['seq_dir'], '{subdir}', '{fname}' + '.bseq')
	output:
		os.path.join(config['algn_dir'], '{subdir}', '{fname}.algn')
	conda:
		"../envs/emboss.yaml"
	threads: 1
	shell:
		'''
		stretcher -asequence {input.aseq} -bsequence {input.bseq} -outfile {output} -auto -aformat markx3
		'''
