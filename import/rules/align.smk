def input_align(wildcards):
	return [x for x in Path(os.path.join(config['cache'], 'fasta')).rglob('*.seq')]

rule align:
	input:
		input_align
	output:
		os.path.join(config['cache'], 'algn', '{subdir}', '{fname}.algn')
	conda:
		"../envs/stretcher.yaml"
	run:
		
