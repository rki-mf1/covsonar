rule diff:
	input:
		os.path.join(config['algn_dir'], '{subdir}', '{fname}.algn')
	output:
		os.path.join(config['diff_dir'], '{subdir}', '{fname}.diff')
	threads: 1
	script:
		"../scripts/markx3_to_diff.py"
