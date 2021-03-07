import re

rule inspect_align:
	input:
		algn = os.path.join(config['algndir'] + "{sbj}.fna")
	output:
		os.path.join(config['algndir'] + "{sbj}.fna.done")
	threads: 1
	run:
		seqhash = os.path.basename(input.algn)[:-4]

		# process alignment file
		seqs = {"ref": [], "sbj": []}
		with open(input.algn, "r") as handle:
			for line in handle:
				if line == ">sonar reference \n":
					key = "ref"
				elif line.startswith(">"):
					key = "sbj"
					seqhash = line[1:].strip()
				else:
					seqs[key].append(line.strip().upper())

		snr.add_alignment("".join(seqs["sbj"]), "".join(seqs["ref"]), seqhash)
		snr.commit()
		open(output[0], "w").close()
