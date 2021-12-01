
rule parse_fasta:
	input:
		config['fasta']
	output:
		dir = directory(TEMPDIR)
	threads: 1
	run:
		with sonarDBManager(db, debug=DEBUG) as dbm:
			for fname in input:
				for record in SeqIO.parse(fname,'fasta'):
					id, header, sequence = record.id, record.description, sonarActions.harmonize(str(record.seq))
					seqhash = sonarActions.hash(sequence)

					# extract metdata from fasta header
					properties = {}
					for match in self.fasta_tag_regex.finditer(header):
						properties[match.group(1)] = match.group(2)

					molecule_accession = molecules[molecules[properties]['molecule']['id']] if "molecule" in properties else config['default_molecule']
					molecule_id = molecules[molecule_accession]['molecule.id']
					molecule_len = molecules[molecule_accession]['molecule.length']
