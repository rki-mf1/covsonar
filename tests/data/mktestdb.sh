#!/usr/bin/env bash

# This is a workaround to the current problem that the import command causes
# pytest to crash. Until that is resolved we can include a small test database
# to do a few downstream tests.

set -x

# setting database up
rm -f test.db
sonar setup --db test.db --debug
sonar add-prop --db test.db --name DATE_DRAW --dtype date --descr descr
sonar add-prop --db test.db --name SEQ_TYPE --dtype text --descr descr
sonar add-prop --db test.db --name SEQ_REASON --dtype text --descr descr
sonar add-prop --db test.db --name SAMPLE_TYPE --dtype text --descr descr
sonar add-prop --db test.db --name OWN_FASTA_ID --dtype text --descr descr
sonar add-prop --db test.db --name RECEIVE_DATE --dtype date --descr descr
sonar add-prop --db test.db --name PROCESSING_DATE --dtype date --descr descr
sonar add-prop --db test.db --name SENDING_LAB_PC --dtype zip --descr descr
sonar add-prop --db test.db --name SEQUENCING_LAB_PC --dtype zip --descr descr
sonar add-prop --db test.db --name GISAID_ACCESSION --dtype text --descr descr
sonar add-prop --db test.db --name LINEAGE --dtype pango --descr descr
sonar update-lineages --db test.db
cp test.db test-with-seqs.db
# Do these two steps together in one go
#sonar import --db test-with-seqs.db --fasta seqs.fasta.gz --cache cache --threads 4
#sonar import --db test-with-seqs.db --tsv meta.tsv --cols sample=IMS_ID --auto-link --cache cache --threads 4 --debug
sonar import --db test-with-seqs.db --fasta seqs.fasta.gz --tsv meta.tsv --cols sample=IMS_ID --auto-link --cache cache --threads 4 --debug
sonar import --db test-with-seqs.db --tsv pango.tsv --cols sample=IMS_ID LINEAGE=lineage --cache cache
rm -rf cache import.log
