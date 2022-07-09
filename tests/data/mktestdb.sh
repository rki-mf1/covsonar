#!/usr/bin/env bash

# This is a workaround to the current problem that the import command causes
# pytest to crash. Until that is resolved we can include a small test database
# to do a few downstream tests.

set -x

# setting database up
rm -f test.db
sonar setup --db test.db
sonar add-prop --db test.db --name SENDING_LAB --dtype integer --descr descr
sonar add-prop --db test.db --name DATE_DRAW --dtype date --descr descr
sonar add-prop --db test.db --name SEQ_TYPE --dtype text --descr descr
sonar add-prop --db test.db --name SEQ_REASON --dtype text --descr descr
sonar add-prop --db test.db --name SAMPLE_TYPE --dtype text --descr descr
sonar add-prop --db test.db --name OWN_FASTA_ID --dtype text --descr descr
sonar add-prop --db test.db --name DOWNLOAD_ID --dtype text --descr descr
sonar add-prop --db test.db --name DEMIS_ID --dtype integer --descr descr
sonar add-prop --db test.db --name RECEIVE_DATE --dtype date --descr descr
sonar add-prop --db test.db --name PROCESSING_DATE --dtype date --descr descr
sonar add-prop --db test.db --name PUBLICATION_STATUS --dtype text --descr descr
sonar add-prop --db test.db --name HASHED_SEQUENCE --dtype text --descr descr
sonar add-prop --db test.db --name TIMESTAMP --dtype text --descr descr
sonar add-prop --db test.db --name STUDY --dtype text --descr descr
sonar add-prop --db test.db --name DOWNLOADING_TIMESTAMP --dtype text --descr descr
sonar add-prop --db test.db --name SENDING_LAB_PC --dtype zip --descr descr
sonar add-prop --db test.db --name DEMIS_ID_PC --dtype zip --descr descr
sonar add-prop --db test.db --name VERSION --dtype integer --descr descr
sonar add-prop --db test.db --name DESH_QC_PASSED --dtype text --descr descr
sonar add-prop --db test.db --name DESH_REJECTION_REASON --dtype text --descr descr
sonar add-prop --db test.db --name DUPLICATE_ID --dtype text --descr descr
sonar add-prop --db test.db --name LINEAGE --dtype text --descr descr
cp test.db test-with-seqs.db
sonar import --db test-with-seqs.db --tsv meta.tsv --fasta seqs.fasta.gz --cols sample=IMS_ID --cache cache --threads 1
rm -rf cache import.log
