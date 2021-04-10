<img src="logo.png"  width="171" height="200" align="right"><br><br><br><br><br><br><br>


# covSonar

covSonar is a database-driven system for handling genomic sequences of SARS-CoV-2 and screening genomic profiles.


## 1. Prerequisites

covSonar has some software-environmental requirements that can most easily be met by building a custom conda environment as described in Section 2.

| software/module       | version  |
|-----------------------|----------|
| python                | 3.9.2    |
| sqlite                | 3.34.0   |
| emboss                | 6.6.0    |
| biopython             | 1.78     |
| numpy                 | 1.20.1   |
| packaging             | 20.9     |
| tqdm                  | 4.59.0   |


## 2. Setup
Proceed as follows to install covSonar:
```sh
# download the repository to the current working directory using git 
git clone https://gitlab.com/s.fuchs/covsonar
# build the custom software environment using conda [recommended]
conda env create -n sonar -f covsonar/sonar.env.yml
# activate the conda evironment if built 
conda activate sonar
# testing
./covsonar/test.sh
```


## 3. Usage

In covSonar there are several tools that can be called via subcommands. 

| subcommand | purpose                                                |
|------------|--------------------------------------------------------|
| add        | to add genome sequences to the database                |
| update     | to import and replace meta information                 |
| match      | to query genome sequences sharing a defined profile    | 
| restore    | to restore genome sequence(s) from the database        |


### 3.1 Adding genomes to the database

Genome sequences of SARS-COV-2 can be added to the database in the form of FASTA files. Intermediate data is stored in a cache directory, which is temporary by default and deleted after import. The SQLite database is stored in a single file that has to be defined. If the defined database file does not exist, a new database is created. 
The import process can be divided into three stages:

1. caching of the sequences to be imported and calculation of sequence hashes.
2. calculating the pairwise alignment and deriving genomic profiles based on it.
3. import of the generated data into the database.

Depending on the number of sequences to be imported and the available system resources, the import may take some time. The import can be accelerated by allocating more CPUs. However, do not underestimate that this may also significantly increase the amount of available RAM. In any case, detailed progress information and time estimates are displayed on the screen during the import.

```sh
# activating conda environment if built (see section 2)
conda activate sonar
# adding all sequences from 'genomes.fasta' to database 'mydb'
# using eight cpus
path/to/covsonar/sonar.py add -f genomes.fasta --db mydb --cpus 8
# as before, but using a permanent cache directory to store 
# intermediate files
path/to/covsonar/sonar.py add -f genomes.fasta --db mydb --cpus 8 --cache mycache
```


### 3.2 Importing meta information

Additional meta-information can be added for each genome sequence, namely lineage information, zip code, collection date, GISAID and ENA identifier. Output files from Pangolin can be used directly to add the appropriate ancestry information to the available genomes. Additional information can be extracted and added from CSV or TSV files. For this, the corresponding column names from the headline have to be defined as follows:

| expression           | description                                        |
|----------------------|----------------------------------------------------|
| accession=_colname1_ | genome accessions are listed in column _colname1_  |
| lineage=_colname2_   | lineage information is listed in column _colname2_ | 
| zip=_colname3_       | zip codes are listed in column _colname3_          |
| date=_colname4_      | sampling dates are listed in column _colname4_     |

```sh
# activating conda environment if built (see section 2)
conda activate sonar
# importing lineage information from pangolin output file 
# to database 'mydb'
path/to/covsonar/sonar.py update --pangolin pangolin.csv --db mydb
# importing zip codes and sampling dates from a custom CSV file
# to database 'mydb'
path/to/covsonar/sonar.py update --csv custom.csv --fields accession=acc zip=zip_codes date=sampling --db mydb
```


### 3.3 Query genome sequences based on genomic profiles 

Genomic profiles can be defined to align genomes. For this purpose, the variants related to the complete genome of the SARS-CoV-2 isolate Wuhan-Hu-1 (NC_045512.2) must be expressed as follows:

| typ       | nucleotide level                                                  | amino acid level              |
|-----------|-------------------------------------------------------------------|-------------------------------|
| SNP       | ref_nuc _followed by_ ref_pos _followed by_ alt_nuc (e.g. A3451T) | protein_symbol:ref_aa _followed by_ ref_pos _followed by_ alt_aa (e.g. S:N501Y) |
| deletion  | del:ref_pos:length_in_bp (e.g. del:3001:8)                        | protein_symbol:del:ref_pos:length_in_aa (e.g. ORF1ab:del:3001:21) | 
| insertion | ref_nuc _followed by_ ref_pos _followed by_ alt_nucs (e.g. A3451TGAT) | protein_symbol:ref_aa _followed by_ ref_pos _followed by_ alt_aas (e.g. N:A34AK)  |  

