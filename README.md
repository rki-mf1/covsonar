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

Each tool provides a help page that can be accessed with the `-h` option.

```sh
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# display help page for adding genomes
path/to/covsonar/sonar.py add -h 
```


### 3.1 Adding genomes to the database

Genome sequences of SARS-COV-2 can be added to the database in the form of FASTA files. Intermediate data is stored in a cache directory, which is temporary by default and deleted after import. The SQLite database is stored in a single file that has to be defined. If the defined database file does not exist, a new database is created. 
The import process can be divided into three stages:

1. caching of the sequences to be imported and calculation of sequence hashes.
2. calculating the pairwise alignment and deriving genomic profiles based on it.
3. import of the generated data into the database.

Each sequence to be added is aligned pairwise to the full genome sequence of the SARS-CoV-2 isolate Wuhan-Hu-1 (NC_045512.2) using EMBOSS Stretcher and a gap-open and gap-extend penalty of 16 and 4, respectively. Gaps are left aligned.

Depending on the number of sequences to be imported and the available system resources, the import may take some time. The import can be accelerated by allocating more CPUs. However, do not underestimate that this may also significantly increase the amount of available RAM. In any case, detailed progress information and time estimates are displayed on the screen during the import.

```sh
# activating conda environment if built and not active yet (see section 2)
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
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# importing lineage information from pangolin output file 
# to database 'mydb'
path/to/covsonar/sonar.py update --pangolin pangolin.csv --db mydb
# importing zip codes and sampling dates from a custom CSV file
# to database 'mydb'
path/to/covsonar/sonar.py update --csv custom.csv --fields accession=acc zip=zip_codes date=sampling --db mydb
```


### 3.3 Query genome sequences based on profiles 

Genomic profiles can be defined to align genomes. For this purpose, the variants related to the complete genome of the SARS-CoV-2 isolate Wuhan-Hu-1 (NC_045512.2) must be expressed as follows:

| typ       | nucleotide level                                                  | amino acid level              |
|-----------|-------------------------------------------------------------------|-------------------------------|
| SNP       | ref_nuc _followed by_ ref_pos _followed by_ alt_nuc (e.g. A3451T) | protein_symbol:ref_aa _followed by_ ref_pos _followed by_ alt_aa (e.g. S:N501Y) |
| deletion  | del:ref_pos:length_in_bp (e.g. del:3001:8)                        | protein_symbol:del:ref_pos:length_in_aa (e.g. ORF1ab:del:3001:21) | 
| insertion | ref_nuc _followed by_ ref_pos _followed by_ alt_nucs (e.g. A3451TGAT) | protein_symbol:ref_aa _followed by_ ref_pos _followed by_ alt_aas (e.g. N:A34AK)  |  

The position specifications refer to the reference in each case and are 1-based. Using the option `-i` multiple variant definitions can be combined into a nucleotide, amino acid, or mixed profile, which means that matching genomes must have all defined variations in common. In contrast, alternative variations can be defined by multiple `-i` options. As an example, `-i S:N501Y S:E484K` matches genomes sharing the _Nelly_ **AND** _Erik_ variation while  `-i S:N501Y -i S:E484K` matches to genomes that share either the _Nelly_ **OR** _Erik_ variation **OR** both. Accordingly, using the option `-e` profiles can be defined that have not to be present in the matched genomes. 

To consider only genomes of a certain lineage, zip code or samplig date, option `--lineage`, `--zip` or `--date` can be used followed by one or more values. To negate a value has to be introduced by ^. As an example, `--lineage B.1.1.7` matches only genomes of the so-called UK variant, while `--lineage B.1.1.7` matches all genomes **NOT** assigned to this lineage. Please consider that zip codes are hierarchically matched, meaning that `--zip 114` includes all zip codes starting with 114. Single dates are formatted as _YYYY-MM-DD_ while date ranges are defined as _from:to_ (_YYYY-MM-DD:YYYY-MM-DD_). 

By default, additional variations are allowed in the matched genomes. Using the `--exclusive` option, genomes with additional variations (including ambiguities such as N) are excluded.

The Output is shown on screen but can be easily rdirected to a file by expanding the command by `> output.csv`. The output contains comma separated values for each matched genome in the following order:

- accession of the matched genome 
- lineage 
- zip code 
- sampling date 
- nucleotide level profile
- amino acid level profile 

By deafult, variations with ambiguities in the alternatie allele (such as N) are not shown which can be changed using the `--ambig` option.

```sh
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# matching B.1.1.7 genomes in DB 'mydb' that share an additional "Erik" mutation 
path/to/covsonar/sonar.py match -i S:E484K --lineage B.1.1.7 --db mydb
# matching genomes in DB 'mydb' sharing the "Nelly" but not the "Erik" mutation
# and that were sampled in 2020
path/to/covsonar/sonar.py match -i S:N501Y -e S:E484K --date 2020-01-01:2020-12-31 --db mydb
# matching genomes in DB 'mydb' sharing the "Nelly" and the "Erik" mutation but not
# belonging to the B.1.1.7 lineage
path/to/covsonar/sonar.py match -i S:N501Y S:E484K --lineage ^B.1.1.7 --db mydb
# as before but redirect the ouptut to a CSV file named out.csv
path/to/covsonar/sonar.py match -i S:N501Y S:E484K --lineage ^B.1.1.7 --db mydb > out.csv
```


### 3.4 Restore genome sequences from the database

Genome sequences can be restored from the database based on their accessions.
The restored sequences are combined with their original FASTA header and  shown on the screen. The screen output can be redirected to a file easily by using `>`.

```sh
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# Restore genome sequences with accessions 'mygenome1' and 'mygenome2' 
# and write to a fasta file named 'restored.fasta'
path/to/covsonar/sonar.py restore --acc mygenome1 mygenome2 > restored.fasta
```

## 4 How to contribute

covSonar has been very carefully programmed and tested, but is still in an early stage of development. You can contribute to this project by reporting problems or writing feature requests to the issue section under https://gitlab.com/s.fuchs/covsonar/-/issues

Your feedback is very welcome!

