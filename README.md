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
| match      | to query genome sequences sharing a defined profile    | 
| restore    | to restore genome sequence(s) from the database        |
| update     | to import and replace meta information                 |


### 3.1 Adding genomes to the database

Genome sequences of SARS-COV-2 can be added to the database in the form of FASTA files. Intermediate data is stored in a cache directory, which is temporary by default and deleted after import. The SQLite database is stored in a single file that has to be defined. If the defined database file does not exist, a new database is created.

```sh
# activate conda environment if built (see section 2)
cond activate sonar
# adding all sequences from 'genomes.fasta' to database 'mydb'
path/to/covsonar/sonar.py add -f genomes.fasta --db mydb
