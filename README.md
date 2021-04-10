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


