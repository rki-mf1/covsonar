<img src="logo.png"  width="134" height="134" align="right"><br><br><br><br><br>

# covSonar

covSonar is a database-driven system for handling genomic sequences of SARS-CoV-2 and screening genomic profiles.

<img src="https://img.shields.io/badge/covSonar-1.1.5-pink" />


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

| subcommand | purpose                                                             |
|------------|---------------------------------------------------------------------|
| add        | to add genome sequences to the database                             |
| update     | to import and replace meta information                              |
| match      | to query genome sequences sharing a defined profile                 | 
| restore    | to restore genome sequence(s) from the database                     |
| info       | show detailed informations about the used sonarversion and database |
| optimize   | optimize the given database                                         |
| db-upgrade | upgrade the database to the latest version.                         |
| update-lineage-info   |  update lineage information (e.g., lib/linage.all.tsv).  |

Each tool provides a help page that can be accessed with the `-h` option.

```sh
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# display help page for adding genomes
path/to/covsonar/sonar.py add -h 
```


### 3.1 Adding genomes to the database

Genome sequences of SARS-COV-2 can be added to the database in the form of FASTA files. Intermediate data is stored in a cache directory, which is temporary by default and deleted after import. The SQLite database is stored in a single file. If the defined database file does not exist, a new database is created. 
The import process can be divided into three stages:

1. caching of the sequences to be imported and calculation of sequence hashes.
2. calculating the pairwise alignment and deriving genomic profiles based on it.
3. import of the generated data into the database.
4. updating meta information of the added genomes.

Each sequence to be added is aligned pairwise to the full genome sequence of the SARS-CoV-2 isolate Wuhan-Hu-1 (NC_045512.2) using EMBOSS Stretcher and a gap-open and gap-extend penalty of 16 and 4, respectively. Gaps are left aligned.

Depending on the number of sequences to be imported and the available system resources, the import may take some time, but can be accelerated by allocating more CPUs. However, more CPUs assigned may also significantly increase the amount of available RAM. Detailed progress information and time estimates are displayed on the screen during the import.

Optionally, selected metadata (data source, data collection or lab) can be imported or updated for all genomes by using the respective option (`--source`, `--collection`, `--lab`).

If you forward the screen output to a file, the progress bar may produce plenty of useless lines. Thus, the progress bar can be disabled with the `--noprogress` option. If necessary, any output can be avoided when adding new genomes with the `--quiet` option.

```sh
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# adding all sequences from 'genomes.fasta' to database 'mydb'
# using eight cpus (the database file will be created if it does not exist)
path/to/covsonar/sonar.py add -f genomes.fasta --db mydb --cpus 8
# as before, but using a permanent cache directory to store 
# intermediate files
path/to/covsonar/sonar.py add -f genomes.fasta --db mydb --cpus 8 --cache mycache
```


### 3.2 Importing meta information

Additional meta-information can be added for each genome sequence, namely lab, data source, data collection, lineage information, zip code, collection date, GISAID and ENA identifier. Output files from Pangolin can be used directly to add the appropriate ancestry information to the available genomes. Additional information can be extracted and added from CSV or TSV files. For this, the corresponding column names from the CSV headline have to be defined by using the `--fields` option followed by the respective column expression(s):

| expression                   | description                                                                     |
|------------------------------|---------------------------------------------------------------------------------|
| accession=_colname1_         | genome accessions are listed in column _colname1_                               |
| lineage=_colname2_           | lineage information is listed in column _colname2_                              | 
| zip=_colname3_               | zip codes are listed in column _colname3_                                       |
| date=_colname4_              | sampling dates are listed in column _colname4_ (needed date format: YYYY-MM-DD) |
| lab=_colname5_               | lab information is listed in column _colname5_                                  | 
| source=_colname6_            | data source is listed in column _colname6_                                      | 
| collection=_colname7_        | data collection is listed in column _colname7_                                  | 
| technology=_colname8_        | used sequencing technology is listed in column _colname8_                       | 
| platform=_colname9_          | used sequencing platform is listed in column _colname9_                         | 
| chemistry=_colname10_        | used sequencing chemistry is listed in column _colname10_                       | 
| software=_colname11_         | software used for genome reconstruction is listed in column _colname11_         | 
| software_version=_colname12_ | software version used for genome reconstruction is listed in column _colname12_ | 
| material=_colname13_         | sampling material is listed in column _colname13_                               | 
| ct=_colname14_               | ct values are listed in column _colname14_                                      | 


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

*More add-on feild to support custom scenario 

| expression                   | description                                                                                                |
|------------------------------|------------------------------------------------------------------------------------------------------------|
| submission_date=_colname_    | This one can be used when the sample is submitted for processing or prcoessing date after sampling date.   | 


### 3.3 Query genome sequences based on profiles 

Genomic profiles can be defined to align genomes. For this purpose, the variants related to the complete genome of the SARS-CoV-2 isolate Wuhan-Hu-1 (NC_045512.2) must be expressed as follows:

| typ       | nucleotide level                                                  | amino acid level              |
|-----------|-------------------------------------------------------------------|-------------------------------|
| SNP       | ref_nuc _followed by_ ref_pos _followed by_ alt_nuc (e.g. A3451T) | protein_symbol:ref_aa _followed by_ ref_pos _followed by_ alt_aa (e.g. S:N501Y) |
| deletion  | del:ref_pos:length_in_bp (e.g. del:3001:8)                        | protein_symbol:del:ref_pos:length_in_aa (e.g. ORF1ab:del:3001:21) | 
| insertion | ref_nuc _followed by_ ref_pos _followed by_ alt_nucs (e.g. A3451TGAT) | protein_symbol:ref_aa _followed by_ ref_pos _followed by_ alt_aas (e.g. N:A34AK)  |  

The positions refer to the reference (first nucleotide in the genome is position 1). Using the option `-i` multiple variant definitions can be combined into a nucleotide, amino acid or mixed profile, which means that matching genomes must have all those variations in common. In contrast, alternative variations can be defined by multiple `-i` options. As an example, `-i S:N501Y S:E484K` matches genomes sharing the _Nelly_ **AND** _Erik_ variation while `-i S:N501Y -i S:E484K` matches to genomes that share either the _Nelly_ **OR** _Erik_ variation **OR** both. Accordingly, using the option `-e` profiles can be defined that have not to be present in the matched genomes. 

To filter genomes based on metadata specific options can be used (see table below). Only genomes linked to the respective metadata are then considered. Metadata values are negated when introduced by ^ (e.g. `--acc ^ID1` matches all genomes accessions but ID1). Metadata filtering is case-insensitive. To see the amount of available metadata in your database use the info tool (see section 3.5). 

| option              | value(s)                                                              | note |
|---------------------|-----------------------------------------------------------------------|------| 
| --acc               | one or more genome accessions (e.g. NC_045512.2)                      |      |
| --lineage           | one or more pangolin lineages (e.g. B.1.1.7)                          |      |
| --zip               | one or more zip codes (e.g. 10627)                                    | zip codes are dynamically extended to the right side, e.g. 033 matches to all zip codes starting with 033|
| --date              | one or more dates or date ranges (e.g. 2021-01-01)                    | single dates are formatted as YYYY-MM-DD while date ranges can be defined by YYYY-MM-DD:YY-MM-DD (from:to) |
| --submission_date   | one or more dates or date ranges (e.g. 2021-01-01)                    | single dates are formatted as YYYY-MM-DD while date ranges can be defined by YYYY-MM-DD:YY-MM-DD (from:to) |
| --lab               | one or more labs (e.g. L1)                                            |      |
| --source            | one or more data sources (e.g. DESH)                                  |      |
| --collection        | one or more data collections (e.g. RANDOM)                            |      |
| --technology        | one or more sequencing technologies (e.g. Illumina)                   |      |
| --platform          | one or more sequencing platforms (e.g. MiSeq)                         |      |
| --chemistry         | one or more sequencing chemistries (e.g. Cleanplex)                   |      |
| --software          | one software tool used for genome reconstruction (e.g. covPipe)       |      |
| --version           | one software tool version used for genome reconstruction (e.g. 3.0.5) | needs --software defined |
| --material          | one or more sample materials (e.g. 'nasal swap')                      |      |
| --min_ct            | minimal ct value (e.g. 20)                                            |      |
| --max_ct            | maximal ct value (e.g. 20)                                            |      |
 

There are additional options to adjust the matching.

| option             | description                                                            |
|--------------------|------------------------------------------------------------------------|
| --count            | count matching genomes only                                            |
| --ambig            | include purely ambiguous variations in the profiles (e.g. N stretches) |
| --only_frameshifts | show only genomes containing frameshift mutations                      |
| --no_frameshifts   | show only genomes containing no frameshift mutations                   |

By default, genome matching produces a comma-separated output (csv). Using the option `--tsv` the output will be tab-delimited (tsv).


```sh
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# matching B.1.1.7 genomes in DB 'mydb' that share an additional "Erik" mutation 
path/to/covsonar/sonar.py match -i S:E484K --lineage B.1.1.7 --db mydb
# as before but matching genomes are counted only
path/to/covsonar/sonar.py match -i S:E484K --lineage B.1.1.7 --count --db mydb
# matching genomes in DB 'mydb' sharing the "Nelly" but not the "Erik" mutation
# and that were sampled in 2020
path/to/covsonar/sonar.py match -i S:N501Y -e S:E484K --date 2020-01-01:2020-12-31 --db mydb
# matching genomes in DB 'mydb' sharing the "Nelly" and the "Erik" mutation but not
# belonging to the B.1.1.7 lineage
path/to/covsonar/sonar.py match -i S:N501Y S:E484K --lineage ^B.1.1.7 --db mydb
# as before but redirect the ouptut to a CSV file named out.csv
path/to/covsonar/sonar.py match -i S:N501Y S:E484K --lineage ^B.1.1.7 --db mydb > out.csv
```

**Wildcard**

CovSonar also supports wildcard query (symbol `%`) to handle more complexity of the query. The operator is used in a **lineage query** to search for a specified pattern in lineage. It can be used in combinations of including and excluding command, for example;

```sh
# edit command at --lineage tag
path/to/covsonar/sonar.py match -i S:N501Y S:E484K --lineage ^B.1.1% AY.4% --db mydb > out.csv
```

This query will include all lineages that start with 'AY.4'.
```
['AY.4', 'AY.4.1', 'AY.4.2', .... , 'AY.46.5', 'AY.46.6', 'AY.47']
```

and exclude all lineages that start with 'B.1.1'.
```
['B.1.1', 'B.1.1.1', 'B.1.1.10', 'B.1.1.121', ... , 'B.1.177.9', 'B.1.179', 'B.1.187']
```
Here are some examples showing different queries with `%`

| query              | description                                       |
|--------------------|---------------------------------------------------|
| AY%                | Finds any lineage that starts with "AY"           |
| %1.1.%             | Finds any lineage that have "1.1." in any position|
| %.1                | Finds any lineage that ends with ".1"             |

**Parent-Child relationship**

If we want to query all sublineages, covSonar offers `--with-sublineage` tag along with match command.

```sh
# We want to get all sublineages of delta variant (B.1.617.2).
# This query will return result from ['B.1.617.2', 'AY.1', 'AY.2', ..., , 'AY.129', 'AY.130']
path/to/covsonar/sonar.py match -i S:N501Y --lineage  B.1.617.2 --with-sublineage --db mydb > out.csv
```
**Note:** covSonar offers lineage-update function. Run `update-lineage-info` flag, it will download the latest version of lineages from https://github.com/cov-lineages/pango-designation/ and install in `lib/lineage.all.tsv`
```sh
# example command
path/to/covsonar/sonar.py update-lineage-info
```
### 3.4 Restore genome sequences from the database

Genome sequences can be restored from the database based on their accessions.
The restored sequences are combined with their original FASTA header and  shown on the screen. The screen output can be redirected to a file easily by using `>`.

```sh
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# Restore genome sequences linked to accessions 'mygenome1' and 'mygenome2' from the 
# database 'mydb' and write these to a fasta file named 'restored.fasta'
path/to/covsonar/sonar.py restore --acc mygenome1 mygenome2 --db mydb > restored.fasta
# as before, but consider all accessions from 'accessions.txt' (the file has to
# contain one accession per line) 
path/to/covsonar/sonar.py restore --file accessions.txt --db mydb > restored.fasta
```


### 3.5 Show infos about the used sonar system and database

Detailed infos about the used sonar system (e.g. version, reference, considered ORFs) and, optionally, a given database (e.g. number of imported genomes, unique sequences, available metadata) can be accessed.

```sh
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# Show infos about the used sonar system
path/to/covsonar/sonar.py info
# Show infos about the used sonar system and database 'mydb'
path/to/covsonar/sonar.py info --db mydb
```

### 3.6 Deleting genomes from the database

If needed, genomes can be removed from the database using the respective genome accession(s)

```sh
# activating conda environment if built and not active yet (see section 2)
conda activate sonar
# removing sequence ACC1 and ACC5 to from 'mydb'
path/to/covsonar/sonar.py remove --acc ACC1 ACC5 --db mydb
# removing accessions listed in file to_delete.txt (one accession per line) from 'mydb'
path/to/covsonar/sonar.py remove --file to_delete.txt --db mydb
```

### 3.7 Export DB to VCF file

covSonar can export accession records in a VCF format using the `var2vcf` command. The output from this feature is a single VCF file that combines all accessions. The output format is in **.gz** form.

```sh
# Export all accessions in the database.
path/to/covsonar/sonar.py var2vcf --db mydb -o merge.vcf
# Just like the option in the match command, we can use  --file, --acc and --date to enable specific accession export.
path/to/covsonar/sonar.py var2vcf --db mydb -f acc.10.txt -o merge.vcf
# To speed up the query, we can use --cpus tag to aid a query performance.
path/to/covsonar/sonar.py var2vcf --db mydb --date 2021-08-01:2021-08-10 -o merge.vcf --cpus 20

# Another solution, we can use --betaV2 tag (x3-5 times faster), 
# The current version is under development, so if you found any bug please report it to us.
path/to/covsonar/sonar.py var2vcf --db mydb --date 2021-08-01:2021-08-10 -o merge.vcf --cpus 20 --betaV2
```
:warning:**Note:** The current performance of this feature still does not perform well (e.g., memory usage and runtime) when trying to export many accessions. However, we are constantly working on improving the performance. :monkey: 

## 4 How to contribute

covSonar has been very carefully programmed and tested, but is still in an early stage of development. You can contribute to this project by reporting problems or writing feature requests to the issue section under https://gitlab.com/s.fuchs/covsonar/-/issues

Your feedback is very welcome!

## 5 FAQ

**Q:** How can I screen for genomes with any SNP or amino acid substitution at a specific position at the genome or gene product?


**A:** covSonar accepts and interpretes the IUPAC nucleotide and amino acid code. Accordingly, you can screen for any nucleotide or amino acid at a certain position using N and X, respectively. Since covSonar stores only sites different from the reference, the reference nucleotide will be not considered when searching for N or X. As an example, use the following command to screen for all genome encoding for any amino acid substitution at position 484 within the Spike protein (reference allele is E at this position).  

```bash 
# screen for all genomes encoding for any amino acid substitution at position 484 within the Spike protein
# using the database 'mydb"
path/to/sonar.py -i S:E484X --db mydb
```


**Q:** Running ´sonar.py optimize´ on my database returns the error "sqlite3.OperationalError: database or disk is full" even I have enough space on the disk where the database is located.


**A:** This happens, when the sqlite temprory directory that might be located on another disk has not enough space to store the intermediate files. You can easily change the temporary directory e.g. to the current working directory using the following Shell command (the changes will only apply to the current Shell session):

```bash 
# changing the sqlite temporary directory to the current working directory 
# replace . by the path to the location you want to use
export SQLITE_TMPDIR="."
```

After executing this command, optimizing your database should work.

**Q:** covSonar returns the error;

```bash
Compatibility error: the given database is not compatible with this version of sonar (Current database version: XXX; Supported database version: XXX)
Please run 'sonar.py  db-upgrade' to upgrade database
```

**A:** This happens, when you use the newest version of covSonar with old database version.

Please use our database upgrade assistant to solve the problem. 
```bash 
# RUN 
python sonar.py db-upgrade --db mydb.db

# Output

Warning: Backup db file before upgrading, Press Enter to continue...

## press Enter
Current version: 3  Upgrade to: 4
Perform the Upgrade: file: mydb.db
Database now version: 4
Success: Database upgrade was successfully completed

```
:warning: Warning: Backup the db file before upgrade.


