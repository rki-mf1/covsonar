<img src="logo.png"  width="134" height="134" align="right"><br><br><br><br><br>

# covSonar2

covSonar is a database-driven system for handling genomic sequences and screening genomic profiles.

What's new in covSonar V.2
* New design
    * Improve workflows
    * Performance improvements
* Exciting new features
	* Support multiple pathogens
	* Flexible in adding meta information 
* New database design
	* New database schema 
    * Retrieval efficiency
	* Significantly smaller than the previous version

## 1. Prerequisites / Setup

covSonar2 has some software-environmental requirements that can most easily be met by building a custom conda environment.

Proceed as follows to install covSonar:
```sh
# download the repository to the current working directory using git 
git clone https://github.com/rki-mf1/covsonar
# build the custom software environment using conda [recommended]
conda env create -n sonar2 -f covsonar/sonar.env.yml
# activate the conda evironment if built 
conda activate sonar2
# testing
./covsonar/test.sh
```

## 3. Usage

In covSonar2 there are several tools that can be called via subcommands.  

| subcommand | purpose                                                             |
|------------|---------------------------------------------------------------------|
| setup      | setup a new database.                             |
| import     | import genome sequences and sample information to the database     |
| list-prop  | view sample properties added to the database             | 
| add-prop    | add a sample property to the database                    |
| delete-prop       | delete a sample property from the database |
| match   |  get mutations profiles for given accessions                                        |
| restore   | restore sequence(s) from the database                                         |
| info   |  show software and database info.                                        |
| optimize   | optimizes the database                                         |
| db-upgrade   | upgrade a database to the latest version                                         |
| update-lineage-info | download latest lineage information                                       |

Each tool provides a help page that can be accessed with the `-h` option.

```sh
# display help page
./sonar.py -h 
# display help page for each tool
./sonar.py import -h 
```

### 3.1 Setup database (setup ⛽)
First, we have to create a new database instance.
```sh
./sonar.py setup --db test.db
```

, or we can create a new database instance with predefined properties.
```sh
./sonar.py setup --db test.db --auto-create
```

By default, the MN908947.3 (SARS-CoV-2) is used as a reference. If we want to set up database for a different pathogen, we can add `--gbk` following with Genbank file.

example;
```sh
./sonar.py setup --db test.db --auto-create --gbk Ebola.gb
```

>  note 📌:  [how to download genbank file](https://ncbiinsights.ncbi.nlm.nih.gov/2017/05/08/genome-data-download-made-easy/)

### 3.2 Property management (list-prop, add-prop, delete-prop)

In covSonar2, users now can arbitrarily add meta information or properties into a database to fit a specific project objective.s

To view the added properties, we can use the `list-prop` command to display all information.
```sh
./sonar.py list-prop --db test.db 
```

To add more properties, we can use the `add-prop` command to add meta information into the database.

The required arguments are listed below when we use `add-prop`
* `--name`, name of sample property
* `--descr`, description of the new property
* `--dtype`, data type of the new property (e.g., 'integer', 'float', 'text', 'date', 'zip')

```sh
./sonar.py add-prop --db test.db  --name SEQ_REASON --dtype text --descr "seq. reason"
```
> tip 🕯️:  `./sonar.py add-prop -h ` to see all available   arguments


The `delete-prop` command is used to delete an unwanted property from the database.

```sh
./sonar.py delete-prop --db test.db  --name SEQ_REASON
```

The program will ask an user for confirmation of the action.
```
Do you really want to delete this property? [YES/no]: YES
```

### 3.3 Adding genomes and meta information to the database (import)


### 3.4 Query genome sequences based on profiles (match)



####  Export DB to VCF file



### 3.5 Show infos about the used sonar system and database (info)

Detailed infos about the used sonar system (e.g. version, reference,  number of imported genomes, unique sequences, available metadata).

```sh
# Show infos about the used sonar system and database 'test.db'
./sonar.py info --db test.db 
```
### 3.6 Restore genome sequences from the database (restore)
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

### 3.7 Database management (db-upgrade, optimize)
Sometimes you might need the `optimize` command to clean the [problems](https://www.sqlite.org/lang_vacuum.html) from database operation (e.g., unused data block or storage overhead ).
```sh
# Show infos about the used sonar system
./sonar.py optimize  --db test.db 
```

When the newest version of covSonar use an old database version, covSonar will return  the following error;
```bash
Compatibility error: the given database is not compatible with this version of sonar (Current database version: XXX; Supported database version: XXX)
Please run 'sonar.py  db-upgrade' to upgrade database
```

We provide our database upgrade assistant to solve the problem. 
```bash 
# RUN 
./sonar.py db-upgrade --db test.db

# Output
Warning: Backup db file before upgrading, Press Enter to continue...

## after pressing the Enter key
Current version: 3  Upgrade to: 4
Perform the Upgrade: file: mydb.db
Database now version: 4
Success: Database upgrade was successfully completed

```
⚠️ Warning: Backup the db file before upgrade.


#### lineage-update function

Run `update-lineage-info` flag, it will download the latest version of lineages from https://github.com/cov-lineages/pango-designation/ and install in `lib/lineage.all.tsv`

```sh
# example command
./sonar.py update-lineage-info
```

## How to contribute 🏗️ 

covSonar has been very carefully programmed and tested, but is still in an early stage of development. You can contribute to this project by reporting problems 🐛 or writing feature requests to the [issue section](https://github.com/rki-mf1/covsonar/issues) 👨‍💻

Your feedback is very welcome 👨‍🔧!

---------------------------------

## Contact

For business inquiries or professional support requests 🍺 please contact [Dr. Stephan Fuchs](https://www.rki.de/SharedDocs/Personen/Mitarbeiter/F/Fuchs_Stephan.html)

