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

[![Tests Status](https://github.com/rki-mf1/covsonar/actions/workflows/tests.yml/badge.svg)](https://github.com/rki-mf1/covsonar/actions/workflows/tests.yml) [![Pre-release/Test Version](https://github.com/rki-mf1/covsonar/actions/workflows/test-pypi.yml/badge.svg)](https://github.com/rki-mf1/covsonar/actions/workflows/test-pypi.yml) [![Stable Version](https://github.com/rki-mf1/covsonar/actions/workflows/release.yml/badge.svg)](https://github.com/rki-mf1/covsonar/actions/workflows/release.yml)

[![Github tag](https://badgen.net/github/release/rki-mf1/covsonar/Strapdown.js)](https://github.com/rki-mf1/covsonar/release/)


## 1. Installation

covSonar2 now can be easily installed by using `pip`, `conda` or `mamba`

### Stable version.🔖
Proceed as follows to install covSonar2

```sh
pip install covsonar
```
or via conda
```sh
conda install covsonar
```

verify installation
```sh
sonar --version
```

### Dev. version.🚧

(Testing or Nightly build or Adding a feature or Fixing a bug)

please check at https://test.pypi.org/project/covsonar/

```sh
# example command
pip install -i https://test.pypi.org/simple/ covsonar==2.0.0a1.dev1657097995
```

## 2. Usage

In covSonar2, the table below shows the several commands that can be called.

| subcommand | purpose                                                             |
|------------|---------------------------------------------------------------------|
| setup      | set up a new database.                                  |
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
sonar -h
# display help page for each tool
sonar import -h
```

### 2.1 Setup database (setup ⛽)

First, we have to create a new database instance.
```sh
sonar setup --db test.db
```

Or we can create a new database instance with predefined properties.
```sh
sonar setup --db test.db --auto-create
```
> TIP 🕯️: We can use [DB Browser](https://sqlitebrowser.org/) to visualise or manipulate the database file.

By default, the MN908947.3 (SARS-CoV-2) is used as a reference. If we want to set up a database for a different pathogen, we can add `--gbk` following with Genbank file. (⚠️WARNING: currently, covSonar2 is tailored for SARS-CoV-2, so some functions might not function with other pathogens.)

Example;
```sh
sonar setup --db test.db --auto-create --gbk Ebola.gb
```

> NOTE 📌:  [how to download genbank file](https://ncbiinsights.ncbi.nlm.nih.gov/2017/05/08/genome-data-download-made-easy/)

### 2.2 Property management (`list-prop`, `add-prop` and `delete-prop`)

In covSonar2, users can now arbitrarily add meta information or properties into a database to fit a specific project objective.

To add properties, we can use the `add-prop` command to add meta information into the database.

The required arguments are listed below when we use `add-prop` command
* `--name`, name of sample property
* `--descr`, description of the new property
* `--dtype`, data type of the new property (e.g., 'integer', 'float', 'text', 'date', 'zip')

```sh
# for example
sonar add-prop --db test.db --name LINEAGE --dtype text --descr "store Lineage"
#
sonar add-prop --db test.db --name AGE --dtype integer --descr "age information"
#
sonar add-prop --db test.db --name DATE_DRAW --dtype date --descr "sampling date"
```
> TIP 🕯️: `sonar add-prop -h ` to see all available arguments.

⚠️ WARNING: We reserve **'sample'** keyword that cannot be used as a property name
(e.g., ⛔❌`--name sample`) because we use this name as the ID in the database schema.⚠️

To view the added properties, we can use the `list-prop` command to display all information.
```sh
sonar list-prop --db test.db
```

The `delete-prop` command is used to delete an unwanted property from the database.

```sh
sonar delete-prop --db test.db  --name SEQ_REASON
```

The program will ask for confirmation of the action.
```
Do you really want to delete this property? [YES/no]: YES
```

### 2.3 Adding genomes and meta information to the database (import)

This example shows how we add sequence along with meta information.

We have sequence file name `valid.fasta` and meta-info file name `day.tsv`.

valid.fasta
```
>IMS-00113
CCAACCAACTTTCGATCTCTTG
```
day.tsv
```
IMS_ID		SAMPLING_DATE	LINEAGE
IMS-00113	2021-02-04		B.1.1.7
```

The required argument for the `import` command are listed as follows;

1. `--fasta` a fasta file containing genome sequences to be imported. A compressed file of fasta is also valid as an input (e.g., `--fasta sample.fasta.gz` or `sample.fasta.xz`).

2. `--tsv` a tab-delimited file containing sample properties to be imported.

3. `--cache` a directory for caching data.

4. `--cols` define column names for sample properties.

So, example
```sh
sonar import --db test.db --fasta valid.fasta --tsv day.tsv --threads 10 --cache tmp_cache  --cols sample=IMS_ID
```
As you can see, we defined `--cols sample=IMS_ID`, in which `IMS_ID` is the column ID that linked the sample name between the fasta file and meta-info file, and `sample` is the reserved word used to link data between tables in the database.

> TIP 🕯️: user might don't need to create an `ID` property because we use the `sample` keyword as the ID to link data in our database schema and also used in the query command, which you will see in the next section.

> TIP 🕯️: use `--threads` to increase the performance.

To update meta information when we add a new property, we can use the same `import` command, but this time, in the `--tsv` tag, we provide a new meta or old file, for example:
```sh
sonar import --db test.db --tsv meta.passed.tsv --threads 200 --cache tmp_cache --cols sample=IMS_ID

```
> NOTE 🤨: please make sure the `--cols sample=IMS_ID` is correctly referenced. If you have a different column name, please change it according to the meta-info file (for example, `--cols sample=IMS_NEW_ID`)

### 2.4 Query genome sequences based on profiles (match)

Genomic profiles can be defined to align genomes. For this purpose, the variants related to the complete genome of the SARS-CoV-2 isolate Wuhan-Hu-1 (NC_045512.2) must be expressed as follows:

| type       | nucleotide level                                                  | amino acid level              |
|-----------|-------------------------------------------------------------------|-------------------------------|
| SNP       | ref_nuc _followed by_ ref_pos _followed by_ alt_nuc (e.g. A3451T) | protein_symbol:ref_aa _followed by_ ref_pos _followed by_ alt_aa (e.g. S:N501Y) |
| deletion  | del:first_NT_deleted-last_NT_deleted (e.g. del:11288-11296)                        | protein_symbol:del:first_AA_deleted-last_AA_deleted (e.g. ORF1ab:del:3001-3004) |
| insertion | ref_nuc _followed by_ ref_pos _followed by_ alt_nucs (e.g. A3451TGAT) | protein_symbol:ref_aa _followed by_ ref_pos _followed by_ alt_aas (e.g. N:A34AK)  |

The positions refer to the reference (first nucleotide in the genome is position 1). Using the option `--profile`, multiple variant definitions can be combined into a nucleotide, amino acid or mixed profile, which means that matching genomes must have all those variations in common. In contrast, alternative variations can be defined by multiple `--profile` options. As an example, `--profile S:N501Y S:E484K` matches genomes sharing the _Nelly_ **AND** _Erik_ variation while `--profile S:N501Y --profile S:E484K` matches to genomes that share either the _Nelly_ **OR** _Erik_ variation **OR** both. Accordingly, using the option **^** profiles can be defined that have not to be present in the matched genomes.


There are additional options to adjust the matching.

| option             | description                                                            |
|--------------------|------------------------------------------------------------------------|
| --count            | count matching genomes only                                            |
| --format {csv,tsv,vcf}| output format (default: tsv) |


> TIP 🕯️: use `sonar match -h ` to see all available arguments.

Example;
```sh
sonar match --profile S:E484K --LINEAGE B.1.1.7 --db test.db

# matching B.1.1.7 genomes in DB 'test.db' that share an additional "Erik" mutation
sonar match --profile S:E484K --LINEAGE B.1.1.7 --db test.db

# --count to count the result
sonar match --profile S:E484K --LINEAGE B.1.1.7 --count --db test.db

# matching genomes in DB 'test.db' sharing the "Nelly" mutation
# and that were sampled on first of January 2020
sonar match --profile S:N501Y  --DATE 2020-01-01 --db test.db

# matching genomes with specific IDs
sonar match --sample ID-001 ID-001 ID-002 --db test.db
```

We use `^` as a **"NOT"** operator. We put it before any conditional statement to negate, exclude or filter the result.
```sh
# matching genomes in DB 'mydb' sharing the "Nelly" and the "Erik" mutation but not
# belonging to the B.1.1.7 lineage
sonar match -profile S:N501Y S:E484K --LINEAGE ^B.1.1.7 --db test.db

```

More example; `--profile` match
```sh
# AA profile OR  NT profile case
sonar match --profile S:del:K418I --profile T418A  --db test.db
# AA profile AND NT profile case
sonar match --profile S:del:K418I del:3677-3677  --db test.db
# exact Match X or N , we use small x for AA and small n for NT
sonar match --profile S:K418x --db test.db
# this will match S:K418X

sonar match --profile A2145n --db test.db
# this will match A2145N

# speacial case, we can combine exact match and any match in alternate postion.
sonar match  --profile A2145nN --db test.db
# this will look in ('NG', 'NB', 'NT', 'NM', 'NS', 'NV', 'NA', 'NH',
# 'ND', 'NY', 'NR', 'NW', 'NK', 'NN', 'NC')

```

More example; property match
```sh
# query with integer type
# by default we use = operator
sonar match  --AGE 25 --db test.db
# however, if we want to query with comparison operators (e.g., >, !=, <, >=, <=)
# , just add " " (double quote) around values.
sonar match  --AGE ">25" --db test.db
sonar match  --AGE ">=25" "<=30" --db test.db # AND Combination: >=25 AND <=30
sonar match  --AGE "!=60" --db test.db

# Range query matches
sonar match  --DEMIS_ID_PC  10641:10658  --db test.db
# 10641, 10642, 10643, .... 10658

# Date
# Sample were sampled in 2020
sonar match  --DATE 2020-01-01:2020-12-31 --db test.db
```
> TIP 🕯️: Don't forget `sonar list-prop --db test.db` to see more details

**Export to CSV/TSV/VCF file**

covSonar can return results in different formats: `--format ["csv", "tsv", "vcf"]`

```sh
# example command
sonar match --profile S:N501Y S:E484K --LINEAGE ^B.1.1.7 --db test.db --format csv -o out.csv

# in vcf format
sonar match -i S:N501Y S:E484K --lineage Q.1 --db test.db --format vcf -o out.vcf

# In case we have a list of ID and it is stored in a file, so we can use --sample-file
# tag to load and query according to the listed ID; example of --sample-file
sonar match --sample-file accessions.txt --db test.db --format vcf -o out.vcf
```

> NOTE 📌: accessions.txt has to contain one ID per line.

<u>Parent-Child relationship</u>

> ⚠️ WARNING: **This function currently works on SARS-CoV-2 only ❗**

First, we have to run `update-lineage-info` command to download the latest version of lineages from https://github.com/cov-lineages/pango-designation/ and install it in the database

```sh
# example command
sonar update-lineage-info
```

We want to search all sublineages with a given lineage, covSonar offers `--with-sublineage PROP_COLUMN` (PROP_COLUMN  means the property name that we added to our database).

In this example; we use `LINEAGE` property to store lineage information.
```sh
sonar match --profile S:E484K --LINEAGE B.1.1.7 --with-sublineage LINEAGE --count --db test.db --debug
```
This query will return results including 'B.1.1.7', 'Q.4', 'Q.5', 'Q.3', 'Q.6', 'Q.1', 'Q.7', 'Q.2', 'Q.8' lineages.

<u>Wildcard search</u>

covSonar also supports wildcard query (symbol `%`) to handle more complexity of the query. The operator is used in a **lineage query** to search for a specified pattern in lineage. It can be used in combinations of including and excluding command, for example;

```sh
# edit command at --lineage tag
sonar match -i S:N501Y S:E484K --LINEAGE ^B.1.1% AY.4% --db test.db --format csv -o out.csv
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

> NOTE 📌: We can also use wildcard with sublineage search. The following example shows when want to exclude all lineage B.1.617.X (X can be 1,2,3 ...) and  we also enable sublineage search, so this result in all lineage and sublineages from B.1.617.X are also filter out.
    ```
    sonar match --LINEAGE ^B.1.617% --with-sublineage LINEAGE --db test.db --count --debug
    ```
    result in
    ```
    NOT IN ('B.1.617.2', 'AY.94', 'AY.43.4', 'AY.39.1', 'AY.122.2', 'AY.122.6', 'AY.9.2.2', 'AY.37', 'AY.112.2', 'AY.40', 'AY.26', 'AY.5.6', 'AY.84', 'AY.28', 'AY.25.2', 'AY.91.1', 'B.1.617.3', 'AY.33', 'AY.95', 'AY.125', 'B.1.617.1', ..... etc.)
    ```

### 2.5 Show infos about the used sonar system and database (info)

Detailed infos about the used sonar system (e.g. version, reference,  number of imported genomes, unique sequences, available metadata).

```sh
# Show infos about the used sonar system and database 'test.db'
sonar info --db test.db
```

### 2.6 Restore genome sequences from the database (restore)
Genome sequences can be restored from the database based on their accessions.
The restored sequences are combined with their original FASTA header and  shown on the screen. The screen output can be redirected to a file easily by using `>`.

```sh
# Restore genome sequences linked to accessions 'mygenome1' and 'mygenome2' from the
# database 'test.db' and write these to a fasta file named 'restored.fasta'
sonar restore --sample mygenome1 mygenome2 --db test.db > restored.fasta
# as before, but consider all accessions from 'accessions.txt' (the file has to
# contain one accession per line)
sonar restore --sample-file accessions.txt --db test.db > restored.fasta
```

### 2.7 Database management (db-upgrade, optimize)
Sometimes you might need the `optimize` command to clean the [problems](https://www.sqlite.org/lang_vacuum.html) from database operation (e.g., unused data block or storage overhead ).
```sh
sonar optimize  --db test.db
```

When the newest version of covSonar use an old database version, covSonar will return  the following error;
```bash
Compatibility error: the given database is not compatible with this version of sonar (Current database version: XXX; Supported database version: XXX)
Please run 'sonar  db-upgrade' to upgrade database
```

We provide our database upgrade assistant to solve the problem.
```bash
# RUN
sonar db-upgrade --db test.db

# Output
Warning: Backup db file before upgrading, Press Enter to continue...

## after pressing the Enter key
Current version: 3  Upgrade to: 4
Perform the Upgrade: file: test.db
Database now version: 4
Success: Database upgrade was successfully completed

```
⚠️ WARNING: Backup the db file before upgrade.

### 2.7 Delete sample (delete)

```sh
sonar delete --db test.db --sample ID_1 ID_2 ID_3
```

## How to contribute 🏗️

covSonar has been very carefully programmed and tested, but is still in an early stage of development. You can contribute to this project by
* reporting problems 🐛
* writing feature requests to the [issue section](https://github.com/rki-mf1/covsonar/issues) 📝
* start hacking -> [read contribution guide]( https://github.com/rki-mf1/covsonar/blob/dev/cov2_mpire/CONTRIBUTING.md)👨‍💻

Please let us know that you plan to contribute before do any coding.

Your feedback is very welcome 👨‍🔧!

With love,

covSonar Team

---------------------------------

## Contact

For business inquiries or professional support requests 🍺 please contact [Dr. Stephan Fuchs](https://www.rki.de/SharedDocs/Personen/Mitarbeiter/F/Fuchs_Stephan.html)
