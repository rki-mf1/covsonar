<img src="logo.png"  width="134" height="134" align="right"><br><br><br><br><br>

# covSonar2

covSonar is a database-driven system for storing, profiling and querying genomic sequences.

What's new in covSonar v.2
* New design
    * Improved workflows
    * Improved database scheme
    * Improved performance
    * Improved alignment (now: optimal semi-global)
* Exciting new features
    * Support of multiple pathogens
    * Support of segmenten genomes
    * Support of user-defined meta information
* Extended query lanquage 
    * more powerful and complex queries
    * greedy and ungreedy deletion handling

[![Tests Status](https://github.com/rki-mf1/covsonar/actions/workflows/tests.yml/badge.svg)](https://github.com/rki-mf1/covsonar/actions/workflows/tests.yml) [![Pre-release/Test Version](https://github.com/rki-mf1/covsonar/actions/workflows/test-pypi.yml/badge.svg)](https://github.com/rki-mf1/covsonar/actions/workflows/test-pypi.yml) [![Stable Version](https://github.com/rki-mf1/covsonar/actions/workflows/release.yml/badge.svg)](https://github.com/rki-mf1/covsonar/actions/workflows/release.yml)

[![Github tag](https://badgen.net/github/release/rki-mf1/covsonar/Strapdown.js)](https://github.com/rki-mf1/covsonar/release/)


## 1. Installation

covSonar2 now can be easily installed by using `pip`, `conda` or `mamba`
 
### Stable version 🔖
Proceed as follows to install covSonar2 via pip

```sh
pip install covsonar
```
or via conda
```sh
conda install -c bioconda covsonar
```
or via mamba
```sh
mamba install -c bioconda covsonar
```
and to verify your installation
```sh
sonar --version
```

### Developer version 🚧

Developer versions are for testing, nightly builds, feature implementation or bug fixing.

For more information, please check at https://test.pypi.org/project/covsonar/

```sh
# example command
pip install -i https://test.pypi.org/simple/ covsonar==2.0.0a1.dev1657097995
```

## 2. Usage

In covSonar2, the table below shows the different sub-commands that can be called.

| subcommand          | purpose                                                                 |
|---------------------|-------------------------------------------------------------------------|
| setup               | set up a new database.                                	                |
| import              | import genome sequences and sample information to the database          |
| list-prop           | list available sample properties for a given  database                  |
| add-prop            | add a new sample property to the database                               |
| delete-prop         | delete an existing sample property from the database                    |
| match               | match stored genomes based on mutation profiles and/or sample properties |
| restore             | restore genome sequence(s) from the stored mutation profiles             |
| info                |  show software and database info.                                       |
| optimize            | optimizes the database                                                  |
| db-upgrade          | upgrade a database to the latest version                                |
| update-lineage-info | download latest lineage information                                     |
| direct-query        | use directly SQLite queries                                             |

Each tool provides a help page that can be accessed with the `-h` option.

```sh
# show general help page on available sub-commands
sonar -h
# show help page for a specific sub-command
sonar import -h
```

### 2.1 Building a database (`setup`)

First, we have to create a new database.
```sh
sonar setup --db test.db
```

To create a new database instance with default sample properties use:
```sh
sonar setup --db test.db --auto-create
```
By default, MN908947.3 (SARS-CoV-2) is used as a reference. If we want to set up a database for a different pathogen, we can add `--gbk` followed by a Genbank file containing both the genome sequence and annotation of the designated reference. Example:
```sh
sonar setup --db test.db --auto-create --gbk Ebola.gb
```
> TIP 🕯️:  [how to download genbank file](https://ncbiinsights.ncbi.nlm.nih.gov/2017/05/08/genome-data-download-made-easy/)

> EXPERT TIP 🕯️: We can use [DB Browser](https://sqlitebrowser.org/) to visualise or manipulate the database file directly.


### 2.2 Property management (`list-prop`, `add-prop` and `delete-prop`)

In covSonar2, users can add properties to store custom meta information in the database. Based on the information to store, the data type (and optionally query type) has to be defined when adding a new property. The table below lists the supported data types with their default query types.

data type      | description                                  | default query type  |
---------------|----------------------------------------------|---------------------|
integer        | property stores integers                     | numeric             |
float           | property stores decimal numbers              | float                |
text           | property stores text                         | text                |
date           | property stores dates in the form YYYY-MM-DD | date                |
zip            | property stores zip codes                    | zip                 |


To add properties, we can use the `add-prop` command to add meta information into the database.

The required arguments are listed below when we use `add-prop` command
* `--name`, name of sample property
* `--descr`, description of the new property
* `--dtype`, data type of the new property (see table above)

```sh
# for example
sonar add-prop --db test.db --name LINEAGE --dtype text --descr "store Lineage"
#
sonar add-prop --db test.db --name AGE --dtype integer --descr "age information"
#
sonar add-prop --db test.db --name DATE_DRAW --dtype date --descr "sampling date"
```
> TIP 🕯️: `sonar add-prop -h ` to see all arguments available.

> NOTE 📌: Property names have to start with letters and have to consist of letters and digits only.

> NOTE 📌: All property names are stored in upper-case letters.

> NOTE 📌: SAMPLE is reserved keyword that cannot be used as property names.

To view added properties, use the `list-prop` command.
```sh
sonar list-prop --db test.db
```
The `delete-prop` command is used to delete an property from the database that is not longer needed any more.

```sh
sonar delete-prop --db test.db  --name SEQ_REASON
```

The program will ask for confirmation of the action, type YES and enter to confirm.
```
Do you really want to delete this property? [YES/no]: YES
```

### 2.3 Adding genomes and meta information to the database (`import`)

This example shows how to import genome sequences along with meta information to a given database.

In this example, sequences to import are stored in a fasta formatted file named `valid.fasta`:
```
>IMS-00113
CCAACCAACTTTCGATCTCTTG
>IMS-00114
CTAACCAACTTTCGATCTCTAG
```
Respective meta information file is stored in a tab-delimited text file named `meta.tsv`:
```
IMS_ID		SAMPLING_DATE	LINEAGE
IMS-00113	2021-02-04		B.1.1.7
IMS-00114	2021-04-21		BA.5
```

The required argument for the `import` command are listed as follows.

* `--fasta`, a fasta file containing genome sequences to be imported. 
* `--tsv`, a tab-delimited file containing sample properties to be imported.
* `--cols`, define column names for sample properties if a file containing meta information is provided

```sh
sonar import --db test.db --fasta valid.fasta --tsv day.tsv --threads 10  --cols sample=IMS_ID
```
> NOTE 📌: When providing meta information, the corresponding file must contain a column containing the sequence IDs to be linked. This column must be defined as sample name column by using `--cols sample=IMS_ID` in our example.
 
> TIP 🕯️: Use `--threads` to define threads to use which increases the performance.

To update meta information only, we can use the same `import` command without referencing any fasta file, for example:
```sh
sonar import --db test.db --tsv meta.passed.tsv --cols sample=IMS_ID
```

### 2.4 Query genome sequences based on mutation profiles and/or meta information (`match`)

Stored genomes can be queried based on their mutation profiles and sample properties. 
We follow the scientific notation for mutations. In the following table *REF* stands for the reference allele, *ALT* for the variant allele, *pos* for the affected reference position (1-based) or *start* and *end* for the affected reference range (1-based). *PROT* stands for symbol of the respective gene product.

| mutation type     | nucleotide level                         | amino acid level                                     |
|-------------------|------------------------------------------|------------------------------------------------------|
| SNP               | *REFposALT* (e.g. A3451T)                | *PROT*:*REFposALT* (e.g. S:N501Y)                    |
| 1bp Deletion      | *del*:*pos* (e.g. del:11288)             | *PROT*:del:*pos* (e.g. ORF1ab:del:3001)              |
| multi-bp Deletion | del:*start*-*end* (e.g. del:11288-11300) | *PROT*:del:*start*-*end* (e.g. ORF1ab:del:3001-3004) |
| Insertion         | *REFposALT* (e.g. A3451TGAT)             | *PROT*:*REFposALT* (e.g. N:A34AK)                    |

> NOTE 📌: Any mutation can be negated by prefixing the notation with an **rooftop character (^)**. For example, *^A1345ATGC* only matches target genomes that do not carry this insertion.

> NOTE 📌: Deletion notations are greedy by default. For example, the notation *del:12* corresponds to any deletion that affects reference position 12. To fix a position, an **equal sign (=)** can be prepended. For instance, the notation *del:=12* matches only 1bp deletions at reference position 12. Likewise, the notation *del:1500-=1505* matches deletions that span at least reference positions 1500 to 1505 and definitely terminate at position 1505.

> NOTE 📌: Any mutation can be negated by prefixing the notation with an **caret (^)**. For example, *^A1345ATGC* only matches target genomes that do not carry this insertion.

> NOTE 📌: By default, the IUPAC characters N (nucleotide level) and X (amino acid level) are interpreted as any polymorphism at the respective reference position. To explicitly search for N or X at specific ones, they use n or x (e.g. *T306n*). 

> NOTE 📌: To include non-informative variant alleles (N resp. X) in the returned genome profiles, use the `--showNX` option. 

Mutation profiles consist of one or multiple mutations and can be queried by `--profile`. Multiple mutations following this argument mean that all notations must be satisfied by the target genome (**AND** logic). Multiple `--profile` arguments can be used to define alternate mutation profiles that are searched simulaatively (**OR** logic). The following examples will illustrate this.

```sh
# querying genomes carrying "Erik" (S:E484K) mutation 
sonar match --profile S:E484K --db test.db

# querying genomes carrying "Erik" (S:E484K) AND "Nelly" (S:N501Y) mutation 
sonar match --profile S:E484K S:N501Y --db test.db

# querying genomes carrying "Erik" (S:E484K) OR "Nelly" (S:N501Y) mutation 
sonar match --profile S:E484K --profile S:N501Y --db test.db

# querying genomes carrying "Erik" (S:E484K) BUT NOT "Nelly" (S:N501Y) mutation 
sonar match --profile S:E484K ^S:N501Y --db test.db

# querying genomes where reference positions 99 to 102 are deleted
sonar match --profile del:99-102 --db test.db

# querying genomes that carry a deletion exactly from reference position 99 to exactly 102
sonar match --profile del:=99-=102 --db test.db
```

Target genomes can also (additionally) be searched based on sample names (sequence IDs). For this, the respective sample names are listed after the `--sample` argument. Alternatively, a text file containing one sample name per line can be specified after `--sample-file`.

Additionally, target genomes can be queried based on linked meta information (properties). For this, the name of the respective property is used as argument (e.g. `--LINEAGE`). Multiple alternate values can follow. Based on the data types of the respective properties, different operators are allowed:

| operator | description                                           |valid data types     |
|----------|-------------------------------------------------------|---------------------|
| >        | larger than (e.g. >1)                                 | integer, float, date |
| <        | larger than (e.g. <1)                                 | integer, float       |
| >=       | larger than or equal to (e.g. >=1)                    | integer, float, date |
| <=       | smaller than or equal to (e.g. <=1)                   | integer, float, date |
| !=       | different than (e.g. !=1)                              | integer, float       |
| :        | range, *from*:*to* (e.g. 2021-01-01:2021-12-31)       | integer, float, date |
| ^        | not the same as (e.g. ^B.1.1.7)                       | text                |
| %        | wildcard standing for any character(s) (e.g. %human%) | text                | 

> NOTE 📌: The zip data type is interpreted hierarchically. This means that the search value 106 matches all zip codes starting with 106.

> NOTE 📌: If a property refers to the lineage classifications of SARS-CoV-2 based on Pango nomenclature (https://cov-lineages.org/resources/pangolin.html), the corresponding sub-lineages can be included when querying specific lineages by using `--with-sublineages` followed by the name of the property to be applied to. Make sure that you use the latest pango lineage definitions by running `sonar update-lineage-info` command before.

```sh
# query genomes collected 2021 and later
sonar match --SAMPLING >=2021-01-01 --db test.db
# query genomes in carrying "Erik" (S:E484K) mutation and NOT classified as B.1.1.7 
sonar match --profile S:E484K --LINEAGE ^B.1.1.7 --db test.db
# query genomes in carrying "Erik" (S:E484K) mutation and classified as B.1.1.7 or any sub-lineage of B.1.1.7 
sonar match --profile S:E484K --LINEAGE B.1.1.7 --db test.db --with-sublineage LINEAGE
# query genomes with any SNP at referebnce genome position 11022
sonar match --profile A11022N --db test.db
# query genomes with N at reference genome position 11022
sonar match --profile A11022n --db test.db
# query genomes resulting from samples with Ct values of 30 and higher
sonar match --CT >=30 --db test.db
# query genomes from patients in the age class 31 to 40
sonar match --AGE 31:40 --db test.db
# query genomes with accession seq01 or seq02
sonar match --sample seq01 seq02
```

By default, mutation profiles of the target genomes (excluding N or X alleles) plus all linked meta information are output as a tab-separated file. The output can be customized by the following options.

| option                   | description                                                                             |
|--------------------------|-----------------------------------------------------------------------------------------|
| --count                  | count target genomes only                                                               |
| --format *{csv,tsv,vcf}* | select output format between csv (comma separated), tsv (tab delimited, default) or vcf |
| --out-cols *(col_names)* | define columns to show (works only for tsv or csv output format)                         |
| -o *(file_name)*          | write the results to the specified file instead of displaying them on the screen          |


> TIP 🕯️: use `sonar match -h` to see all available query and output options.

### 2.5 Native database queries (`direct-query`)

Advanced users familiar with the covSonar database scheme, can use the SQLite query syntax to directly query the database as demonstrated by the following example.

```sh
sonar direct-query --sql "SELECT COUNT(*) FROM sequences" --db test.db
```

> Note 📌: For added security, `direct-query` allows read-only access to the database, not write access.

### 2.6 Show infos about the used sonar system and database (`info`)

The `info` sub-command lists details about the used sonar software and database (e.g. version,  number of imported genome sequences and unique sequences).

```sh
sonar info --db test.db
```

### 2.7 Restore genome sequences from the database (`restore`)
Genome sequences can be retrieved from the database based on their accessions. Optionally, the sequences can be displayed in quasi-aligned form to the standard reference, with insertions indicated by lowercase letters and deletions indicated by minus (-) signs. The argument "-o" can be used to forward the output to a file.

```sh
sonar restore --sample mygenome1 mygenome2 --db test.db -o restored.fasta
```

### 2.8 Database management (`db-upgrade`, `optimize`)
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
> WARNING ⚠️: Backup the db file before upgrade. 

> WARNING ⚠️: Due to the change in the alignment algorithm used, it is not recommended to upgrade databases from covsonar v.1.x to v.2.x.

### 2.7 Delete sample (`delete`)

As shown by the follwing example, genome sequences can be deleted from the databse based on their sample name (sequence ID).

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
