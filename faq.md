# Frequently Asked Questions

## What are the differences between `covSonar` and `covSonar2`?

In comparing `covSonar` and `covSonar2`, several differences become evident, particularly regarding the alignment and representation of deletions in genomic profiles.

**Q: Why do some deletions appear shifted when comparing genomic profiles of `covSonar` and `covSonar2`?**

**A:** The observed shifts in deletion positions result from a deliberate design decision in `covSonar2`. In contrast to `covSonar`, where deletions were right-aligned, `covSonar2` adopts a left-aligned approach for displaying deletions.
While there isn't a universal standard for the alignment direction of deletions in DNA sequences, the shift to left-aligned gaps offers advantages such as preserving gene reading frames and ensuring a more consistent representation of deletions.
This shift can lead to the repositioning of deletions by one or more bases due to changes in alignment direction, as illustrated below:

**Example:**
```plaintext
Exemplary reference sequence:   ATGCCATGATTAGGAATTCTGA
Right-aligned (covSonar):       ATGCCATGATT-------CTGA
Left-aligned (covSonar2):       ATGCCATGA-------TTCTGA
```

In situations where a deletion is adjactend to a stretch of "N" characters, indicating ambiguous or unknown bases, the deletion itself might undergo additional shifts, potentially spanning multiple bases. This behavior is demonstrated in the following examples:

**Example:**
```plaintext
Exemplary reference sequence:   ATGCCATGATTAGGAATTCTGA
Right-aligned (covSonar):       ATGCCATNNNNNN---TTCTGA
Left-aligned (covSonar2):       ATGCCAT---NNNNNNTTCTGA
```

**Q: Why has the formatting of deletions changed between `covSonar` and `covSonar2`?**

**A:** The way deletions are formatted has undergone a significant transformation from `del:1021:9` in `covSonar` to the more intuitive `del:1021-1030` in `covSonar2` where the deletion length is replaced by the end coordinate.
This shift represents a move towards a clearer and more informative style that explicitly includes both the start and end positions to meet standard DNA sequence representation conventions and to minimize any confusion.

**Q: Why are insertions following a SNP notated differently in `covSonar2` compared to `covSonar`?**

**A:** In `covSonar2`, the notation for insertions following a single nucleotide polymorphism (SNP) has been refined to accurately reflect their position within the genomic coordinate system, providing improved clarity and specificity. Unlike `covSonar`, which previously used a condensed mutation notation, `covSonar2` employs a more explicit approach.

In the case of insertions the base directly before the insertion, known as the anchor base, is used to reference the insertion's location. Let's consider the following example:

**Example:**
```plaintext
Reference Sequence:      ATGCCATGATT-----AGGAATTCTGA
Sequence with Insertion: ATGCCATGATTGTTTAAGGAATTCTGA
```

Both `covSonar` and `covSonar2` would notate this insertion as `T11TGTTTA`. However, it's important to note that the `T` at position 11 is not part of the insertion itself; it serves as an anchor base to indicate the insertion's position.

In rare cases, the anchor base itself might be mutated, leading to scenarios like:

**Example:**
```plaintex
Reference Sequence:                    ATGCCATGATT-----AGGAATTCTGA
Sequence with SNP preceding Insertion: ATGCCATGATAGTTTAAGGAATTCTGA
```

In such cases, `covSonar2` addresses the complexity more explicitly. While covSonar used to express this as a single mutation (`T11AGTTTA`), which might not accurately represent both independent mutation events, `covSonar2` dissects the mutations into separate notations. This approach references both independent mutational events as `T11A` and `T11TGTTTA` in the genomic profile. This distinction  also ensures accurate matching and interpretation of the mutations independently from each other.

**Q: Why are some deletions very close to one of the sequence termini missing in the genomic profile of `covSonar2`?**

**A:** The omission of certain deletions located near both the start and end of the sequence in the `covSonar2` genomic profile is a result of transitioning from global to semi-global alignment. Semi-global alignment prioritizes aligning subsequences without introducing penalties for gaps at the beginning or end of the reference sequence. This alignment approach aligns well with amplicon sequencing limitations, as sequence ends often have no coverage due to experimental constraints. Deletions located near the sequence's beginning or end, when globally aligned, can undergo shifts towards the sequence termini in a semi-global alignment (see example below). To ensure alignment accuracy, `covSonar2`, by default, excludes these deletions from the profile. Consider the following scenario:

**Example:**

```plaintex
Reference Sequence:             ATGCCATGATTAGGTGAAATTCTGA
Globally aligned Sequence:      ATGCCATGATTAGG--------TGA
Semi-globally aligned Sequence: ATGCCATGATTAGGTGA--------
```

**Q: Why does `covSonar2`, by default, display inserted N (stretches) in the genomic profile?**

**A:** In `covSonar2`, the default inclusion of inserted N base or N stretch in the genomic profile stems from the recognition that these stretches hold informative value. Unlike uniformative mutations (N SNPs), inserted N (stretches) provide information about the inserted sequence's specific length. This rationale sets the stage for a distinction from `covSonar`'s approach, where only N SNPs are hidden from the default profile due to their uninformative nature.

## Adapting covSonar 1 commands to covSonar 2

|Description                                                                               |covSonar 1 command                                                                                                      |covSonar 2 command                                                                                                   |
|---                                                                                       |---                                                                                                                     |---                                                                                                                  |
|Create a new database using default SC2 reference                                         |This is done automatically when you `./sonar.py add` sequences to a non-existent database                               |`sonar setup --db test.db`                                                                                           |
|Add sequences to a database, forcing existing sequnces to be updated. Use 8 parallel jobs.|`./sonar.py add --noprogress --force --db test.db --file seqs.fasta --cpus 8`                                           |`sonar import --db test.db --fasta seqs.fasta --threads 8 --no-progress` (updating of existing seqs is the default)  |
|Add/update sequence metadata from a TSV format table                                      |`./sonar.py update --fields accession=covv_accession_id software=covv_assembly_method [...] --tsv metadata.tsv --db $db`|`sonar import --db test.db --tsv metadata.tsv --cols accession=covv_accession_id software=covv_assembly_method [...]`|
|Do the above two steps at once (add sequnces and metadata simultaneously)                 |(not possible)                                                                                                          |`sonar import --db test.db --fasta seqs.fasta --tsv metadata.tsv --cols prop=col prop2=col2 [...]`                   |
|Update a single metadata column (e.g. for updating lineage information)                   |`./sonar.py update --db test.db --fields accession=taxon lineage=lineage --csv pangolin.csv`                            |`sonar import --db test.db --csv pangolin.csv --cols SAMPLE_NAME=taxon LINEAGE=lineage`                              |
|Update Pangolin lineage parent-child relationships                                        |`./sonar.py update-lineage-info --db test.db`                                                                           |`sonar update-lineage-info --db test.db`                                                                             |
|Optimize sqlite database                                                                  |`./sonar.py optimize --db test.db`                                                                                      |`sonar optimize --db test.db`                                                                                        |
|Calculate database statistics                                                             |`./sonar.py info --db ${db}`                                                                                            |`sonar info --db test.db`                                                                                            |
|Remove a set of sequences by ID                                                           |`./sonar.py remove --db test.db  --file sequences-to-delete.ids`                                                        |`sonar delete --db test.db --sample-file sequences-to-delete.ids`                                                    |
|Identify sequences with a given lineage or any of its sublineages                         |`./sonar.py match --lineage B.1.1.7 --with-sublineage --db test.db --tsv`                                               |`sonar match --db test.db --LINEAGE B.1.1.7 --with-sublineage LINEAGE`                                               |
|Output all sequences, including purely ambiguous variations in profiles (e.g. N stretches)|`./sonar.py match --ambig --tsv --db test.db`                                                                           |`sonar match --db test.db --showNX`                                                                                  |
|Output all sequences that contain frameshift mutations                                    |`./sonar.py match --only_frameshifts --db test.db`                                                                      |`sonar match --db test.db --frameshifts-only`                                                                        |
|Output all sequences that do not belong to a given lineage, in a date range               |`./sonar.py match --lineage ^B.1.1.529 --lineage ^BA.% --date 2021-12-01:2021-12-31  --db test.db --tsv`                |`sonar match --db test.db --LINEAGE ^B.1.1.529 --LINEAGE ^BA% --date 2021-12-01:2021-12-31`                          |

## Troubleshooting

**Q:  `covSonar2` crashes with: `sqlite3.OperationalError: attempt to write a readonly database`**

**A**: You're encountering the SQLite error attempt to write a readonly database, even though your file permissions appear to be correct.

SQLite uses locking to prevent concurrent access to the database by multiple processes. Typically, these locks are released once a transaction is completed. However, under certain circumstances, a lock might not be properly released. This can occur for a variety of reasons, including a process crashing before it can release its lock, or multiple threads within a process trying to acquire a lock concurrently. If another process or thread tries to access the database while it's locked, SQLite will return an error.

In `covSonar2`, the database establishes a write connection only during the data import process. To minimize issues such as database locks, concurrent write connections are explicitly disallowed. However, if a system crash occurs during the data import, while the integrity of the data is preserved, the lock on the database may persist. In such cases, a persistent journal file would likely be present in the same directory as the database, contributing to the `sqlite3.OperationalError: attempt to write a readonly database` error.

SQLite creates a journal file as part of its ACID compliance to ensure data consistency, even in the event of an interruption during a write operation. If a journal file isn't cleaned up properly after the write operation, it can cause persistent database lock and other issues with subsequent attempts to write to the database.

Here's what you can do:

* **Check for the Journal File:** Look for a file in the same directory as your database file with the same name and a `-journal` extension. For example, if your database is `gisaid.db`, the journal file will be `gisaid.db-journal`.

* **Safe Removal of the Journal File:** If the `-journal` file exists, it could be locking your database. You might want to delete this file, but proceed with caution. Only remove it if you're certain that no other processes are accessing the database, as deleting the journal file while a process is still using it can lead to data corruption. Ensure that you have a backup of your database before deleting the journal file.

* **Try executing the sonar `optimize` command:** This command rebuilds the database, repacking it into a minimal amount of disk space. Note that this operation can take a while on large databases and locks the database file for the duration.

Please exercise caution when handling database and journal files. Always maintain backups and ensure your actions won't interrupt an active process or lead to data loss.
