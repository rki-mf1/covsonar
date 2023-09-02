# Frequently Asked Questions

## What are the differences between covSonar and covSonar2?

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

In situations where a deletion is immediately adjactend to a stretch of "N" characters, indicating ambiguous or unknown bases, the deletion itself might undergo additional shifts, potentially spanning multiple bases. This behavior is demonstrated in the following examples:

**Example:**
```plaintext
Exemplary reference sequence:   ATGCCATGATTAGGAATTCTGA
Right-aligned (covSonar):       ATGCCATNNNNNN---TTCTGA
Left-aligned (covSonar2):       ATGCCAT---NNNNNNTTCTGA
```

**Q:** Why has the formatting of deletions changed between `covSonar` and `covSonar2`?

**A:** The way deletions are formatted has undergone a significant transformation from `del:1021:9` in `covSonar` to the more intuitive `del:1021-1030` in `covSonar2` where the deletion length is replaced by the end coordinate. 
This shift represents a move towards a clearer and more informative style that explicitly includes both the start and end positions to meet standard DNA sequence representation conventions and to minimize any confusion.

**Q: Why are insertions following a SNP notated differently in covSonar2 compared to covSonar?**

**A:** In covSonar2, the notation for insertions following a single nucleotide polymorphism (SNP) has been refined to accurately reflect their position within the genomic coordinate system, providing improved clarity and specificity. Unlike covSonar, which previously used a condensed mutation notation, covSonar2 employs a more explicit approach.

**Example:**

In the case of insertions the base directly before the insertion, known as the anchor base, is used to reference the insertion's location. Let's consider the following example:

**Example:**
```plaintext
Reference Sequence:      ATGCCATGATT-----AGGAATTCTGA
Sequence with Insertion: ATGCCATGATTGTTTAAGGAATTCTGA
```

Both covSonar and covSonar2 would notate this insertion as `T11TGTTTA`. However, it's important to note that the `T` at position 11 is not part of the insertion itself; it serves as an anchor base to indicate the insertion's position.

In rare cases, the anchor base itself might be mutated, leading to scenarios like:

**Example:**
```plaintex
Reference Sequence:      ATGCCATGATT-----AGGAATTCTGA
Sequence with Insertion: ATGCCATGATAGTTTAAGGAATTCTGA
```

In such cases, covSonar2 addresses the complexity more explicitly. While covSonar used to express this as a single mutation (`T11AGTTTA`), which might not accurately represent both independent mutation events, covSonar2 dissects the mutations into separate notations. This approach references both independent mutational events as `T11A` and `T11TGTTTA` in the genomic profile. This distinction  also ensures accurate matching and interpretation of the mutations independently from each other.
