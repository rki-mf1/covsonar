BEGIN;

-- 0. Upgrade scheme version
PRAGMA user_version = 4;

-- 1. Update the genome table and add index
ALTER TABLE genome ADD submission_date TEXT;
CREATE INDEX idx_meta_submission_date ON genome (submission_date);
-- 2. Remove old VIEW
DROP VIEW dna_view;
DROP VIEW essence;
DROP VIEW prot_view;
-- 3.Create new VIEW
CREATE VIEW IF NOT EXISTS essence
AS
SELECT
	accession,description,lab,source,collection,technology,platform,
	chemistry,material,ct,software,software_version,gisaid,ena,
	zip,date,submission_date,lineage,genome.seqhash,dna_profile,aa_profile,fs_profile
FROM
	genome
LEFT JOIN sequence USING (seqhash)
LEFT JOIN profile USING (seqhash);

CREATE VIEW IF NOT EXISTS dna_view
AS
SELECT
	accession,description,lab,source,collection,technology,platform,chemistry,
	material,ct,software,software_version,gisaid,ena,zip,date,submission_date,
	lineage,genome.seqhash,start,end,ref,alt
FROM
	genome
LEFT JOIN sequence USING (seqhash)
LEFT JOIN sequence2dna USING (seqhash)
LEFT JOIN dna USING (varid);

CREATE VIEW IF NOT EXISTS prot_view
AS
SELECT
	accession,description,lab,source,collection,technology,platform,chemistry,material,ct,software,
	software_version,gisaid,ena,zip,date,submission_date,lineage,
	genome.seqhash,protein,locus,start,end,ref,alt
FROM
	genome
LEFT JOIN sequence USING (seqhash)
LEFT JOIN sequence2prot USING (seqhash)
LEFT JOIN prot USING (varid);

COMMIT;