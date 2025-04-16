
-- Bash commands
  --sudo apt install mysql-server
  --sudo mysql -u root -p < iedb_public.sql 
  --sudo mysql -u root -p

-- Select iedb database
USE iedb_public;

-- Select id for species Homo sapiens
SELECT
  organism_id
FROM
  organism_names
WHERE
  LOWER(name_txt) = 'homo sapiens';

-- Select epitope id and object id
SELECT
  curated_epitope_id, e_object_id
FROM
  curated_epitope
INTO
  OUTFILE '/var/lib/mysql-files/epitope_object.tsv'
  FIELDS TERMINATED BY '\t'
  ENCLOSED BY ''
  LINES TERMINATED BY '\n';

-- Select epitope id, disease id and mhc type from mhc
SELECT
  curated_epitope_id,
  iv1_disease_id,
  h_mhc_types_present,
  mhc_allele_name
FROM
  mhc_elution
WHERE
  h_organism_id = '9606'
  AND iv1_disease_id IS NOT NULL
INTO
  OUTFILE '/var/lib/mysql-files/mhc_epitope.tsv'
  FIELDS TERMINATED BY '\t'
  ENCLOSED BY ''
  LINES TERMINATED BY '\n';

-- Select epitope id and disease id from tcell
SELECT
  curated_epitope_id,
  h_mhc_types_present,
  mhc_allele_name,
  iv1_disease_id,
  iv2_disease_id,
  adt_iv_disease_id
FROM
  tcell
WHERE
  h_organism_id = '9606'
  AND (
    iv1_disease_id IS NOT NULL
    OR iv2_disease_id IS NOT NULL
    OR adt_iv_disease_id IS NOT NULL
  )
INTO
  OUTFILE '/var/lib/mysql-files/t_cell.tsv'
  FIELDS TERMINATED BY '\t'
  ENCLOSED BY ''
  LINES TERMINATED BY '\n';

-- Select epitope id and disease id from bcell
SELECT
  curated_epitope_id,
  h_mhc_types_present,
  iv1_disease_id,
  iv2_disease_id,
  adt_iv_disease_id
FROM
  bcell
WHERE
  h_organism_id = '9606'
  AND (
    iv1_disease_id IS NOT NULL
    OR iv2_disease_id IS NOT NULL
    OR adt_iv_disease_id IS NOT NULL
  )
INTO
  OUTFILE '/var/lib/mysql-files/b_cell.tsv'
  FIELDS TERMINATED BY '\t'
  ENCLOSED BY ''
  LINES TERMINATED BY '\n';

-- Select disease id from disease
SELECT
  disease_id, disease_name
FROM
  disease
WHERE
  LOWER(disease_name) REGEXP 'cancer|carcinoma|leukemia|lymphoma|sarcoma|tumor|neoplasm|malignancy|metastasis|melanoma|myeloma|adenoma|lipoma|glioma|meningioma|blastoma|
grade|cytoma|atoma'
INTO
  OUTFILE '/var/lib/mysql-files/disease.tsv'
  FIELDS TERMINATED BY '\t'
  ENCLOSED BY ''
  LINES TERMINATED BY '\n';

-- Select object id and sequence from human
SELECT
  object_id, mol1_seq
FROM
  object
WHERE
  organism_id = '9606' OR organism2_id = '9606'
INTO
  OUTFILE '/var/lib/mysql-files/object_seq.tsv'
  FIELDS TERMINATED BY '\t'
  ENCLOSED BY ''
  LINES TERMINATED BY '\n';

-- Bash commands
  -- sudo mv /var/lib/mysql-files/*tsv ../results/db_tables
