CREATE TABLE snpmap(
  name,
  chromosome,
  position 
); 

CREATE TABLE snps(
  snp VARCHAR(50),
  animal VARCHAR(50), 
  allele1 VARCHAR(1), 
  allele2 VARCHAR(1), 
  x FLOAT, 
  y FLOAT, 
  gcscore FLOAT 
);

CREATE INDEX snp_idx ON snps(snp);
CREATE INDEX animal_idx ON snps(animal);   

CREATE INDEX chromosome_idx ON snpmap(chromosome);
