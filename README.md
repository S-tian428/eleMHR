# eleMHR
eleMHR is a computational method for discovering mutation hotspot regions in genomic elements (CDS, 3UTR, 5UTR, promoters, splice site, etc.) to identify cancer driver genes genome-wide.eleMHR quantifies the enrichment of somatic mutations and detects mutation hotspot regions in genomic elements based on mutation set enrichment analysis. 

Two major files are required by eleMHR: a set of genomic regions of interest element, and a set of somatic mutations in a tumor cohort, they are all in tab-separated text files. The element files we provide have two reference genomes containing hg19 and hg38. Users can make their selection based on the reference genome of somatic mutations.

The match between somatic mutations and element region is completed by R code match_mutations_elements.R and generates element_mutations files. The eleMHR code has both R and Java versions for users.
