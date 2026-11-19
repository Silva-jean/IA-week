This release includes inbred genotypic data from 2014 to 2023 of the G x E Project within the Genomes to Fields (G2F) initiative.

Variant calls for the 2014-2023 the G x E Project materials were obtained using the Practical Haplotype Graph (PHG, Bradbury et al. 2022). The PHG aligns sequencing reads against a database populated with genome assemblies’ sequences, and imputes genotypes based on the haplotypes stored in the database. 

The Maize 2.1 PHG database’s haplotypes were created from 86 genome assemblies aligned using anchorwave (Song et al. 2022) to the B73 v5 assembly (Hufford et al. 2021). 35 assemblies were downloaded from MaizeGDB (Woodhouse et al.2021), one from (Yang et al.2019), 5 Stiff-stalk assemblies from (Bornowski et al.2021), and 46 were private assemblies shared by collaborators. Nodes for the graph were defined by dividing the B73 genome into genic and intergenic reference ranges based on the Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3. Reference range end points were selected in regions where at least 23 of the NAM genomes have 10 or more conserved base pairs. Genic regions with 0.0001 or less divergence and intergenic regions with 0.001 or less divergence were collapsed to a single consensus haplotype.

Reads for each inbred were aligned against the PHG Pangenome, and haplotypes matching each read were identified. The imputed path through the graph was then used to identify variants for the positions of the high density 600 k SNP genotyping array (Unterseer et al. 2014). Array positions were uplifted to v5 coordinate equivalents through Crossmap. Positions that could not be uplifted, missing in more than 21 assemblies, or monomorphic for all 86 assemblies were filtered out, for a final dataset with 437,214 variant positions. 

2193 inbreds represented in the 2014-2023 G2F hybrids were sequenced with different technologies:
• 2014-2017 germplasm was sequenced with GBS (genotyping by sequencing) following the Elshire et al 2011 protocol, and ApeKI as restriction enzyme. Raw data is available under BioProject IDs PRJNA385022, SRP021921, SRP009896, and SRP004282.1.
• 2018-2019 germplasm was skim sequenced at ~5x coverage with Nextera libraries in a HiSeq X Ten or NovaSeq 6000. Raw data available under BioProject ID PRJNA530187.
• 2020-2021 germplasm was sequenced using Exome capture, with a few samples sequenced using GBS with PstI-MspI as restriction enzymes. Raw data available under BioProject IDs PRJNA781986 (GBS), and PRJNA781987 (Exome capture).
• 2022-2023 germplasm was sequenced using GBS with PstI-MspI as restriction enzymes with the NovaSeq 6000 sequencer. Raw data available under BioProject ID PRJNA956551.

Please note that there are errors in this dataset, and it is minimally filtered. Some inbreds are not present in the data set because they are either publicly unavailable or did not pass quality control filters.

For general information about the Genomes to Fields (G2F) initiative, the G x E Project, and questions regarding the released data set, please contact g2f@wisc.edu.



Inventory:

/readme.txt (this file).

/key_inbreds_G2F_2014-2023.txt
	   Contains Cultivar (genotype name), Dataset (technology type), SourceName, Bioproject (when available), BioSample (when available), Alternative name (for inbred that used more than one name, for example, DK78010 and 78010), and Comments.

/inbreds_G2F_2014-2023_437k.vcf
	Inbreds genotypic data, 2193 unique genotypes and 437,214 variant sites.  



References:
1. Bradbury et al. The Practical Haplotype Graph, a platform for storing and using pangenomes for imputation. Bioinformatics. 2022 38:3698–702 https://doi.org/10.1093/bioinformatics/btac410
2. Woodhouse et al. A pan-genomic approach to genome databases using maize as a model system. BMC Plant Biol 2021, 385 https://doi.org/10.1186/s12870-021-03173-5
3. Yang et al. Genome assembly of a tropical maize inbred line provides insights into structural variation and crop improvement. Nat Genet 2019 51 https://doi.org/10.1038/s41588-019-0427-6
4. Bornowski et al. Genomic variation within the maize stiff-stalk heterotic germplasm pool. Plant Genome. 2021, 14: e20114 https://doi.org/10.1002/tpg2.20114
5. Hufford et al. De novo assembly, annotation, and comparative analysis of 26 diverse maize genomes. Science. 2021 373:655-662 https://doi.org/10.1126/science.abg5289
6. Unterseer et al. A powerful tool for genome analysis in maize: development and evaluation of the high density 600 k SNP genotyping array. BMC Genomics. 2014 15: 823 https://doi.org/10.1186/1471-2164-15-823
7. Bradbury et al. TASSEL: software for association mapping of complex traits in diverse samples, Bioinformatics.2007 23: 2633-2635 https://doi.org/10.1093/bioinformatics/btm308