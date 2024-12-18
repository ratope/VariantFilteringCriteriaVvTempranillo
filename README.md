# Filtering criteria of the variants called in 53 samples of *Vitis vinifera* cv. Tempranillo

This is the list of the criteria applied to the raw list of variants obtained from the variant calling with `freebayes`. The list of criteria depicted here and the number of variants produced in each step is also provided in the file **Stats_vcf_versions_v71_v71c_20221005.xlsx**.  

The following filtering steps were done:

- **Only SNP and INDELS were considered**:  
  `bcftools-1.11 view` and `grep -P TYPE=(snp(,snp)?|ins|del)` were used.

- **New tags were calculated and included in the VCF info**:  
  A Perl script provided new VCF tags to filter by:  
  - `N11`: Number of homozygous 1/1 samples  
  - `N01`: Number of heterozygous 0/1 samples  
  - `AFS`: Alt allele(s) frequency  
  - `AODP`: Quotient AO/DP  
  - `QD`: Quotient QUAL/DP  
  - `NMS`: Number of samples with missing genotype  
  - `NHM`: Number of samples with homozygous genotype  
  - `NHT`: Number of samples with heterozygous genotype  

- **Filtering according to alternative counts (AO), depth (DP), quality (QD) and homozygosity/heterozygosity**:  
  A Perl script removed variants if any sample triggers any of the following filters at that position:  
  - `HIGH_AFS_TRIALLELIC`: AFS of the second allele > 0.025.  
  - `MISSING_SAMPLE_S`: Remove if there is a/some missing sample(s).  
  - `LOW_AO_AODP_HM`: For any homozygous sample, remove if AO < 9 or AODP < 0.975.  
  - `LOW_AO_AODP_HT`: For any heterozygous sample, remove if AO < 5 or AODP < 0.25.  
  - `LOW_QD`: QD < 2 in any sample.  
  - `ALL_HOMOZ`: Remove if all samples are homozygous.  
  - `ALL_HETEROZ`: Remove if all samples are heterozygous.  

- **Removing variants present in the Tempranillo RJ51 parental**:  
  RJ51 sample FASTQ files provided by the user (`RJ51_root_PCRfree.R1.fastq.gz` and `RJ51_root_PCRfree.R2.fastq.gz`) were processed using:  
  - BWA-MEM aligner  
  - Optical removal with Picard  
  - Variant calling with `freebayes` (`-p 2`, genome reference: `benedicto_v1.2_scaffolds.fasta`)  
  A raw (unfiltered) list of 11,082,405 parental variants was obtained.  
  Using `bcftools isec -C <prefiltered_multi_vcf> <parental_RJ51_variants_vcf> -w1`, parental variants were subtracted from the multisample VCF.

- **Removing variants present in the cv. Albillo haplotype**:  
  Using `bcftools isec -C <prefiltered_multi_vcf> <albillo_haplotype_variants_vcf> -w1`, the variants of the haplotype from cv. Albillo (user-provided file `albillo_v1.2_scaffolds.chrs.nuc.benedicto_v1.2_scaffolds.chrs.m.i88.l100.2syri.renameIDs.SNP_INS-DEL.vcf.gz`) were subtracted from the multisample VCF.

- **Strand Bias filtering**:  
  1. `bcftools query -f '%CHROM:%POS\n'` was used to obtain the list of variant positions.  
  2. `bcftools mpileup -a FORMAT/AD, FORMAT/ADF, FORMAT/ADR, FORMAT/SP` added relevant tags to the file (multisample VCF):  
     - `AD`: Allelic depths  
     - `ADF`: Allelic depths on the forward strand  
     - `ADR`: Allelic depths on the reverse strand  
     - `SP`: Phred-scaled strand bias P-value  
  3. Perl scripts removed variants with:  
     - `SP > 201`  
     - `ADF < 1` or `ADR < 1`

- **Remove SNPs located within 100 base pairs around an INDEL**:  
  `bcftools filter --SnpGap 100:indel`

- **Filter clusters of indels separated by 100 or fewer bp allowing only one to pass**:  
  `bcftools filter --IndelGap 100`

- **Remove SNPs located too close (50 bp) one from another**:  
  `vcftools --gzvcf --stdout --thin 50 --remove-filtered-all --recode --recode-INFO-all`

> *The higher the SP of the variant, the more the probability of being strand bias.*

- **Remove variants present in homopolymeric regions**:  
  1. Using Kevin Blighe's AWK solution ([Biostars thread](https://www.biostars.org/p/379454/)), a "masked bed file" was generated containing regions of the reference genome that are homopolymeric (6 bp or more) plus adjacent positions.  
  2. Variants located within any homopolymeric region were removed using:  
     `bedtools intersect -a <prefiltered_multi_vcf> -b <masked_bed_file> -wa -v -header`

- **Remove variants located in regions with high mappability**:  
  User provided BED files of high mappability regions:  
  - **v71 (GTEp0.5)**: `P-value >= 0.5`  
    File: `benedicto_v1.2_scaffolds.genmap.GTEp0.5.sorted.merged.bed`.  
  - **v71c (GTp0.5)**: `P-value > 0.5`  
    File: `benedicto_v1.2_scaffolds.genmap.GTp0.5.sorted.merged.bed`.  

  Variants were filtered using:  
  `bedtools intersect -a <prefiltered_multi_vcf> -b <HighMappabiltiy_BED> -wa -header`.  
  Statistics were obtained using `RTG-Tools (vcfstats)` and `bcftools stats`.

---

## Samples Subsettings and Second Filtering

Using `bcftools view` and `grep`, SNPs were extracted from the multisample VCF file.  

**Subsets** of samples (30 total) were extracted, and a new set of filters was applied via Perl scripts to generate **26 final VCF files** of filtered SNPs:

- `ALL_HOMOZ`: Remove SNP if all samples are homozygous.  
- `ALL_HETEROZ`: Remove SNP if all samples are heterozygous.  
- `MISSING_SAMPLE_S`: Remove SNP if there is a/some missing sample(s).  
- `LOW_QD`: Remove SNP if QD < 2 in any sample.  
- `LOW_DP`: Remove SNP if DP < 10 in any sample.  
- `HIGH_DP`: Remove SNP if DP > 100 in any sample.  
- `LOW_MQMx`: Remove SNP if MQM < 10 or MQMR < 103.  
- `HIGH_AODP_00`: Remove SNP if, having also 1/1- or 0/1-genotyped samples, there are any 0/0 samples with `AO/DP` quotient >= 0.025.
