# VariantFilteringCriteriaVvTempranillo
List of the filtering criteria of the variants called in 53 samples of Vitis vinifera cv. Tempranillo

A list of the filtering criteria described in this and the following sections
and the number of variants produced in each step is provided in the file
Stats_vcf_versions_v71_v71c_20221005.xlsx.
The following filtering steps were done:
• Only SNP and INDELS were considered: bcftools-1.11 view and grep
-P TYPE=(snp(,snp)?|ins|del) were used.
• New tags were calculated and included in the VCF info: A Perl script
provided new VCF tags to filter by:
– N11: Number of homozygous 1/1 samples
– N01: Number of heterozygous 0/1 samples
– AFS: Alt allele(s) frequency
– AODP: Quotient AO/DP
– QD: Quotient QUAL/DP
– NMS: Number of samples with missing genotype
– NHM: Number of samples with homozygous genotype
– NHT: Number of samples with heterozygous genotype
• Fitering according to alternative counts (AO), depth (DP), quality
(QD) and homozygosity/heterozygosity: A Perl script removed variants
if any sample triggers any of the following filters at that position:
– HIGH_AFS_TRIALLELIC: AFS of the second allele > 0.025.
– MISSING_SAMPLE_S: Remove if there is a/some missing sample(
s).
– LOW_AO_AODP_HM: For any homozygous sample, remove if
AO < 9 or AODP < 0.975.
– LOW_AO_AODP_HT: For any heterozygous sample, remove if
AO < 5 or AODP < 0.25.
– LOW_QD: QD < 2 in any sample.
– ALL_HOMOZ: Remove if all samples are homozygous.
– ALL_HETEROZ: Remove if all samples are heterozygous.
• Removing variants present in the Tempranillo RJ51 parental: RJ51
sample FASTQ files provided by the user (files RJ51_root_PCRfree
.R1.fastq.gz and RJ51_root_PCRfree.R2.fastq.gz) were processed
using BWA-MEM aligner, optical removal with Picard then variant
calling with freebayes (-p 2 and genome reference: "benedicto_v1.2
_scaffolds.fasta"). A raw (unfiltered) list of 11,082,405 parental
variants was obtained.
By using bcftools isec -C <prefiltered_multi_vcf> <parental_RJ51
_variants_vcf> -w1, parental variants were substracted from the
multisample VCF.
• Removing variants present in the cv. Albillo haplotype: By using
bcftools isec -C <prefiltered_multi_vcf> <albillo_haplotype_variants
_vcf> -w1, the variants of the haplotype from cv. albillo (user-provided
file "albillo_v1.2_scaffolds.chrs.nuc.benedicto_v1.2_scaffolds.chrs
.m.i88.l100.2syri.renameIDs.SNP_INS-DEL.vcf.gz" were substracted
from the multisample VCF.
• Strand Bias filtering: 1) bcftools query -f ’%CHROM:%POS\n’ was
used to obtain the list of variant positions at this point of the filtering;
2) bcftools mpileup -a FORMAT/AD, FORMAT/ADF, FORMAT/ADR,
FORMAT/SP was used to add AD (Allelic depths), ADF (Allelic depths
on the forward strand), ADR (Allelic depths on the reverse strand)
and SP (Phred-scaled strand bias P-value) tags to the file (multisample
VCF) of variants. 3) Perl scripts remove those variants having
SP > 201 or ADF < 1 or ADR < 1.
• Remove SNPs located within 100 base pairs around an INDEL: bcftools
filter –SnpGap 100:indel.
• Filter clusters of indels separated by 100 or fewer bp allowing only
one to pass: bcftools filter –IndelGap 100.
• Remove SNPs located too close (50 bp) one from another: vcftools[12]
–gzvcf –stdout –thin 50 –remove-filtered-all –recode –recode-INFOall.
1The higher the SP of the variant, the more the probability of being an strand bias
6
RTP20220218_ICVV_SNP_DIGEVIDA (2/2)
• Remove variants present in homopolymeric regions: 1) By using the
Kevin Blighe’s AWK solution described in https://www.biostars.
org/p/379454/, we generated a "masked bed file" containing the
regions of the reference genome that are homopolymeric (sized
6 bp or more) plus one up- and downstream adjacent position; 2)
Variants located withing any homopolymeric region is removed: bedtools
intersect -a <prefiltered_multi_vcf> -b <masked_bed_file> -
wa -v -header.
• Remove variants located in regions with high mappability: User
provided several BED files containing regions of the Primitivo reference
genome having high mappability, at different levels of confidence.
Two of these mappability files were used to generate ‘two
parallel versions of results requested by the user:
– v71 (GTEp0.5): P-value >= 0.5, corresponding to mappability
file: benedicto_v1.2_scaffolds.genmap.GTEp0.5.sorted.merged
.bed.
– v71c (GTp0.5): P-value > 0.5, corresponding to mappability
file: benedicto_v1.2_scaffolds.genmap.GTp0.5.sorted.merged
.bed.
We filtered out those variants located in high-mappable regions of
each of the mappability files by using bedtools intersect -a <prefiltered_
multi_vcf> -b <HighMappabiltiy_BED> -wa -header.
We obtained statistics of the variants remaining after this filtering
process by using RTG-Tools (vcfstats)[13] and bcftools stats.

Samples subsettings and second filtering
By using bcftools view and grep we extracted the SNPs from the variants
multisample VCF file resulting in the previous section2
We extracted diverse subsets of samples from this SNP multisample
file.

Then, we applied a new set of filters on each of those 30 subsets by
using a Perl script to generate the final set of 26 VCF files of filtered
SNPs:
• ALL_HOMOZ: Remove SNP if all samples are homozygous.
• ALL_HETEROZ: Remove SNP if all samples are heterozygous.
• MISSING_SAMPLE_S: Remove SNP if there is a/some missing sample(
s).
• LOW_QD: Remove SNP is QD < 2 in any sample.
• LOW_DP: Remove SNP if DP < 10 in any sample.
• HIGH_DP: Remove SNP if DP > 100 in any sample.
• LOW_MQMx: Remove SNP if MQM < 10 or MQMR < 103.
• HIGH_AODP_00: Remove SNP if, having also 1/1- or 0/1-genotyped
samples, there are any 0/0 sample wiith AO/DP quotient >= 0.025.

References
