#!/bin/bash -       
#title           :bioinformatic_pipeline
#description     :This script contains all command used for the HyRAD eDNA Frogs paper
#date            :01.06.2024

#==============================================================================


#Demultiplexing and cleaning

awk '{print $3}' list_barcodes | sort | uniq > list_lib

while read a
	do
	grep "$a" list_barcodes | awk '{print ">"$1"#^"$2}' | tr "#" "\n" > temp_barcodes.fasta
	cutadapt -e 0.17 --minimum-length 30 -q 10 --no-indels -g file:temp_barcodes.fasta -o demux-{name}_R1_.fastq.gz -p demux-{name}_R2_.fastq.gz "$a"*_R1_001.fastq.gz "$a"*_R2_001.fastq.gz 
	done < list_lib

for i in `ls demux-*_R1_.fastq.gz`
	do
	cutadapt -q 10 -m 30 -a AGATCGGAAGAGC -o clean_"$i" "$i"
	cutadapt -u -5 -q 10 -m 30 -o clean2_"$i" clean_"$i"
	done

for i in `ls demux-*_R2_.fastq.gz`
	do
	cutadapt -u 5 -q 10 -m 30 -a AGATCGGAAGAGC -o clean_"$i" "$i"
	cutadapt -u -6 -q 10 -m 30 -o clean2_"$i" clean_"$i"
	done

for i in `ls clean2_demux-*_R1_.fastq.gz`
		do
		name=`echo $i | sed -e 's/_R1_.fastq.gz//g'`
		gzip -d "$i"
		gzip -d "$name"_R2_.fastq.gz
		fastq_pair "$name"_R1_.fastq "$name"_R2_.fastq
		done

#Mapping on Reference genome

for i in `ls clean2_demux-*_R1_.fastq.paired.fq`
	do
	sample=`echo $i | sed -e 's/clean2_demux-//g' -e 's/_R1_.fastq.paired.fq//g'`
	bwa mem GCF_905171775.1_aRanTem1.1_genomic_cut.fna clean2_demux-"$sample"_R1_.fastq.paired.fq clean2_demux-"$sample"_R2_.fastq.paired.fq > "$sample"_on_ref.sam
	done

for i in `ls *_on_ref.sam`
		do
		sample=`echo $i | sed -e 's/_on_ref.sam//g'`
		samtools sort "$i" -o temp_1_sorted.bam
		samtools view -bF 4 -q 20 temp_1_sorted.bam > "$sample"_keep.bam
		java -jar -mx128G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar CreateSequenceDictionary R=GCF_905171775.1_aRanTem1.1_genomic_cut.fna O=GCF_905171775.1_aRanTem1.1_genomic_cut.dict
		java -jar -mx512G /home/jeremy/local/envpicard/share/picard-2.20.2-0/picard.jar AddOrReplaceReadGroups I="$sample"_keep.bam O=temp_1_sorted_keep_rg.bam ID=["$sample"] RGLB=[id] PL=[pl] PU=[pu] SM=["$sample"]
		samtools index temp_1_sorted_keep_rg.bam
		GenomeAnalysisTK -T RealignerTargetCreator -I temp_1_sorted_keep_rg.bam -R GCF_905171775.1_aRanTem1.1_genomic_cut.fna -o temp.intervals
		java -jar -Xmx512G /home/jeremy/local/envgatk/opt/gatk-3.8/GenomeAnalysisTK.jar -T IndelRealigner -I temp_1_sorted_keep_rg.bam -R GCF_905171775.1_aRanTem1.1_genomic_cut.fna -targetIntervals temp.intervals -o "$sample"_sorted_keep_rg_realign.bam
		rm temp*
		done

#Example merging
samtools merge -o 1M_sorted_keep_rg_realign.bam SPY_201_298_sorted_keep_rg_realign.bam SPY_211_197_sorted_keep_rg_realign.bam SPY_211_203_sorted_keep_rg_realign.bam 

#Coverage

for i in `ls *_rg_realign.bam`
	do
	genomeCoverageBed -ibam "$i" -d | awk '{if ($3 >= 1) print $0}' > "$i"_bedcov_sup1
	grep -c "." "$i"_bedcov_sup1 > count_cov1
	awk '{if ($3 >= 20) print $0}' "$i"_bedcov_sup1 | grep -c "." > count_cov20
	done


#Venn diagrams
awk '{print $1}' list_samples | sort | uniq > list_condi
while read a 
	do
	awk '{if ($1 == "'$a'") print $2}' list_samples > temp_list_s
	S1=`cat temp_list_s | tr '\n' '\t' | awk '{print $1}'`
	S2=`cat temp_list_s | tr '\n' '\t' | awk '{print $2}'`
	S3=`cat temp_list_s | tr '\n' '\t' | awk '{print $3}'`
	echo $a
	Rscript ./Venn_plot_simple.R "$S1" "$S2" "$S3"
	mv venn_diagramm.png "$a"_venn_diagramm.png
	done < list_condi

#Rscript
library(VennDiagram)

args <- commandArgs(trailingOnly = TRUE)

A<-read.table(args[1],h=F,sep="\t")
B<-read.table(args[2],h=F,sep="\t")
C<-read.table(args[3],h=F,sep="\t")

venn.diagram(x = list(A$V1, B$V1, C$V1), category.names = c("Set 1" , "Set 2 " , "Set 3"), filename = 'venn_diagramm.png', output=TRUE)



#SNP calling

samples=""

for data in `ls *_realign.bam`
	do
	samples=$samples" -I "$data
	done

gatk HaplotypeCaller -R GCF_905171775.1_aRanTem1.1_genomic_cut.fna -O gatk4_calling.vcf $samples


#SNP filters and Structure format

##eDNA
vcftools --vcf gatk4_calling.vcf --keep list_eDNA --maf 0.2 --max-missing 0.6 --recode --out eDNA_maf02_miss06_gatk4_calling

##Individuals
vcftools --vcf gatk4_calling.vcf --keep list_indiv --maf 0.15 --max-missing 0.6 --recode --out Indiv_maf015_miss06_gatk4_calling

##Conversion in .str
vcf2structure_gn.sh eDNA_maf02_miss06_gatk4_calling


#PCoA in R

library(adegenet)
library(ape)
library(ggplot2)
library(dartR)
library(vcfR)
library(dplyr)
library(corrplot)

## Import vcf
vcf.pond <- read.vcfR("16721snps_4edna.recode.vcf")

## Transform vcf to genind
genlight.pond <- vcfR2genlight(vcf.pond)
### Add pop info
pop.pond <- read.table("edna-pond-pooled.txt", header=FALSE)
colnames(pop.pond) <- c("SAMPLE","TYPE","POP")
genlight.pond@pop <- as.factor(pop.pond$POP)

### Make euclidean distance
source("pcoa_function.R")
pcoa.pond.graph <- pcoa.function(genlight.pond, pop.pond)
pcoa.pond.graph

## Import vcf
vcf.aqua <- read.vcfR("207975snps_4aqua.recode.vcf")

## Transform vcf to genind
genlight.aqua <- vcfR2genlight(vcf.aqua)
### Add pop info
pop.aqua <- read.table("edna-aqua-pooled.txt", header=TRUE)
colnames(pop.aqua) <- c("SAMPLE","TYPE","POP")
genlight.aqua@pop <- as.factor(pop.aqua$POP)

### Make euclidean distance
pcoa.aqua.graph <- pcoa.function(genlight.aqua,pop.aqua)
pcoa.aqua.graph

## Import vcf
vcf.ind <- read.vcfR("250739snps_76ind.recode.vcf")
## Transform vcf to genind
genlight.I <- vcfR2genlight(vcf.ind)

### Add pop info
pop.ind <- read.table("dna-76samples.txt", header=FALSE)
colnames(pop.ind) <- c("SAMPLE","TYPE","POP")
genlight.I@pop <- as.factor(pop.ind$POP)

### Make PCoA
pcoa.ind.graph <- pcoa.function(genlight.I,pop.ind)
pcoa.ind.graph

### Save the PCoA
cowplot::plot_grid(pcoa.pond.graph,pcoa.aqua.graph, pcoa.ind.graph, nrow=1, labels=c("A","B","C"))

### Remove the population 4
genlight.I <- gl.compliance.check(genlight.I)
genlight.I.3pop <- gl.drop.pop(genlight.I, pop.list = "4", mono.rm = TRUE)

### Make PCoA on the 3 sampling location
pcoa.ind.graph.3pop <- pcoa.function(genlight.I.3pop,pop.ind[1:57,])
pcoa.ind.graph.3pop


#Venn
library(VennDiagram)
library(RColorBrewer)

### Check the common SNPs based on list of list extract from the vcf
set1 <- read.table("16721snps.txt")
set2 <- read.table("207975snps.txt")
set3 <- read.table("250739snps.txt")

pos1 <- as.vector(paste(set1$V1, set1$V2, sep="_"))
pos2 <- as.vector(paste(set2$V1, set2$V2, sep="_"))
pos3 <- as.vector(paste(set3$V1, set3$V2, sep="_"))

### Set colors that match the figure
myCol <- c("green3", "blue4", "grey")

# Chart
venn.diagram.hyrad <- venn.diagram(
  x = list(pos1, pos2, pos3),
  category.names = c("Ponds" , "Aquarium " , "Individuals"),
  filename = '#14_venn_diagramm.png',
  output=TRUE,
  
  # Output features
  imagetype="png" ,
  height = 480 , 
  width = 480 , 
  resolution = 300,
  compression = "lzw",
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol,
  
  # Numbers
  cex = .4,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 0.4,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)
overlap <- calculate.overlap(x=list("pond"=pos1,"aquarium" = pos2,"ind"=pos3))
snps.all.samples <- as.data.frame(overlap[1])
write.table(snps.all.samples, "1556snps.txt", quote=FALSE, row.names=FALSE)

#### STEP 03: SUBSET ONLY SNPS IN COMMON ####

### Set up the list of SNPs to keep for all samples
loc.to.keep <- as.vector(snps.all.samples$a5)
genlight_1556snps_pond <- gl.keep.loc(genlight.pond, loc.list = loc.to.keep)
genlight_1556snps_aqua <- gl.keep.loc(genlight.aqua, loc.list = loc.to.keep)
genlight_1556snps_ind <- gl.keep.loc(genlight.I, loc.list = loc.to.keep)

### Make the PCoA on each genlight
pcoa.pond.subset <- pcoa.function(genlight_1556snps_pond, pop.pond)
pcoa.aqua.subset <- pcoa.function(genlight_1556snps_aqua, pop.aqua)
pcoa.ind.subset <- pcoa.function(genlight_1556snps_ind,pop.ind)

### Save the PCoA
cowplot::plot_grid(pcoa.pond.subset,pcoa.aqua.subset, pcoa.ind.subset, nrow=1, labels=c("A","B","C"))

ggsave("PCoA-all-subset.pdf", width=10, height=5)




