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

#ddRAD loci reconstruction

for i in `ls *_realign.bam`
	do
	bedtools bamtobed -i "$i" > "$i".bed
	bedtools cluster -d 500 -i "$i".bed > "$i".clustered
	echo $i >> output_loci 
	cut -f 4 "$i".clustered | sort | uniq | wc -l >> output_loci 
	done


#SNP calling

samples=""

for data in `ls *_realign.bam`
	do
	samples=$samples" -I "$data
	done

gatk HaplotypeCaller -R GCF_905171775.1_aRanTem1.1_genomic_cut.fna -O gatk4_calling.vcf $samples


#SNP filters and Structure format

vcftools --vcf gatk4_calling.vcf --remove-indels --min-alleles 2 --max-alleles 2 --keep list_indiv --maf 0.05 --max-missing 0.8 --recode --out snp_bi_miss08_maf005_gatk4

##Conversion in .str
vcf2structure_gn.sh snp_bi_miss08_maf005_gatk4.recode.vcf


#PCA in R
data<-read.structure("snp_bi_miss08_maf005_gatk4.recode.str",n.ind=108,n.loc=17617,onerowperind=FALSE,col.lab=1,col.pop=NULL,col.others=NULL,row.marknames=NULL,NA.char=-9,ask=F) 
is.genind(data)
sansna<-scaleGen(data,scale=F, NA.method="mean")
pca1 <- dudi.pca(sansna,cent=F,scale=F,scannf=FALSE,nf=4)
barplot(pca1$eig[1:20],main="PCA eigenvalues", col=heat.colors(20))
s.class(pca1$li,xax=1,yax=2,pch=19)
(pca1$eig/sum(pca1$eig))*100

pca<-data.frame(pca1$l1)
ggplot(data,aes(x=PC1, y=PC2)) +
  geom_point(aes(color=pop, shape=type,size=1.1)) +
  scale_shape_manual(values = c(22, 15, 24, 17, 20)) +
  scale_color_manual(values = c("red","dodgerblue1","purple","orange"))+
  theme_classic()


#Poolseq analyses
samtools mpileup -B 1M_sorted_keep_rg_realign_rg.bam 2M_sorted_keep_rg_realign_rg.bam 3M_sorted_keep_rg_realign_rg.bam 4M_sorted_keep_rg_realign_rg.bam > 1M2M3M4M.mpileup
perl mpileup2sync.pl --input 1M2M3M4M.mpileup --output 1M2M3M4M.sync --fastq-type sanger --min-qual 1
perl fst-sliding.pl --input 1M2M3M4M.sync --output 1M2M3M4M.fst --suppress-noninformative --min-count 2 --min-coverage 20 --max-coverage 2000 --window-size 1 --step-size 1 --pool-size 40

ggplot(data, aes(x = ddRAD , y = eDNA_Pond)) +
  geom_point(size=3) +
  stat_smooth(data=data, color = "grey", alpha=0.1, method= "lm") +
  theme_classic() +
  xlim(c(0,0.05)) +
  ylim(c(0,0.05)) +
  xlab("FST individual samples") + 
  ylab("FST pond eDNA")
  theme(axis.text=element_text(size=12))
  
