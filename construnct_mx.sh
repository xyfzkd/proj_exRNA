path0=/home/xieyufeng/proj_exRNA

## calculate raw counts/rpkm/rpm
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
l=1;
for k in rRNA miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp;
do echo $k;
samtools view -h $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.bam > $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.sorted.sam
let "l++";
done;
l=12;
samtools view -h $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.bam > $path0/output/02.mapping/$l.hg38other/bam/$j.hg38other.sorted.sam
done;

mkdir $path0/output/04.counts
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j;
mkdir $path0/output/04.counts/$j
# HTSeq
#htseq-count -m intersection-strict --idattr=ID --type=miRNA_primary_transcript $path0/02.mapping/$j/miRNA/$j.miRNA.rsem.clean.sorted.sam /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gff > $path0/04.counts/$j/$j.miRNA.htseq.ct
# featureCounts
#featureCounts -t miRNA_primary_transcript -g ID -a /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gff -o $path0/04.counts/$j/$j.miRNA.featureCounts.counts $path0/02.mapping/$j/miRNA/$j.miRNA.rsem.clean.sorted.bam
l=1;
for RNA_type in rRNA miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp;
do echo $RNA_type;
# featurecounts
featureCounts  -t exon -g transcript_id -s 1 -a $path0/anno/gtf/$RNA_type.* -o $path0/output/04.counts/$j/$j.$RNA_type.featureCounts.counts $path0/output/02.mapping/$l.no_$RNA_type/rsem_bam/$j.$RNA_type.rsem.clean.bam
let "l++";
done
# calculate rpm/cpm
#analyzeRepeats.pl /BioII/lulab_b/shared/genomes/human_hg38/gtf/miRNA.gtf hg38 -count exons -d $path0/03.tags/$j/merged/ -gid -norm 1e7 > $path0/04.counts/$j/$j.miRNA.homer.rpm
done;



path0=/home/xieyufeng/proj_exRNA
## merge featurecounts expression matrix for each RNA type
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;

cut -f 7 $path0/output/04.counts/${j}/${j}.${k}.featureCounts.counts | sed -e "1i ${j}" | sed 's/-/_/g' > $path0/output/tmp/${j}.${k}.featureCounts.counts.tmp
sed -i '1,3d' $path0/output/tmp/${j}.${k}.featureCounts.counts.tmp
done;
done;

for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;

cut -f 1 $path0/output/04.counts/Sample_N13/Sample_N13.${k}.featureCounts.counts > $path0/output/tmp/$k.featureCounts.counts.header
sed -i '1,2d' $path0/output/tmp/$k.featureCounts.counts.header
done;
mkdir $path0/output/05.matrix
for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;
paste $path0/output/tmp/$k.featureCounts.counts.header $path0/output/tmp/*.${k}.featureCounts.counts.tmp > $path0/output/05.matrix/proj_exRNA.featureCounts.counts.$k.mx
#sed -i "1s/......\t//1" $path0/output/05.matrix/proj_exRNA.featureCounts.counts.$k.mx
done;

## merge all RNA-type
for i in `ls $path0/raw_data | egrep -v 'zip|html'`;
do j=${i%.*};
echo $j
for k in  miRNA piRNA Y_RNA snRNA snoRNA srpRNA tRNA lncRNA mRNA tucp;
do echo $k;
cut -f 7 $path0/output/04.counts/${j}/${j}.${k}.featureCounts.counts | sed -e "1i ${j}" | sed 's/-/_/g'> $path0/output/tmp/${j}.${k}.featureCounts.counts.tmp
sed -i '1,3d' $path0/output/tmp/${j}.${k}.featureCounts.counts.tmp
done;
cat $path0/output/tmp/${j}.miRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.piRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.Y_RNA.featureCounts.counts.tmp $path0/output/tmp/${j}.snRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.snoRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.srpRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.tRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.lncRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.mRNA.featureCounts.counts.tmp $path0/output/tmp/${j}.tucp.featureCounts.counts.tmp | sed -e "1i ${j}" | sed 's/-/_/g' > $path0/output/tmp/${j}.merged.featureCounts.counts.tmp
cat $path0/output/tmp/miRNA.featureCounts.counts.header $path0/output/tmp/piRNA.featureCounts.counts.header $path0/output/tmp/Y_RNA.featureCounts.counts.header $path0/output/tmp/snRNA.featureCounts.counts.header $path0/output/tmp/snoRNA.featureCounts.counts.header $path0/output/tmp/srpRNA.featureCounts.counts.header $path0/output/tmp/tRNA.featureCounts.counts.header $path0/output/tmp/lncRNA.featureCounts.counts.header $path0/output/tmp/mRNA.featureCounts.counts.header $path0/output/tmp/tucp.featureCounts.counts.header | grep -v "geneID" | sed -e "1i geneID" > $path0/output/tmp/merged.featureCounts.counts.header
done;

paste $path0/output/tmp/merged.featureCounts.counts.header $path0/output/tmp/*.merged.featureCounts.counts.tmp > $path0/output/05.matrix/proj_exRNA.featureCounts.counts.merged.mx

#sed -i "1s/......\t//1" $path0/output/05.matrix/proj_exRNA.featureCounts.counts.merged.mx
#sed -i '2,3d' $path0/output/05.matrix/proj_exRNA.featureCounts.counts.merged.mx
#sed -e "s/^/${k}_/g" | sed -e "1i geneID" $path0/output/05.matrix/proj_exRNA.featureCounts.counts.merged.mx

