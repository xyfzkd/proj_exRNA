## summary the results and statistics
## number of mapped reads for each RNA types
path0=<path to>/proj_exRNA;
for i in `ls $path0/data/raw_data`;
do j=${i%.*};
echo $j;
# total library size
libSizeN=`echo $(cat $path0/data/raw_data/${j}.fastq | wc -l)/4 | bc`
echo -e "$j\tpreprocess\tlibSizeN\t$libSizeN" >> $path0/stat/$j.readsN.stat.tsv
# clean reads
cleanN=`echo $(cat $path0/output/01.trim/cutadapt/$j.cutadapt.fastq | wc -l)/4 | bc`
echo -e "$j\tpreprocess\tcleanN\t$cleanN" >> $path0/stat/$j.readsN.stat.tsv
# rRNA mapped reads
rRNA_N=`samtools flagstat $path0/output/02.mapping/1.no_rRNA/sam/$j.rRNA.sam | awk 'NR==5{print $1}'`
echo -e "$j\tpreprocess\trRNA_N\t$rRNA_N" >> $path0/stat/$j.readsN.stat.tsv
# kept reads
keptN=`echo $(cat $path0/output/02.mapping/1.no_rRNA/fastq/$j.no_rRNA.fq | wc -l)/4 | bc`
echo -e "$j\tpreprocess\tkeptN\t$keptN" >> $path0/stat/$j.readsN.stat.tsv
# map to different RNA types
l=2;
for k in miRNA piRNA Y_RNA srpRNA tRNA snRNA snoRNA lncRNA mRNA tucp;
do echo $k;
RNAtypes_N=`samtools flagstat $path0/output/02.mapping/$l.no_$k/rsem_bam/$j.$k.rsem.clean.bam | awk 'NR==5{print $1}'`
echo -e "$j\tsequentialMap\t$k\t$RNAtypes_N" >> $path0/stat/$j.readsN.stat.tsv
let "l++";
done;
# map to hg38 other region
hg38other_N=`samtools flagstat $path0/output/02.mapping/12.no_hg38other/sam/$j.hg38other.sam | awk 'NR==5{print $1}'`
echo -e "$j\tmap2hg38other\thg38other\t$hg38other_N" >> $path0/stat/$j.readsN.stat.tsv
# non-human
nonHuman_N=`echo $(cat $path0/output/02.mapping/12.no_hg38other/fastq/$j.hg38other.unAligned.fastq | wc -l)/4 | bc`
echo -e "$j\tmap2hg38other\tnonHuman_N\t$nonHuman_N" >> $path0/stat/$j.readsN.stat.tsv
done;

cut -f 2,3 $path0/stat/$j.readsN.stat.tsv | sed -e "1i sample" > $path0/stat/readsN.stat.header
for i in `ls $path0/data/raw_data`;
do j=${i%.*};
echo $j;
cut -f 4 $path0/stat/$j.readsN.stat.tsv | sed -e "1i $j" > $path0/stat/$j.readsN.stat.foo
done;
paste $path0/stat/readsN.stat.header $path0/stat/*.readsN.stat.foo > $path0/stat/readsN.stat.tsv
