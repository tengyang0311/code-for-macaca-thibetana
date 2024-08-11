sample=$1
  
date
bwa mem -t 17 -k 32 -M -R "@RG\tID:$sample\tLB:$sample\tSM:$sample" sortlenthofgenome.fa "$sample"_R1.fq.gz "$sample"_R2.fq.gz |samtools view -bS -t sortlenthofgenome.fa.fai -> $sample.bam
date

perl stat_bwa_mem_v5.txx.pl -bam $sample.bam -out ../bam |samtools view -bS - > $sample.filted.bam

date
samtools sort -@16 -m 4000000000 $sample.filted.bam -o $sample.sort.bam
date
samtools rmdup $sample.sort.bam $sample.rmdup.bam
date
samtools index $sample.rmdup.bam
date
date
mkdir $sample
date
perl depth_v2.txx.pl -l 2842224279 $sample.rmdup.bam $sample
echo $sample is done

gatk --java-options "-Xmx20g -XX:ParallelGCThreads=2" \
    GenotypeGVCFs \
    -R sortlenthofgenome.fa \
    -V combined.1.g.vcf \
    -O raw_variants.1.vcf

