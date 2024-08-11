
samtools=/mnt/qjw/software/samtools-1.3.1/samtools
bcftools=/mnt/qjw/software/yes/envs/atac/bin/bcftools
vcfutils=/mnt/qjw/software/yes/envs/atac/bin/vcfutils.pl

sample=$1



/mnt/qjw/software/samtools-1.3.1/samtools mpileup -C 50 -q 1 -m 2 -F 0.002 -uf /media/primates/ty/ref/sortlenthofgenome.fa /media/primates/ty/02bam/${sample}.rmdup.bam | /mnt/qjw/software/yes/envs/atac/bin/bcftools call -c - | /mnt/qjw/software/yes/envs/atac/bin/vcfutils.pl vcf2fq -d 10 -D 100 | gzip > ${sample}.fq.gz

/mnt/qjw/software/psmc-master/utils/fq2psmcfa -q 20 ${sample}.fq.gz > ./${sample}/${sample}.psmcfa

/mnt/qjw/software/psmc-master/psmc -N 25 -t 15 -r 5 -p "4+25*2+4+6" -o ./${sample}/${sample}.all.psmc ./${sample}/${sample}.psmcfa
perl /mnt/qjw/software/psmc-master/utils/psmc_plot.pl -g 11 -u 4.67e-09 -p -w 3 -m 100 -R ./${sample}/${sample}.all ./${sample}/${sample}.all.psmc


