/mnt/qjw/software/plink/plink --vcf ../07passvcf/raw.20scaffold.q20.biallel.vcf --make-bed --autosome-num 20 --out th
king -b thi.bed --kinship
king -b thi.bed --kinship --degree 1
/home/primates/miniconda3/bin/Rscript picture.R kinship.pdf
