#!/bin/sh

rsync -av -f"+ */" -f"- *" "/data/dataset/" "./"

for sample in /data/dataset/**/*
do
  parent=$(basename $(dirname "$sample"))
  fileName=$(basename "$sample")
  fileNameWithoutExt=$(basename "$sample" | cut -d. -f1)
  if [ ! -f "$parent"/"$fileNameWithoutExt"_htseq_count.txt ]; then
    samtools sort -n -O sam -T "$parent" /data/picardtools/output/"$parent"/"$fileNameWithoutExt"_sorted.bam | htseq-count -m intersection-nonempty -s yes -t exon -i gene_name - /data/reference-genome/gencode.v7.annotation/gencode.v7.annotation_goodContig.gtf > "$parent"/"$fileNameWithoutExt"_htseq_count.txt
  fi
done

perl /data/scripts/htseq/HTSeq_combine_output.pl -c -o Casp2KO.GFP-Cas9.GFP_HTSeq_Expression_Matrix.txt -r Casp2KO.GFP-Cas9.GFP_HTSeq_report.txt /data/scripts/htseq/Casp2-KO-GFP_Cas9-GFP.txt
mkdir /data/edgeR
mkdir /data/edgeR/output
Rscript /data/scripts/edgeR/r_script.R