#!/bin/bash

module load cellranger
cd /data/Project/FASTQ/

cellranger count --id=132 \
        --transcriptome=/data/Genomes/Human/10XCellRanger_GRCh38_2020-A/GRCh38 \
        --fastqs=/data/Project/FASTQ/ \
        --sample=HNCL7DSX2_19144504,HTLYHDSX2_19144504,HNCL7DSX2_19144526,HTLYHDSX2_19144526,HNCL7DSX2_19144546,HTLYHDSX2_19144546,HNCL7DSX2_19144566,HTLYHDSX2_19144566,H7LTGDSX3_19144504,H7LTGDSX3_19144526,H7LTGDSX3_19144546,H7LTGDSX3_19144566 \
        --chemistry=SC3Pv3 \
        --r1-length=28 \
        --r2-length=91
