#!/bin/bash
  
module load cellranger

cd /data/Project/FASTQ/

cellranger aggr --id=HV4_2 --csv=cellranger_aggr_HV_revised.csv
cellranger aggr --id=HV2 --csv=aggr_HV2.csv
#cellranger aggr --id=HV1 --csv=aggr_HV1.csv
#cellranger aggr --id=MD_Leigh2 --csv=aggr_MD2.csv
#cellranger aggr --id=MD_Leigh1 --csv=aggr_MD1.csv
#cellranger aggr --id=MD_Leigh3 --csv=aggr_MD3.csv
#cellranger aggr --id=MD_Leigh4 --csv=aggr_MD4.csv
#cellranger aggr --id=MD_Leigh5 --csv=aggr_MD1.csv
#cellranger aggr --id=MD_Leigh --csv=aggr_MD.csv
#cellranger aggr --id=HV4 --csv=cellranger_aggr_HV.csv
#cellranger aggr --id=HV4_1 --csv=cellranger_aggr_HV.csv --normalize=none
# HV4 crashed becaues the sequncing depth differentiation
#cellranger aggr --id=HV5 --csv=aggr_HV1.csv
