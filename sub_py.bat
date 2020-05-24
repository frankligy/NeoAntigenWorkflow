#BSUB -W 10:00
#BSUB -M 250000
#BSUB -n 10
#BSUB -o /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.out
#BSUB -e /data/salomonis2/LabFiles/Frank-Li/job_dump/%J.err


# you need to make sure necessary packages has been pre-installed in the following virtual enviroment
cd $(pwd)
module load anaconda3
conda init --all bash
source activate python3
# proxy_on 
python3 mhcPresent.py -i /data/salomonis2/LabFiles/Frank-Li/breast_cancer_run/all_events_checkGTEx/Hs_RNASeq_top_alt_junctions-PSI_EventAnnotation.txt -t breast_all_GTEx -o /data/salomonis2/LabFiles/Frank-Li/breast_cancer_run/all_events_checkGTEx -d /data/salomonis2/LabFiles/Frank-Li/python3/data -k 8 -H HLA-A01:01,HLA-A03:01,HLA-B07:02,HLA-B27:05,HLA-B58:01 -s /data/salomonis2/LabFiles/Frank-Li/python3/netMHCpan-4.1/netMHCpan -M MHCI -m singleSample -C 10 -c True
