# batch_GDrun: Will make both 'html' and 'tsv' comparision table from gd files in "gdfiles" folder
#!/bin/sh
eval "$(conda shell.bash hook)"
conda activate breseq # activating environment with breseq
cd gdfiles
pwd

# edit reference file name
gdtools COMPARE -o ../AllComparison.html  -r ../../MG1655_NC_00913.gbk -f HTML *.gd 
gdtools COMPARE -o ../ALLGD.tsv -r ../../MG1655_NC_00913.gbk -f TSV *.gd

cd ..

conda deactivate # deactivating enviroment with breseq 
