# batch_processData: Will get the coverage of each basepair and convert index.html to csv 

#!/bin/sh
eval "$(conda shell.bash hook)"

# Getting coverage of each amplification
conda activate breseqAnalysis # has html2csv installed 
cd bamfiles
python ./plotCoverageAllBAM.py
cd ..

# Getting CSV from html file
cd htmlfiles 

# Loop through files in the directory
for file in *; do
     if [ -f "$file" ]; then  # Check if it's a regular file
         # Check if the file is an HTML file
        if [[ "$file" == *.html ]]; then
            # Convert HTML to CSV using html2csv
            html2csv "$file" -o "${file%.html}.csv"
        fi
    fi
done

cd ..  # Return to parent directory
#conda deactivate
