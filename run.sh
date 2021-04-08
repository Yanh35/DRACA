#/bin/bash
if [ $# -lt 1 ]; then
    echo "
usage: bash run.sh lncRNA_symbol
"
    exit 1
fi

number=$(grep -wrn $1 data/lncrna-id.txt | wc -l)
if [ $number  -gt 0 ];then
    matlab -nodisplay -nosplash -nodesktop -r "run('DRACA/DRACA_breast.m');exit;"
    matlab -nodisplay -nosplash -nodesktop -r "run('DRACA/DRACA_lung.m');exit;"
    matlab -nodisplay -nosplash -nodesktop -r "run('DRACA/DRACA_colorectal.m');exit;"
    python DRACA/judge.py $1 breast_score.csv lung_score.csv colorectal_score.csv
else
    echo $1 'is not in the list of predict' 
fi
