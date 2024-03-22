#!/usr/bin/env bash

set -e -o pipefail

main() {
    
    export PYSPARK_PYTHON=python3
    export PYTHONPATH=/opt/:$PYTHONPATH

    dx-spark-submit \
        --collect-log \
        --log-level INFO \
        --conf spark.executorEnv.PYTHONPATH="$PYTHONPATH" \
        /opt/is-gwas-spark.py "${case_cohort}" "${control_cohort}" "${geno_tb}" "${min_maf}" "${firth}" "${test}"

    hdfs dfs -copyToLocal hdfs:///tmp/csv /tmp
    mkdir -p $HOME/out/csv
    mkdir -p $HOME/out/manhattan

    # Concatenate all part files into one and add header
    echo "Variant,FreqCases,Ncases,FreqControls,Ncontrols,Beta,SE,log10p" > $HOME/out/csv/"${out}".csv
    cat /tmp/csv/part* >> $HOME/out/csv/"${out}".csv

    # Produce interactive Manhattan plot
    case_id=$(dx-jobutil-parse-link "${case_cohort}")
    control_id=$(dx-jobutil-parse-link "${control_cohort}")
    case_name=$(dx describe "${case_id}" --name)
    control_name=$(dx describe "${control_id}" --name)
    python3 /opt/manhattan.py $HOME/out/csv/"${out}".csv --out "${out}.html" --title "${case_name} vs ${control_name}"
    mv "${out}.html" $HOME/out/manhattan/
    dx-upload-all-outputs

}

