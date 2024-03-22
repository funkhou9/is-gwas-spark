import argparse
import json
import pyspark
from pyspark.sql.functions import lit
from calculate_af import get_geno_database, get_cohort_ids, compute_af
from is_gwas import isGWAS_udf

sc = pyspark.SparkContext()
spark = pyspark.sql.SparkSession(sc)
spark.sql("SET spark.sql.shuffle.partitions=64")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("case_cohort_link")
    parser.add_argument("control_cohort_link")
    parser.add_argument("geno_tb")
    parser.add_argument("min_maf", type = float)
    parser.add_argument("firth")
    parser.add_argument("test")
    args = parser.parse_args()
    firth_bool = json.loads(args.firth)
    anno_tb = args.geno_tb.replace('genotype', 'annotation')
    
    # Get database dxlink
    case_dxlink = json.loads(args.case_cohort_link)
    control_dxlink = json.loads(args.control_cohort_link)
    database_dxlink = get_geno_database(case_dxlink, control_dxlink)

    # Get cohort IDs
    case_ids = get_cohort_ids(case_dxlink)
    control_ids = get_cohort_ids(control_dxlink)

    # Compute cohort-specific allele frequencies
    df = compute_af(db_dxlink = database_dxlink,
                    case_ids = case_ids,
                    control_ids = control_ids,
                    tb = args.geno_tb,
                    anno = anno_tb,
                    spark = spark)

    # Apply is-GWAS
    df_results = df.withColumn("isgwas", isGWAS_udf('Ncases',
                                                    'Ncontrols',
                                                    'FreqCases',
                                                    'FreqControls',
                                                    lit(args.min_maf),
                                                    lit(False),
                                                    lit(firth_bool),
                                                    lit(args.test),
                                                    lit(1e-10),
                                                    lit(20)))
    df_results = df_results.selectExpr(
        "*",
        "isgwas.Beta as Beta",
        "isgwas.SE as SE",
        "isgwas.log10p as log10p"
    ).drop("isgwas")

    df_results.write.csv(f'hdfs:///tmp/csv', header = False)