from pyspark.sql.dataframe import DataFrame
from pyspark.sql import SparkSession
import dxpy
import subprocess
from typing import List

def get_geno_database(case_dxlink: dict, control_dxlink: dict) -> dict:
    '''
    Returns dxlink to shared database between two cohorts
    '''
    case_databases = dxpy.DXRecord(case_dxlink).get_details()['databases']
    control_databases = dxpy.DXRecord(control_dxlink).get_details()['databases']
    if case_databases != control_databases:
        raise ValueError("Case and Control cohorts must be derived from the same database")
    database_dxlink = list(case_databases[0].values())[0]
    return database_dxlink

def get_id_field(dxlink: dict) -> str:
    '''
    Returns the <entity_name>,<field_name> corresponding to participant IDs
    '''
    details = dxpy.DXRecord(dxlink).get_details()
    containers = details['dashboardConfig']['cohort_browser']['containers']
    ted_container = next((i for i in containers if i['id'] == 'ted_container'), None)
    patient = ted_container['tiles'][0]
    entity, field = list(patient['dataQuery']['fields'])[0].split("$")
    return f'{entity}.{field}'

def get_cohort_ids(dxlink: dict) -> list:
    '''
    Returns list of cohort IDs
    '''
    dxrecord = f'{dxpy.PROJECT_CONTEXT_ID}:{dxpy.DXRecord(dxlink).get_id()}'
    id_field = get_id_field(dxlink)
    cmd = ['dx', 'extract_dataset', dxrecord, '--fields', id_field, '-o', '-']
    print(f'Running "{" ".join(cmd)}"')
    data = subprocess.run(cmd, stdout = subprocess.PIPE).stdout.decode().splitlines()
    data.pop(0)
    return data

def compute_af(db_dxlink: dict,
               case_ids: List[str],
               control_ids: List[str],
               tb: str,
               anno: str,
               spark: SparkSession) -> DataFrame:
    '''
    Computes allele frequencies for case and control cohorts given cohort IDs and genotype table
    '''
    db = dxpy.DXDatabase(db_dxlink).describe()['uniqueDatabaseName']
    case = ",".join(case_ids)
    control = ",".join(control_ids)
    Ncase = 2 * len(case_ids)
    Ncontrol = 2 * len(control_ids)
    sql_cmd = f'''
    SELECT ac_cases.a_id AS Variant,
           ( ac_cases.ac / CAST({Ncase} AS float) ) AS FreqCases,
           {len(case_ids)}l AS Ncases,
           ( ac_controls.ac / CAST({Ncontrol} AS float) ) AS FreqControls,
           {len(control_ids)}l AS Ncontrols
    FROM
        (SELECT a_id, Sum(sample_allele_count) AS ac
         FROM (SELECT *,
              ( CASE
                  WHEN type = 'hom' THEN 2
                  ELSE 1
               END ) AS sample_allele_count
               FROM   {db}.{tb}
               WHERE  {tb}.a_id IS NOT NULL
               AND  {tb}.sample_id IN ({case}))
         GROUP BY a_id) AS ac_cases
    INNER JOIN
        (SELECT a_id, Sum(sample_allele_count) AS ac
         FROM (SELECT *,
              ( CASE
                  WHEN type = 'hom' THEN 2
                  ELSE 1
               END ) AS sample_allele_count
               FROM   {db}.{tb}
               WHERE  {tb}.a_id IS NOT NULL
               AND  {tb}.sample_id IN ({control}))
         GROUP BY a_id) AS ac_controls
    ON ac_cases.a_id = ac_controls.a_id
    '''
    df = spark.sql(sql_cmd)
    return df