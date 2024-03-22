import argparse
import dash_bio as dashbio
import pandas as pd

def parse_variant_col(df: pd.DataFrame) -> pd.DataFrame:
    '''
    Split variant column into CHR, BP, REF, and ALT
    '''
    df[['CHR', 'BP', 'REF', 'ALT']] = df['Variant'].str.split('_', expand=True)
    return df

def remove_sex_chr(df: pd.DataFrame) -> pd.DataFrame:
    return df[~df['CHR'].isin(['X', 'Y'])]

def thin_results(df: pd.DataFrame, thin: int = 10000) -> pd.DataFrame:
    '''
    For ease of plotting, thin results
    '''
    df.sort_values(['CHR', 'BP'], inplace=True, ignore_index=True)
    tmp = df[df['log10p'] <= 1.3]
    tmp = df.sample(n = thin, random_state = 5131989)
    df = df[df['log10p'] > 1.3]
    df = pd.concat([df, tmp])
    df['CHR'] = df['CHR'].astype(int)
    df['BP'] = df['BP'].astype(int)
    df.sort_values(['CHR', 'BP'], inplace=True, ignore_index=True)
    return df

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("csv")
    parser.add_argument("--out", default="manhattan.html")
    parser.add_argument("--title", default="Manhattan Plot")
    args = parser.parse_args()

    df = pd.read_csv(args.csv)
    df.dropna(subset=['log10p'], inplace=True, ignore_index=True)
    df = parse_variant_col(df)
    df = remove_sex_chr(df)
    df = thin_results(df)
    fig = dashbio.ManhattanPlot(df,
                                p='log10p',
                                snp='Variant',
                                gene=None,
                                logp=False,
                                title=args.title)
    
    # Workaround for apparent bug where manhattan plot shrinks each time window is resized
    fig.layout.template.layout.xaxis.automargin = False
    
    fig.write_html(args.out)

