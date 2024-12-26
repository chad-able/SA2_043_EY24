import pandas as pd
import re

filename = "E:/codes/RC24/PWMapping/usgs_newts_data.csv"
df = pd.read_csv(filename, encoding='ISO-8859-1', low_memory = False)

plays = ['Permian']
pattern = '|'.join(re.escape(play) for play in plays)

df = df[df[['BASIN', 'FIELD', 'FORMATION']].apply(
    lambda row: row.str.contains(pattern, case=False, na=False).any(), axis=1
)]

print(df['Li'].median())