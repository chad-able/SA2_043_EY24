import pandas as pd
import re

filename = "E:/codes/RC24/PWMapping/usgs_newts_data.csv"
df = pd.read_csv(filename, encoding='ISO-8859-1', low_memory = False)

basin_medians = df.groupby('BASIN')['Li'].median().reset_index()
field_medians = df.groupby('FIELD')['Li'].median().reset_index()
formation_medians = df.groupby('FORMATION')['Li'].median().reset_index()

basin_medians.columns = ['name', 'Li_concentration']
field_medians.columns = ['name', 'Li_concentration']
formation_medians.columns = ['name', 'Li_concentration']

basin_medians['type'] = 'basin'
formation_medians['type'] = 'formation'
field_medians['type'] = 'field'

all_medians = pd.concat([basin_medians, field_medians, formation_medians], ignore_index=True)

all_medians.columns = ['name', 'Li_concentration', 'type']

# all_medians.to_csv('Li_medians.csv', index=False)

plays = ['Permian']
pattern = '|'.join(re.escape(play) for play in plays)

df = df[df[['BASIN', 'FIELD', 'FORMATION']].apply(
    lambda row: row.str.contains(pattern, case=False, na=False).any(), axis=1
)]

print(df['Li'].median())