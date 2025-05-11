import pandas as pd
import numpy as np

file_path='or131-rwl-noaa.txt'
df=pd.read_csv(file_path,comment='#',sep='\t',na_values='NA')
print(df.shape)
print(df.describe())
print(df.info())

#outliers
#separate year and treering data
year_col='age_CE'
df_cols = [col for col in df.columns if col != year_col]

for col in df_cols:
    series = df[col]
    z = (series - series.mean())/series.std()
    df.loc[z.abs()>3, col] = np.nan

# missing values
df.fillna(df.median(),inplace=True)
print(df.isnull().sum())

# check format of col names: good
print(df.columns.values)

# check for duplicates: no duplicate
print(df.duplicated().sum())

# standardize treering
tree_cols = [col for col in df.columns if col != year_col]
df[tree_cols] = df[tree_cols].apply(pd.to_numeric,errors='coerce')
df[tree_cols]=df[tree_cols].apply(lambda x: (x-x.mean())/x.std() if x.std()!=0 else x)
print(df.describe())

#df.to_csv('cleaned_data.csv',index=False,float_format='%.4f')
