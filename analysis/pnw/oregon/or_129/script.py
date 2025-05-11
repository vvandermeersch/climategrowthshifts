import pandas as pd
import numpy as np

# locate the dataset
file_path = 'or129-rwl-noaa.txt'
df= pd.read_csv(file_path,comment='#',sep='\t',na_values="NA")
print(df.shape)
print(df.info())


# outlier
# separate year and treering data
year_col="age_CE"
df_cols= [col for col in df.columns if col != year_col]

for col in df_cols:
    series = df[col]
    z = (series-series.mean())/series.std()
    # replace outlier with NAN
    df.loc[z.abs()>3,col]=np.nan
    
# missing values: use median number to replace NAN
df.fillna(df.median(),inplace=True)
print(df.isnull().sum())

#check the format of col names: good format
col=df.columns.values
print(col)

# check for duplicate: no duplicate
print(df.duplicated().sum())

#standardize treering
tree_cols = [col for col in df.columns if col != year_col]
df[tree_cols]=df[tree_cols].apply(lambda x: (x-x.mean())/x.std())
print(df.describe())

#df.to_csv('cleaned_data.csv',index=False)