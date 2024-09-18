#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 10:06:44 2024

@author: heitor
"""

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re
from io import StringIO



df =  pd.read_csv('stars_data.csv')



# Extract star name from the filename using a regex
df['Star'] = df['star_name'].apply(lambda x: re.match(r'^[^_]+', x).group(0))

# Group by the star name
grouped = df.groupby('Star')

# Plotting
plt.figure(figsize=(12, 8))
sns.set(style="whitegrid")

# Scatter plot of RVs for each star group
for name, group in grouped:
    sns.scatterplot(x=[name] * len(group), y=group['rv'], label=name)

plt.xlabel('Star Name')
plt.ylabel('Radial Velocity (RV)')
plt.title('Radial Velocities Grouped by Star Name with Outliers')
plt.xticks(rotation=45)
plt.legend(title='Star Name', bbox_to_anchor=(1.05, 1), loc='upper left')
plt.show()






























































#