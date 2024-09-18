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



df =  pd.read_csv('stars_data-v1.csv')




df['Star'] = df['star_name'].apply(lambda x: re.match(r'^[^_]+', x).group(0))

# Create a table with all RVs and corresponding file names
table = df.groupby('Star').agg({
    'rv': lambda x: list(x), 
    'star_name': lambda x: list(x)
}).reset_index()

# Flatten the RV and star_name lists in the table
table = table.explode(['rv', 'star_name']).reset_index(drop=True)

# Print the table
print(table)

# Plotting box plots for each star group
plt.figure(figsize=(12, 8))
sns.set(style="whitegrid")

# Boxplot for RVs grouped by Star
sns.boxplot(x='Star', y='rv', data=df)

plt.xlabel('Star Name')
plt.ylabel('Radial Velocity (RV)')
plt.title('Box Plot of Radial Velocities Grouped by Star Name')
plt.grid()
plt.xticks(rotation=90)
plt.show()


table.to_csv('rv-check-corrected.csv')  






















































#