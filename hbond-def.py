import pandas as pd
import matplotlib.pyplot as plt
import re
import numpy as np

# Set font sizes
plt.rcParams.update({'font.size': 15, 'axes.titlesize': 20, 'axes.labelsize': 20, 'xtick.labelsize': 16, 'ytick.labelsize': 16})

""" Load all solvent_avg
    Filter and get only Frac>0.1 into final_df
    Merge the df of all the reps and normalize
"""

def process_df(df):
    # Group by unique combinations of "Acceptor" and "Donor" columns and sum the "Frac" values
    grouped = df.groupby(['#Acceptor', 'Donor'])['Frac'].sum()

    # Create a list to store DataFrames for each unique combination
    dfs = []
    
    # Loop over unique combinations
    for (acceptor, donor), frac_sum in grouped.items():
        # Filter the original DataFrame for the current combination
        filtered_df = df[(df['#Acceptor'] == acceptor) & (df['Donor'] == donor)].copy()
        # Add a new column with the summed "Frac" values
        filtered_df.loc[:, 'Frac_Sum'] = frac_sum
        # Filter rows with "Frac" > 0.1
        filtered_df = filtered_df[filtered_df['Frac'] > 0.1]
        filtered_df.drop(columns=['Frac', 'AvgDist', 'DonorH', 'AvgAng'], inplace=True)        
        # Append the filtered DataFrame to the list
        dfs.append(filtered_df)

    # Concatenate all DataFrames in the list to get the final DataFrame
    final_df = pd.concat(dfs, ignore_index=True)

    final_df.drop_duplicates(subset='Frac_Sum', inplace=True)
    return final_df

def modify_tick_label(label):
    match = re.search(r'_(\d+)@', label)
    if match:
        numeric_part = int(match.group(1)) - 4
        modified_label = re.sub(r'_(\d+)@', '_' + str(numeric_part) + '@', label)
        return modified_label
    return label

# Read the files into DataFrames and process them
#df1 = pd.read_csv("solvent_avg.dat", delim_whitespace=True)
#df2 = pd.read_csv("solvent_avg-r1.dat", delim_whitespace=True)
#df3 = pd.read_csv("solvent_avg-r2.dat", delim_whitespace=True)

df1 = pd.read_csv("solute_avg.dat", delim_whitespace=True)
df2 = pd.read_csv("solute_avg-r1.dat", delim_whitespace=True)
df3 = pd.read_csv("solute_avg-r2.dat", delim_whitespace=True)
pdf1 = process_df(df1)
pdf2 = process_df(df2)
pdf3 = process_df(df3)

# Concatenate and group by unique pairs of "#Acceptor" and "Donor", then take the mean of 'Frac_Sum'
merged_df = pd.concat([pdf1, pdf2,pdf3], ignore_index=True)
merged_df = merged_df.groupby(['#Acceptor', 'Donor'])['Frac_Sum'].mean().reset_index()
merged_df = merged_df.sort_values(by='Frac_Sum', ascending=False).reset_index()

#Replace the COV_285 to ASP-Glc_285 for paper clarity
# Define the replacement function
def replace_cov_285(text):
    return re.sub(r'COV_289', 'ASP-Glc_289', text)

merged_df = merged_df.map(lambda x: replace_cov_285(x) if isinstance(x, str) else x)
print(merged_df)

# Pivot the DataFrame for heatmap
heatmap_data = merged_df.pivot(index='#Acceptor', columns='Donor', values='Frac_Sum')

# Modify tick labels
heatmap_data.index = [modify_tick_label(label) for label in heatmap_data.index]
heatmap_data.columns = [modify_tick_label(label) for label in heatmap_data.columns]

# Create the heatmap using imshow
plt.figure(figsize=(20, 15))  # Increase figure size for better readability
plt.imshow(heatmap_data, cmap='GnBu', aspect='auto')
plt.colorbar()#label='Frac_Sum')
#plt.title('WT')
#plt.xlabel('Donor')
#plt.ylabel('Acceptor')
plt.xticks(np.arange(len(heatmap_data.columns)), heatmap_data.columns, rotation=45, ha='right')
plt.yticks(np.arange(len(heatmap_data.index)), heatmap_data.index)

# Annotate each cell with its value, excluding NaN values
for i in range(len(heatmap_data.index)):
    for j in range(len(heatmap_data.columns)):
        value = heatmap_data.iloc[i, j]
        if not pd.isnull(value):
            plt.text(j, i, '{:.2f}'.format(value), ha='center', va='center', color='black')

plt.tight_layout()
plt.savefig("./aWT-hbond-solute-TESt-PAPER-QUALITY2.pdf")
plt.show()

