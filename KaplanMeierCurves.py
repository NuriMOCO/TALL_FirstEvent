import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test

# Load the DataFrame
FirstEvent = pd.read_csv('/Users/nuricornejo/Documents/LABORATORIOS/INMEGEN/TALL/BinaryFirstEvent242.csv') 
BinaryClusters = pd.read_csv('/Users/nuricornejo/Documents/LABORATORIOS/INMEGEN/TALL/final/BinaryClusters1825.csv') # Metrics dichotomized by Youden's index

# Merge the dataframes
merged = pd.concat([FirstEvent, BinaryClusters], axis=1)


# Ensure 'EventFreeSurvDays' and 'BinaryFirstEvent1825' are numeric
merged['EventFreeSurvDays1825'] = pd.to_numeric(merged['EventFreeSurvDays1825'], errors='coerce')
merged['BinaryFirstEvent1825'] = pd.to_numeric(merged['BinaryFirstEvent1825'], errors='coerce')

# Rename columns for better readability
new_column_names = {
    'num_communities_youden': 'Number of Communities',
    'diameter_youden': 'Diameter',
    'num_nodes_youden': 'Number of Nodes',
    'modularity_score_youden': 'Modularity Score',
    'global_clustering_youden': 'Global Clustering',
    'efficiency_youden': 'Efficiency',
    'mean_CommunitySize_youden': 'Mean Community Size',
    'max_community_size_youden': 'Maximum Community Size'
}

# Rename the columns
BinaryClusters.rename(columns=new_column_names, inplace=True)

# Merge again to ensure the renamed columns are in the merged DataFrame
merged = pd.concat([FirstEvent, BinaryClusters], axis=1)

# Initialize the figure and subplots
fig, axes = plt.subplots(2, 4, figsize=(12, 7))
fig.subplots_adjust(hspace=0.1, wspace=0.3, top=0.7, bottom=0.1)  # Adjusted bottom to lower the chart

# Create Kaplan-Meier plots for each metric
kmf = KaplanMeierFitter()
metrics = BinaryClusters.columns[0:8]  # Assuming the first column is patient ID and the next 8 are metrics
labels = ['A)', 'B)', 'C)', 'D)', 'E)', 'F)', 'G)', 'H)']

for i, (metric, label) in enumerate(zip(metrics, labels)):
    ax = axes[i // 4, i % 4]
    group1 = merged[merged[metric] == 1]
    group2 = merged[merged[metric] == 2]
    
    kmf.fit(group1['EventFreeSurvDays1825'], event_observed=group1['BinaryFirstEvent1825'], label="Group 1")
    kmf.plot_survival_function(ax=ax, ci_show=False, color="tab:orange")
    
    kmf.fit(group2['EventFreeSurvDays1825'], event_observed=group2['BinaryFirstEvent1825'], label="Group 2")
    kmf.plot_survival_function(ax=ax, ci_show=False, color="tab:purple")
    
    # Calculate the p-value
    result = logrank_test(group1['EventFreeSurvDays1825'], group2['EventFreeSurvDays1825'],
                          event_observed_A=group1['BinaryFirstEvent1825'],
                          event_observed_B=group2['BinaryFirstEvent1825'])
    p_value = result.p_value
    
    # Add p-value text inside the plot
    ax.text(0.2, 0.1, f'p = {p_value:.3f}', transform=ax.transAxes, horizontalalignment='center')
    
    ax.set_title(f'{metric}', fontsize=10)
    ax.set_xlabel('Event-Free Time in Days', fontsize=8)
    ax.set_ylabel('Cumulative Event-Free Probability', fontsize=8, labelpad=10)  # Adjusted labelpad to move y-axis label to the left
    ax.tick_params(axis='both', which='major', labelsize=8)
    ax.set_xlim(left=0)  # Adjust the width of the plots
    ax.set_ylim(0.6, 1.005) # Adjust the height of the plots .9 for 365, .85 for 730, .6 for 1825
    ax.set_yticks([ 0.6, 0.7, 0.8, 0.90,  1.0])
    ax.yaxis.set_major_formatter(plt.FuncFormatter(lambda y, _: f'{y:.2f}'))
    
    # Remove the legend from individual plots
    ax.get_legend().remove()
    
    # Add risk table
    n_patients_group1 = len(group1)
    n_events_group1 = group1['BinaryFirstEvent1825'].sum()
    n_patients_group2 = len(group2)
    n_events_group2 = group2['BinaryFirstEvent1825'].sum()
    
    total_patients = n_patients_group1 + n_patients_group2
    perc_group1 = (n_patients_group1 / total_patients) * 100
    perc_group2 = (n_patients_group2 / total_patients) * 100
    
    risk_table_data = [
        ['', 'Patients', 'Events'],
        ['Group 1', f'{n_patients_group1} ({perc_group1:.1f}%)', n_events_group1],
        ['Group 2', f'{n_patients_group2} ({perc_group2:.1f}%)', n_events_group2]
    ]
    
    table = ax.table(cellText=risk_table_data, colLabels=None, cellLoc='center', loc='bottom', bbox=[0, -0.4, 1, 0.2])
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.2)
    table.auto_set_column_width([0, 1, 2])
    for key, cell in table.get_celld().items():
        cell.set_linewidth(0)
        cell.set_height(0.01)  # Reduce the space between lines

# Add a single legend outside the plots
handles, labels = ax.get_legend_handles_labels()
fig.legend(handles, ['Group 1 (> Cutoff point)', 'Group 2 (< Cutoff point)'], loc='upper right', ncol=2)
fig.text(0.05, 0.95, 'Year 9+', ha='left', fontsize=12, fontweight='bold')
# Add legends on the left and right
#fig.legend(handles, ['Group 1 (> Cutoff point)', 'Group 2 (< Cutoff point)'], loc='upper left', ncol=1)
#fig.legend(handles, ['First Event Occurrence Analysis Year 1'], loc='upper right', ncol=1)

plt.tight_layout(rect=[0, 0, 1, 0.93])  # Adjust layout to make room for the legend
plt.savefig('/Users/nuricornejo/Documents/LABORATORIOS/INMEGEN/TALL/KaplanMeierPlots9.png', dpi=300)  # Adjust the DPI value as needed
plt.show()
