#------------- Explanation -------------#
"""
Author : Rick Vonk

This script is used to make a 3D histogram plot to assist in analysis of combination prediction
aswel as a txt file with analytics on the 3d plot data

Input:
    - Total experiment folder
    - metrics wanted for plotting
    - highlighted cases for plotting
    - view angle for 3d plot

Output:
    - csv-file with total plot data
    - txt-file with analysis
    - HTML-file for plotting 1d
    - 1d plot with all cases and the selected highlights
    - 3d plot with all cases
    - Heatmap with total data for all cases
    - weighted Heatmap with data of all cases

"""
#------------ Import -------------#
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import numpy as np
import pandas as pd
import glob
import json
import os
import csv

#----------------- inputs ------------------#
working_dir = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments'
plot_view_angle= [60,-30] # standard is [30,-60]

# total possibilities : ['AF3_confidence', 'total_energy', 'elec_energy','deltaG', 'haddock_score']
metrics = ['AF3_confidence', 'haddock_score', 'deltaG', 'elec_energy']

highlight = ['a1_b12', 'a1_b8', 'a10_b8', 'a10_b11', 'a11_b12', 'a11_b4', 'a11_b6', 'a11_b10', 'a2_b7', 'a2_b12', 'a2_b11', 'a2_b1', 'a2_b2', 'a3_b11', 'a3_b8', 'a3_b7', 'a3_b1', 'a3_b12', 'a4_b8', 'a4_b7', 'a4_b12', 'a5_b10', 'a5_b4', 'a6_b8', 'a7_b12', 'a7_b8', 'a7_b10', 'a8_b11', 'a9_b11', 'a9_b8'] 


# ----------------- code ------------------ #
def Extraction_data(Exp_dir):
    df_list = []
    for experiment in glob.glob(os.path.join(Exp_dir, "tcr_alpha_*")):
        experiment_name = os.path.basename(experiment)
        csv_data_file = os.path.join(experiment,'Result',f'{experiment_name}.csv')
        csv_df = pd.read_csv(glob.glob(csv_data_file)[0])

        cols_to_keep = [col for col in csv_df.columns if not col.startswith("norm_")]
        csv_df = csv_df[cols_to_keep]

        df_list.append(csv_df)


    df_summed = pd.concat(df_list, ignore_index=True)

    return df_summed

def Normalisation(df_summed):
    norm_cols = ["total_energy", "elec_energy", "confidence", "deltaG", "haddock_score"]
    for col in norm_cols:
        min_val = df_summed[col].min()
        max_val = df_summed[col].max()
        range_val = (max_val - min_val) if max_val != min_val else 1

        norm_series = (df_summed[col] - min_val) / range_val

        if col in ["elec_energy", "total_energy", "deltaG", "haddock_score"]:
            norm_series = 1 - norm_series

        if col == "confidence":
            out_col = "norm_AF3_confidence"
        else:
            out_col = f"norm_{col}"

        df_summed[out_col] = norm_series

    return df_summed

def Calculation_totals(df_summed, selection, output_dir):
    df_list = []
    df_summed['base_ID'] = df_summed['experiment'].str.extract(r'^(.*)_rs')
    for base, group in df_summed.groupby('base_ID'):
        avg_df = group[selection].mean().to_frame().T
        avg_df.insert(0, 'base_ID', base)
        df_list.append(avg_df)
    
    df_combined = pd.concat(df_list, ignore_index=True)

    columns = selection + ['experiment']
    csv_file = os.path.join(output_dir, '3D_plot_info.csv')

    df_csv = df_combined.copy()
    df_csv['Total'] = df_csv.select_dtypes(include='number').sum(axis=1).round(2)
    df_csv.to_csv(csv_file, index=False)

    return df_combined, csv_file

def Coordinate_collection(df_combined):
    coordinate_xyz = {}
    for index, row in df_combined.iterrows():
        base_id = row['base_ID']
        parts = base_id.split('_')
        x_pos = parts[0]
        y_pos = parts[1]
        z_pos = round(row.drop('base_ID').sum(),2)
        coordinate_xyz[base_id] = [x_pos, y_pos, z_pos]

    return coordinate_xyz

def Plotting_3d(data_dict, output_dir, plot_view_angle=None):
    fig = plt.figure(figsize=(10, 8))
    ax = fig.add_subplot(projection='3d')
    xs_labels = sorted(list(set(val[0] for val in data_dict.values())), key=lambda x: int(x[1:]))
    ys_labels = sorted(list(set(val[1] for val in data_dict.values())), key=lambda x: int(x[1:]))

    x_map = {label: i for i, label in enumerate(xs_labels)}
    y_map = {label: i for i, label in enumerate(ys_labels)}

    xpos = []
    ypos = []
    dz = []

    for key, val in data_dict.items():
        x_label, y_label, height = val
        xpos.append(x_map[x_label])
        ypos.append(y_map[y_label])
        dz.append(height)

    zpos = np.zeros(len(dz))
    dx = dy = 0.4 

    ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='skyblue', edgecolor='black', linewidth=0.5, alpha=0.8, zsort='average')
    ax.set_xticks([i + 0.25 for i in range(len(xs_labels))])
    ax.set_xticklabels(xs_labels)
    ax.set_yticks([i + 0.25 for i in range(len(ys_labels))])
    ax.set_yticklabels(ys_labels)

    ax.set_xlabel('A-chain')
    ax.set_ylabel('B-Chain')
    ax.set_zlabel('normalised Totals')

    if not plot_view_angle:
        plot_view_angle = [30, -60]
    ax.view_init(elev=plot_view_angle[0],azim=plot_view_angle[1])

    plt.savefig(os.path.join(output_dir,'3d_plot'), dpi=300)
    plt.close()
    return

def Plotting_1d(df_combined, output_dir, highlight):
    df = df_combined.copy()
    numeric_cols = df.columns.drop('base_ID')
    df['Total'] = df[numeric_cols].sum(axis=1)
    
    df = df.sort_values(by='Total', ascending=False).reset_index(drop=True)
    ids = df['base_ID'].tolist()
    x = np.arange(len(ids))
    bottom = np.zeros(len(ids))

    edge_colors = ['black' if i in highlight else 'none' for i in ids]
    color_map = {
        'norm_AF3_confidence': '#1f77b4',
        'norm_total_energy': '#ff7f0e',
        'norm_elec_energy': '#2ca02c',
        'norm_deltaG' : '#e8e337',
        'norm_haddock_score' : '#FC9483'
    }

    plt.figure(figsize=(12, 6))

    for dataset_name in numeric_cols:
        values = df[dataset_name].values
        plt.bar(x, values, bottom=bottom, edgecolor=edge_colors, 
                linewidth=1.5, color=color_map.get(dataset_name, 'gray'), label=dataset_name)
        
        bottom += values
    
    plt.xticks(x, ids, rotation=45, ha='right', fontsize=6)
    plt.xlabel('alpha/beta chain combination')
    plt.ylabel('Normalized Combined Score')
    plt.title('1d plot of all cases')
    plt.legend()
    plt.tight_layout()

    plot_out = os.path.join(output_dir,'1d_plot')
    plt.savefig(plot_out, dpi=300)
    plt.close()
    return

def Plotting_HTML(df_combined, output_dir):
    df = df_combined.copy()
    numeric_cols = [col for col in df.columns if col.startswith('norm_')]

    for col in numeric_cols:
        df[col] = pd.to_numeric(df[col], errors='coerce')
        
    df['Total'] = df[numeric_cols].sum(axis=1)
    
    df = df.sort_values(by='Total', ascending=False).reset_index(drop=True)
    ids = df['base_ID'].tolist()

    color_map = {
        'norm_AF3_confidence': '#1f77b4',
        'norm_total_energy': '#ff7f0e',
        'norm_elec_energy': '#2ca02c',
        'norm_deltaG' : '#e8e337',
        'norm_haddock_score' : '#FC9483'
    }

    fig = go.Figure()

    for dataset_name in numeric_cols:
        fig.add_trace(go.Bar(
            x=ids,
            y=df[dataset_name],
            name=dataset_name,
            marker=dict(
                color=color_map.get(dataset_name, 'gray'),
                line=dict(color='rgba(0,0,0,0)', width=0) 
            ),
            hovertemplate=f"<b>%{{x}}</b><br>{dataset_name}: %{{y:.3f}}<extra></extra>"
        ))

    fig.update_layout(
        title='1D Interactive Plot of All Cases (Sorted by Combined Score)',
        xaxis=dict(
            title='Alpha/Beta Chain Combination',
            tickangle=45,
            tickfont=dict(size=9),
            type='category' 
        ),
        yaxis=dict(title='Normalized Combined Score'),
        barmode='stack',
        hovermode='x unified', 
        legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
        template="plotly_white"
    )

    plotly_html_div = fig.to_html(full_html=False, include_plotlyjs='cdn', div_id='plotly-graph')

    presets = {
        "top_5": ['a1_b12', 'a1_b8', 'a10_b8', 'a10_b11', 'a10_b12', 'a11_b12', 'a11_b4', 'a11_b6', 'a11_b10', 'a2_b7', 'a2_b11', 'a2_b12', 'a2_b1', 'a2_b2', 'a3_b1', 'a3_b8', 'a3_b11', 'a3_b7', 'a3_b12', 'a4_b8', 'a4_b7', 'a4_b12', 'a5_b10', 'a5_b4', 'a5_b8', 'a6_b8', 'a7_b12', 'a7_b8', 'a7_b10', 'a8_b11', 'a9_b11', 'a9_b8'],
        "top_5_weighted": ['a10_b8', 'a10_b11', 'a11_b4', 'a11_b6', 'a11_b12', 'a2_b7', 'a2_b2', 'a2_b12', 'a2_b1', 'a2_b11', 'a3_b3', 'a3_b1', 'a3_b7', 'a3_b9', 'a3_b2', 'a4_b7', 'a4_b8', 'a5_b10', 'a5_b4', 'a7_b12', 'a8_b11', 'a9_b11'],
        "bottom_5": ['a1_b7', 'a1_b1', 'a1_b3', 'a1_b9', 'a1_b2', 'a10_b6', 'a10_b9', 'a10_b10', 'a10_b5', 'a10_b3', 'a11_b1', 'a11_b7', 'a11_b11', 'a11_b5', 'a11_b3', 'a2_b3', 'a2_b5', 'a3_b5', 'a3_b6', 'a4_b4', 'a4_b1', 'a4_b6', 'a4_b2', 'a4_b11', 'a5_b12', 'a5_b7', 'a5_b1', 'a5_b9', 'a5_b3', 'a6_b1', 'a6_b2', 'a6_b9', 'a6_b3', 'a6_b5', 'a7_b6', 'a7_b1', 'a7_b3', 'a7_b9', 'a7_b11', 'a8_b2', 'a8_b3', 'a8_b7', 'a8_b4', 'a8_b5', 'a9_b7', 'a9_b2', 'a9_b9', 'a9_b5', 'a9_b3'] 
    }

    full_html_content = fr"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Interactive Plot with Search & Presets</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            .controls-container {{ 
                display: flex; 
                flex-wrap: wrap; 
                align-items: center; 
                gap: 20px; 
                margin-bottom: 15px; 
            }}
            .search-container {{ font-size: 16px; }}
            #graphSearchInput {{ padding: 8px; width: 400px; border: 1px solid #ccc; border-radius: 4px; }}
            
            .preset-container {{ display: flex; gap: 8px; }}
            .preset-btn {{
                padding: 8px 12px;
                font-size: 14px;
                font-weight: bold;
                border: 1px solid #007bff;
                background-color: #ffffff;
                color: #007bff;
                border-radius: 4px;
                cursor: pointer;
                transition: background-color 0.2s, color 0.2s;
            }}
            .preset-btn:hover {{
                background-color: #007bff;
                color: white;
            }}
            .preset-btn.active {{
                background-color: #0056b3;
                color: white;
                border-color: #0056b3;
            }}
            .clear-btn {{
                border-color: #dc3545;
                color: #dc3545;
            }}
            .clear-btn:hover {{
                background-color: #dc3545;
                color: white;
            }}
        </style>
    </head>
    <body>

        <div class="controls-container">
            <div class="search-container">
                <label for="graphSearchInput"><b>Search:</b> </label>
                <input type="text" id="graphSearchInput" placeholder="e.g. a1_b12, a7 or (a2_b7, a2_b2)">
            </div>

            <div class="preset-container">
                <button class="preset-btn" onclick="applyPreset('top_5', this)">Top 5</button>
                <button class="preset-btn" onclick="applyPreset('top_5_weighted', this)">Top 5 Weighted</button>
                <button class="preset-btn" onclick="applyPreset('bottom_5', this)">Bottom 5</button>
                <button class="preset-btn clear-btn" onclick="clearHighlights()">Clear All</button>
            </div>
        </div>

        {plotly_html_div}

        <script>
            // Safely inject Python data structures into JS via JSON
            const originalIds = {json.dumps(ids)};
            const presets = {json.dumps(presets)};
            
            // Centralized function to manage updating the Plotly DOM borders
            function updateGraphHighlights(matchConditionFn) {{
                const graphDiv = document.getElementById('plotly-graph');
                const newColors = [];
                const newWidths = [];

                originalIds.forEach((id) => {{
                    if (matchConditionFn(id)) {{
                        newColors.push('black');
                        newWidths.push(4.0);
                    }} else {{
                        newColors.push('rgba(0,0,0,0)');
                        newWidths.push(0.0);
                    }}
                }});

                const update = {{
                    'marker.line.color': [newColors],
                    'marker.line.width': [newWidths]
                }};

                Plotly.restyle(graphDiv, update);
            }}

            // 1. Text Search Input Event (Updated to support list inputs)
            document.getElementById('graphSearchInput').addEventListener('input', function(e) {{
                clearActiveButtonStyles(); // Clear button highlights if user types
                
                // Remove brackets or parentheses if the user pastes them like (x1,x2) or [x1,x2]
                let rawInput = e.target.value.replace(/[\[\]\(\)]/g, '');
                
                if (rawInput.trim() === "") {{
                    updateGraphHighlights(() => false); // Reset if search is empty
                    return;
                }}

                // Split by commas, trim extra spaces, and filter out any empty strings
                const searchTerms = rawInput.split(',')
                                            .map(term => term.trim().toLowerCase())
                                            .filter(term => term !== "");

                // Highlight an ID if it partially matches *any* of the terms in the user's list
                updateGraphHighlights((id) => {{
                    const currentId = String(id).toLowerCase();
                    return searchTerms.some(term => currentId.includes(term));
                }});
            }});

            // 2. Preset Button Logic
            function applyPreset(presetKey, buttonElement) {{
                document.getElementById('graphSearchInput').value = ""; // Clear text search
                clearActiveButtonStyles();
                buttonElement.classList.add('active'); // Style clicked button

                const targetList = presets[presetKey].map(id => String(id).toLowerCase());
                
                // Highlight if the ID strictly exists in our preset array
                updateGraphHighlights((id) => targetList.includes(String(id).toLowerCase()));
            }}

            // 3. Helper: Clear Everything
            function clearHighlights() {{
                document.getElementById('graphSearchInput').value = "";
                clearActiveButtonStyles();
                updateGraphHighlights(() => false);
            }}

            // Helper to reset button color states
            function clearActiveButtonStyles() {{
                document.querySelectorAll('.preset-btn').forEach(btn => btn.classList.remove('active'));
            }}
        </script>

    </body>
    </html>
    """

    plot_out = os.path.join(output_dir, 'HMTL_1d_plot_search_bar.html')
    with open(plot_out, 'w', encoding='utf-8') as f:
        f.write(full_html_content)
        
    return

def Heatmap_alpha_beta_total_comparison(df_combined, output_dir):
    df_combined['a_group'] = df_combined['base_ID'].str.extract(r'(a\d+)')
    df_combined['b_group'] = df_combined['base_ID'].str.extract(r'(b\d+)')
    df_combined['Total'] = df_combined.select_dtypes(include='number').sum(axis=1).round(2)
    
    success_matrix = df_combined.pivot_table(
    index="a_group", columns="b_group", values="Total", aggfunc="mean"
    )

    sorted_b = sorted(success_matrix.columns, key=lambda x: int(x.replace("b", "")))
    sorted_a = sorted(success_matrix.index, key=lambda x: int(x.replace("a", "")))
    success_matrix = success_matrix.reindex(index=sorted_a, columns=sorted_b)

    vmin_val = np.nanmin(success_matrix.values)
    vmax_val = np.nanmax(success_matrix.values)

    color_threshold = vmin_val + (vmax_val - vmin_val) * 0.4
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(success_matrix.values, cmap="YlGnBu", vmin=vmin_val, vmax=vmax_val, origin="upper")

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("Combined normalised score", rotation=-90, va="bottom")

    ax.set_xticks(np.arange(len(success_matrix.columns)))
    ax.set_yticks(np.arange(len(success_matrix.index)))
    ax.set_xticklabels(success_matrix.columns)
    ax.set_yticklabels(success_matrix.index)

    for i in range(len(success_matrix.index)):
        for j in range(len(success_matrix.columns)):
            val = success_matrix.values[i, j]
            if pd.isna(val):
                continue
            text_color = "white" if val > color_threshold else "black"
            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                color=text_color,
                fontweight="bold",
            )

    ax.set_title("Totals score per alpha/beta chain combination")
    ax.set_ylabel("A chain")
    ax.set_xlabel("B chain")
    plt.tight_layout()

    plot_out = os.path.join(output_dir, 'Heatmap_chain_and_totals_comparison.png')
    plt.savefig(plot_out, dpi=300, bbox_inches='tight')
    plt.close()
    return

def heatmap_weighted_alpha_chain(df_weighted_heatmap, output_dir, threshold_val=0):
    success_matrix = df_weighted_heatmap.pivot_table(
    index="a_group", columns="b_group", values="Adjusted_Total",
    )

    sorted_b = sorted(success_matrix.columns, key=lambda x: int(x.replace("b", "")))
    sorted_a = sorted(success_matrix.index, key=lambda x: int(x.replace("a", "")))
    success_matrix = success_matrix.reindex(index=sorted_a, columns=sorted_b)

    vmin_val = np.nanmin(success_matrix.values)
    vmax_val = np.nanmax(success_matrix.values)

    color_threshold = vmin_val + (vmax_val - vmin_val) * 0.6
    fig, ax = plt.subplots(figsize=(10, 6))
    im = ax.imshow(success_matrix.values, cmap="YlGnBu", vmin=vmin_val, vmax=vmax_val, origin="upper")

    cbar = ax.figure.colorbar(im, ax=ax)
    cbar.ax.set_ylabel("weighted normalised score", rotation=-90, va="bottom")

    ax.set_xticks(np.arange(len(success_matrix.columns)))
    ax.set_yticks(np.arange(len(success_matrix.index)))
    ax.set_xticklabels(success_matrix.columns)
    ax.set_yticklabels(success_matrix.index)

    for i in range(len(success_matrix.index)):
        for j in range(len(success_matrix.columns)):
            val = success_matrix.values[i, j]
            if pd.isna(val):
                continue
            if val < threshold_val:
                continue
            text_color = "white" if val > color_threshold else "black"
            ax.text(
                j,
                i,
                f"{val:.2f}",
                ha="center",
                va="center",
                color=text_color,
                fontweight="bold",
            )

    ax.set_title(f"weighted score per alpha/beta chain combination, with threshold {threshold_val}")
    ax.set_ylabel("A chain")
    ax.set_xlabel("B chain")
    plt.tight_layout()

    plot_out = os.path.join(output_dir, 'Heatmap_weighted_alpha_chain.png')
    plt.savefig(plot_out, dpi=300, bbox_inches='tight')
    plt.close()
    return

def Analytics(df_combined, coordinates, output_dir, selection):
    with open (os.path.join(output_dir, 'Total_plot_analytics.txt'), 'w', newline='') as f:
        
        # Best match per A-chain with top 5 B-chain's
        print('#-------------------- Best match per A-chain, top 5 --------------------#', file=f)
        
        df_match = pd.DataFrame.from_dict(coordinates, orient='index', 
                                         columns=['A_chain', 'B_chain', 'Total'])
        df_match = df_match.sort_values(['A_chain', 'Total'], ascending=[True, False])
        df_match['Rank'] = df_match.groupby('A_chain').cumcount() + 1

        top_5 = df_match[df_match['Rank'] <= 5].copy()
        top_5['Entry'] = top_5['B_chain'] + " (" + top_5['Total'].astype(str) + ")"

        side_by_side = top_5.pivot(index='Rank', columns='A_chain', values='Entry')
        print(side_by_side.to_string(), file=f)

        # Alpha / Beta scoring comparison
        # Here for each alpah chain i look at all beta chain values and compare them all to each other and make a map of comparison
        df = df_combined.copy()

        b_group_means = df.groupby('b_group')['Total'].transform('mean')
        df['Adjusted_Total'] = df['Total'] - b_group_means
        df_present = df[['a_group', 'b_group','Adjusted_Total']]
        
        df_sorted = df_present.sort_values(by=['a_group', 'Adjusted_Total'], ascending=[True, False]).copy()
        df_sorted['Formatted_Value'] = df_sorted['b_group'] + ' (' + df_sorted['Adjusted_Total'].round(2).astype(str) + ')'
        df_sorted['Rank'] = df_sorted.groupby('a_group').cumcount() + 1
        df_ranked = df_sorted.pivot(index='Rank', columns='a_group', values='Formatted_Value')
        df_ranked.columns.name = 'A_chain'

        print('\n#-------------------- results per alpha chain, weighted to beta chains --------------------#', file=f)
        print(df_ranked.head(5), file=f)


        # Make a print list for future plotting for highlights

        # top 5 with threshold per alpha chain
        df_print = df[['a_group', 'b_group', 'base_ID', 'Total']]
        threshold = 1.3
        df_top5_per_group = (
            df_print[df_print['Total'] > threshold]
            .sort_values(by=['a_group', 'Total'], ascending=[True, False])
            .groupby('a_group')
            .head(5)
        )

        # bottom 5 with threshold per alpha chain 
        df_bottom5_per_group = (
            df_print[df_print['Total'] < threshold]
            .sort_values(by=['a_group', 'Total'], ascending=[True, False])
            .groupby('a_group')
            .tail(5)
        )


        # Weighted alpha based on threshold
        df_weighted = df[['a_group', 'b_group', 'Adjusted_Total', 'Total']].copy()

        filtered_df = df_weighted[(df_weighted["Adjusted_Total"] > 0.1) & (df_weighted['Total'] > threshold)]
        filtered_df = filtered_df.sort_values(by=['a_group', 'Adjusted_Total'], ascending=[True, False])
        
        top_5_per_group = filtered_df.groupby("a_group").head(5)
        print_list_weighted = [
            f"{row['a_group']}_{row['b_group']}"
            for _, row in top_5_per_group.iterrows()
        ]


        print('\n#-------------------- suggested highlight list --------------------#', file=f)
        
        print(f'The following list is the top 5 alpha/beta combination per alpha chain, threshold = {threshold}', file=f)
        print(f"Highlight : {df_top5_per_group['base_ID'].tolist()} \n", file=f)

        print(f'The following list is the bottom 5 alpha/beta combination per alpha chain, threshold = {threshold}', file=f)
        print(f"Highlight : {df_bottom5_per_group['base_ID'].tolist()} \n", file=f)

        print('The following list is from the weighted alpha chains, with a threshold of > 0.1', file=f)
        print(f'Highlight : {print_list_weighted} \n', file=f)



    return df_present



#----------------- Activation ------------------#

# Some non function calcs
selection = []
for metric in metrics:
    selection.append(f'norm_{metric}')

Result_dir = os.path.join(working_dir, 'Combined_results')
os.makedirs(Result_dir, exist_ok=True)

if __name__ == "__main__":
    inf_data_df = Extraction_data(working_dir)
    Norm_data = Normalisation(inf_data_df)
    combined_data_df, csv_file = Calculation_totals(Norm_data, selection, Result_dir)
    coordinates = Coordinate_collection(combined_data_df)

    Plotting_3d(coordinates, Result_dir, plot_view_angle)
    Plotting_1d(combined_data_df, Result_dir, highlight)
    Plotting_HTML(combined_data_df, Result_dir)
    Heatmap_alpha_beta_total_comparison(combined_data_df, Result_dir)

    weighted_heatmap_data = Analytics(combined_data_df, coordinates, Result_dir, selection)
    heatmap_weighted_alpha_chain(weighted_heatmap_data, Result_dir, threshold_val=0.1)


