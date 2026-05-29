csv_file = '/projects/0/prjs1135/report_Rick/3_Jolada_data/Experiments/Combined_results/combined_orientation_file.csv'

import pandas as pd




df = pd.read_csv(csv_file)
df_sorted = df.sort_values(by=['proper_docking','AF3_confidence_score'], ascending=[True, False])
print(df_sorted.head(5))