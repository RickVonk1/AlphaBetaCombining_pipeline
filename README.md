 

# Sequence generation (pre-pipeline_script) 

The start of the pipeline consists of a script that generates sequences from gene locations.
This script makes use of Stitchr to generate the whole sequence with a provided CDR3 sequence (J. Heather, 2026).
Stitchr return the sequence in parts, of these parts the leader sequence is removed, and the constant domain is replaced with a provided murine constant domain.
The result is a csv file containing the full sequences, aswel as an updated version of ‘Alpha_Beta_single_chains.csv’ in the 1_pre-AF3 section. 

# AlphaFold preparation (1_Pre-AF3) 

The next step is preparation of directories for AlphaFold.
This script makes a csv file within the newly made directories for AlphaFold's use.
This csv file contains the full sequence of the TCR:pMHC complex with the appropriate chain ID, the MHC sequence was obtained by use of PANDORA (Rademaker et al., 2025).
The format of the directories made is based on previous work of a AlphaFold TCR pipeline by Dario Marzella for the RadboudUMC (D. Marzella, 2026). 

# AlphaFold structure generation (2_AF3) 

This step of the pipeline is a reduced version of the AlphaFold pipeline for the TCR:pMHC complex by D. Marzella.
This reduced version does not use the 3rd and final step in the standard pipeline but only generates the structures and renumbers them to IMGT numbering. 

# Geometric Filtering and Orientation Scoring (3_csv-data) 

To filter out poorly generated structures, the two quality identifier checks are applied.
The first script filters out docking orientation through the use of SwiftTCR utils code.
This script also extracts the AlphaFold confidence score per model.
The second script filters out models that are reversely docked.
Additionally, this script provides small analysis on the quality of the AlphaFold run and tars the AlphaFold directory to reduce project space. 

# Energy and Prodigy calculation (4_haddock) 

The next scripts are the creation of Haddock files, activation of haddock files and analysis of haddock files.
First is the generation of haddock files by combining TCR alpha and beta into a single chain ID and combining alpha-MHC and b2m into a single MHC chain.
The haddock files are combined into a runs file and activated, for memory efficiency runs are combined to groups of 50 files.
Post haddock activation several metrics are extracted from haddock and the AlphaFold output and combined into a single csv file.
This final script additionally generated a plot with the results of this specific run as well as an analytics txt file containing minor analysis.
Lastly this final script also tars the haddock folder to reduce project space. 

# Comparative Analytics and Visualizations (5_Plotting_Analysis) 

With the completion of all alpha chain runs, generated results of individual runs can be combined and compared.
This is done for both the model quality output of AlphaFold and the discussed metrics desired form AlphaFold and Haddock.
The experiment wide analysis is then returned into several plots and an analytics files containing a: top 5, weighted top5 and bottom 5 lists calculated per alpha and per beta chain. 

# Utils 

This contains scripts either used multiple times within other sections or miscelanious scripts not required for the proper functioning of the pipeline. 

# Sources
#### D. Marzella, AF3_TCRpMHC_snellius, (2026), GitHub repository, X-lab-3D/AF3_TCRpMHC_snellius: SPed up TCRpMHC AlphaFold3 pipeline for snellius
#### J. Heather, stitchr, (12-03-2026), GitHub repository, https://github.com/JamieHeather/stitchr
#### Rademaker, D. T., Parizi, F. M., Van Vreeswijk, M., Eerden, S., Marzella, D. F., & Xue, L. C. (2025). Predicting reverse-bound peptide conformations in MHC Class II with PANDORA. Frontiers in Immunology, 16, 1525576. https://doi.org/10.3389/fimmu.2025.1525576  


