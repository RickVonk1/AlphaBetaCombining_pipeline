# TCR:pMHC Modeling Pipeline

This pipeline describes the generation, structural modeling, filtering, and thermodynamic analysis of T-Cell Receptor to peptide-Major Histocompatibility Complex (TCR:pMHC) structures for the identification of potential binding pairs.

---

## Script Architecture

Every Python script within this pipeline follows a standard internal structure:
*   **Explanation Sections:** Quickly explain the exact function of the script and outline which input variables are modifiable.
*   **Input Sections:** Clearly isolate the variables that must be adjusted by the user before executing the script.
*   **Usage:** python /path/to/script/
---

## Pipeline Stages

### 0. Sequence Generation (`0_pre-pipeline_script`)
The pipeline begins by generating full-length sequences from gene locations.
*   Uses **Stitchr** (Heather, 2026) to generate the full sequence using a provided CDR3 sequence.
*   Modifies the raw Stitchr output by removing the leader sequence and replacing the native constant domain with a provided murine constant domain.
*   **Output:** A CSV file containing the finalized full sequences, alongside an updated version of `Alpha_Beta_single_chains.csv` deposited into the `1_pre-AF3` directory.

### 1. AlphaFold Preparation (`1_Pre-AF3`)
Prepares the directory structures and input files required for structural prediction.
*   Generates a targeted CSV file within newly created directories for AlphaFold's use.
*   Integrates the full sequence of the TCR:pMHC complex mapped to appropriate chain IDs. The pMHC sequence is obtained via **PANDORA** (Rademaker et al., 2025). Directory formatting is modeled after the AlphaFold TCR pipeline developed by Dario Marzella for RadboudUMC (Marzella, 2026).

### 2. AlphaFold Structure Generation (`2_AF3`)
Executes a streamlined version of the AlphaFold3 TCRpMHC pipeline.
*   This is a shortend version of the standard D. Marzella pipeline that skips the third final refinement step, focusing purely on generating the structural models and renumbering them according to IMGT unique numbering system.

### 3. Geometric Filtering & Orientation Scoring (`3_csv-data`)
Applies structural quality filters to weed out poorly generated models.
*   **Script 1:** Filters out invalid docking orientations using modified `SwiftTCR-utils` logic and extracts AlphaFold confidence scores per model.
*   **Script 2:** Flags and filters out models displaying reversed docking orientations. It also generates a brief run-quality summary and compresses (`.tar`) the massive AlphaFold directories to optimize storage space.

### 4. Energy & Prodigy Calculation (`4_haddock`)
Automates the setup, execution, and analysis of HADDOCK refinement runs.
*   **Preparation:** Generates HADDOCK-compatible files by merging the TCR alpha/beta chains into a single chain ID, and merging the alpha-MHC/b2m components into a single MHC chain ID.
*   **Execution:** Combines HADDOCK files into unified run scripts. To optimize memory efficiency, runs are batched into groups of 50.
*   **Analysis:** Extracts HADDOCK thermodynamic metrics alongside AlphaFold outputs into a centralized CSV. The final script generates a run-specific performance plot, compiles an analytics summary text file, and archives the HADDOCK folder into a `.tar` file.

### 5. Comparative Analytics & Visualizations (`5_Plotting_Analysis`)
Final script of the pipeline, aggregates individual alpha chain runs into a comprehensive experiment-wide analysis.
*   Compares overall AlphaFold model quality alongside combined metrics from AlphaFold and HADDOCK.
*   **Output:** Returns multiple comparative plots and an analytics summary file detailing the top 5, weighted top 5, and bottom 5 structural models, calculated individually per alpha and per beta chain.

### Utils
Contains shared utility scripts utilized across multiple pipeline sections, alongside miscellaneous helper scripts not strictly required for the core pipeline workflow.

---

## Sources & Citations

*   **D. Marzella.** *AF3_TCRpMHC_snellius*, (2026), GitHub repository, [X-lab-3D/AF3_TCRpMHC_snellius](https://github.com/X-lab-3D/AF3_TCRpMHC_snellius) (Accelerated TCRpMHC AlphaFold3 pipeline for Snellius).
*   **J. Heather.** *stitchr*, (2026), GitHub repository, [JamieHeather/stitchr](https://github.com/JamieHeather/stitchr).
*   **Rademaker, D. T., Parizi, F. M., Van Vreeswijk, M., Eerden, S., Marzella, D. F., & Xue, L. C.** (2025). *Predicting reverse-bound peptide conformations in MHC Class II with PANDORA*. Frontiers in Immunology, 16, 1525576. https://doi.org/10.3389/fimmu.2025.1525576