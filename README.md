# HEROES-Proteomics

## Workflow:
1. Breakthrough_preprocess.Rmd
    * Produces `breakthrough_sample_soma.rds` and `breakthrough_proteins.rds`
2. PCLR_breakthrough.slurm
    * Uses PCLR_breakthrough.R
    * Writes into Breakthrough_PCLR_min and Breakthrough_PCLR_1SE
        * These directories must be created ahead of time
3. Breakthrough_postproc.Rmd
    * Produces `Breakthrough_PCLR_min-proc.csv` and `Breakthrough_PCLR_1SE-proc.csv`
4. Overrepresentation Analysis (ORA).
    * I'll update this when the reinfection data comes around. It's more complicated.
  
