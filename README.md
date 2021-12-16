# MicrogliaMegaAtlas

The biggest microglia single-cell RNAseq dataset ever(!)

There are two main aims to this project:

1. Generate a single-cell atlas to map the entire transcriptional repertoire of microglia and their response to all possible perturbations (*e.g.* neurodegeneration, traumatic injury, drug exposure, etc.)
2. Test the hypothesis that the microglial response to neurodegenerative injury is similar to the microglial response to traumatic injury.

Two sub-aims of this project are:

1. Use *in silico* cell lineage-tracing (*i.e.* trajectory inference) to identify transcriptional determinants of microglial proliferation and characterize microglia progeny - are they neuroprotective/neurotoxic?
2. Identify compounds/drugs that can alter microglial trajectories towards neuroprotective transcriptional phenotypes.


![UMAP of microglia across tissue and disease](results/Integration_Microglia/tissue_study_umap.png)

## Summary Statistics

Study
Hammond2018 | Somebang2021 | Milich2021 | Hamel2020 | Safaiyan2021 | Lee2021 | Witcher2021
--- | --- | --- | --- | --- | --- | ---
79937 | 77358 | 21221 | 16598 | 10626 | 8941 | 8910

Tissue
cortex | whole_brain | spinal_cord | hippocampus | GreyM | WhiteM
--- | --- | --- | --- | --- | ---
86268 | 79937 | 37819 | 8941 | 6036 | 4590

Injury_model
TBI | Healthy | SCI | Aging | mouseAD
--- | --- | --- | --- | ---
86268 | 79937 | 37819 | 10626 | 8941

Sex
M | F
--- | ---
57864 | 31014

Age
P720 | P600
--- | ---
10626 | 8941

