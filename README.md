# pCSMC-AS
Implementation of the pCSMC‑AS algorithms from the examples in Diz‑Lois Palomares et al. (2026) for both parameter and latent state inference.

The repository is organized into two main folders:

- Synthetic_data: Which contains the examples on synthetic data. One where only the $\theta$ parameters are estimated (i.e. $R_t$ is assumed known). And a more demanding one where both $\theta$ and $R_t$ are estimated.
- Covid_example: Where we apply the pCSMC-AS algorithm to estimate the reproductive number and infectivity profile from  daily hospital covid-19 admissions during the arrival of the Alpha SARS-CoV-2 variant to Norway between February $1^{st}$ and March $15^{th}$ 2021. Data was originally gathered by the Beredskapsregisteret for COVID-19 (Beredt C19) and made it available for use. 
