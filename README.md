# pCSMC-AS
Implementation of the pCSMC‑AS algorithms from the examples in Diz‑Lois Palomares et al. (2026) for both parameter and latent state inference.

The repository is organized into two main folders:

- Synthetic_data: Contains the examples on synthetic data. One example estimates only the $\theta$ parameters (i.e., $R_t$​ is assumed known), and a more demanding example estimates both $\theta$ and $R_t$.​. In both cases, each algorithm has its own main script, which assumes that the synthetic data are available as syn_data.RDS. This file can be generated using the SIR_model.R script.
- Covid_example: Where we apply the pCSMC-AS algorithm to estimate the reproductive number and infectivity profile from  daily hospital covid-19 admissions during the arrival of the Alpha SARS-CoV-2 variant to Norway between February $1^{st}$ and March $15^{th}$ 2021. Data was originally gathered by the Beredskapsregisteret for COVID-19 (Beredt C19) and made it available for use.

  
