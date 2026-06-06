# Train_Timetabling_with_Operational_Strategies

This repository contains the supplementary materials for our paper, "Train Timetabling with Rolling Stock Assignment, Short-Turning and Skip-Stop Strategy for a Bidirectional Metro Line." It includes the datasets, model parameters, and Python scripts required to reproduce the main computational results reported in the manuscript.

The repository provides the network topology, static OD matrices, and simulated time-dependent OD demand based on the 30–50–20 station demand profile, together with all key model parameters (e.g., (h_{\min}), (h_{\max}), (C), (\kappa), (s_k), station demand profiles (\phi_i[u]), and OD shares).

In addition to the original scripts used to reproduce the computational experiments, the repository now includes a Framework directory containing a modular implementation of the proposed methodology. The framework is organized into reusable components for data processing, optimization, queue-feedback analysis, and performance evaluation. A data adapter module is also provided to facilitate the integration of alternative network and demand datasets, enabling users to apply the framework to other metro systems with minimal code modifications.

For reproducibility, a script to regenerate Table 7 (Row 1) is available in the data/Reproducibility_script/ folder, and step-by-step instructions for running the experiments are provided in the main README file.
