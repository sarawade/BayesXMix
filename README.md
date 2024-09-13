This repo contains the code and data accompanying the experiements in the paper "Bayesian dependent mixtures: a predictive comparison and survey" by Wade and Inacio (2024) under reveision for Statistical Science:
https://arxiv.org/abs/2307.16298

Data generation for three examples of the paper can be found in the 'data_simulation' folder:
- Example 1: drawbacks of the joint (generative) approach
- Example 2: drawbacks of the conditional (discriminative) approach with depedent atoms and single weights
- Example 3: drawbacks of the conditional (discriminative) approach with depedent weights
  
In addition, four extra examples are explored in the Supplementary Material:
- Example 1a: extending Example 1 with student-t errors
- Example 1b: extending Example 1 with dependent skewed normal errors
- Example 1c: extending Example 1 with dependent bimodal errors
- Implicit example: response and covariates are implicitly defined based on a latent variable

We compare the predictive performance of seven models:
- Joint approach:
  - Joint Dirichlet process mixture ('Joint' folder)
  - Joint Enriched Dirichlet process mixture ('Joint' folder)
- Conditional approach with depedent atoms and single weights
  - Linear dependent Dirichlet process ('single_weights_ddp' folder)
  - Linear dependent Dirichlet process with B-splines basis expansion ('single_weights_ddp' folder)
- Conditional approach with depedent weights
  - Logisitic stick-breaking process with linear dependence in the stick-breaking proportions ('LSBP' folder)
  - Logisitic stick-breaking process with nonlinear dependence (natural spline expansion) in the stick-breaking proportions ('LSBP' folder)
  - Normalized weights ('NormalizedWeights' folder) 
