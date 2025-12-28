# Speciale
Repository for master's thesis written by David Senderovitz and Kasper Niss. 

Consider a time series $(X_t, Y_t)$ where $X_t$ is $d_X$-dimensional and $Y$ is $d_Y$-dimensional. We are interested in the hypothesis 

$$ H_0\ : \quad  X_{t+1} \perp Y_t \mid X_t $$

under the assumption that $(X_t, Y_t)$ is a Markov process. Our test is built upon a doubly robust estimating equation. We calculate the following test statistic

$$ \hat S = \max_{b \in [B]}\sqrt{T - L}\max \\{|\Re\hat\Gamma(b)|, |\Im\hat\Gamma(b)|\\} $$ 

using `lightbgm` to estimate conditional characteristic functions. We then calculate a covariance matrix from $\hat \Gamma(b)$ for each $b$ and use that to construct critical values for the test, $\hat c_\alpha$, for the desired significance level $\alpha$. We then reject the $H_0$ at level $\alpha$ if 

$$  \hat S > \hat c_{1-\alpha}.$$

The code implemented in this repo is designed for testing that specific hypothesis. For further details than those above, we refer to our thesis: _working on link..._.

# File Overview
Relevant files the testing procedure:
- `estimate_test.R`: This R script contains all the relevant code for computing $\hat S$ and the covariance estimate.
- `sim_crit_value.R`: This R script contains the code for calculating the critical value given a covariance estimate.
- `test_hypothesis.R`: This R script contains the code for testing the hypothesis based on a single dataset.
- `sim_rej_rate.R`: This R script contains the code for calculating the rejection based on multiple datasets. This is used for all the generated calibration and power plots.

Relevant files for the data generating procedure:
- `simulate_parameters`: This R script generates parameters for VAR(1) and CIR processes with certain properties.
- `simulate_AR_process.R`: This R script generates VAR(1) processes.
- `CIR_drift_diffusion`: This R script generates the drift and diffusion functions for the CIR process given a parameter $\theta$.
- `simulate_SDE.R`: This R script generates CIR processes. Both a 'full' version and a 'discretized' version.
