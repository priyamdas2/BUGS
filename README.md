# 🚀 BUGS: Bayesian Univariate-Guided Shrinkage

A scalable Bayesian framework for high-dimensional regression that integrates **marginal signal guidance** with **global–local shrinkage** to achieve:

- strong control of false discoveries  
- near-perfect signal recovery  
- scalability to ultra-high dimensions (\(p \sim 10^6\))  

---

## 🔍 Motivation

Classical shrinkage priors such as the regularized horseshoe (RHS) treat all variables symmetrically. However, in many high-dimensional settings:

- marginal signals contain useful information  
- ignoring them leads to inflated false discoveries  
- computational cost becomes prohibitive  

**BUGS addresses these challenges by guiding shrinkage using marginal evidence while preserving Bayesian rigor.**

---

## 🧠 Core Idea

BUGS augments the regularized horseshoe prior with a **guidance mechanism** based on marginal scores:

\[
\tilde{z}_j^{*} \quad \longrightarrow \quad \exp(\eta \tilde{z}_j^{*})
\]

which adaptively modulates shrinkage.

---

## ⚙️ Model Structure

### Hierarchical formulation

\[
\beta_j \mid \lambda_j, \tau, c, \eta, \sigma^2
\sim \mathcal{N}(0, \sigma^2 \tilde{\kappa}_j^2)
\]

where

\[
\tilde{\kappa}_j^2 =
\frac{c^2 \tau^2 \lambda_j^2 \exp(\eta \tilde{z}_j^{*})}
{c^2 + \tau^2 \lambda_j^2 \exp(\eta \tilde{z}_j^{*})}
\]

and

- \( \lambda_j \): local scale  
- \( \tau \): global scale  
- \( c^2 \): slab parameter (RHS)  
- \( \eta \): guidance strength  

---

### Intuition

- Large \( \tilde{z}_j^{*} \) → **less shrinkage**  
- Small \( \tilde{z}_j^{*} \) → **strong shrinkage**  
- \( \eta = 0 \) → reduces to standard RHS  

---

## 📌 Model Diagram

![BUGS structure](images/bugs_structure.png)

---

## 📈 Effect of Guidance

![Guidance effect](images/guidance_effect.png)

- Guidance selectively inflates variance for promising variables  
- Produces stronger separation between signal and noise  

---

## 📊 Posterior Behavior

![Posterior densities](images/posterior_densities.png)

- Signals: well-separated posterior mass  
- Noise: sharply concentrated at zero  

---

## 🔁 MCMC Diagnostics

![Trace plots](images/trace_plots.png)

- Stable mixing across parameters  
- No pathological behavior in high dimensions  

---

## 🧪 Signal vs Noise Separation

![Marginal + local scales + exceedance](images/marginal_local_exceedance.png)

- True variables: high posterior support  
- Noise variables: aggressively shrunk  

---

## ⚡ BUGS-Active: Scalable Inference

To scale to ultra-high dimensions, we introduce **BUGS-Active**:

- restrict updates to a data-adaptive active set  
- reduce complexity from \(O(p)\) to \(O(|A_n|)\)  

---

### Performance with different active sizes

![Active 20](images/bugs_active_20.png)  
![Active 50](images/bugs_active_50.png)  
![Active 100](images/bugs_active_100.png)

---

### Full BUGS vs Active Approximation

![Full BUGS](images/bugs_full.png)

- Maintains high accuracy  
- Dramatically reduces computation time  

---

## 🧬 Real Data Application: DNA Methylation

### Posterior effects (top CpGs)

![CpG effects](images/cpg_effects.png)

---

### Posterior inclusion probabilities

![CpG probabilities](images/cpg_probabilities.png)

---

### Predictive performance

![Predicted vs true](images/predicted_vs_true.png)

---

### Distribution across age groups

![Age violin](images/age_violin.png)

---

### Threshold stability

![Threshold curves](images/threshold_curves.png)

---

## 🏆 Key Advantages

### Statistical
- Near-perfect signal recovery  
- Strong FDR control  
- Adaptive shrinkage driven by data  

### Computational
- Scales to \(p \sim 10^6\)  
- Active-set MCMC  
- Efficient posterior updates  

### Practical
- No manual threshold tuning  
- Interpretable posterior summaries  
- Robust to weak signals  

---

## 📌 Summary

BUGS provides a **principled and scalable Bayesian solution** for high-dimensional variable selection by:

- integrating marginal information into shrinkage  
- preserving global–local structure  
- enabling reliable inference at massive scale  

---

## 📄 Reference

(Coming soon — manuscript under submission)

---

## ⭐ If you find this useful

Please consider starring the repository and citing the work once available.
