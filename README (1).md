# Genomic Prediction of Flowering Time in *Arabidopsis thaliana*

### High-Dimensional Regression Using Screening, Dimension Reduction and Genomic Kernel Methods

**Course:** Regression Techniques  
**Project type:** Statistical modelling and genomic prediction  
**Authors:** Abhishek Kumar, Bhagavath Chukka, Tanish Kothari  
**Institution:** Indian Statistical Institute

---

## 1. Project Overview

This project investigates the prediction of flowering time in *Arabidopsis thaliana* using high-dimensional genomic data.

The central statistical challenge is that we have a relatively small number of individuals but an extremely large number of genetic predictors.

Our dataset contains approximately 1,003 Arabidopsis accessions and more than 1.7 million SNP markers. Thus, the number of predictors is much larger than the number of observations:

```math
p \gg n
```


In this setting, ordinary least squares cannot uniquely estimate unrestricted coefficients for all SNPs.

We therefore investigate several approaches to high-dimensional regression, including variable screening, dimension reduction, penalised regression, Bayesian regression, random projections and genomic kernel methods.

The primary objective is to compare their predictive performance and computational requirements.

## 2. Research Objectives

The project addresses the following questions:

1. How can we perform regression when the number of genomic predictors greatly exceeds the number of individuals?
2. Can screening methods identify informative SNPs before model fitting?
3. Can PCA, PLS and random projections reduce dimensionality while retaining useful predictive information?
4. Does aggregating multiple random projections improve prediction stability?
5. Can randomized SVD approximate exact PCA at lower computational cost?
6. Does using the complete genomic relationship structure provide better predictions than screened regression methods?
7. Can targeted genomic kernels improve upon ordinary GBLUP?
8. How do population structure and linkage disequilibrium affect SNP screening?

## 3. Dataset and Response Variables

The analysis uses genomic and flowering-time information for Arabidopsis accessions.

| Description | Number |
|---|---:|
| FT10 accessions | 1,003 |
| FT16 observed accessions | 970 |
| FT16 missing observations | 33 |
| SNPs after quality control | 1,712,884 |
| Non-constant SNPs used for the GRM | 1,712,865 |

Two flowering-time measurements are considered:

- **FT10:** Flowering time measured at 10°C. This is the primary response variable.
- **FT16:** Flowering time measured at 16°C. This is used for secondary analysis.

FT10 was selected for the main comparison because it was available for all 1,003 accessions. FT16 was used to examine whether the relative performance patterns persisted under a second temperature condition.

## 4. Data Cleaning and Preprocessing

The genotype data underwent quality control using PLINK 2.

The main filtering steps were:

| Step | SNPs remaining |
|---|---:|
| Raw VCF | 12,883,854 |
| Retain biallelic sites | 11,769,920 |
| Remove SNPs with more than 10% missing calls | 6,414,265 |
| Retain SNPs with MAF ≥ 0.01 | 1,712,884 |

The genotype information was coded using dosage values 0 and 2. Missing genotype values were imputed using the mean dosage of the corresponding SNP (`mean2` imputation).

Approximately 4.97% of genotype calls were missing before imputation. Phenotype values were not imputed.

The genotype matrix was subsequently centred and scaled for the genomic relationship matrix and PCA calculations.

## 5. Statistical Methods

The project compares several families of high-dimensional regression techniques.

### 5.1 Principal Component Analysis (PCA)

PCA is used to identify major directions of genomic variation.

Starting from the standardised genotype matrix $Z$, we construct the accession-level Gram matrix:

```math
G=ZZ^T
```


We then perform eigendecomposition:

```math
G=U\Lambda U^T
```


The principal components summarise genomic variation and are also used to investigate population structure.

### 5.2 Partial Least Squares (PLS)

Unlike PCA, PLS uses the response variable when constructing components.

Its objective is to find directions in the genotype matrix that have high covariance with flowering time.

This allows PLS to construct response-guided summaries of the original SNP predictors.

### 5.3 Screening Methods

We investigate three screening approaches:

- **SIS:** Sure Independence Screening ranks SNPs according to their marginal association with flowering time and retains the top candidates.
- **RIS:** Randomized Independence Screening selects SNPs probabilistically using association-based inclusion probabilities.
- **ISIS:** Iterative SIS uses fitted-model residuals to screen for additional signals.

We also investigate structure-adjusted SIS, which first adjusts flowering time for the first five genomic PCs before performing SNP screening.

### 5.4 TARP: Targeted Random Projection

TARP combines phenotype-informed screening with sparse random projection.

The main steps are:

1. Screen SNPs using RIS.
2. Project the retained SNP matrix into a lower-dimensional random feature space.
3. Fit Bayesian regression using a Normal–Inverse-Gamma framework.
4. Repeat the procedure with multiple random projections.
5. Aggregate predictions across projection draws.

The aggregated prediction is:

```math
\widehat y_{\mathrm{TARP}}
=
\frac{1}{B}\sum_{b=1}^{B}\widehat y^{(b)}
```


The main experiment uses $B=100$ projections.

### 5.5 Randomized Singular Value Decomposition (RSVD)

RSVD approximates the leading singular-vector subspace of a large matrix using random sketches and power iterations.

It is investigated as a computationally efficient alternative to exact PCA.

### 5.6 Penalised Regression

The screening methods are followed by regularised regression approaches.

**Ridge regression:**

```math
\widehat\beta
=
\arg\min_\beta
\left\{
\frac{\|y-X\beta\|^2}{2n}
+
\frac{\lambda}{2}\|\beta\|_2^2
\right\}
```


**Lasso regression:**

```math
\widehat\beta
=
\arg\min_\beta
\left\{
\frac{\|y-X\beta\|^2}{2n}
+
\lambda\|\beta\|_1
\right\}
```


Ridge shrinks coefficients continuously, while Lasso can set some coefficients exactly to zero.

### 5.7 Genomic Best Linear Unbiased Prediction (GBLUP)

GBLUP models accession-level genetic effects using a genomic relationship matrix.

The mixed-model formulation is:

```math
y=\mathbf{1}\mu+g+\epsilon
```


with:

```math
g\sim N(0,\sigma_g^2K)
```


```math
\epsilon\sim N(0,\sigma_e^2I)
```


The genomic relationship matrix is constructed as:

```math
K=\frac{ZZ^T}{p^*}
```


where $p^*$ is the number of non-constant SNPs.

GBLUP uses genome-wide genetic similarity to predict the genetic effects of unseen accessions.

It is also mathematically connected to Kernel Ridge Regression, allowing prediction through an accession-level kernel system rather than a SNP-level regression system.

### 5.8 Targeted GBLUP

We further investigate whether screening can improve the full-genome GBLUP kernel.

A second kernel is constructed using screened SNPs:

```math
K_{\mathrm{screened}}
=
\frac{Z_SZ_S^T}{p_\gamma}
```


The combined kernel is:

```math
K_w=K_{\mathrm{full}}+wK_{\mathrm{screened}}
```


Here, $w$ controls the additional weight assigned to the screened SNPs.

This approach retains genome-wide information while giving additional emphasis to potentially informative markers.

An additional experiment applies correlation-based pruning to the SIS-selected candidates before constructing the targeted kernel.

## 6. Model Validation and Evaluation

The primary comparison uses ten shared outer 80/20 train-test splits.

Approximately:

- 80% of accessions are used for training.
- 20% are reserved for testing.

Screening, preprocessing and hyperparameter tuning are performed using training data only. Inner cross-validation is used for model selection.

The principal prediction metric is Mean Squared Prediction Error:

```math
MSPE
=
\frac{1}{n_{\mathrm{test}}}
\sum_{i\in T_{\mathrm{test}}}
(y_i-\widehat y_i)^2
```


Lower MSPE indicates better predictive performance.

We also use test $R^2$ relative to a mean-only training baseline:

```math
R^2_{\mathrm{test}}
=
1-
\frac{MSPE_{\mathrm{model}}}
{MSPE_{\mathrm{mean-only}}}
```


The reported mean-only reference MSPE for FT10 is 334.80.

## 7. Main Results

### 7.1 FT10 Prediction Results

The following table summarises important results from the primary FT10 experiment.

| Method | Mean test MSPE |
|---|---:|
| Full-genome GBLUP | 107.03 |
| RIS-targeted GBLUP | 106.459 |
| RIS + Ridge | 113.04 |
| PLS + RIS | 115.18 |
| TARP × 100 | 117.53 |
| Exact PCA | 118.25 |
| Randomized SVD (q = 2) | 119.79 |
| ISIS + Lasso | 119.62 |

**Additional targeted-kernel experiment:**

SIS initially selected approximately 2,000 candidate SNPs. After pruning highly correlated markers using the threshold $|r|>0.9$, approximately 535 SNPs remained.

The resulting SIS-pruned targeted GBLUP achieved a mean MSPE of 104.481 and performed better than ordinary GBLUP in 6 out of 10 splits.

This is an exploratory result rather than evidence of reliable superiority. The main overall leaderboard reports RIS-targeted GBLUP.

### 7.2 TARP Aggregation

A single random projection produced:

```math
MSPE=188.05
```


Aggregating 100 random projections reduced the mean MSPE to:

```math
MSPE=117.53
```


This demonstrates the importance of aggregation in reducing the variability associated with an individual random projection.

### 7.3 PCA and Randomized SVD

Randomized SVD with $q=2$ power iterations captured approximately 99.0% of the exact-PCA energy.

The reported runtimes were:

| Method | Runtime |
|---|---:|
| Exact PCA | 14.335 seconds |
| Randomized SVD, q = 2 | 3.645 seconds |

The randomized approximation was computationally faster in this experiment, with a small difference in predictive MSPE.

### 7.4 Population Structure and Screening

The first five genomic PCs were used in structure-adjusted screening.

The project also investigated how ancestry can influence marginal SNP associations.

A semi-synthetic experiment was conducted using real genotypes and artificially generated responses, for which the causal SNP set was known.

In the ancestry-confounded simulation setting:

| Metric | Ordinary SIS | Structure-adjusted SIS |
|---|---:|---:|
| Exact causal recovery | 0 | 0.346 |
| Locus recovery | 0 | 0.400 |
| Ancestry tagging | 0.766 | 0.068 |

This controlled experiment demonstrates that population structure can mislead marginal screening and that PC-based adjustment can improve recovery under the simulated conditions.

It does not identify causal SNPs for the real flowering-time phenotype.

### 7.5 Secondary Analysis: FT16

FT16 was used to investigate whether the relative method-performance pattern persisted at a second temperature.

| Method | Mean test MSPE |
|---|---:|
| GBLUP | 223.53 |
| Targeted GBLUP | 219.55 |
| TARP × 100 | 239.31 |

A broadly similar ordering was observed.

Since FT16 uses largely overlapping accessions, this is replication across measurement conditions rather than independent-population validation.

## 8. Conclusions and Key Findings

The project leads to several important conclusions.

**1. Dense genomic prediction is highly competitive.**

Full-genome GBLUP achieved a mean FT10 MSPE of 107.03, outperforming the main screening and projection-based alternatives.

This suggests that the aggregate information contained in genome-wide genetic relationships is highly useful for predicting flowering time.

**2. Screening is useful, but does not automatically outperform genome-wide modelling.**

RIS, SIS and ISIS reduce the number of predictors and make subsequent regression more manageable. However, the strongest screening-based pipelines did not outperform the full-genome GBLUP benchmark in the main FT10 comparison.

**3. Aggregation matters in random projection methods.**

The substantial reduction in MSPE from a single TARP projection to the 100-projection aggregate highlights the importance of averaging over multiple random projections.

**4. Randomized SVD offers a computational advantage.**

The q = 2 randomized SVD approximation captured 99% of exact-PCA energy while requiring considerably less runtime in the reported experiment.

**5. Population structure can mislead marginal screening.**

The semi-synthetic experiment demonstrated that ancestry-related associations can cause ordinary SIS to miss genuinely causal simulated SNPs. Structure-adjusted SIS improved recovery in the confounded setting.

**6. Targeted kernels are promising but require further validation.**

Giving additional weight to screened SNPs produced modest improvements in the observed mean MSPE. However, the RIS-targeted improvement was small relative to its paired variability, and the SIS-pruned result was exploratory.

Overall, the results support GBLUP as a strong predictive benchmark for this Arabidopsis accession panel, while screening and projection methods offer alternative approaches when computational efficiency or variable reduction is important.

## 9. Limitations

- The evaluation uses ten overlapping train-test splits, not ten independent populations.
- FT16 is measured on largely the same accession panel as FT10.
- The semi-synthetic experiment demonstrates behaviour under controlled simulation settings and does not establish actual causal FT10 SNPs.
- Screening can be affected by linkage disequilibrium and population structure.
- Improvements from targeted GBLUP are modest and require independent validation.
- Runtime measurements depend on hardware, implementation and experimental settings.
- Predictive associations should not be interpreted as evidence of biological causation.

**Prediction is not the same as association, causal discovery or external validation.**

## 10. References

Selected methodological references:

1. 1001 Genomes Consortium (2016). 1,135 genomes reveal the global pattern of polymorphism in *Arabidopsis thaliana*. *Cell*, 166, 481–491.
2. Fan, J. and Lv, J. (2008). Sure independence screening for ultrahigh dimensional feature space. *Journal of the Royal Statistical Society: Series B*, 70, 849–911.
3. Halko, N., Martinsson, P. G. and Tropp, J. A. (2011). Finding structure with randomness: Probabilistic algorithms for constructing approximate matrix decompositions. *SIAM Review*, 53, 217–288.
4. Meuwissen, T. H. E., Hayes, B. J. and Goddard, M. E. (2001). Prediction of total genetic value using genome-wide dense marker maps. *Genetics*, 157, 1819–1829.
5. Mukhopadhyay, M. and Dunson, D. B. (2020). Targeted random projection for prediction with large-scale data. *Journal of the American Statistical Association*, 115, 1998–2010.
6. VanRaden, P. M. (2008). Efficient methods to compute genomic predictions. *Journal of Dairy Science*, 91, 4414–4423.
7. Wold, S., Sjöström, M. and Eriksson, L. (2001). PLS-regression: A basic tool of chemometrics. *Chemometrics and Intelligent Laboratory Systems*, 58, 109–130.

---

**Project takeaway:** This study compares screening, dimension reduction, Bayesian regression and genomic kernel approaches for high-dimensional genomic prediction. For the analysed Arabidopsis panel, genome-wide GBLUP provides a strong predictive benchmark, while targeted kernels, aggregation and population-structure adjustment offer useful directions for further investigation.
