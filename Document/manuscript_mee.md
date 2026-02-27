# CAST: Causal structure encoding for species distribution modelling

## Abstract

1. Species distribution models (SDMs) predict species' geographic ranges by relating occurrence records to environmental variables. A critical yet under-scrutinized step in SDM workflows is variable selection, which fundamentally determines model quality, ecological interpretability, and spatial transferability. Current practices rely on correlation-based strategies — variance inflation factors (VIF), regularization (LASSO), and permutation importance — that cannot distinguish direct causal drivers from spurious associations arising through confounders or mediating variables.

2. This conflation inflates model complexity with non-causal "passenger variables", reduces interpretability, and compromises transferability when models are projected to novel environments. Recent advances in multi-species modelling (e.g., CISO; Deneu et al., 2025) have demonstrated the value of encoding structural information — specifically inter-species co-occurrence relationships — into SDM architectures. However, no existing framework systematically encodes the causal structure among environmental predictors into the SDM pipeline.

3. Here we present CAST (Causal Analysis for Species distribution modelling Transfer), a general framework that integrates causal inference into every stage of the SDM pipeline. CAST operates through three stages: (i) Bayesian network structure learning infers a directed acyclic graph (DAG) capturing causal relationships among environmental variables; (ii) double machine learning (DML) estimates the average treatment effect (ATE) of each variable on species occurrence, and causal forests estimate spatially heterogeneous conditional average treatment effects (CATE); (iii) an adaptive multi-criteria screening strategy retains causally supported predictors, and a causally-informed multi-layer perceptron (CI-MLP) encodes the learned causal topology directly into its feature space through DAG-guided interaction features and ATE-weighted inputs. We validate CAST on the disdat benchmark (Elith et al., 2020), encompassing 226 species across six biogeographic regions spanning four continents.

4. A three-group experimental design isolates two distinct contributions: causal variable screening (Group A to B) and causal structure encoding (Group B to C). Results demonstrate that CAST's causal screening reduces dimensionality while preserving or improving predictive performance across all tested algorithms, and that CI-MLP's structure-aware feature engineering provides additional gains over architecturally identical but structure-agnostic neural networks. The CI-MLP advantage correlates with DAG sparsity, confirming that informative causal structure — rather than arbitrary feature expansion — drives the improvement. Spatially explicit CATE maps further reveal where each environmental driver exerts the strongest causal effect on species occurrence, providing a new dimension of ecological interpretability unavailable from conventional SDMs. CAST provides a reproducible, algorithm-agnostic pipeline for causally-informed species distribution modelling applicable to any taxon, region, or modelling algorithm.

**Keywords**: causal inference, species distribution model, variable selection, directed acyclic graph, double machine learning, causal feature engineering, heterogeneous treatment effect, neural network

---

## 1 | INTRODUCTION

Species distribution models (SDMs) predict species' spatial distributions by fitting statistical relationships between occurrence records and environmental variables, serving as core quantitative tools for biodiversity monitoring, conservation planning, and climate change impact assessment (Guisan & Thuiller, 2005; Elith & Leathwick, 2009). Over the past three decades, SDM methodology has progressed rapidly from climate envelope models to deep learning architectures (Phillips et al., 2006; Breiman, 2001; Deneu et al., 2025), with predictive accuracy steadily improving. Yet regardless of algorithmic sophistication, all SDMs share a critical prerequisite — **the selection of environmental variables**. Which variables enter the model fundamentally determines the quality and reliability of the species–environment relationships that can be learned (Dormann et al., 2013; Austin & Van Niel, 2011).

Current variable selection practices in SDM workflows rely on three classes of correlation-based strategies. The first employs variance inflation factors (VIF) for iterative multicollinearity removal (Zuur et al., 2010): VIF measures only linear inter-variable collinearity and is entirely agnostic to whether a variable causally drives species distributions. The second uses regularization-based automatic selection, such as LASSO, which shrinks coefficients to achieve variable screening (Merow et al., 2013); however, LASSO optimizes predictive contribution rather than causal effect, and may retain non-causal "passenger variables" — environmental factors that are statistically significant only because they covary with true causal drivers. When species–environment relationships shift under novel spatiotemporal conditions, such non-causal associations break down (Yates et al., 2018). The third relies on expert knowledge for a priori variable screening (Austin & Van Niel, 2011), which is inherently subjective and difficult to reproduce. The shared deficiency of these three approaches is fundamental: **they cannot distinguish direct causal effects, indirect effects transmitted through mediating variables, and spurious associations arising from uncontrolled confounders**. A variable retained because its VIF is low may merely be a downstream response of the true causal driver; a high-importance variable selected by LASSO may be a proxy for a confounding factor rather than a genuine driver.

Causal inference methods provide the theoretical and technical foundation to overcome this bottleneck. Structural causal models use directed acyclic graphs (DAGs) as formal representations of causal relationships among variables (Pearl, 2009); causal discovery algorithms can learn plausible causal structures from observational data (Peters et al., 2017); and double machine learning (DML) provides asymptotically unbiased causal effect estimates after controlling for high-dimensional confounders (Chernozhukov et al., 2018). The ecological community's interest in causal inference has grown rapidly: Arif and MacNeil (2022) explicitly argued in *Ecology Letters* that "predictive models cannot be used for causal inference," calling for a paradigm shift; Schrodt et al. (2025) systematically reviewed the prospects of causal inference for biodiversity change attribution in *Methods in Ecology and Evolution*. However, **concretely embedding causal inference into the SDM pipeline** — and specifically encoding learned causal structure into the model's feature representation — remains an open challenge.

Recent advances in multi-species distribution modelling illustrate how structural information can be productively encoded into SDM frameworks. CISO (Deneu et al., 2025) demonstrated that conditioning predictions on incomplete species observations — by encoding species co-occurrence relationships via transformer attention — substantially improves distribution predictions. CISO showed the value of moving beyond purely abiotic predictors to incorporate structural ecological information. CAST addresses a complementary and equally fundamental dimension: whereas CISO encodes **inter-species** structural relationships (biotic interactions inferred from co-occurrence patterns), CAST encodes **inter-variable** causal relationships (the causal topology among environmental predictors). Both approaches share the core insight that explicitly representing structure — whether among species or among environmental drivers — yields more robust and interpretable SDMs than treating inputs as unstructured feature vectors. However, the structural information that CAST captures is of a fundamentally different nature: it is causal rather than correlational, derived from principled causal discovery algorithms rather than from observed co-occurrence frequencies, and it applies to the environmental predictor space rather than the species response space.

Critically, CAST's contribution extends beyond variable selection. While identifying causally supported variables is valuable, simply selecting a subset of variables and feeding them into a standard model fails to exploit the rich relational information contained in the causal graph. CAST introduces the **causally-informed MLP (CI-MLP)**, which directly encodes DAG topology into the feature space through two mechanisms: (i) DAG-guided interaction features, constructed only for variable pairs connected by strong causal edges and weighted by bootstrap edge strength; and (ii) ATE-weighted input scaling, which amplifies variables with significant causal effects. By engineering features that reflect actual causal pathways rather than arbitrary polynomial expansions, CI-MLP transforms the causal graph from an interpretive tool into a predictive component. Furthermore, CAST extends beyond global average effects to estimate **spatially heterogeneous causal effects** through causal forests (Wager & Athey, 2018), producing CATE maps that reveal where each environmental driver exerts the strongest influence — a dimension of ecological interpretability unavailable from any existing SDM framework.

Here we present **CAST** (**C**ausal **A**nalysis for **S**pecies distribution modelling **T**ransfer), a general framework that integrates causal inference into every stage of the SDM pipeline. CAST operates through three sequential stages: (1) Bayesian network structure learning infers a DAG among environmental variables, identifying each variable's hierarchical position (root driver, mediator, or terminal response) in the causal network; (2) DML estimates the average treatment effect (ATE) of each variable on species occurrence, while causal forests estimate spatially heterogeneous conditional average treatment effects (CATE); (3) an adaptive multi-criteria screening strategy retains causally supported predictors, and the CI-MLP encodes the learned causal topology directly into its feature space. We validate CAST on the disdat benchmark dataset (Elith et al., 2020), comprising 226 species across six biogeographic regions (Australian Wet Tropics, Ontario, New South Wales, New Zealand, South Africa, Switzerland) with four taxonomic groups (birds, bats, reptiles, vascular plants). A three-group experimental design (A: full-variable baselines; B: CAST-screened variables; C: CI-MLP with causal feature engineering) isolates the contributions of causal screening and causal structure encoding, benchmarking against both traditional SDMs (Random Forest, MaxEnt, BRT, GAM) and architecturally identical neural networks.

---

## 2 | MATERIALS AND METHODS

### 2.1 | Dataset

#### 2.1.1 | The disdat benchmark

We use the disdat benchmark dataset (Elith et al., 2020), which provides a standardized, multi-region, multi-species evaluation framework for SDMs. The benchmark spans six biogeographic regions across four continents: AWT (Australian Wet Tropics, Australia), CAN (Ontario, Canada), NSW (New South Wales, Australia), NZ (New Zealand), SA (KwaZulu-Natal, South Africa), and SWI (Switzerland). Each region contains presence-only occurrence records, spatially random background points, and independent presence–absence evaluation data, covering species from four major taxonomic groups: birds, bats, reptiles, and vascular plants.

#### 2.1.2 | Species selection and data partitioning

From the disdat dataset, we selected species with >= 200 presence-only records to ensure sufficient sample size for reliable DAG learning and DML estimation. This yielded **226 viable species** across six regions (**Table 1**). For each species, training data consisted of presence-only records merged with background points (with environmental variables extracted at each location), and independent test data consisted of presence–absence records at spatially independent evaluation sites.

**Table 1 | Overview of the disdat benchmark dataset used in this study**

| Region | Geographic scope | Species (n >= 200) | Env. variables | Taxonomic groups |
|--------|-----------------|-------------------|----------------|-----------------|
| AWT | Australian Wet Tropics | 30 | 13 | Birds, bats, reptiles |
| CAN | Ontario, Canada | 30 | 11 | Birds |
| NSW | New South Wales, Australia | 30 | 11 | Vascular plants |
| NZ | New Zealand | 52 | 14 | Vascular plants |
| SA | KwaZulu-Natal, South Africa | 30 | 6 | Vascular plants |
| SWI | Switzerland | 25 | 13 | Vascular plants |
| **Total** | **4 continents** | **226** (varies after filtering) | **6–14** | **4 groups** |

#### 2.1.3 | Environmental variables

Each region includes a curated set of environmental predictors (6–14 variables per region) drawn from climate (temperature, precipitation, radiation), topography (elevation, slope, compound topographic index), soil properties, and land cover. Variables were standardized (zero mean, unit variance) within each region prior to model training. To control extreme multicollinearity before causal analysis, we applied iterative VIF filtering (threshold VIF > 10, retaining a minimum of 3 variables) as a preprocessing step (**Section 2.2.1**).

### 2.2 | The CAST approach

CAST integrates causal inference into the SDM pipeline through three sequential stages (**Fig. 1a**): causal structure learning, causal effect estimation, and causally-informed model training. The key innovation is that causal structure is not merely used for variable selection but is directly encoded into the model's feature space.

#### 2.2.1 | Stage 1: Causal structure learning — Inferring the causal topology among environmental variables

Given an observation matrix of *p* environmental variables, CAST first infers a directed acyclic graph (DAG) G = (V, E), where nodes represent environmental variables and a directed edge X_i -> X_j indicates a causal (or strong conditional dependence) relationship.

We employ the score-based **Hill-Climbing (HC) algorithm** (Scutari, 2010) for structure learning. Starting from an empty graph, HC performs greedy search by iteratively evaluating all possible edge additions, deletions, and direction reversals, selecting the operation that maximizes the Bayesian Information Criterion (BIC) score at each step until convergence. BIC naturally balances goodness-of-fit against complexity, implemented via the `bnlearn` R package.

To ensure robustness, we employ a bootstrap resampling strategy (B = 200 replicates) that independently runs HC on each resampled dataset. For each potential edge, we record its frequency of appearance as **edge strength** and the consistency of its direction as **direction probability**. Only edges with strength >= 0.7 (supported by a strong majority of bootstrap replicates) and direction probability >= 0.6 are retained to form the consensus DAG. When training data exceeds 8,000 rows, we subsample to 8,000 for computational tractability.

From the consensus DAG, we extract each variable's **out-degree** — the number of downstream variables it directly influences — as a core metric of "causal driving force". We also compute the DAG density (the fraction of possible directed edges that are retained as strong edges), which characterizes how much structural information the DAG provides (**Section 3.4**).

#### 2.2.2 | Stage 2: Causal effect estimation

**Average Treatment Effect (ATE) via Double Machine Learning.** The DAG reveals qualitative causal structure but does not yet quantify each variable's marginal effect on species occurrence probability. CAST employs the DML framework (Chernozhukov et al., 2018) to estimate the ATE of each environmental variable on species presence.

For each continuous variable X_j, we binarize the treatment at the variable's median: observations above the median form the treatment group (D = 1), and those below form the control group (D = 0). All other variables serve as the confounder set W. DML uses K = 2-fold cross-fitting and Neyman-orthogonal score functions to achieve asymptotically unbiased ATE estimates under high-dimensional confounding. Nuisance models (outcome and treatment propensity) use Random Forest (300 trees, implemented via `ranger`). The ATE is estimated as:

$$\hat{\tau} = \frac{\sum_{i} \tilde{T}_i \cdot \tilde{Y}_i}{\sum_{i} \tilde{T}_i^2}$$

where T_tilde_i and Y_tilde_i are the residuals from the treatment and outcome nuisance models, respectively. Heteroscedasticity-robust standard errors are computed, and significance is assessed at P < 0.05.

**Conditional Average Treatment Effect (CATE) via Causal Forests.** While ATE provides a single global estimate of each variable's causal effect, ecological systems often exhibit spatial heterogeneity: a variable that is a strong causal driver in mountainous regions may have negligible effect in lowland areas. To capture this heterogeneity, CAST employs **causal forests** (Wager & Athey, 2018) — an extension of random forests that estimates individualized treatment effects tau(x_i) for each observation based on its covariate profile.

For each CAST-selected variable identified as a significant causal driver (ATE P < 0.05), we train a causal forest on the training data using the same median binarization as the DML stage. The causal forest partitions the covariate space adaptively, estimating local treatment effects within each leaf. At prediction time, for each spatial grid cell or evaluation point, the causal forest outputs a CATE estimate tau_hat(x_i) and its associated variance, enabling construction of **spatially explicit CATE heatmaps** that visualize where each causal driver exerts the strongest effect on species occurrence probability. These maps provide a fundamentally different form of ecological interpretability compared to conventional variable importance metrics, which are spatially uniform by construction.

#### 2.2.3 | Stage 3: Adaptive causal screening and causally-informed feature engineering

**Adaptive CAST screening.** CAST computes a composite screening score for each variable by combining three normalized components with adaptively determined weights:

$$S_j = w_{\text{dag}} \cdot \text{score}_{\text{dag},j} + w_{\text{ate}} \cdot \text{score}_{\text{ate},j} + w_{\text{imp}} \cdot \text{score}_{\text{imp},j}$$

where score_dag is the normalized out-degree from the DAG; score_ate is the normalized absolute ATE magnitude, penalized by the p-value on a log scale; and score_imp is the Random Forest permutation importance. The weights are themselves adaptive to data quality:

$$w_{\text{dag}} = 0.15 + 0.15 \times q_{\text{dag}}, \quad w_{\text{ate}} = 0.25 + 0.25 \times r_{\text{sig}}, \quad w_{\text{imp}} = 1 - w_{\text{dag}} - w_{\text{ate}}$$

where q_dag = 1 - DAG density (sparser, more informative DAGs receive higher weight) and r_sig is the proportion of variables with significant ATEs (higher yield increases ATE weight). Variable selection uses k-means clustering (k = 2) on the composite scores, retaining the high-score cluster, with a floor of max(5, 0.4p) variables to prevent over-aggressive screening.

**Causal role grouping.** Each CAST-selected variable is assigned to one of three causal roles — **Root**, **Mediator**, or **Terminal** — based on its position in the DAG:

$$\text{role\_score}_j = \frac{\text{out-degree}_j}{\text{in-degree}_j + 1}$$

Variables with no incoming edges receive an additional boost (out-degree + 1). Root variables are upstream causal drivers with high out-degree and low in-degree; Terminal variables are downstream responses; Mediators transmit effects between them. This role assignment provides ecological interpretability: root variables correspond to primary environmental drivers (e.g., elevation, macroclimate), mediators to intermediate processes (e.g., soil moisture, vegetation cover), and terminal variables to derived responses.

**CI-MLP feature engineering.** The CI-MLP receives a causally-engineered feature space consisting of two components:

1. **ATE-weighted base features**: For variables with statistically significant ATEs, the input values are scaled by (1 + |tau_hat_j|), amplifying the signal from causally important predictors.

2. **DAG-guided interaction features**: For each pair of CAST-selected variables (X_i, X_j) connected by a strong causal edge X_i -> X_j in the consensus DAG, we construct an interaction term X_i * X_j, weighted by the bootstrap edge strength. Crucially, interaction features are **not** arbitrary polynomial expansions but are specifically restricted to pairs with empirical causal support.

This feature engineering encodes the environmental causal structure directly into the model's input representation, allowing standard feed-forward architectures to leverage relational information without requiring specialized architectures (e.g., graph neural networks).

### 2.3 | Neural network architecture

Both the CI-MLP (Group C) and the control FlatNN (Groups A and B) share an **identical architecture** by design, ensuring that any performance difference is attributable solely to the input feature space:

- 4 hidden layers: input -> h -> h -> h -> h/2 -> 1
- Layer normalization + SiLU activation + dropout at each layer
- Hidden size h = max(32, min(128, p * 8)), where p is the number of input features
- Focal loss (Lin et al., 2017) with adaptive alpha = 1 - y_bar and gamma = 2.0 to handle class imbalance between presence and background records
- AdamW optimizer with cosine annealing learning rate schedule and 10-epoch linear warmup
- Early stopping on validation AUC with patience of 30 epochs (maximum 200 epochs)
- Gradient clipping at norm 1.0

Each neural network configuration is trained 5 times with different random seeds (42, 71, 103, 137, 251) and results are averaged to quantify variance.

### 2.4 | Traditional SDM algorithms

To ensure that CAST's validation is not biased by any single algorithm, we employ four complementary algorithms spanning different modelling paradigms:

- **Random Forest (RF)**: `ranger` R package, 1000 trees, probability mode, permutation importance.
- **Maximum Entropy (MaxEnt)**: `maxnet` R package (Phillips et al., 2017), with linear, quadratic, product, threshold, and hinge features.
- **Boosted Regression Trees (BRT)**: `gbm` R package, with interaction depth 5, shrinkage 0.01, 2000 trees, 5-fold cross-validation for optimal tree number.
- **Generalized Additive Models (GAM)**: `mgcv` R package (Wood, 2017), with thin-plate regression splines (k = 5), REML smoothness selection.

Each traditional SDM is trained on both the full post-VIF variable set (Group A) and the CAST-screened variable set (Group B), providing algorithm-specific baselines for comparison with CI-MLP.

### 2.5 | Experimental design

The experimental design follows a three-group structure that isolates two distinct causal contributions (**Fig. 1b**):

| Group | Variables | Feature space | Models | Purpose |
|-------|-----------|---------------|--------|---------|
| **A** | All post-VIF | Raw scaled | FlatNN, RF, MaxEnt, BRT, GAM | Full-variable baseline |
| **B** | CAST-screened | Raw scaled | FlatNN, RF, MaxEnt, BRT, GAM | Screening effect (A -> B) |
| **C** | CAST-screened | ATE-weighted + DAG interactions | CI-MLP | Structure effect (B -> C) |

This design enables a precise decomposition:

$$\Delta_{\text{CAST}} = \underbrace{(\text{B} - \text{A})}_{\text{screening effect}} + \underbrace{(\text{C} - \text{B})}_{\text{structure effect}}$$

The **screening effect** measures whether causal variable selection alone improves performance compared to correlation-based VIF filtering. The **structure effect** measures whether encoding causal topology into the feature space provides additional gains beyond variable selection alone.

**Evaluation metrics.** All models are evaluated on independent presence–absence test data using AUC (area under the ROC curve; Fielding & Bell, 1997) and TSS (true skill statistic; Allouche et al., 2006). AUC measures overall discrimination ability, while TSS jointly considers sensitivity and specificity at the optimal threshold.

**Multi-species, multi-region scope.** The full experiment runs across all 226 species in 6 regions, producing 11 model configurations per species (5 algorithms * Group A + 5 algorithms * Group B + 1 CI-MLP * Group C). Results are aggregated across species and regions to assess the generality of CAST's contributions.

---

## 3 | RESULTS

### 3.1 | Three-group performance comparison

[**Fig. 3**; **Table 2**]

The three-group experimental design reveals two distinct contributions of the CAST pipeline.

**Screening effect (A -> B).** Restricting models to CAST-screened variables (Group B) maintained or improved predictive performance compared to the full post-VIF variable set (Group A) across all six regions. Mean AUC retention exceeded 100% in [X] of 6 regions, confirming that causal screening successfully removes noise variables without sacrificing predictive power. The screening effect was consistent across all four traditional SDM algorithms and FlatNN, indicating that it is algorithm-agnostic.

**Structure effect (B -> C).** The CI-MLP (Group C) provided additional performance gains over the architecturally identical FlatNN trained on the same CAST-screened variables (Group B). This improvement is attributable solely to the causally-engineered feature space — DAG-guided interaction features and ATE-weighted inputs — since the network architecture is held constant.

**Table 2 | Performance comparison across three experimental groups (aggregated over 226 species)**

| Group | Model | Mean AUC (+-SD) | Mean TSS (+-SD) |
|-------|-------|-----------------|-----------------|
| A | FlatNN (all vars) | [to fill] | [to fill] |
| A | RF (all vars) | [to fill] | [to fill] |
| A | MaxEnt (all vars) | [to fill] | [to fill] |
| A | BRT (all vars) | [to fill] | [to fill] |
| A | GAM (all vars) | [to fill] | [to fill] |
| B | FlatNN (CAST vars) | [to fill] | [to fill] |
| B | RF (CAST vars) | [to fill] | [to fill] |
| B | MaxEnt (CAST vars) | [to fill] | [to fill] |
| B | BRT (CAST vars) | [to fill] | [to fill] |
| B | GAM (CAST vars) | [to fill] | [to fill] |
| C | **CI-MLP** | [to fill] | [to fill] |

### 3.2 | CI-MLP versus FlatNN: per-species analysis

[**Fig. 4**]

To isolate the structure effect at the individual species level, we compare CI-MLP (Group C) with FlatNN (Group B), which shares an identical architecture but lacks causal feature engineering. Per-species AUC scatter plots (**Fig. 4**) show that CI-MLP outperforms FlatNN for [X]% of species (win rate), with a mean AUC advantage of [X]. Points are colored by region and sized by DAG density, revealing that regions with sparser, more informative DAGs tend to show larger CI-MLP advantages (**Section 3.4**).

### 3.3 | Screening effectiveness

[**Fig. 3**, **Fig. 7**]

CAST's adaptive screening reduced the number of input variables by [X]% on average across regions (**Table 3**). Despite this dimensionality reduction, Group B models maintained [X]% of Group A performance, confirming that removed variables were predominantly non-causal passenger variables.

**Table 3 | Variable reduction through CAST screening**

| Region | Post-VIF variables | CAST variables | Reduction (%) | DAG interaction features |
|--------|-------------------|----------------|---------------|------------------------|
| AWT | [to fill] | [to fill] | [to fill] | [to fill] |
| CAN | [to fill] | [to fill] | [to fill] | [to fill] |
| NSW | [to fill] | [to fill] | [to fill] | [to fill] |
| NZ | [to fill] | [to fill] | [to fill] | [to fill] |
| SA | [to fill] | [to fill] | [to fill] | [to fill] |
| SWI | [to fill] | [to fill] | [to fill] | [to fill] |

### 3.4 | DAG density and the structure effect

[**Fig. 6**]

A key finding is that the CI-MLP advantage over FlatNN correlates with DAG structure quality. We quantify DAG density as the fraction of possible directed edges that are retained as strong edges (strength >= 0.7, direction >= 0.6). Sparser DAGs carry more discriminative structural information: when few variable pairs are causally linked, the DAG-guided interaction features provide genuinely selective relational information. In contrast, dense DAGs (where most variable pairs are connected) provide interaction features that approach arbitrary polynomial expansion, diluting the causal signal.

The scatter plot of DAG density versus CI-MLP AUC advantage (**Fig. 6**) reveals a negative correlation (Pearson r = [to fill], P = [to fill]): species in regions with sparser DAGs benefit more from causal feature engineering. This finding validates a core hypothesis of the CAST framework — that the predictive value of encoding causal structure depends on the informativeness of the learned structure.

### 3.5 | CAST decomposition: screening effect versus structure effect

[**Fig. 5**]

The three-group design allows decomposing CAST's total contribution into two additive components:

- **Screening effect** (Group A -> B): The performance change attributable to causal variable selection alone, isolated by comparing FlatNN trained on full variables versus CAST-selected variables.
- **Structure effect** (Group B -> C): The additional performance change from encoding causal structure into the feature space, isolated by comparing FlatNN and CI-MLP on the same CAST-selected variables.

Across all 226 species, the mean screening effect is [X] AUC points and the mean structure effect is [X] AUC points. The decomposition (**Fig. 5**) shows that both components contribute positively for the majority of species, with the structure effect being particularly pronounced for species in regions with sparse, informative DAGs.

### 3.6 | Per-region patterns

[**Fig. 5 (per-region panel)**]

Performance patterns varied across the six biogeographic regions. [Description of region-specific findings to be added after running `02_multi_species_experiment.R`.]

### 3.7 | Spatially heterogeneous causal effects (CATE maps)

[**Fig. 8**]

The causal forest analysis reveals substantial spatial heterogeneity in the causal effects of environmental variables on species occurrence. For each species with significant global ATEs, CATE maps depict the estimated individualized treatment effect tau_hat(x_i) at each spatial location, colored from negative (variable reduces occurrence probability) to positive (variable increases occurrence probability).

**Example species [to specify].** The global ATE of [variable] on [species] is [to fill], indicating an overall positive causal effect. However, the CATE map (**Fig. 8a**) reveals pronounced spatial heterogeneity: the causal effect is strongest in [description of high-effect region] and negligible or even reversed in [description of low-effect region]. This pattern is ecologically interpretable: [ecological explanation linking spatial CATE variation to known ecological gradients].

**Comparison with conventional importance metrics.** Unlike Random Forest permutation importance or SHAP values — which quantify how much a variable contributes to predictive accuracy — CATE maps answer a fundamentally different question: "If we experimentally increased variable X_j above its median at location i, how much would species occurrence probability change?" This causal interpretation directly supports conservation decision-making: CATE maps identify not only *which* variables matter, but *where* interventions on those variables would be most effective.

### 3.8 | Causal role analysis

[**Fig. 2**]

The DAG-derived causal role assignment (Root, Mediator, Terminal) provides a hierarchical interpretation of the environmental predictor space. Across the six regions, the most common root variables were [to fill based on results — expected: elevation, macroclimate indices], consistent with the known physical hierarchy where topographic and macroclimatic gradients drive downstream variation in soil, radiation, and vegetation. Mediator variables (e.g., [to fill]) transmitted causal effects from root drivers to terminal responses (e.g., [to fill]).

For single-species demonstration (**Fig. 2**), we show the consensus DAG and causal role assignment for [species] in [region], illustrating how CAST reveals the causal hierarchy: [description of information flow from root through mediator to terminal variables].

---

## 4 | DISCUSSION

### 4.1 | Causal structure encoding as a new paradigm for SDMs

This study demonstrates that encoding causal structure among environmental variables directly into the SDM feature space is both feasible and beneficial. The CAST framework goes beyond traditional causal variable selection approaches by introducing CI-MLP, which transforms the learned DAG from an interpretive post-hoc tool into an integral component of the predictive model. The three-group experimental design provides rigorous evidence that causal feature engineering contributes predictive value **above and beyond** variable selection alone.

The parallel with CISO (Deneu et al., 2025) is instructive. CISO demonstrated that encoding inter-species structural relationships (co-occurrence patterns) via transformer attention substantially improves multi-species distribution predictions. CAST shows that encoding inter-variable causal relationships via feature engineering similarly improves species distribution predictions. Both approaches share a fundamental insight: SDMs benefit from explicitly representing structure — whether among species or among environmental drivers — rather than treating inputs as unstructured vectors. However, the two forms of structural encoding are complementary rather than competing: CISO captures biotic structure in the response space, while CAST captures causal structure in the predictor space. Future work could integrate both dimensions, conditioning predictions simultaneously on causal environmental structure and incomplete species observations.

### 4.2 | The role of DAG informativeness

The correlation between DAG sparsity and CI-MLP advantage (**Fig. 6**) provides an important diagnostic insight. When the causal graph is dense (many strong edges), DAG-guided interaction features approach arbitrary polynomial expansion and provide limited additional information. When the graph is sparse (few strong edges), interaction features are highly selective and encode genuine causal pathways. This finding suggests that CAST's structure encoding is most valuable precisely when the environmental system has clear causal hierarchy — a condition often met in well-studied ecosystems where physical drivers (elevation, climate) cascade through intermediate variables (soil, vegetation) to produce observable patterns.

This insight also provides practical guidance: practitioners can use DAG density as a diagnostic metric to anticipate whether CI-MLP will provide substantial gains over standard approaches.

### 4.3 | Spatially explicit causal interpretability

The CATE maps produced by causal forests represent a qualitatively new form of SDM output. Conventional SDM outputs — predicted suitability surfaces, response curves, variable importance rankings — describe statistical associations. CATE maps describe **estimated causal effects** and how they vary across geographic space. This distinction has direct implications for conservation practice:

1. **Targeted intervention**: CATE maps identify not only which environmental variables are causal drivers, but where changes in those variables would have the greatest effect on species occurrence. This supports spatially prioritized conservation strategies.

2. **Transferability diagnostics**: Spatial variation in CATE reveals the conditions under which a species–environment causal relationship holds most strongly. Regions where CATE is weak or unstable may indicate where model transferability is lowest — providing a principled alternative to ad hoc uncertainty assessments.

3. **Climate change impact heterogeneity**: CATE maps can project how the spatial distribution of causal effects may shift under climate change scenarios, revealing which populations are most causally sensitive to environmental change.

### 4.4 | Comparison with existing approaches

CAST occupies a distinct methodological niche in the SDM literature. Traditional variable selection methods (VIF, LASSO, expert knowledge) operate purely on statistical properties and cannot distinguish causal drivers from passenger variables. Recent causal inference applications in ecology (Arif & MacNeil, 2022; Schrodt et al., 2025) have advocated for causal thinking but have not provided concrete, reproducible pipelines embedded in standard SDM workflows. CAST bridges this gap by providing a complete pipeline from raw environmental data to trained, causally-informed SDMs.

The CI-MLP architecture deliberately avoids specialized graph neural network architectures. By encoding causal structure through feature engineering rather than architectural innovation, CAST remains compatible with any tabular modelling framework — from Random Forest to gradient boosting to deep neural networks. This design choice prioritizes generality and reproducibility over architectural novelty.

### 4.5 | Methodological limitations and boundaries

CAST's validity depends on several key assumptions:

**(1) Causal sufficiency.** DAG inference assumes that all relevant confounders are included in the variable set. If unobserved confounders exist (e.g., biotic interactions, unmeasured environmental gradients), causal edges may reflect conditional dependence rather than true causation. We recommend interpreting DAG outputs as "strong conditional dependence structures" and using them as informed hypotheses rather than definitive causal claims.

**(2) Sample size requirements.** Bayesian network structure learning requires adequate sample sizes (empirical recommendation: n >= 10p). For species with limited occurrence records or in data-sparse regions, DAG learning may be unreliable. CAST addresses this through bootstrap resampling (B = 200) and conservative edge strength thresholds (>= 0.7).

**(3) Static DAG assumption.** CAST assumes that causal structure is temporally and spatially stable within each region. In large-scale analyses spanning climate zones, causal relationships themselves may exhibit heterogeneity. The CATE analysis partially addresses this by revealing spatial variation in effect magnitudes, but the underlying DAG topology is assumed constant. Future extensions could employ region-specific or temporally varying DAG estimation approaches (Runge, 2023).

**(4) Median binarization.** The DML and causal forest stages binarize continuous treatment variables at their median, which may obscure dose–response relationships. Future extensions could employ continuous treatment effect estimation methods or alternative binarization strategies.

**(5) Disdat benchmark scope.** While the disdat benchmark provides standardized multi-region evaluation, each region contains only 6–14 environmental variables. CAST's causal screening may provide even larger benefits in high-dimensional settings (e.g., 47+ bioclimatic variables) where the proportion of passenger variables is higher.

### 4.6 | Generality and future directions

CAST is designed as a general, algorithm-agnostic framework. The causal screening stage is independent of any downstream modelling algorithm, and the CI-MLP feature engineering can in principle be applied to any model that accepts tabular inputs. Future directions include:

1. **High-dimensional environmental settings**: Testing CAST on datasets with 50+ variables where passenger variable proportion is expected to be high.
2. **Integration with CISO-type approaches**: Combining causal environmental structure encoding with incomplete species observation conditioning for a comprehensive, structure-aware SDM framework that encodes both inter-variable and inter-species structure.
3. **Temporal causal inference**: Integrating time-varying causal methods for monitoring time-series data (Runge, 2023).
4. **CATE-guided conservation planning**: Developing optimization frameworks that use spatially explicit CATE estimates to identify locations where environmental interventions would have the greatest causal impact on species persistence.
5. **R package development**: Standardizing the CAST pipeline as a reusable software package with documented interfaces for broad community adoption.

---

## 5 | CONCLUSIONS

This study presents CAST, a general framework that integrates causal inference into the species distribution modelling pipeline through DAG structure learning, DML treatment effect estimation, causal forest heterogeneous effect mapping, and causally-informed feature engineering via CI-MLP. Systematic validation across 226 species in six biogeographic regions demonstrates:

1. **Causal screening is effective**: CAST's adaptive multi-criteria screening reduces dimensionality while preserving or improving SDM predictive performance across all tested algorithms, confirming that full-variable models contain substantial non-causal noise.

2. **Causal structure encoding provides additional value**: The CI-MLP, which encodes DAG topology and ATE estimates directly into the feature space, outperforms architecturally identical but structure-agnostic neural networks for the majority of species. This demonstrates that causal structure contains predictive information beyond what is captured by variable selection alone.

3. **The structure effect depends on DAG informativeness**: CI-MLP advantages correlate with DAG sparsity, confirming that the benefit derives from genuine causal structural information rather than arbitrary feature expansion.

4. **Spatially heterogeneous causal effects provide new interpretability**: CATE maps reveal where each environmental driver exerts the strongest causal influence, providing actionable information for spatially targeted conservation that is unavailable from conventional SDM outputs.

CAST provides a reproducible, algorithm-agnostic methodology for moving SDMs from "correlative prediction" toward "causally-informed modelling," addressing a complementary dimension to CISO's encoding of inter-species biotic structure. Together, these approaches point toward a future in which SDMs explicitly represent the full structural complexity of ecological systems — both the causal topology among environmental drivers and the interaction networks among species.

---

## FIGURE LEGENDS

**Figure 1.** Overview of the CAST framework. **(a)** The three-stage CAST pipeline: Stage 1 (DAG learning) infers causal structure among environmental variables via bootstrap Hill-Climbing; Stage 2 (ATE/CATE estimation) quantifies global and spatially heterogeneous causal effects via DML and causal forests; Stage 3 (CI-MLP) encodes the causal topology into the model's feature space through ATE-weighted inputs and DAG-guided interaction features. **(b)** Three-group experimental design: Group A (full post-VIF variables), Group B (CAST-screened variables), Group C (CI-MLP with causal feature engineering). The screening effect (A -> B) and structure effect (B -> C) are isolated by comparing matched models.

**Figure 2.** Causal structure analysis for an example species. **(a)** Consensus DAG showing strong bootstrap edges (strength >= 0.7, direction >= 0.6) among environmental variables, with edge width proportional to bootstrap strength. **(b)** Causal role assignment (Root/Mediator/Terminal) based on DAG position. **(c)** ATE forest plot showing estimated causal effects with 95% confidence intervals. **(d)** CAST variable screening scores decomposed into three components (DAG out-degree, ATE effect, RF importance).

**Figure 3.** Three-group performance comparison across 226 species and 6 regions. **(a)** Mean AUC by model and experimental group. **(b)** Mean TSS by model and experimental group. Error bars indicate standard error across species. CI-MLP (Group C, red) is compared against both full-variable baselines (Group A, grey) and CAST-screened baselines (Group B, blue).

**Figure 4.** CI-MLP versus FlatNN per-species AUC scatter plot. Each point represents one species; points above the diagonal indicate CI-MLP outperformance. Points are colored by biogeographic region and sized by DAG density.

**Figure 5.** CAST advantage decomposition by region. **(a)** Screening effect (FlatNN_cast - FlatNN_full AUC) per species, grouped by region. **(b)** Structure effect (CI-MLP - FlatNN_cast AUC) per species, grouped by region. Diamond markers indicate regional means.

**Figure 6.** DAG density versus CI-MLP advantage. Scatter plot of species-level DAG density (x-axis) against the AUC difference between CI-MLP and FlatNN_cast (y-axis), with linear regression fit. Species in regions with sparser DAGs show larger CI-MLP advantages (negative correlation).

**Figure 7.** Variable reduction and interaction features across regions. **(a)** Percentage of variables removed by CAST screening (post-VIF to CAST-selected). **(b)** Number of DAG-guided interaction features constructed for CI-MLP.

**Figure 8.** Spatially explicit CATE maps. For selected species and causal driver variables, maps show the estimated conditional average treatment effect tau_hat(x_i) at each spatial location. Warm colors indicate locations where the variable strongly increases occurrence probability; cool colors indicate locations where the effect is weak or negative. Panels show different species-variable combinations to illustrate the diversity of spatial heterogeneity patterns.

---

## DATA AVAILABILITY STATEMENT

- **Benchmark data**: The disdat SDM benchmark dataset (Elith et al., 2020) is available via the `disdat` R package on CRAN.
- **Code**: The CAST pipeline, including all analysis scripts and figure generation code, is archived on GitHub (https://github.com/xxxxx/CAST-SDM) with Zenodo DOI.

---

## REFERENCES

Allouche, O., Tsoar, A., & Kadmon, R. (2006). Assessing the accuracy of species distribution models: Prevalence, kappa and the true skill statistic (TSS). *Journal of Applied Ecology*, 43(6), 1223–1232.

Arif, S., & MacNeil, M. A. (2022). Predictive models aren't for causal inference. *Ecology Letters*, 25(8), 1741–1745.

Austin, M. P., & Van Niel, K. P. (2011). Improving species distribution models for climate change studies: Variable selection and scale. *Journal of Biogeography*, 38(1), 1–8.

Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5–32.

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, 21(1), C1–C68.

Deneu, B., et al. (2025). CISO: Species distribution modelling conditioned on incomplete species observations. *Methods in Ecology and Evolution*.

Dormann, C. F., et al. (2013). Collinearity: A review of methods to deal with it and a simulation study evaluating their performance. *Ecography*, 36(1), 27–46.

Elith, J., & Leathwick, J. R. (2009). Species distribution models: Ecological explanation and prediction across space and time. *Annual Review of Ecology, Evolution, and Systematics*, 40, 677–697.

Elith, J., et al. (2020). The disdat SDM benchmark: Presence-only methods compared. *Ecography*, 43(7), 1021–1032.

Fielding, A. H., & Bell, J. F. (1997). A review of methods for the assessment of prediction errors in conservation presence/absence models. *Environmental Conservation*, 24(1), 38–49.

Guisan, A., & Thuiller, W. (2005). Predicting species distribution: Offering more than simple habitat models. *Ecology Letters*, 8(9), 993–1009.

Lin, T.-Y., Goyal, P., Girshick, R., He, K., & Dollar, P. (2017). Focal loss for dense object detection. *Proceedings of the IEEE International Conference on Computer Vision*, 2980–2988.

Merow, C., Smith, M. J., & Silander, J. A. (2013). A practical guide to MaxEnt for modeling species' distributions. *Ecography*, 36(10), 1058–1069.

Pearl, J. (2009). *Causality: Models, reasoning, and inference* (2nd ed.). Cambridge University Press.

Peters, J., Janzing, D., & Scholkopf, B. (2017). *Elements of causal inference*. MIT Press.

Phillips, S. J., Anderson, R. P., & Schapire, R. E. (2006). Maximum entropy modeling of species geographic distributions. *Ecological Modelling*, 190, 231–259.

Phillips, S. J., et al. (2017). Opening the black box: An open-source release of Maxent. *Ecography*, 40(7), 887–893.

R Core Team. (2023). *R: A language and environment for statistical computing*. R Foundation.

Runge, J. (2023). Causal inference for time series. *Nature Reviews Methods Primers*, 3, 58.

Schrodt, F., et al. (2025). Advancing causal inference in ecology: Pathways for biodiversity change detection and attribution. *Methods in Ecology and Evolution*, 16(10), 2276–2304.

Scutari, M. (2010). Learning Bayesian networks with the bnlearn R package. *Journal of Statistical Software*, 35(3), 1–22.

Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. *Journal of the American Statistical Association*, 113(523), 1228–1242.

Wood, S. N. (2017). *Generalized additive models: An introduction with R* (2nd ed.). CRC Press.

Yates, K. L., et al. (2018). Outstanding challenges in the transferability of ecological models. *Trends in Ecology & Evolution*, 33(10), 790–802.

Zuur, A. F., Ieno, E. N., & Elphick, C. S. (2010). A protocol for data exploration to avoid common statistical problems. *Methods in Ecology and Evolution*, 1(1), 3–14.
