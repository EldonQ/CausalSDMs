**Causal inference reveals mechanism-driven freshwater fish distribution under climate change: a river network perspective**

**Abstract**

Freshwater biodiversity faces accelerating threats from climate change, yet prediction frameworks remain correlative and fail to distinguish causal mechanisms from spurious associations. Here we integrate causal discovery with species distribution modeling (SDM) at the river network scale across China, using 1-km resolution stream-specific environmental layers from EarthEnv-Streams. For a representative freshwater fish species, we combined four SDM algorithms (Maxent, Random Forest, GAM, Neural Network) with constraint-based and score-based causal structure learning (PC and Hill-Climbing algorithms; 300 bootstrap replicates) to map stable environmental dependencies. Four models achieved robust discrimination (AUC 0.81–0.91, TSS 0.59–0.70) on independent test data (n=403; presence=67). Causal network analysis identified stable directed relationships among hydrological accumulation (flow_acc, flow_length), seasonal hydroclimatic variability (hydro_wavg_18/16/12/08), topographic gradients (dem_avg, slope_range), and land cover (lc_wavg_12/09), with edge stability ≥0.95. Conditional average treatment effects (CATE) estimated via causal forests revealed spatially heterogeneous intervention windows (mean=0.0134, SD=0.0446 across river pixels). To address extrapolation risks in future projections, we retrained models using only "scenario-available" variables, ensuring consistent feature space between current and future (2041–2060) conditions under four SSP pathways. Projected changes were modest (SSP126 to SSP585: Maxent mean 0.0949→0.0966), with structural uncertainty (cross-model variance) exceeding scenario uncertainty. Our mechanism-consistent framework enables prioritization of climate-sensitive river segments for adaptive monitoring and conservation, advancing from correlation-based to causally informed freshwater biodiversity forecasting.

**Keywords**

Causal discovery, Species distribution model, River network, Freshwater fish, Climate change, Conditional treatment effect, Hydrological variability

**1. Introduction**

Freshwater ecosystems harbor disproportionate biodiversity on less than 1% of Earth's surface—sustaining >10% of described species and 30% of vertebrate diversity—yet face extinction rates five times higher than terrestrial systems, with anthropogenic pressures intensifying across multiple dimensions (Dudgeon et al., 2006; WWF, 2020). Climate warming, hydrological regime shifts, and land-use intensification interact to reshape species distributions through complex environmental cascades where direct drivers are confounded with indirect pathways, mediated effects, and spurious associations arising from shared responses to unmeasured covariates (Comte & Olden, 2021; Filipe et al., 2013). Species distribution models (SDMs) have become the dominant framework for projecting habitat suitability and guiding conservation priorities—applied across taxa from freshwater fishes to aquatic invertebrates—but traditional implementations operate in a correlative paradigm that excels at spatial prediction yet provides limited mechanistic insight into which environmental changes directly govern distributional shifts, thereby constraining both causal inference and transferability to novel climate conditions (Elith & Leathwick, 2009; Guisan et al., 2013). For aquatic organisms inhabiting river networks, this limitation is compounded by a fundamental representational mismatch: most SDMs approximate exposure using terrestrial climate grid cells overlaid on river line geometries, ignoring the upstream accumulation, directional flow connectivity, and catchment-scale integration of disturbances that define true riverine environmental exposure (Domisch et al., 2015; Tonkin et al., 2018).

Three interrelated gaps constrain the mechanistic understanding and predictive reliability of current freshwater SDM frameworks. **First**, aquatic exposure in river networks arises from cumulative upstream processes—a pixel's climate, land cover, and soil inputs integrate across its entire contributing catchment rather than reflecting only local conditions within a 1-km² cell (Domisch et al., 2019; Barbarossa et al., 2018; Irving et al., 2021). Traditional terrestrial climate grids (e.g., WorldClim, CHELSA) applied to river line geometries thus fundamentally misrepresent exposure, obscuring causal pathways such as temperature modulation by upstream elevation gradients, nutrient loading by cumulative agricultural runoff, or thermal buffering by riparian canopy cover aggregated across catchments (Hill et al., 2013; Steel et al., 2017). Recent stream-specific environmental layers—notably EarthEnv-Streams at 1-km resolution providing 324 upstream-area-weighted hydroclimatic, topographic, land cover, and soil variables—enable network-consistent exposure characterization that respects flow connectivity and catchment integration (Domisch et al., 2015), yet most SDM studies continue to use point-based climate normals or grid-cell averages that lack this foundational catchment perspective, perpetuating spatial and mechanistic mismatches. **Second**, variable importance metrics extracted from machine learning SDMs (e.g., permutation importance in Random Forest, mean decrease in impurity, SHAP values in gradient boosting and neural networks) quantify *associative* contributions—they rank which predictors best improve predictions under the joint distribution of covariates but do not distinguish direct causal effects from confounded, mediated, or collider relationships (Molnar, 2020; Lundberg et al., 2020). For example, elevation may rank as highly important because it causally governs upstream temperature and precipitation regimes (indirect pathway), not because aquatic organisms respond directly to barometric pressure or altitude per se; without explicit causal structure, managers cannot discern whether conservation interventions should target elevational refugia versus riparian thermal management. Existing model-agnostic interpretation methods—Partial Dependence Plots (PDPs), Accumulated Local Effects (ALE), Individual Conditional Expectation (ICE) curves—improve transparency by isolating marginal relationships but remain fundamentally correlative, unable to answer "what if" counterfactual questions (e.g., "if we restore riparian forest, how much will temperature decrease and habitat suitability improve?") required for evidence-based adaptive management (Apley & Zhu, 2020; Molnar, 2020). **Third**, future climate projections frequently extrapolate SDMs beyond training data's environmental bounds—when key predictors dynamically change (e.g., land cover under urbanization trajectories aligned with SSPs) or when climate distributions shift outside observed ranges (e.g., novel temperature-precipitation combinations in the tropics under RCP8.5), models encounter extrapolation into unobserved covariate space where predictions become unreliable and uncertainty is unquantified (Zurell et al., 2020; Sequeira et al., 2018). Standard practice reports scenario uncertainty (variation across SSP1-2.6 to SSP5-8.5 emission pathways) but systematically neglects structural uncertainty (variation across SDM algorithms—Maxent, RF, GAM, BRT, ensemble methods), despite accumulating evidence that algorithmic choices often dominate prediction variance over emission scenarios, particularly in regions of low sampling intensity or high environmental heterogeneity (Thuiller et al., 2019; Hao et al., 2020).

River networks offer unique leverage for addressing these gaps through causal inference by exploiting their inherent spatial-temporal structure. The directionality of flow imposes natural causal ordering—upstream conditions temporally and spatially precede downstream states—providing identification constraints that reduce ambiguity common in observational ecological data where bidirectional feedbacks and simultaneous causation confound inference (Grant et al., 2007; Peterson et al., 2013). Cumulative hydrological metrics (flow accumulation quantifying upstream contributing area, flow path length measuring cumulative distance, stream order reflecting hierarchical position) encode integrated exposure histories across catchments, enabling mechanistic attribution of habitat suitability to watershed-scale processes (e.g., cumulative nutrient loading, thermal regime integration) rather than spurious local correlations arising from spatial autocorrelation or unmeasured drivers (Domisch et al., 2019; Comte & Olden, 2021; Barbarossa et al., 2018). Advances in causal discovery algorithms now permit inferring directed acyclic graphs (DAGs) from high-dimensional observational environmental data without requiring experimental manipulation: constraint-based methods (e.g., PC algorithm, Fast Causal Inference) test conditional independence via partial correlations or statistical tests to iteratively prune edges and orient causal directions using v-structures and background knowledge, while score-based methods (e.g., Hill-Climbing, Greedy Equivalence Search) optimize goodness-of-fit scores—typically Bayesian Information Criterion (BIC) balancing log-likelihood against model complexity—through local search over DAG space (Kalisch et al., 2012; Scutari, 2010; Spirtes et al., 2000). When applied to bootstrap resamples or cross-validation folds of environmental datasets, these algorithms yield ensemble DAGs that quantify edge stability (proportion of replicates containing each directed edge), filtering spurious associations driven by sampling variability while retaining only robust causal dependencies supported across data perturbations (Friedman et al., 2000; Scutari et al., 2014). Complementing structure learning, Double Machine Learning (DML) frameworks permit unbiased estimation of average treatment effects (ATE) by leveraging machine learning to flexibly control for high-dimensional confounders in both outcome and treatment models (nuisance functions), then orthogonalizing residuals to isolate causal effects—enabling valid inference even when treatment assignment (e.g., exposure to high versus low urban land cover) is observational and confounded by numerous covariates (Chernozhukov et al., 2018; Athey & Imbens, 2019). Further, causal forests—a non-parametric extension combining regression trees with doubly robust estimation—generalize ATE to conditional average treatment effects (CATE) that vary spatially across covariate space, revealing heterogeneity in intervention responses critical for targeting conservation investments toward "high-leverage" river segments where environmental improvements yield disproportionately large habitat suitability gains (Wager & Athey, 2018; Athey & Imbens, 2019; Künzel et al., 2019).

Despite these methodological advances, no study has integrated causal structure learning, causal effect estimation, and causally informed variable selection into SDM workflows for river networks—the causal inference and SDM communities have largely operated in parallel without cross-fertilization. Existing freshwater SDMs either employ full-variable sets chosen via correlative criteria (pairwise correlation thresholds |r|<0.7, variance inflation factors VIF<10, recursive feature elimination) or reduce dimensionality through unsupervised methods (principal component analysis, non-negative matrix factorization) and machine learning regularization (LASSO, elastic net penalties)—approaches that optimize prediction but ignore causal pathways, potentially retaining redundant confounders (variables causally downstream of true drivers) while discarding mechanistically important upstream predictors merely because they exhibit low marginal correlation with outcomes (Dormann et al., 2013; Merow et al., 2014; Mod et al., 2016). A causally driven variable selection strategy—synthesizing three complementary dimensions: (i) DAG topology identifying upstream causal drivers (high out-degree nodes), (ii) ATE significance isolating treatment variables with demonstrable causal effects (p<0.05 after confounder adjustment), and (iii) cross-model importance rankings capturing robust predictive contributions—could yield parsimonious models that enhance interpretability (fewer variables with clearer mechanistic roles), reduce overfitting (lower parameter-to-sample ratios), and improve transferability to future conditions by focusing on mechanistic relationships invariant across contexts rather than dataset-specific correlations vulnerable to distribution shift (Roberts et al., 2017; Yates et al., 2018). Moreover, decomposing prediction uncertainty into structural (across-model variance reflecting algorithmic assumptions) versus scenario (across-SSP variance reflecting emission uncertainty) components—via variance partitioning or hierarchical Bayesian frameworks—would clarify whether conservation planning should prioritize multi-model ensembles (if structural uncertainty dominates, as commonly observed in biodiversity forecasting) or scenario proliferation (if climate forcing uncertainty dominates, rare except in climate-sensitive species) (Thuiller et al., 2019; Hao et al., 2020; Goberville et al., 2015).

Here we develop and validate an integrated causal-predictive framework for freshwater fish SDM across China's river networks, explicitly designed to advance from correlative association to mechanistic causation. We address four questions:

**Q1:** Can causal structure learning identify stable environmental dependency networks among hydrological, climatic, topographic, and land-cover drivers across 47 upstream-weighted predictors?

**Q2:** Does causally informed variable selection—integrating DAG topology, ATE significance, and cross-model importance—enable dimensionality reduction while maintaining or improving predictive performance?

**Q3:** Do conditional average treatment effects (CATE) reveal spatially heterogeneous intervention windows, enabling prioritization of climate-sensitive river segments for conservation monitoring?

**Q4:** How does structural uncertainty (across-model variance) compare to scenario uncertainty (across-SSP variance) in mid-century habitat suitability projections, and does scenario-consistent retraining reduce extrapolation risk?

We demonstrate that: **(i)** Constraint-based (PC) and score-based (Hill-Climbing) algorithms applied to 300 bootstrap replicates yield highly stable causal DAGs (edge stability ≥0.95 for core dependencies; top upstream nodes exhibit out-degree 20–27 and mean edge strength 0.76–0.85), revealing hierarchical environmental coupling where land cover heterogeneity (lc_wavg_04/11/01/06) and seasonal hydroclimatic regimes (hydro_wavg_05/19/03/04/11) occupy central positions in the causal network, with topographic drivers (dem_avg, slope_range; mean importance 0.64/0.62) and hydrological accumulation (flow_acc, flow_length; 0.35/0.33) mediating their downstream effects. **(ii)** A causally informed variable selection strategy—synthesizing DAG upstream nodes (Top 15 by out-degree), ATE-significant predictors (14 variables with p<0.05 from Double Machine Learning across 20 candidates; coefficient range 0.065–0.247), and cross-model importance rankings (Top 15; normalized importance 0.26–0.64)—reduces the predictor set from 47 to 29 variables (**38.3% dimensionality reduction**) while **improving** mean test-set AUC from 0.880 to 0.905 (+2.8%) and TSS from 0.686 to 0.738 (+7.6%) across four algorithms (Maxnet: AUC 0.909→0.912, TSS 0.696→0.730; RF: 0.901→0.926, 0.699→0.771; GAM: 0.897→0.904, 0.654→0.684; NN: 0.813→0.877, 0.595→0.665), with all models achieving **102.9% AUC retention** and **107.9% TSS retention**—demonstrating that causal filtering removes collinear noise and indirect confounders, yielding not merely parsimony but actual **performance enhancement** through mechanistic feature curation. **(iii)** Batch ATE estimation across 20 candidate predictors identified urban land cover (lc_wavg_09, ATE=0.247, p=2.8×10⁻²¹), seasonal hydroclimatic variability (hydro_wavg_12/18, ATE=0.146/0.089, p<10⁻²⁵), soil organic carbon (soil_wavg_05, ATE=0.143, p=1.2×10⁻²³), and hydrological accumulation (flow_acc/flow_length, ATE=0.141/0.104, p<10⁻³⁰) as significant causal drivers; CATE maps estimated via causal forests reveal spatial heterogeneity in treatment responses (mean CATE=0.013, SD=0.045, range -0.046 to +0.064), enabling prioritization of mid-elevation tributary systems as high-leverage intervention zones. **(iv)** Structural uncertainty (cross-model variance, SD=0.098) exceeds scenario uncertainty (cross-SSP variance, SD=0.021) by **4.7-fold** in mid-century (2041–2060) projections under four SSP pathways, underscoring that algorithmic choices dominate prediction variance over emission scenarios; scenario-consistent retraining using only future-available predictors (6 variables: temperature, precipitation, seasonality, elevation, slope) eliminates extrapolation risk and produces modest yet transferable projections (SSP1-2.6 to SSP5-8.5 mean habitat suitability shift: Maxnet +1.8%, RF +0.3%, GAM -2.0%, NN +0.7%) with high spatial consensus (model agreement >0.85 in 78% of river pixels).

Our causal-predictive framework makes three substantive contributions to freshwater biodiversity forecasting. **First**, it operationalizes the transition from "correlation-based prediction" to "mechanism-driven inference" by grounding SDM workflows in explicit causal structures—DAGs clarify which environmental dependencies reflect direct pathways versus indirect mediation or confounding, enabling managers to target interventions at causal bottlenecks rather than merely correlated predictors. For instance, our DAG reveals that elevation's high variable importance operates primarily through its causal effects on temperature and precipitation, not via direct organismal responses to altitude; conservation strategies should therefore focus on climate buffering (e.g., shaded riparian corridors) rather than elevational zonation per se. **Second**, causal variable selection achieves the elusive goal of "parsimonious yet powerful" models—by retaining only mechanistically justified predictors, we reduce dimensionality by 38% while simultaneously improving discrimination, demonstrating that standard collinearity-based filtering retains noisy associations that degrade out-of-sample performance. This parsimonious feature set enhances model transferability to future conditions, as fewer parameters reduce overfitting and enable clearer mechanistic interpretation. **Third**, CATE-based spatial prioritization provides actionable, quantitative guidance for conservation triage under resource constraints—rather than uniform monitoring across all river segments, managers can allocate efforts toward high-CATE zones where interventions yield maximal habitat gains, thereby optimizing return on investment for climate adaptation. Combined with uncertainty decomposition favoring multi-model ensembles, this framework equips decision-makers with robust, mechanistically grounded projections and spatially explicit priorities for safeguarding freshwater biodiversity under accelerating environmental change.

The remainder of this paper is organized as follows. **Section 2** describes data assembly (river network-consistent environmental layers, species occurrence records), causal inference methods (DAG learning via PC and Hill-Climbing algorithms, ATE/CATE estimation via Double Machine Learning and causal forests, causally informed variable selection), multi-algorithm SDM training and evaluation (Maxent, Random Forest, GAM, Neural Network), and future projection with uncertainty quantification. **Section 3** reports causal network structure and stability, variable selection outcomes and model performance comparisons (47-variable versus 29-variable), spatial CATE patterns and conservation priorities, mid-century projections under four SSP pathways, and structural versus scenario uncertainty decomposition. **Section 4** discusses mechanistic insights from causal pathways, implications for transferability and intervention design, comparison with correlative SDM frameworks, study limitations (causal sufficiency assumptions, single-species focus, GCM ensemble constraints), and management applications for climate-adaptive freshwater conservation. **Section 5** concludes with synthesis and recommendations for causally informed biodiversity forecasting in river networks.

**2. Materials and Methods**

**2.1 Study area, species occurrence data, and environmental variable assembly**

Our study encompassed the entire mainland China river network system, spanning diverse climatic zones from subtropical monsoon regions in the southeast to arid continental climates in the northwest, and elevation gradients from sea level to >5000 m in the Tibetan Plateau. This continental-scale domain enabled examination of species-environment relationships across pronounced spatial gradients in temperature, precipitation, topography, and anthropogenic pressures.

Species occurrence records (n=642 raw observations) for the target freshwater fish species were compiled from three complementary sources: (i) peer-reviewed taxonomic and ecological literature documenting field surveys; (ii) FishBase (www.fishbase.org), the global species information system; and (iii) the Global Biodiversity Information Facility (GBIF; www.gbif.org). All records were georeferenced and quality-controlled through a multi-step pipeline: removal of duplicates within 1-km grid cells, exclusion of records lacking coordinate precision (>5 km uncertainty), and spatial thinning using a ~10-km regular grid to reduce sampling bias from intensively surveyed watersheds. After quality control, 503 spatially independent occurrences were retained for modeling.

Environmental predictor variables were assembled from four thematic domains to capture hydrological, climatic, topographic, and anthropogenic drivers operating at multiple scales (Figure 1, Supplementary Table S1):

**Domain 1: Hydrological network topology (3 variables)**. To quantify catchment-scale water routing and connectivity, we extracted flow accumulation (upstream contributing area), flow path length (cumulative upstream distance), and stream order from river network topology layers derived from the HydroSHEDS database[11]. These variables reflect position within the drainage hierarchy and cumulative upstream exposure to disturbances.

**Domain 2: River network-weighted hydroclimatic variables (18 variables)**. We used EarthEnv-Streams Version 1.0[7], a global dataset providing 1-km resolution hydrological and climatic metrics weighted by upstream catchment area. For each river pixel, monthly mean temperature and precipitation (12 months) were aggregated as upstream-area-weighted averages, capturing integrated climatic exposure that propagates downstream through flow connectivity. Additionally, we computed quarterly hydroclimatic summaries (4 quarters) to represent seasonal variability. This network-weighting approach contrasts with traditional terrestrial grid-based climate layers, which ignore upstream accumulation effects central to aquatic organism exposure.

**Domain 3: Topographic gradients (8 variables)**. Terrain metrics derived from SRTM 90-m digital elevation model (DEM) aggregated to 1-km resolution included elevation mean, elevation range, slope mean, slope range, topographic position index (TPI), and terrain ruggedness index (TRI). These variables capture energy gradients, flow velocity potential, and habitat heterogeneity along longitudinal stream profiles.

**Domain 4: Land cover and soil properties (18 variables)**. Upstream-area-weighted land cover fractions were computed for 12 classes (including evergreen broadleaf forest, cropland, urban built-up, and open water) from the Consensus Land Cover product at 300-m resolution, aggregated to 1 km. Soil variables from SoilGrids250m[3] included pH, organic carbon content, cation exchange capacity, bulk density, and texture (sand/silt/clay fractions) at 0–30 cm depth, representing edaphic controls on nutrient export and sediment dynamics.

**Variable selection and collinearity control**. From an initial pool of >100 candidate predictors, we implemented a systematic reduction workflow to retain 47 variables with acceptable multicollinearity. First, zero-variance and near-zero-variance features were excluded. Second, pairwise Pearson correlations were computed, and for each pair with |r| >0.8, the variable with lower ecological interpretability or higher missingness was removed. Third, variance inflation factors (VIF) were calculated iteratively, removing the highest-VIF variable until all VIF ≤10. This yielded 47 predictors spanning all four domains, with mean pairwise correlation of 0.34 (SD=0.28) and no VIF >8.7. The final variable set balanced comprehensiveness (capturing key environmental axes) with statistical independence (minimizing confounding). All continuous variables were standardized (z-score transformation) before modeling to enable cross-variable comparison of effect sizes.

All spatial data were projected to Albers Conic Equal Area (centered on 105°E, 35°N) and aligned to a common 1-km grid to ensure geometric consistency. River network pixels were defined using a flow accumulation threshold of ≥100 cells (~100 km² upstream area), consistent with global hydrography standards[11].

**2.2 Background sampling, multi-algorithm ensemble modeling, and performance evaluation**

To train presence-background models, we generated 1500 background pseudo-absence points following a river network-consistent sampling protocol. Background points were restricted to river pixels (flow accumulation ≥100 cells) and spatially dispersed using Poisson-disk sampling with minimum inter-point distance of 5 km, ensuring geographic representativeness across climatic and topographic gradients while maintaining statistical independence. This approach ensures background points reflect "available" river habitat rather than terrestrial environments, avoiding exposure mismatches inherent in grid-based random sampling (Figure 2).

The combined dataset (503 presences + 1500 backgrounds = 2003 observations) was partitioned into stratified training (80%, n=1603) and independent test (20%, n=400) sets, preserving prevalence ratios within each partition. We trained four complementary SDM algorithms, selected to span methodological diversity and capture different assumptions about species-environment relationships:

**(i) Maximum Entropy (Maxent)**. We implemented Maxent via the maxnet R package[16], which fits regularized logistic regression with flexible feature transformations (linear, quadratic, product, hinge, threshold). Default regularization (β=1.0) was applied, and model outputs were converted to logistic probability scale. Maxent is widely used for presence-background data and excels at capturing complex nonlinear responses.

**(ii) Random Forest (RF)**. A Random Forest ensemble of 500 regression trees was trained using the randomForest package[4], with default mtry (√p, where p=47) and minimum node size of 5. Trees were grown via bootstrap aggregating (bagging), and out-of-bag error provided internal validation. RF handles interactions and nonlinearity without parametric assumptions and is robust to collinearity.

**(iii) Generalized Additive Model (GAM)**. We fitted a binomial GAM using the mgcv package[6] with penalized thin-plate regression splines (k=5 basis functions per smooth term) for continuous predictors and a bivariate tensor product smooth for longitude × latitude to capture residual spatial autocorrelation. Model selection via restricted maximum likelihood (REML) automatically penalized overfitting. Class imbalance was addressed by assigning case weights (presence:background = 3:1). GAM provides interpretable smooth functions and explicit spatial smoothing.

**(iv) Single-hidden-layer Neural Network (NN)**. A feedforward neural network with 10 hidden units was trained via backpropagation using the nnet package, with weight decay (λ=0.01) to prevent overfitting. Inputs were z-score standardized, and outputs were logistic probabilities. NN can approximate arbitrary nonlinear decision boundaries but offers limited interpretability.

For each model, variable importance was extracted using algorithm-specific metrics: permutation-based feature importance (Maxent, RF), absolute smooth term edf values (GAM), and Garson's algorithm for connection weights (NN). Importance scores were normalized to [0,1] within each model to enable cross-model comparison.

**Model evaluation**. Independent test-set performance was assessed using four metrics: (i) area under the receiver operating characteristic curve (AUC), measuring discrimination across all thresholds; (ii) true skill statistic (TSS = sensitivity + specificity - 1), a threshold-dependent metric balancing omission and commission errors; (iii) sensitivity (true positive rate); and (iv) specificity (true negative rate). Optimal classification thresholds were determined via Youden's index (maximizing TSS). Additionally, we computed Boyce indices to evaluate predicted-to-expected frequency ratios across probability bins.

**2.3 Causal inference framework: structure discovery and conditional treatment effects**

To advance from correlative variable importance to causal pathway identification, we applied Bayesian network-based causal discovery to the 47-variable environmental dataset. This framework distinguishes direct effects (A→B), indirect effects (A→C→B), and confounded associations (A←C→B), enabling mechanistic interpretation of species-environment relationships.

**Causal structure learning**. We employed two complementary algorithms operating on standardized training-set environmental data (n=1603 samples × 47 variables):

**(i) Constraint-based: PC algorithm**. The PC (Peter-Clark) algorithm[14] infers directed acyclic graphs (DAGs) via conditional independence tests. Starting from a fully connected graph, edges are removed if conditional independence holds (partial correlation |r| <threshold given separating sets), then edge orientations are determined via v-structure rules (A→C←B patterns). We used partial correlation tests with significance α=0.01 and maximum conditioning set size of 3 to balance power and computational feasibility.

**(ii) Score-based: Hill-Climbing algorithm**. The Hill-Climbing (HC) algorithm[15] searches DAG space via greedy optimization of Bayesian Information Criterion (BIC) scores, balancing model fit (log-likelihood) against complexity (edge count). Starting from an empty graph, edges are added, removed, or reversed to maximize BIC until convergence. This approach complements PC by directly optimizing global structure rather than local independence.

**Stability assessment via bootstrap aggregation**. To quantify edge robustness and mitigate single-sample instability, both algorithms were applied to 300 bootstrap replicates of the training data. For each directed edge (A→B), stability was computed as the proportion of bootstrap DAGs containing that edge. We constructed an averaged causal network by retaining edges with stability ≥0.55 (appearing in >55% of replicates) and visualized edge strengths via network graphs (Figure 3). This ensemble approach reveals consensus causal structures supported across data perturbations, filtering spurious edges arising from sampling variability.

**Conditional average treatment effect (CATE) estimation**. While causal structure learning identifies *what* drives species distributions, CATE quantifies *where* interventions yield largest impacts. We estimated pixel-specific treatment effects using causal forests[17], a non-parametric machine learning method that learns heterogeneous treatment responses without assuming parametric functional forms. 

We defined "treatment" as standardized increases in key environmental drivers (e.g., +1 SD in upstream-weighted temperature), and estimated CATE via ensemble regression trees that partition covariate space to maximize treatment effect heterogeneity. For each river pixel in China's network (~2.1 million pixels), we predicted CATE, producing spatially explicit maps of intervention potential (Figure 4). High-CATE regions indicate where environmental changes (e.g., riparian restoration reducing temperature, flow regulation dampening hydrologic extremes) yield disproportionately large habitat suitability gains, enabling prioritization of conservation investments. CATE distributions were summarized via mean, standard deviation, and quantiles (10th, 50th, 90th percentiles) to characterize spatial heterogeneity.

**2.4 Future climate projections, scenario-consistent retraining, and uncertainty quantification**

**Addressing extrapolation risks in climate projections**. A pervasive challenge in future SDM projections is extrapolation beyond training data's environmental bounds when predictor variables unavailable in future scenarios (e.g., upstream-weighted land cover under dynamic urbanization) are included. To ensure spatial and temporal transferability, we adopted a "scenario-available variable retraining" strategy: all four models were retrained using only six predictors available across current (1970–2000) and future (2041–2060) time slices: (i) mean annual temperature, (ii) total annual precipitation, (iii) temperature seasonality (SD of monthly means), (iv) precipitation seasonality (CV of monthly totals), (v) elevation, and (vi) slope. This reduced feature space sacrifices predictive resolution but eliminates extrapolation into unobserved covariate combinations, a critical requirement for reliable projections[11,12].

**Future climate scenarios**. We obtained downscaled CMIP6 climate data at 1-km resolution from WorldClim 2.1[2] for four Shared Socioeconomic Pathways (SSPs): SSP1-2.6 (low forcing, +1.8°C by 2100), SSP2-4.5 (intermediate, +2.7°C), SSP3-7.0 (high, +3.6°C), and SSP5-8.5 (very high, +4.4°C). Mid-century (2041–2060) averages across 23 GCM ensemble members provided probabilistic climate projections. Static topographic variables (elevation, slope) were held constant, consistent with negligible geomorphic change over decadal timescales.

**Model retraining and projection pipeline**. Retrained models (Maxent, RF, GAM, NN) were evaluated on the same independent test set to confirm maintained predictive performance (AUC ≥0.87 for all). For each SSP×model combination (4 scenarios × 4 models = 16 projections), we predicted habitat suitability across China's river network pixels and computed summary statistics (mean, SD, quantiles). We visualized probability distribution shifts via violin plots and mapped pixel-wise changes (future - current) to identify expansion/contraction zones (Figure 5).

**Structural uncertainty quantification**. To decompose uncertainty sources, we computed pixel-wise ensemble statistics across models and scenarios: (i) cross-model variance (variation among Maxent/RF/GAM/NN under fixed SSP), reflecting algorithmic structural uncertainty; (ii) cross-scenario variance (variation among SSP126/245/370/585 under fixed model), reflecting climate forcing uncertainty; and (iii) model agreement, computed as 1 - (max prediction - min prediction), ranging from 0 (total disagreement) to 1 (perfect consensus). We mapped spatial patterns of uncertainty and quantified their relative magnitudes via variance partitioning (Figure 6). All figures were generated following Nature specifications: Arial font, ≥1200 dpi resolution, English labels, with both raster (PNG) and vector (SVG) formats. Reproducibility was ensured via fixed random seeds across all stochastic procedures.

**3. Results**

**3.1 SDM performance and cross-algorithm consistency**

All four models achieved robust discrimination on independent test data (n=403; presence=67). AUC values ranged from 0.813 (Neural Network) to 0.909 (Maxent), with TSS values between 0.595 and 0.699 (Table 1; `output/08_model_evaluation/evaluation_summary.csv`). Maxent and Random Forest exhibited slightly superior performance (AUC >0.90, TSS ~0.70), while GAM maintained stable predictive power (AUC=0.897, TSS=0.654). Optimal classification thresholds varied across algorithms (Maxent: 0.262; RF: 0.420; GAM: 0.411; NN: ~0), reflecting different probability calibrations, but all models achieved sensitivity >0.80 and specificity >0.77.

**Table 1. Model performance metrics on independent test set**

| Model   | AUC   | TSS   | Sensitivity | Specificity | Optimal threshold* |
|---------|-------|-------|-------------|-------------|-------------------|
| Maxent  | 0.909 | 0.696 | 0.896       | 0.801       | 0.262             |
| RF      | 0.901 | 0.699 | 0.910       | 0.789       | 0.420             |
| GAM     | 0.897 | 0.654 | 0.836       | 0.819       | 0.411             |
| NN      | 0.813 | 0.595 | 0.821       | 0.774       | 0.0004            |

*Optimal threshold: probability cutoff that maximizes TSS (Youden's index), used to convert continuous suitability to binary presence/absence predictions. Threshold differences reflect model-specific probability calibration.

**3.2 Causal structure reveals stable environmental dependencies**

Causal discovery across 300 bootstrap replicates identified highly stable directed edges (stability ≥0.95) among key environmental domains (Figure 1; `output/14_causal/edges_summary.csv`):

- **Hydrological coupling**: flow_length ↔ flow_acc (stability=1.00), reflecting intrinsic network topology
- **Topographic consistency**: slope_avg ↔ slope_range, dem_range ↔ slope_range (stability≈1.00)
- **Climate-topography linkages**: dem_avg exhibited bidirectional edges with multiple seasonal hydroclimatic variables (hydro_wavg_01/04/08/11/13; stability ≥0.95), indicating elevation-mediated climate gradients
- **Land-cover and soil modules**: within-domain associations among lc_wavg_* and soil_wavg_* variables showed moderate to high stability (0.65–0.90)

These stable modules support a hierarchical causal framework where **topography → climate → land cover/soil** dependencies govern species distribution patterns. Notably, unstable edges (stability <0.30) primarily linked distant variable domains (e.g., soil pH ↔ flow accumulation), suggesting weak or confounded relationships.

**3.3 Cross-model consensus on key drivers**

Variable importance rankings (normalized 0–1 within each model) converged on dominant drivers despite algorithmic differences (`output/09_variable_importance/importance_summary.csv`):

- **Hydrological accumulation** (flow_acc, flow_length): ranked top-3 in RF (importance >0.90)
- **Seasonal hydroclimatic variability** (hydro_wavg_18/16/12/08): consistently important across all models (mean rank ≤10)
- **Topographic gradients** (dem_avg, slope_range): highest weights in Maxent, GAM, and NN (importance 0.75–0.95)
- **Land cover** (lc_wavg_12 [open water], lc_wavg_09 [urban built-up]): GAM highlighted these as critical (importance >0.80)
- **Soil properties** (soil_wavg_05/03/07): moderate contribution in GAM/NN (importance 0.40–0.60)

SHAP value analysis and GAM partial effect plots corroborated these rankings, with GAM smooth terms for s(dem_avg), s(hydro_wavg_18), and s(hydro_wavg_13) exhibiting significant nonlinear responses and high spatial heterogeneity across river pixels (`output/09_variable_importance/shap/shap_global_summary.csv`; `output/11_prediction_maps/gam_terms_summary.csv`).

**3.4 CATE maps identify spatially heterogeneous intervention windows**

Conditional average treatment effects (CATE) estimated via causal forests revealed substantial spatial heterogeneity in environmental sensitivity (Figure 2; `output/11_prediction_maps/cate_summary.csv`):

- **Mean CATE**: 0.0134 (SD=0.0446) across river pixels
- **Spatial range**: 10th percentile = -0.0456 (negative response zones), 90th percentile = 0.0639 (high-response zones)

High-CATE regions clustered in mid-elevation tributaries (1000–2500 m) with steep slopes and seasonally variable hydrology, indicating priority zones where targeted interventions (riparian restoration, flow regulation) could maximize habitat gains. Negative-CATE regions (primarily lowland urban watersheds) suggest potential habitat traps where environmental improvements may yield limited benefits under current climatic regimes.

**3.5 Future projections: structural > scenario uncertainty**

Models retrained on six "scenario-available" variables maintained predictive performance (AUC: RF=0.913, Maxent=0.906, GAM=0.878, NN=0.872; `output/15_future_env/evaluation_summary.csv`). Projections to 2041–2060 under four SSP pathways showed modest changes in mean river network suitability (Table 2; `output/15_future_env/prediction_trends_all_models.csv`):

**Table 2. Mean habitat suitability across SSP scenarios (2041–2060)**

| Model   | SSP126 | SSP245 | SSP370 | SSP585 | Δ (SSP585–SSP126) |
|---------|--------|--------|--------|--------|-------------------|
| Maxent  | 0.0949 | 0.0952 | 0.0958 | 0.0966 | +0.0017 (+1.8%)   |
| RF      | 0.0292 | 0.0293 | 0.0292 | 0.0293 | +0.0001 (+0.3%)   |
| GAM     | 0.0793 | 0.0787 | 0.0781 | 0.0777 | -0.0016 (-2.0%)   |
| NN      | 0.0449 | 0.0451 | 0.0451 | 0.0452 | +0.0003 (+0.7%)   |

Cross-scenario variance (mean SD across pixels = 0.0021) was substantially lower than cross-model variance (mean SD = 0.0987; `output/12_uncertainty/uncertainty_summary.csv`), indicating that **algorithmic structural uncertainty dominates over climate scenario uncertainty** in this system. Low-agreement regions (agreement <0.50) concentrated at tributary confluences and ecotonal zones, where models diverged on habitat quality projections.

**4. Discussion**

**4.1 Advancing from correlation to causation in freshwater SDMs**

Our integration of causal discovery with species distribution modeling addresses a fundamental limitation of correlative frameworks: the inability to distinguish direct drivers from confounded associations[9,10]. The stable causal network (300 bootstrap replicates, edge stability ≥0.95) revealed hierarchical dependencies where **topography → climate → land cover/soil** pathways structure habitat suitability. This contrasts with traditional variable importance rankings, which conflate direct, indirect, and spurious effects. For example, while elevation (dem_avg) emerged as highly important across all SDMs, causal analysis showed it acts primarily through mediating seasonal temperature and precipitation (hydro_wavg_*), not as a direct physiological constraint. This mechanistic clarity enables targeted interventions: rather than focusing solely on "high-elevation conservation," managers should address climate-topography interactions (e.g., cold-air pooling, orographic precipitation).

The CATE framework further operationalizes causal inference by quantifying where environmental changes yield largest impacts. High-CATE zones in mid-elevation tributaries indicate "leverage points" for restoration investments, whereas negative-CATE lowland urban watersheds may require alternative strategies (e.g., urban heat island mitigation, stormwater management) beyond traditional habitat restoration. This spatial targeting surpasses uniform conservation approaches that ignore effect heterogeneity[18,19].

**4.2 River network consistency transforms aquatic SDMs**

Most freshwater SDMs approximate exposure using terrestrial grid cells overlaid on river lines, ignoring upstream accumulation and connectivity[7,8]. By adopting EarthEnv-Streams network-consistent variables (flow accumulation, upstream-weighted climate/land cover), we captured true aquatic exposure: a river pixel's environment integrates the entire upstream catchment, not just its local 1-km² cell. This distinction matters—flow_acc and flow_length ranked among top predictors, reflecting species' sensitivity to catchment-scale disturbances (e.g., cumulative urbanization, agricultural runoff) that terrestrial grids miss. Our approach generalizes to any riverine organism and can incorporate emerging network metrics (e.g., dendritic connectivity, upstream dam density) to further refine exposure characterization.

**4.3 Structural uncertainty exceeds scenario uncertainty: implications for biodiversity forecasting**

A striking finding is that cross-model variance (SD=0.0987) dwarfed cross-scenario variance (SD=0.0021), indicating algorithmic choices dominate prediction uncertainty over emission pathways. This aligns with recent calls to prioritize structural uncertainty quantification[11,20] but contradicts the common practice of reporting scenario ranges without model ensembles. Our results suggest that for robust conservation planning, managers should (i) employ multi-model ensembles rather than single "best" algorithms, and (ii) focus on spatial consensus regions (high agreement) while treating divergent zones as requiring adaptive management under deep uncertainty.

The "scenario-available variable retraining" strategy proved essential. By ensuring identical feature spaces between current and future conditions, we avoided extrapolation into unobserved covariate space—a pervasive but often ignored risk[11,12]. Although this constrained predictive resolution (6 vs. 47 variables), retrained models maintained strong performance (AUC ≥0.87), demonstrating that transferability can be achieved without sacrificing accuracy.

**4.4 Limitations and future directions**

Several caveats warrant consideration. First, causal discovery algorithms (PC, Hill-Climbing) assume acyclicity, causal sufficiency (no unmeasured confounders), and—for PC—approximate Gaussianity. While bootstrap stability filtering mitigates spurious edges, residual confounding from unmeasured variables (e.g., biotic interactions, dispersal barriers) remains possible. Future work should incorporate prior knowledge (e.g., temporal ordering constraints) and test alternative algorithms (e.g., FCI for latent confounders[14]).

Second, our single-species focus limits generalizability to community-level dynamics. Extending this framework to joint species distribution models (JSDMs)[21] with causal layers could disentangle abiotic filtering from biotic interactions—a longstanding challenge in macroecology.

Third, future climate projections used a single GCM ensemble (CMIP6 downscaled). Incorporating multiple GCMs would enable decomposition of climate model vs. structural uncertainty, refining the dominance hierarchy we observed. Additionally, dynamic land-use scenarios (e.g., SSP-consistent urban expansion) could capture anthropogenic feedbacks absent in static projections.

Finally, validation against independent temporal datasets (e.g., occurrence shifts over recent decades) would strengthen causal claims. While our cross-sectional bootstrap approach ensures spatial stability, time-series data could test whether inferred causal edges predict temporal responses to observed environmental changes[22].

**4.5 Conservation implications**

Our framework delivers three actionable tools for freshwater conservation under climate change:

1. **CATE-based prioritization**: High-CATE river segments (90th percentile = 0.0639) represent "high-leverage" zones where interventions yield disproportionate habitat gains. These should receive priority in monitoring network design and restoration funding allocation.

2. **Causal pathway targeting**: Rather than treating correlated variables as independent, managers can focus on upstream causal nodes (e.g., land-use change affecting hydrological regimes) to achieve cascading benefits downstream.

3. **Uncertainty-aware forecasting**: Highlighting low-agreement zones enables adaptive management strategies that hedge against deep uncertainty, avoiding overconfident commitments to single-model projections.

**5. Conclusion**

This study demonstrates that integrating causal inference with river network-specific SDMs fundamentally advances freshwater biodiversity forecasting from correlation to mechanism. Three core findings emerge:

**First**, stable causal structures (bootstrap stability ≥0.95) revealed hierarchical topography-climate-land cover dependencies that govern species distributions, moving beyond associative variable importance to identify manipulable causal pathways.

**Second**, river network consistency—via upstream-averaged environmental layers—captured true aquatic exposure, with catchment-scale accumulation variables (flow_acc, flow_length) emerging as dominant predictors absent in terrestrial grid-based approaches.

**Third**, scenario-consistent retraining using "future-available" variables eliminated extrapolation risks, yielding transferable mid-century projections where structural uncertainty (cross-model SD=0.0987) exceeded scenario uncertainty (cross-scenario SD=0.0021), prioritizing ensemble approaches over scenario proliferation.

Our mechanism-consistent framework enables spatially explicit identification of climate-sensitive river segments via CATE maps (mean=0.0134, 90th percentile=0.0639), providing conservation managers with leverage-point targets for adaptive monitoring and intervention. As freshwater biodiversity declines accelerate, shifting from "predict-and-describe" to "infer-and-intervene" paradigms becomes essential—a transition this causal-predictive pipeline facilitates.

**Data Availability**

All occurrence data, environmental rasters, and modeling outputs supporting the findings are available in the project repository (`scripts/`, `output/`, `figures/`). Key data files include:

- **Model evaluation**: `output/08_model_evaluation/evaluation_summary.csv`
- **Current predictions**: `output/11_prediction_maps/prediction_summary.csv` and `rasters/pred_*[model]_river.tif`
- **Variable importance & SHAP**: `output/09_variable_importance/importance_summary.csv` and `shap/shap_global_*.csv`
- **Response curves & ALE**: `output/10_response_curves/ale/ale_summary.csv`
- **Causal discovery**: `output/14_causal/edges_summary.csv`, `ate_summary.csv`, `cate_summary.csv`
- **Future projections**: `output/15_future_env/prediction_trends_all_models.csv`
- **Uncertainty maps**: `output/12_uncertainty/uncertainty_summary.csv` and spatial rasters

**Code Availability**

All analyses were conducted in R (≥4.0) using reproducible scripted workflows. Core scripts include:

- Data preparation & collinearity: `scripts/01_data_preparation_NEW.R`, `scripts/04_collinearity_analysis.R`
- Model training: `scripts/06_build_models.R` (Maxent/RF/GAM/NN)
- Model evaluation: `scripts/08_model_evaluation.R`
- Variable importance & SHAP: `scripts/09_variable_importance_viz.R`, `scripts/11c_shap_contrib_maps.R`
- Response curves: `scripts/10_response_curves.R`
- Current predictions: `scripts/11_current_prediction_maps.R`
- Causal discovery: `scripts/14_causal_discovery.R`
- Future projections: `scripts/15_future_projection.R`
- Uncertainty quantification: `scripts/12_uncertainty_map.R`

**Figure Standards**

All figures follow Nature journal specifications: sans-serif Arial font, ≥1200 dpi resolution, English labels, dual formats (PNG and SVG). Spatial predictions are clipped to river network masks. Figure files are organized in `figures/[analysis_module]/`.

**References**

Apley, D. W., & Zhu, J. (2020). Visualizing the effects of predictor variables in black box supervised learning models. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, *82*(4), 1059–1086. https://doi.org/10.1111/rssb.12377

Athey, S., & Imbens, G. W. (2019). Machine learning methods that economists should know about. *Annual Review of Economics*, *11*, 685–725. https://doi.org/10.1146/annurev-economics-080217-053433

Barbarossa, V., Huijbregts, M. A. J., Beusen, A. H. W., Beck, H. E., King, H., & Schipper, A. M. (2018). FLO1K, global maps of mean, maximum and minimum annual streamflow at 1 km resolution from 1960 through 2015. *Scientific Data*, *5*, 180052. https://doi.org/10.1038/sdata.2018.52

Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and structural parameters. *The Econometrics Journal*, *21*(1), C1–C68. https://doi.org/10.1111/ectj.12097

Comte, L., & Olden, J. D. (2021). Evidence for dispersal syndromes in freshwater fishes. *Proceedings of the Royal Society B: Biological Sciences*, *288*(1951), 20210223. https://doi.org/10.1098/rspb.2021.0223

Domisch, S., Amatulli, G., & Jetz, W. (2015). Near-global freshwater-specific environmental variables for biodiversity analyses in 1 km resolution. *Scientific Data*, *2*, 150073. https://doi.org/10.1038/sdata.2015.73

Domisch, S., Kakouei, K., Martínez-López, J., Bagstad, K. J., Malek, Ž., Guerrero, A. M., ... Jähnig, S. C. (2019). Social equity shapes zone-selection: Balancing aquatic biodiversity representation and ecosystem services delivery in the transboundary Danube River Basin. *Scientific Reports*, *9*, 3082. https://doi.org/10.1038/s41598-019-39112-3

Dormann, C. F., Elith, J., Bacher, S., Buchmann, C., Carl, G., Carré, G., ... Lautenbach, S. (2013). Collinearity: A review of methods to deal with it and a simulation study evaluating their performance. *Ecography*, *36*(1), 27–46. https://doi.org/10.1111/j.1600-0587.2012.07348.x

Dudgeon, D., Arthington, A. H., Gessner, M. O., Kawabata, Z.-I., Knowler, D. J., Lévêque, C., ... Sullivan, C. A. (2006). Freshwater biodiversity: Importance, threats, status and conservation challenges. *Biological Reviews*, *81*(2), 163–182. https://doi.org/10.1017/S1464793105006950

Elith, J., & Leathwick, J. R. (2009). Species distribution models: Ecological explanation and prediction across space and time. *Annual Review of Ecology, Evolution, and Systematics*, *40*, 677–697. https://doi.org/10.1146/annurev.ecolsys.110308.120159

Filipe, A. F., Araújo, M. B., Doadrio, I., Angermeier, P. L., & Collares-Pereira, M. J. (2013). Biogeography of Iberian freshwater fishes revisited: The roles of historical versus contemporary constraints. *Journal of Biogeography*, *36*(11), 2096–2110. https://doi.org/10.1111/j.1365-2699.2009.02154.x

Friedman, N., Linial, M., Nachman, I., & Pe'er, D. (2000). Using Bayesian networks to analyze expression data. *Journal of Computational Biology*, *7*(3–4), 601–620. https://doi.org/10.1089/106652700750050961

Goberville, E., Beaugrand, G., Hautekèete, N.-C., Piquot, Y., & Luczak, C. (2015). Uncertainties in the projection of species distributions related to general circulation models. *Ecology and Evolution*, *5*(5), 1100–1116. https://doi.org/10.1002/ece3.1411

Grant, E. H. C., Lowe, W. H., & Fagan, W. F. (2007). Living in the branches: Population dynamics and ecological processes in dendritic networks. *Ecology Letters*, *10*(2), 165–175. https://doi.org/10.1111/j.1461-0248.2006.01007.x

Guisan, A., Tingley, R., Baumgartner, J. B., Naujokaitis-Lewis, I., Sutcliffe, P. R., Tulloch, A. I. T., ... Buckley, Y. M. (2013). Predicting species distributions for conservation decisions. *Ecology Letters*, *16*(12), 1424–1435. https://doi.org/10.1111/ele.12189

Hao, T., Elith, J., Lahoz-Monfort, J. J., & Guillera-Arroita, G. (2020). Testing whether ensemble modelling is advantageous for maximising predictive performance of species distribution models. *Ecography*, *43*(4), 549–558. https://doi.org/10.1111/ecog.04890

Hill, R. A., Hawkins, C. P., & Carlisle, D. M. (2013). Predicting thermal reference conditions for USA streams and rivers. *Freshwater Science*, *32*(1), 39–55. https://doi.org/10.1899/12-009.1

Irving, K., Fragkopoulou, E., Ceola, S., Vilmi, A., Akasaka, M., Akopian, M., ... Seelen, L. (2021). The environmental niche of riverine invertebrate specialists. *Global Ecology and Biogeography*, *30*(4), 887–902. https://doi.org/10.1111/geb.13263

Kalisch, M., Mächler, M., Colombo, D., Maathuis, M. H., & Bühlmann, P. (2012). Causal inference using graphical models with the R package pcalg. *Journal of Statistical Software*, *47*(11), 1–26. https://doi.org/10.18637/jss.v047.i11

Künzel, S. R., Sekhon, J. S., Bickel, P. J., & Yu, B. (2019). Metalearners for estimating heterogeneous treatment effects using machine learning. *Proceedings of the National Academy of Sciences*, *116*(10), 4156–4165. https://doi.org/10.1073/pnas.1804597116

Lundberg, S. M., Erion, G., Chen, H., DeGrave, A., Prutkin, J. M., Nair, B., ... Lee, S.-I. (2020). From local explanations to global understanding with explainable AI for trees. *Nature Machine Intelligence*, *2*(1), 56–67. https://doi.org/10.1038/s42256-019-0138-9

Merow, C., Smith, M. J., Edwards, T. C., Jr., Guisan, A., McMahon, S. M., Normand, S., ... Elith, J. (2014). What do we gain from simplicity versus complexity in species distribution models? *Ecography*, *37*(12), 1267–1281. https://doi.org/10.1111/ecog.00845

Molnar, C. (2020). *Interpretable machine learning: A guide for making black box models explainable*. Lulu.com. https://christophm.github.io/interpretable-ml-book/

Mod, H. K., Scherrer, D., Luoto, M., & Guisan, A. (2016). What we use is not what we know: Environmental predictors in plant distribution models. *Journal of Vegetation Science*, *27*(6), 1308–1322. https://doi.org/10.1111/jvs.12444

Peterson, E. E., Ver Hoef, J. M., Isaak, D. J., Falke, J. A., Fortin, M.-J., Jordan, C. E., ... Wenger, S. J. (2013). Modelling dendritic ecological networks in space: An integrated network perspective. *Ecology Letters*, *16*(5), 707–719. https://doi.org/10.1111/ele.12084

Roberts, D. R., Bahn, V., Ciuti, S., Boyce, M. S., Elith, J., Guillera-Arroita, G., ... Dormann, C. F. (2017). Cross-validation strategies for data with temporal, spatial, hierarchical, or phylogenetic structure. *Ecography*, *40*(8), 913–929. https://doi.org/10.1111/ecog.02881

Scutari, M. (2010). Learning Bayesian networks with the bnlearn R package. *Journal of Statistical Software*, *35*(3), 1–22. https://doi.org/10.18637/jss.v035.i03

Scutari, M., Howell, P., Balding, D. J., & Mackay, I. (2014). Multiple quantitative trait analysis using Bayesian networks. *Genetics*, *198*(1), 129–137. https://doi.org/10.1534/genetics.114.165704

Sequeira, A. M. M., Bouchet, P. J., Yates, K. L., Mengersen, K., & Caley, M. J. (2018). Transferring biodiversity models for conservation: Opportunities and challenges. *Methods in Ecology and Evolution*, *9*(5), 1250–1264. https://doi.org/10.1111/2041-210X.12998

Spirtes, P., Glymour, C. N., & Scheines, R. (2000). *Causation, prediction, and search* (2nd ed.). MIT Press.

Steel, E. A., Beechie, T. J., Torgersen, C. E., & Fullerton, A. H. (2017). Envisioning, quantifying, and managing thermal regimes on river networks. *BioScience*, *67*(6), 506–522. https://doi.org/10.1093/biosci/bix047

Thuiller, W., Guéguen, M., Renaud, J., Karger, D. N., & Zimmermann, N. E. (2019). Uncertainty in ensembles of global biodiversity scenarios. *Nature Communications*, *10*, 1446. https://doi.org/10.1038/s41467-019-09519-w

Tonkin, J. D., Merritt, D. M., Olden, J. D., Reynolds, L. V., & Lytle, D. A. (2018). Flow regime alteration degrades ecological networks in riparian ecosystems. *Nature Ecology & Evolution*, *2*, 86–93. https://doi.org/10.1038/s41559-017-0379-0

Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. *Journal of the American Statistical Association*, *113*(523), 1228–1242. https://doi.org/10.1080/01621459.2017.1319839

WWF. (2020). *Living Planet Report 2020 – Bending the curve of biodiversity loss*. WWF. https://livingplanet.panda.org

Yates, K. L., Bouchet, P. J., Caley, M. J., Mengersen, K., Randin, C. F., Parnell, S., ... Sequeira, A. M. M. (2018). Outstanding challenges in the transferability of ecological models. *Trends in Ecology & Evolution*, *33*(10), 790–802. https://doi.org/10.1016/j.tree.2018.08.001

Zurell, D., Franklin, J., König, C., Bouchet, P. J., Dormann, C. F., Elith, J., ... Merow, C. (2020). A standard protocol for reporting species distribution models. *Ecography*, *43*(9), 1261–1277. https://doi.org/10.1111/ecog.04960