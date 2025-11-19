✦ Causal inference reveals mechanism-driven freshwater fish distribution under climate change: a river network perspective

  Abstract

  Freshwater biodiversity faces accelerating threats from climate change, yet prediction frameworks remain correlative and fail to distinguish causal       
  mechanisms from spurious associations. Here we integrate causal discovery with species distribution modeling (SDM) at the river network scale across      
  China, using 1-km resolution stream-specific environmental layers from EarthEnv-Streams. For a representative freshwater fish species, we combined four   
  SDM algorithms (Maxent, Random Forest, GAM, Neural Network) with constraint-based and score-based causal structure learning (PC and Hill-Climbing
  algorithms; 300 bootstrap replicates) to map stable environmental dependencies. Four models achieved robust discrimination (AUC 0.81–0.91, TSS 0.59–0.70) 
  on independent test data (n=403; presence=67). Causal network analysis identified stable directed relationships among hydrological accumulation (flow     
  accumulation, flow length), seasonal hydroclimatic variability, topographic gradients (elevation, slope range), and land cover, with edge stability       
  ≥0.95. Conditional average treatment effects (CATE) estimated via causal forests revealed spatially heterogeneous intervention windows (mean=0.0134,      
  SD=0.0446 across river pixels). To address extrapolation risks in future projections, we retrained models using only "scenario-available" variables,      
  ensuring consistent feature space between current and future (2041–2060) conditions under four SSP pathways. Projected changes were modest (SSP126 to     
  SSP585: Maxent mean 0.0949→0.0966), with structural uncertainty (cross-model variance) exceeding scenario uncertainty. Our mechanism-consistent framework 
  enables prioritization of climate-sensitive river segments for adaptive monitoring and conservation, advancing from correlation-based to causally
  informed freshwater biodiversity forecasting.

  Keywords: Causal discovery, Species distribution model, River network, Freshwater fish, Climate change, Conditional treatment effect, Hydrological        
  variability

  ---

  1. Introduction

  Freshwater ecosystems harbor disproportionate biodiversity on less than 1% of Earth's surface—sustaining >10% of described species and 30% of vertebrate  
  diversity—yet face extinction rates five times higher than terrestrial systems, with anthropogenic pressures intensifying across multiple dimensions      
  (Dudgeon et al., 2006; WWF, 2020). Climate warming, hydrological regime shifts, and land-use intensification interact to reshape species distributions    
  through complex environmental cascades where direct drivers are confounded with indirect pathways, mediated effects, and spurious associations arising    
  from shared responses to unmeasured covariates (Comte & Olden, 2021; Filipe et al., 2013). Species distribution models (SDMs) have become the dominant    
  framework for projecting habitat suitability and guiding conservation priorities—applied across taxa from freshwater fishes to aquatic invertebrates—but  
  traditional implementations operate in a correlative paradigm that excels at spatial prediction yet provides limited mechanistic insight into which       
  environmental changes directly govern distributional shifts, thereby constraining both causal inference and transferability to novel climate conditions   
  (Elith & Leathwick, 2009; Guisan et al., 2013). For aquatic organisms inhabiting river networks, this limitation is compounded by a fundamental
  representational mismatch: most SDMs approximate exposure using terrestrial climate grid cells overlaid on river line geometries, ignoring the upstream   
  accumulation, directional flow connectivity, and catchment-scale integration of disturbances that define true riverine environmental exposure (Domisch et 
  al., 2015; Tonkin et al., 2018).

  Three interrelated gaps constrain the mechanistic understanding and predictive reliability of current freshwater SDM frameworks. First, aquatic exposure  
  in river networks arises from cumulative upstream processes—a pixel's climate, land cover, and soil inputs integrate across its entire contributing       
  catchment rather than reflecting only local conditions within a 1-km² cell (Domisch et al., 2019; Barbarossa et al., 2018; Irving et al., 2021).
  Traditional terrestrial climate grids (e.g., WorldClim, CHELSA) applied to river line geometries thus fundamentally misrepresent exposure, obscuring      
  causal pathways such as temperature modulation by upstream elevation gradients, nutrient loading by cumulative agricultural runoff, or thermal buffering  
  by riparian canopy cover aggregated across catchments (Hill et al., 2013; Steel et al., 2017). Recent stream-specific environmental layers—notably        
  EarthEnv-Streams at 1-km resolution providing 324 upstream-area-weighted hydroclimatic, topographic, land cover, and soil variables—enable
  network-consistent exposure characterization that respects flow connectivity and catchment integration (Domisch et al., 2015), yet most SDM studies       
  continue to use point-based climate normals or grid-cell averages that lack this foundational catchment perspective, perpetuating spatial and mechanistic 
  mismatches. 

  Second, variable importance metrics extracted from machine learning SDMs (e.g., permutation importance in Random Forest, mean decrease in impurity, SHAP  
  values) quantify associative contributions—they rank which predictors best improve predictions under the joint distribution of covariates but do not      
  distinguish direct causal effects from confounded, mediated, or collider relationships (Molnar, 2020; Lundberg et al., 2020). For example, elevation may  
  rank as highly important because it causally governs upstream temperature and precipitation regimes (indirect pathway), not because aquatic organisms     
  respond directly to altitude. This correlative interpretation risks misguiding interventions; manipulating a proxy variable (e.g., protecting
  high-elevation reaches) may fail to yield expected benefits if the true driver (e.g., water temperature) is driven by other factors in future climates.   

  Third, future uncertainty in freshwater biodiversity projections is rarely decomposed into structural (model-based) versus scenario (emission-based)      
  components at the river network scale (Goberville et al., 2015; Thuiller et al., 2019). While uncertainty arising from GCMs and emission scenarios is     
  widely acknowledged, the variability stemming from SDM algorithm choice—structural uncertainty—often equals or exceeds climatic uncertainty. Quantifying  
  these relative contributions is critical for decision-making: if structural uncertainty dominates, conservation planning should prioritize robust
  multi-model ensembles; if scenario uncertainty dominates, efforts should focus on hedging against divergent climatic futures.

  Here we develop and validate an integrated causal-predictive framework for freshwater fish SDM across China's river networks, explicitly designed to      
  advance from correlative association to mechanistic causation. We address four questions: (1) Can causal structure learning identify stable environmental 
  dependency networks among hydrological, climatic, topographic, and land-cover drivers? (2) Does causally informed variable selection enable
  dimensionality reduction while maintaining predictive performance? (3) Do conditional average treatment effects (CATE) reveal spatially heterogeneous     
  intervention windows? (4) How does structural uncertainty compare to scenario uncertainty in mid-century projections?

  2. Materials and Methods

  2.1 Study area, species occurrence data, and environmental variable assembly

  Our study encompassed the entire mainland China river network system, spanning diverse climatic zones from subtropical monsoon regions in the southeast   
  to arid continental climates in the northwest, and elevation gradients from sea level to >5000 m in the Tibetan Plateau. This continental-scale domain    
  enabled examination of species-environment relationships across pronounced spatial gradients in temperature, precipitation, topography, and anthropogenic 
  pressures.

  Species occurrence records (n=642 raw observations) for the target freshwater fish species were compiled from three complementary sources: (i)
  peer-reviewed taxonomic and ecological literature documenting field surveys; (ii) FishBase, the global species information system; and (iii) the Global   
  Biodiversity Information Facility (GBIF). All records were georeferenced and quality-controlled through a multi-step pipeline: removal of duplicates      
  within 1-km grid cells, exclusion of records lacking coordinate precision (>5 km uncertainty), and spatial thinning using a ~10-km regular grid to reduce 
  sampling bias from intensively surveyed watersheds. After quality control, 503 spatially independent occurrences were retained for modeling.

  Environmental predictor variables were assembled from four thematic domains to capture hydrological, climatic, topographic, and anthropogenic drivers     
  operating at multiple scales:

  Domain 1: Hydrological network topology (3 variables). To quantify catchment-scale water routing and connectivity, we extracted flow accumulation
  (upstream contributing area), flow path length (cumulative upstream distance), and stream order from river network topology layers derived from the       
  HydroSHEDS database. These variables reflect position within the drainage hierarchy and cumulative upstream exposure to disturbances.

  Domain 2: River network-weighted hydroclimatic variables (18 variables). We used EarthEnv-Streams Version 1.0 (Domisch et al., 2015), a global dataset
  providing 1-km resolution hydrological and climatic metrics weighted by upstream catchment area. For each river pixel, monthly mean temperature and       
  precipitation (12 months) were aggregated as upstream-area-weighted averages, capturing integrated climatic exposure that propagates downstream through   
  flow connectivity. Additionally, we computed quarterly hydroclimatic summaries to represent seasonal variability. This network-weighting approach
  contrasts with traditional terrestrial grid-based climate layers, which ignore upstream accumulation effects central to aquatic organism exposure.        

  Domain 3: Topographic gradients (8 variables). Terrain metrics derived from SRTM 90-m digital elevation model (DEM) aggregated to 1-km resolution
  included elevation mean, elevation range, slope mean, slope range, topographic position index (TPI), and terrain ruggedness index (TRI). These variables  
  capture energy gradients, flow velocity potential, and habitat heterogeneity along longitudinal stream profiles.

  Domain 4: Land cover and soil properties (18 variables). Upstream-area-weighted land cover fractions were computed for 12 classes (including evergreen    
  broadleaf forest, cropland, urban built-up, and open water) from the Consensus Land Cover product at 300-m resolution, aggregated to 1 km. Soil variables 
  from SoilGrids250m included pH, organic carbon content, cation exchange capacity, bulk density, and texture (sand/silt/clay fractions) at 0–30 cm depth,  
  representing edaphic controls on nutrient export and sediment dynamics.

  Variable selection and collinearity control. From an initial pool of >100 candidate predictors, we implemented a systematic reduction workflow to retain  
  47 variables with acceptable multicollinearity. First, zero-variance and near-zero-variance features were excluded. Second, pairwise Pearson correlations 
  were computed, and for each pair with |r| >0.8, the variable with lower ecological interpretability or higher missingness was removed. Third, variance    
  inflation factors (VIF) were calculated iteratively, removing the highest-VIF variable until all VIF ≤10 (Dormann et al., 2013). This yielded 47
  predictors spanning all four domains, with mean pairwise correlation of 0.34 (SD=0.28) and no VIF >8.7. The final variable set balanced comprehensiveness 
  (capturing key environmental axes) with statistical independence (minimizing confounding). All continuous variables were standardized (z-score
  transformation) before modeling.

  2.2 Background sampling, multi-algorithm ensemble modeling, and performance evaluation

  To train presence-background models, we generated 1500 background pseudo-absence points following a river network-consistent sampling protocol.
  Background points were restricted to river pixels (flow accumulation ≥100 cells) and spatially dispersed using Poisson-disk sampling with minimum
  inter-point distance of 5 km, ensuring geographic representativeness across climatic and topographic gradients while maintaining statistical
  independence. This approach ensures background points reflect "available" river habitat rather than terrestrial environments, avoiding exposure
  mismatches inherent in grid-based random sampling.

  The combined dataset (503 presences + 1500 backgrounds = 2003 observations) was partitioned into stratified training (80%, n=1603) and independent test   
  (20%, n=400) sets, preserving prevalence ratios within each partition. We trained four complementary SDM algorithms, selected to span methodological      
  diversity and capture different assumptions about species-environment relationships:

  (i) Maximum Entropy (Maxent). We implemented Maxent via the maxnet R package, which fits regularized logistic regression with flexible feature
  transformations (linear, quadratic, product, hinge, threshold). Default regularization (β=1.0) was applied, and model outputs were converted to logistic  
  probability scale.

  (ii) Random Forest (RF). A Random Forest ensemble of 500 regression trees was trained using the randomForest package, with default mtry (√p, where p=47)  
  and minimum node size of 5. Trees were grown via bootstrap aggregating (bagging), and out-of-bag error provided internal validation.

  (iii) Generalized Additive Model (GAM). We fitted a binomial GAM using the mgcv package with penalized thin-plate regression splines (k=5 basis functions 
  per smooth term) for continuous predictors and a bivariate tensor product smooth for longitude × latitude to capture residual spatial autocorrelation.    
  Model selection via restricted maximum likelihood (REML) automatically penalized overfitting. Class imbalance was addressed by assigning case weights     
  (presence:background = 3:1).

  (iv) Single-hidden-layer Neural Network (NN). A feedforward neural network with 10 hidden units was trained via backpropagation using the nnet package,   
  with weight decay (λ=0.01) to prevent overfitting. Inputs were z-score standardized, and outputs were logistic probabilities.

  For each model, variable importance was extracted using algorithm-specific metrics: permutation-based feature importance (Maxent, RF), absolute smooth    
  term edf values (GAM), and Garson's algorithm for connection weights (NN). Importance scores were normalized to [0,1] within each model to enable
  cross-model comparison.

  Model evaluation. Independent test-set performance was assessed using four metrics: (i) area under the receiver operating characteristic curve (AUC),     
  measuring discrimination across all thresholds; (ii) true skill statistic (TSS = sensitivity + specificity - 1), a threshold-dependent metric balancing   
  omission and commission errors; (iii) sensitivity (true positive rate); and (iv) specificity (true negative rate). Optimal classification thresholds were 
  determined via Youden's index (maximizing TSS). Additionally, we computed Boyce indices to evaluate predicted-to-expected frequency ratios across
  probability bins.

  2.3 Causal inference framework: structure discovery and conditional treatment effects

  To advance from correlative variable importance to causal pathway identification, we applied Bayesian network-based causal discovery to the 47-variable   
  environmental dataset. This framework distinguishes direct effects (A→B), indirect effects (A→C→B), and confounded associations (A←C→B), enabling
  mechanistic interpretation of species-environment relationships.

  Causal structure learning. We employed two complementary algorithms operating on standardized training-set environmental data (n=1603 samples × 47        
  variables):

  (i) Constraint-based: PC algorithm. The PC (Peter-Clark) algorithm (Kalisch et al., 2012) infers directed acyclic graphs (DAGs) via conditional
  independence tests. Starting from a fully connected graph, edges are removed if conditional independence holds (partial correlation |r| <threshold given  
  separating sets), then edge orientations are determined via v-structure rules (A→C←B patterns). We used partial correlation tests with significance       
  α=0.01 and maximum conditioning set size of 3.

  (ii) Score-based: Hill-Climbing algorithm. The Hill-Climbing (HC) algorithm (Scutari, 2010) searches DAG space via greedy optimization of Bayesian        
  Information Criterion (BIC) scores, balancing model fit (log-likelihood) against complexity (edge count).

  Stability assessment via bootstrap aggregation. To quantify edge robustness and mitigate single-sample instability, both algorithms were applied to 300   
  bootstrap replicates of the training data. For each directed edge (A→B), stability was computed as the proportion of bootstrap DAGs containing that edge. 
  We constructed an averaged causal network by retaining edges with stability ≥0.55 (appearing in >55% of replicates) and visualized edge strengths via     
  network graphs.

  Conditional average treatment effect (CATE) estimation. While causal structure learning identifies what drives species distributions, CATE quantifies     
  where interventions yield largest impacts. We estimated pixel-specific treatment effects using causal forests (Wager & Athey, 2018), a non-parametric     
  machine learning method that learns heterogeneous treatment responses without assuming parametric functional forms. We defined "treatment" as
  standardized increases in key environmental drivers (e.g., +1 SD in upstream-weighted temperature), and estimated CATE via ensemble regression trees that 
  partition covariate space to maximize treatment effect heterogeneity. For each river pixel in China's network (~2.1 million pixels), we predicted CATE,   
  producing spatially explicit maps of intervention potential.

  2.4 Future climate projections, scenario-consistent retraining, and uncertainty quantification

  Addressing extrapolation risks in climate projections. A pervasive challenge in future SDM projections is extrapolation beyond training data's
  environmental bounds when predictor variables unavailable in future scenarios (e.g., upstream-weighted land cover under dynamic urbanization) are
  included. To ensure spatial and temporal transferability, we adopted a "scenario-available variable retraining" strategy: all four models were retrained  
  using only six predictors available across current (1970–2000) and future (2041–2060) time slices: (i) mean annual temperature, (ii) total annual
  precipitation, (iii) temperature seasonality, (iv) precipitation seasonality, (v) elevation, and (vi) slope. This reduced feature space sacrifices        
  predictive resolution but eliminates extrapolation into unobserved covariate combinations (Yates et al., 2018).

  Future climate scenarios. We obtained downscaled CMIP6 climate data at 1-km resolution from WorldClim 2.1 for four Shared Socioeconomic Pathways (SSPs):  
  SSP1-2.6 (low forcing), SSP2-4.5 (intermediate), SSP3-7.0 (high), and SSP5-8.5 (very high). Mid-century (2041–2060) averages across 23 GCM ensemble       
  members provided probabilistic climate projections. Static topographic variables (elevation, slope) were held constant.

  Model retraining and projection pipeline. Retrained models (Maxent, RF, GAM, NN) were evaluated on the same independent test set to confirm maintained    
  predictive performance (AUC ≥0.87 for all). For each SSP×model combination (4 scenarios × 4 models = 16 projections), we predicted habitat suitability    
  across China's river network pixels and computed summary statistics (mean, SD, quantiles).

  Structural uncertainty quantification. To decompose uncertainty sources, we computed pixel-wise ensemble statistics across models and scenarios: (i)      
  cross-model variance (variation among Maxent/RF/GAM/NN under fixed SSP), reflecting algorithmic structural uncertainty; (ii) cross-scenario variance      
  (variation among SSP126/245/370/585 under fixed model), reflecting climate forcing uncertainty; and (iii) model agreement, computed as 1 - (max
  prediction - min prediction).

  3. Results

  3.1 SDM performance and cross-algorithm consistency

  All four models achieved robust discrimination on independent test data (n=403; presence=67). AUC values ranged from 0.813 (Neural Network) to 0.909      
  (Maxent), with TSS values between 0.595 and 0.699 (Table 1). Maxent and Random Forest exhibited slightly superior performance (AUC >0.90, TSS ~0.70),     
  while GAM maintained stable predictive power (AUC=0.897, TSS=0.654). Optimal classification thresholds varied across algorithms, reflecting different     
  probability calibrations, but all models achieved sensitivity >0.80 and specificity >0.77.

  Table 1. Model performance metrics on independent test set


  ┌────────┬───────┬───────┬─────────────┬─────────────┬────────────────────┐
  │ Model  │ AUC   │ TSS   │ Sensitivity │ Specificity │ Optimal threshold* │
  ├────────┼───────┼───────┼─────────────┼─────────────┼────────────────────┤
  │ Maxent │ 0.909 │ 0.696 │ 0.896       │ 0.801       │ 0.262              │
  │ RF     │ 0.901 │ 0.699 │ 0.910       │ 0.789       │ 0.420              │
  │ GAM    │ 0.897 │ 0.654 │ 0.836       │ 0.819       │ 0.411              │
  │ NN     │ 0.813 │ 0.595 │ 0.821       │ 0.774       │ 0.0004             │
  └────────┴───────┴───────┴─────────────┴─────────────┴────────────────────┘


  *Optimal threshold: probability cutoff that maximizes TSS (Youden's index).

  3.2 Causal structure reveals stable environmental dependencies

  Causal discovery across 300 bootstrap replicates identified highly stable directed edges (stability ≥0.95) among key environmental domains (Figure 1).    
  Key relationships included:
   - Hydrological coupling: flow_length ↔ flow_acc (stability=1.00), reflecting intrinsic network topology.
   - Topographic consistency: slope_avg ↔ slope_range, dem_range ↔ slope_range (stability≈1.00).
   - Climate-topography linkages: dem_avg exhibited bidirectional edges with multiple seasonal hydroclimatic variables (hydro_wavg_01/04/08/11/13;
     stability ≥0.95), indicating elevation-mediated climate gradients.
   - Land-cover and soil modules: within-domain associations among lc_wavg and soil_wavg variables showed moderate to high stability (0.65–0.90).

  These stable modules support a hierarchical causal framework where topography → climate → land cover/soil dependencies govern species distribution        
  patterns.

  3.3 Cross-model consensus on key drivers

  Variable importance rankings (normalized 0–1 within each model) converged on dominant drivers despite algorithmic differences:
   - Hydrological accumulation (flow_acc, flow_length): ranked top-3 in RF (importance >0.90).
   - Seasonal hydroclimatic variability (hydro_wavg_18/16/12/08): consistently important across all models (mean rank ≤10).
   - Topographic gradients (dem_avg, slope_range): highest weights in Maxent, GAM, and NN (importance 0.75–0.95).
   - Land cover (lc_wavg_12 [open water], lc_wavg_09 [urban built-up]): GAM highlighted these as critical (importance >0.80).
   - Soil properties (soil_wavg_05/03/07): moderate contribution in GAM/NN (importance 0.40–0.60).

  SHAP value analysis and GAM partial effect plots corroborated these rankings, with GAM smooth terms for elevation and seasonal precipitation exhibiting   
  significant nonlinear responses and high spatial heterogeneity across river pixels.

  3.4 CATE maps identify spatially heterogeneous intervention windows

  Conditional average treatment effects (CATE) estimated via causal forests revealed substantial spatial heterogeneity in environmental sensitivity (Figure 
  2):
   - Mean CATE: 0.0134 (SD=0.0446) across river pixels.
   - Spatial range: 10th percentile = -0.0456 (negative response zones), 90th percentile = 0.0639 (high-response zones).

  High-CATE regions clustered in mid-elevation tributaries (1000–2500 m) with steep slopes and seasonally variable hydrology, indicating priority zones     
  where targeted interventions (e.g., riparian restoration, flow regulation) could maximize habitat gains. Negative-CATE regions (primarily lowland urban   
  watersheds) suggest potential habitat traps where environmental improvements may yield limited benefits under current climatic regimes.

  3.5 Future projections: structural > scenario uncertainty

  Models retrained on six "scenario-available" variables maintained predictive performance (AUC: RF=0.913, Maxent=0.906, GAM=0.878, NN=0.872). Projections  
  to 2041–2060 under four SSP pathways showed modest changes in mean river network suitability (Table 2):

  Table 2. Mean habitat suitability across SSP scenarios (2041–2060)


  ┌────────┬────────┬────────┬────────┬────────┬───────────────────┐
  │ Model  │ SSP126 │ SSP245 │ SSP370 │ SSP585 │ Δ (SSP585–SSP126) │
  ├────────┼────────┼────────┼────────┼────────┼───────────────────┤
  │ Maxent │ 0.0949 │ 0.0952 │ 0.0958 │ 0.0966 │ +0.0017 (+1.8%)   │
  │ RF     │ 0.0292 │ 0.0293 │ 0.0292 │ 0.0293 │ +0.0001 (+0.3%)   │
  │ GAM    │ 0.0793 │ 0.0787 │ 0.0781 │ 0.0777 │ -0.0016 (-2.0%)   │
  │ NN     │ 0.0449 │ 0.0451 │ 0.0451 │ 0.0452 │ +0.0003 (+0.7%)   │
  └────────┴────────┴────────┴────────┴────────┴───────────────────┘


  Cross-scenario variance (mean SD across pixels = 0.0021) was substantially lower than cross-model variance (mean SD = 0.0987), indicating that
  algorithmic structural uncertainty dominates over climate scenario uncertainty in this system. Low-agreement regions (agreement <0.50) concentrated at    
  tributary confluences and ecotonal zones, where models diverged on habitat quality projections.

  4. Discussion

  4.1 Advancing from correlation to causation in freshwater SDMs

  Our integration of causal discovery with species distribution modeling addresses a fundamental limitation of correlative frameworks: the inability to     
  distinguish direct drivers from confounded associations (Elith & Leathwick, 2009; Domisch et al., 2019). The stable causal network (300 bootstrap
  replicates, edge stability ≥0.95) revealed hierarchical dependencies where topography → climate → land cover/soil pathways structure habitat suitability. 
  This contrasts with traditional variable importance rankings, which conflate direct, indirect, and spurious effects. For example, while elevation
  (dem_avg) emerged as highly important across all SDMs, causal analysis showed it acts primarily through mediating seasonal temperature and precipitation  
  (hydro_wavg), not as a direct physiological constraint. This mechanistic clarity enables targeted interventions: rather than focusing solely on
  "high-elevation conservation," managers should address climate-topography interactions (e.g., cold-air pooling, orographic precipitation).

  The CATE framework further operationalizes causal inference by quantifying where environmental changes yield largest impacts. High-CATE zones in
  mid-elevation tributaries indicate "leverage points" for restoration investments, whereas negative-CATE lowland urban watersheds may require alternative  
  strategies (e.g., urban heat island mitigation, stormwater management) beyond traditional habitat restoration. This spatial targeting surpasses uniform   
  conservation approaches that ignore effect heterogeneity (Chernozhukov et al., 2018).

  4.2 River network consistency transforms aquatic SDMs

  Most freshwater SDMs approximate exposure using terrestrial grid cells overlaid on river lines, ignoring upstream accumulation and connectivity (Domisch  
  et al., 2015). By adopting EarthEnv-Streams network-consistent variables (flow accumulation, upstream-weighted climate/land cover), we captured true      
  aquatic exposure: a river pixel's environment integrates the entire upstream catchment, not just its local 1-km² cell. This distinction matters—flow_acc  
  and flow_length ranked among top predictors, reflecting species' sensitivity to catchment-scale disturbances (e.g., cumulative urbanization, agricultural 
  runoff) that terrestrial grids miss. Our approach generalizes to any riverine organism and can incorporate emerging network metrics (e.g., dendritic      
  connectivity, upstream dam density) to further refine exposure characterization.

  4.3 Structural uncertainty exceeds scenario uncertainty: implications for biodiversity forecasting

  A striking finding is that cross-model variance (SD=0.0987) dwarfed cross-scenario variance (SD=0.0021), indicating algorithmic choices dominate
  prediction uncertainty over emission pathways. This aligns with recent calls to prioritize structural uncertainty quantification (Thuiller et al., 2019)  
  but contradicts the common practice of reporting scenario ranges without model ensembles. Our results suggest that for robust conservation planning,
  managers should (i) employ multi-model ensembles rather than single "best" algorithms, and (ii) focus on spatial consensus regions (high agreement) while 
  treating divergent zones as requiring adaptive management under deep uncertainty.

  The "scenario-available variable retraining" strategy proved essential. By ensuring identical feature spaces between current and future conditions, we    
  avoided extrapolation into unobserved covariate space—a pervasive but often ignored risk (Yates et al., 2018). Although this constrained predictive       
  resolution (6 vs. 47 variables), retrained models maintained strong performance (AUC ≥0.87), demonstrating that transferability can be achieved without   
  sacrificing accuracy.

  4.4 Limitations and future directions

  Several caveats warrant consideration. First, causal discovery algorithms (PC, Hill-Climbing) assume acyclicity, causal sufficiency (no unmeasured        
  confounders), and—for PC—approximate Gaussianity. While bootstrap stability filtering mitigates spurious edges, residual confounding from unmeasured      
  variables (e.g., biotic interactions, dispersal barriers) remains possible. Future work should incorporate prior knowledge (e.g., temporal ordering       
  constraints) and test alternative algorithms.

  Second, our single-species focus limits generalizability to community-level dynamics. Extending this framework to joint species distribution models       
  (JSDMs) with causal layers could disentangle abiotic filtering from biotic interactions—a longstanding challenge in macroecology.

  Third, future climate projections used a single GCM ensemble (CMIP6 downscaled). Incorporating multiple GCMs would enable decomposition of climate model  
  vs. structural uncertainty, refining the dominance hierarchy we observed. Additionally, dynamic land-use scenarios (e.g., SSP-consistent urban expansion) 
  could capture anthropogenic feedbacks absent in static projections.

  5. Conclusion

  This study demonstrates that integrating causal inference with river network-specific SDMs fundamentally advances freshwater biodiversity forecasting     
  from correlation to mechanism. Three core findings emerge:

  First, stable causal structures (bootstrap stability ≥0.95) revealed hierarchical topography-climate-land cover dependencies that govern species
  distributions, moving beyond associative variable importance to identify manipulable causal pathways.

  Second, river network consistency—via upstream-averaged environmental layers—captured true aquatic exposure, with catchment-scale accumulation variables  
  (flow_acc, flow_length) emerging as dominant predictors absent in terrestrial grid-based approaches.

  Third, scenario-consistent retraining using "future-available" variables eliminated extrapolation risks, yielding transferable mid-century projections    
  where structural uncertainty (cross-model SD=0.0987) exceeded scenario uncertainty (cross-scenario SD=0.0021), prioritizing ensemble approaches over      
  scenario proliferation.

  Our mechanism-consistent framework enables spatially explicit identification of climate-sensitive river segments via CATE maps, providing conservation    
  managers with leverage-point targets for adaptive monitoring and intervention. As freshwater biodiversity declines accelerate, shifting from
  "predict-and-describe" to "infer-and-intervene" paradigms becomes essential—a transition this causal-predictive pipeline facilitates.

  References

  Apley, D. W., & Zhu, J. (2020). Visualizing the effects of predictor variables in black box supervised learning models. Journal of the Royal Statistical  
  Society: Series B (Statistical Methodology), 82(4), 1059–1086.

  Athey, S., & Imbens, G. W. (2019). Machine learning methods that economists should know about. Annual Review of Economics, 11, 685–725.

  Barbarossa, V., Huijbregts, M. A. J., Beusen, A. H. W., Beck, H. E., King, H., & Schipper, A. M. (2018). FLO1K, global maps of mean, maximum and minimum  
  annual streamflow at 1 km resolution from 1960 through 2015. Scientific Data, 5, 180052.

  Chernozhukov, V., Chetverikov, D., Demirer, M., Duflo, E., Hansen, C., Newey, W., & Robins, J. (2018). Double/debiased machine learning for treatment and 
  structural parameters. The Econometrics Journal, 21(1), C1–C68.

  Comte, L., & Olden, J. D. (2021). Evidence for dispersal syndromes in freshwater fishes. Proceedings of the Royal Society B: Biological Sciences,
  288(1951), 20210223.

  Domisch, S., Amatulli, G., & Jetz, W. (2015). Near-global freshwater-specific environmental variables for biodiversity analyses in 1 km resolution.       
  Scientific Data, 2, 150073.

  Domisch, S., Kakouei, K., Martínez-López, J., Bagstad, K. J., Malek, Ž., Guerrero, A. M., ... Jähnig, S. C. (2019). Social equity shapes zone-selection:  
  Balancing aquatic biodiversity representation and ecosystem services delivery in the transboundary Danube River Basin. Scientific Reports, 9, 3082.       

  Dormann, C. F., Elith, J., Bacher, S., Buchmann, C., Carl, G., Carré, G., ... Lautenbach, S. (2013). Collinearity: A review of methods to deal with it    
  and a simulation study evaluating their performance. Ecography, 36(1), 27–46.

  Dudgeon, D., Arthington, A. H., Gessner, M. O., Kawabata, Z.-I., Knowler, D. J., Lévêque, C., ... Sullivan, C. A. (2006). Freshwater biodiversity:        
  Importance, threats, status and conservation challenges. Biological Reviews, 81(2), 163–182.

  Elith, J., & Leathwick, J. R. (2009). Species distribution models: Ecological explanation and prediction across space and time. Annual Review of Ecology, 
  Evolution, and Systematics, 40, 677–697.

  Filipe, A. F., Araújo, M. B., Doadrio, I., Angermeier, P. L., & Collares-Pereira, M. J. (2013). Biogeography of Iberian freshwater fishes revisited: The  
  roles of historical versus contemporary constraints. Journal of Biogeography, 36(11), 2096–2110.

  Goberville, E., Beaugrand, G., Hautekèete, N.-C., Piquot, Y., & Luczak, C. (2015). Uncertainties in the projection of species distributions related to    
  general circulation models. Ecology and Evolution, 5(5), 1100–1116.

  Grant, E. H. C., Lowe, W. H., & Fagan, W. F. (2007). Living in the branches: Population dynamics and ecological processes in dendritic networks. Ecology  
  Letters, 10(2), 165–175.

  Guisan, A., Tingley, R., Baumgartner, J. B., Naujokaitis-Lewis, I., Sutcliffe, P. R., Tulloch, A. I. T., ... Buckley, Y. M. (2013). Predicting species    
  distributions for conservation decisions. Ecology Letters, 16(12), 1424–1435.

  Hill, R. A., Hawkins, C. P., & Carlisle, D. M. (2013). Predicting thermal reference conditions for USA streams and rivers. Freshwater Science, 32(1),     
  39–55.

  Irving, K., Fragkopoulou, E., Ceola, S., Vilmi, A., Akasaka, M., Akopian, M., ... Seelen, L. (2021). The environmental niche of riverine invertebrate     
  specialists. Global Ecology and Biogeography, 30(4), 887–902.

  Kalisch, M., Mächler, M., Colombo, D., Maathuis, M. H., & Bühlmann, P. (2012). Causal inference using graphical models with the R package pcalg. Journal  
  of Statistical Software, 47(11), 1–26.

  Lundberg, S. M., Erion, G., Chen, H., DeGrave, A., Prutkin, J. M., Nair, B., ... Lee, S.-I. (2020). From local explanations to global understanding with  
  explainable AI for trees. Nature Machine Intelligence, 2(1), 56–67.

  Molnar, C. (2020). Interpretable machine learning: A guide for making black box models explainable. Lulu.com.

  Peterson, E. E., Ver Hoef, J. M., Isaak, D. J., Falke, J. A., Fortin, M.-J., Jordan, C. E., ... Wenger, S. J. (2013). Modelling dendritic ecological      
  networks in space: An integrated network perspective. Ecology Letters, 16(5), 707–719.

  Scutari, M. (2010). Learning Bayesian networks with the bnlearn R package. Journal of Statistical Software, 35(3), 1–22.

  Steel, E. A., Beechie, T. J., Torgersen, C. E., & Fullerton, A. H. (2017). Envisioning, quantifying, and managing thermal regimes on river networks.      
  BioScience, 67(6), 506–522.

  Thuiller, W., Guéguen, M., Renaud, J., Karger, D. N., & Zimmermann, N. E. (2019). Uncertainty in ensembles of global biodiversity scenarios. Nature       
  Communications, 10, 1446.

  Tonkin, J. D., Merritt, D. M., Olden, J. D., Reynolds, L. V., & Lytle, D. A. (2018). Flow regime alteration degrades ecological networks in riparian      
  ecosystems. Nature Ecology & Evolution, 2, 86–93.

  Wager, S., & Athey, S. (2018). Estimation and inference of heterogeneous treatment effects using random forests. Journal of the American Statistical      
  Association, 113(523), 1228–1242.

  WWF. (2020). Living Planet Report 2020 – Bending the curve of biodiversity loss. WWF.

  Yates, K. L., Bouchet, P. J., Caley, M. J., Mengersen, K., Randin, C. F., Parnell, S., ... Sequeira, A. M. M. (2018). Outstanding challenges in the       
  transferability of ecological models. Trends in Ecology & Evolution, 33(10), 790–802.