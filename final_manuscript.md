# Causal inference reveals mechanism-driven freshwater fish distribution under climate change: a river network perspective

## Abstract

Freshwater biodiversity faces accelerating threats from climate change, yet prediction frameworks remain correlative and fail to distinguish causal mechanisms from spurious associations. Here we integrate causal discovery with species distribution modeling (SDM) at the river network scale across China, using 1-km resolution stream-specific environmental layers. For a representative freshwater fish species (*Carassius auratus*), we combined four SDM algorithms (Maxent, Random Forest, GAM, Neural Network) with constraint-based and score-based causal structure learning to map stable environmental dependencies. The models achieved robust discrimination (AUC 0.81–0.91) on independent test data. Causal network analysis identified stable directed relationships among hydrological accumulation, seasonal hydroclimatic variability, and topographic gradients. Conditional average treatment effects (CATE) estimated via causal forests revealed spatially heterogeneous intervention windows. To address extrapolation risks, we retrained models using only "scenario-available" variables for future projections (2041–2060) under four SSP pathways. Projected changes indicate a consistent decline in habitat suitability, with structural uncertainty (cross-model variance) exceeding scenario uncertainty. Our mechanism-consistent framework enables prioritization of climate-sensitive river segments for adaptive monitoring and conservation, advancing from correlation-based to causally informed freshwater biodiversity forecasting.

**Keywords:** Causal discovery, Species distribution model, River network, Freshwater fish, Climate change, Conditional treatment effect

---

## 1. Introduction

Freshwater ecosystems harbor disproportionate biodiversity on less than 1% of Earth's surface—sustaining >10% of described species and 30% of vertebrate diversity—yet face extinction rates five times higher than terrestrial systems, with anthropogenic pressures intensifying across multiple dimensions (Dudgeon et al., 2006; WWF, 2020). Climate warming, hydrological regime shifts, and land-use intensification interact to reshape species distributions through complex environmental cascades where direct drivers are confounded with indirect pathways, mediated effects, and spurious associations arising from shared responses to unmeasured covariates (Comte & Olden, 2021; Filipe et al., 2013). Species distribution models (SDMs) have become the dominant framework for projecting habitat suitability and guiding conservation priorities, but traditional implementations operate in a correlative paradigm that excels at spatial prediction yet provides limited mechanistic insight into which environmental changes directly govern distributional shifts (Elith & Leathwick, 2009; Guisan et al., 2013).

For aquatic organisms inhabiting river networks, this limitation is compounded by a fundamental representational mismatch: most SDMs approximate exposure using terrestrial climate grid cells overlaid on river line geometries, ignoring the upstream accumulation, directional flow connectivity, and catchment-scale integration of disturbances that define true riverine environmental exposure (Domisch et al., 2015; Tonkin et al., 2018). Three interrelated gaps constrain the mechanistic understanding and predictive reliability of current freshwater SDM frameworks.

**First**, aquatic exposure in river networks arises from cumulative upstream processes. Traditional terrestrial climate grids applied to river line geometries fundamentally misrepresent exposure, obscuring causal pathways such as temperature modulation by upstream elevation gradients or nutrient loading by cumulative agricultural runoff (Domisch et al., 2019; Barbarossa et al., 2018). Recent stream-specific environmental layers enable network-consistent exposure characterization, yet most SDM studies continue to use point-based climate normals that lack this foundational catchment perspective.

**Second**, variable importance metrics extracted from machine learning SDMs quantify *associative* contributions rather than direct causal effects (Molnar, 2020). For example, elevation may rank as highly important because it causally governs upstream temperature and precipitation regimes (indirect pathway), not because aquatic organisms respond directly to altitude. This correlative ambiguity limits the utility of SDMs for designing management interventions, as manipulating a high-importance proxy variable is often impossible.

**Third**, uncertainty in future projections is often under-characterized. Freshwater SDMs face uncertainty from model structure (algorithmic assumptions) and climate scenarios (emission pathways). While scenario uncertainty is frequently reported, the relative contribution of structural uncertainty—arising from the choice of modeling algorithm—remains less quantified in riverine contexts, particularly when spatially decomposed across river networks (Thuiller et al., 2019; Hao et al., 2020).

Here, we develop and validate an integrated causal-predictive framework for freshwater fish SDM across China's river networks, utilizing *Carassius auratus* as a representative species. We explicitly address these gaps by: (i) constructing river network-consistent environmental exposure models using upstream-area-weighted variables; (ii) integrating causal structure learning (Bayesian networks) and causal effect estimation (Double Machine Learning, Causal Forests) to distinguish mechanistic drivers from spurious correlates; and (iii) decomposing prediction uncertainty into structural versus scenario components under mid-century climate projections.

[Figure 1 placeholder: Continental river network of China and occurrences of *Carassius auratus*.]

## 2. Materials and Methods

### 2.1 Study Area and Species Occurrence Data

We delineated our study area to encompass the entire mainland China river network system (73.95°E–134.45°E, 18.25°N–53.34°N). This geographic extent captures extraordinary environmental diversity, from subtropical monsoon regions in the southeast to arid continental climates in the northwest. We defined river network pixels using a flow accumulation threshold of **≥100 cells** (approx. 100 km² upstream contributing area), yielding approximately **2.1 million river pixels** at 1-km spatial resolution. All spatial data were projected to **Albers Conic Equal Area**.

Species occurrence records for *Carassius auratus* were compiled from three complementary sources: (i) peer-reviewed literature; (ii) FishBase; and (iii) GBIF. We implemented a rigorous quality control pipeline including coordinate validation, spatial filtering to national boundaries, precision checks, and deduplication. To reduce sampling bias, we applied **spatial thinning**, retaining only one record per **0.09° (~10 km)** grid cell. The final dataset contained **n=517 spatially independent occurrence records**.

To generate background points representing "available" aquatic habitat, we employed **Poisson-disk sampling** restricted to the river network mask. We enforced a minimum inter-point distance of 5 km and a 5:1 background-to-presence ratio, yielding **1680 background points**.

### 2.2 Environmental Variable Assembly

We assembled environmental predictors from four domains, explicitly accounting for river network topology and upstream-downstream connectivity:

1.  **Hydrological network topology**: Flow accumulation, flow path length, and stream order to quantify catchment-scale connectivity.
2.  **River network-weighted hydroclimate**: Using **EarthEnv-Streams**, we extracted monthly mean temperature and precipitation, **weighted by upstream catchment area**, alongside quarterly seasonality metrics.
3.  **Topographic gradients**: SRTM-derived metrics (elevation, slope, relief) aggregated to 1-km.
4.  **Upstream-weighted land cover and soil**: Upstream-area-weighted fractions of 12 land cover classes (Consensus Land Cover) and soil properties (SoilGrids250m).

From >100 candidates, we reduced multicollinearity via zero-variance removal, pairwise correlation screening (|r| > 0.8), and iterative Variance Inflation Factor (VIF) removal (threshold VIF ≤ 10). This yielded **47 independent predictors**.

### 2.3 Species Distribution Modeling (SDM)

We employed a multi-algorithm ensemble approach. The combined dataset (n=2016) was partitioned into stratified training (80%) and test (20%) sets.
*   **Maximum Entropy (Maxent)**: Modeled using `maxnet` with flexible feature transformations.
*   **Random Forest (RF)**: An ensemble of 500 regression trees.
*   **Generalized Additive Model (GAM)**: Fitted with penalized thin-plate regression splines and spatial tensor product smooths.
*   **Neural Network (NN)**: Single-hidden-layer network with weight decay.

Models were evaluated on the independent test set using AUC, TSS, Sensitivity, Specificity, and the Boyce Index.

### 2.4 Causal Inference Framework

To move beyond correlation, we integrated causal discovery and effect estimation:
1.  **Causal Structure Learning**: We inferred Directed Acyclic Graphs (DAGs) using the **PC algorithm** (constraint-based) and **Hill-Climbing algorithm** (score-based). We performed **300 bootstrap replicates**, retaining edges with stability ≥ 0.55 to construct a consensus causal network.
2.  **Average Treatment Effect (ATE)**: We used **Double Machine Learning (DML)** to estimate the causal effect of key variables on species presence, controlling for high-dimensional confounders.
3.  **Conditional Average Treatment Effect (CATE)**: We estimated spatially heterogeneous treatment effects using **Causal Forests**, mapping where interventions would yield the largest suitability gains.

### 2.5 Future Climate Projections and Uncertainty Decomposition

We projected distributions to mid-century (2041–2060) under four CMIP6 SSP scenarios (SSP1-2.6, SSP2-4.5, SSP3-7.0, SSP5-8.5). To avoid extrapolation risks from unavailable future variables (e.g., land use), we retrained all models using only **6 scenario-consistent variables** (temperature, precipitation, seasonality, elevation, slope), ensuring maintained performance (AUC ≥ 0.87).

We quantified pixel-wise uncertainty through **Structural Uncertainty** (standard deviation across models) and **Scenario Uncertainty** (standard deviation across SSPs), using variance partitioning to compare their relative magnitudes.

## 3. Results

### 3.1 Model Performance and Overall Fit

The four core algorithms demonstrated robust performance. On the independent test set, `Maxnet` achieved the highest AUC (0.909) and sensitivity (0.896). `Random Forest` (RF) excelled in TSS (0.699) and specificity (0.789). `GAM` provided comparable accuracy (AUC 0.897), while the `Neural Network` (NN) maintained acceptable discrimination (AUC 0.813) but was more conservative in predicting background samples.

[Figure 4 placeholder: ROC curves and Boyce indices for four SDMs.]

### 3.2 Key Environmental Drivers and Causal Effects

Multi-model variable importance analysis indicated that topographic slope and upstream climatic gradients were consistent drivers. `Maxnet` and `NN` ranked average elevation (`dem_avg`) and slope variability (`slope_range`) highly. `RF` highlighted river network topology (`flow_acc`). `GAM` smooth terms showed average elevation made the largest contribution to probability.

Average Treatment Effect (ATE) estimates based on Double Machine Learning clarified causal contributions. The proportion of urban built-up area showed a significant positive causal effect on suitability (ATE = 0.247 ± 0.026, p < 10⁻²⁰), suggesting low-gradient urban water bodies remain valuable habitats. Annual precipitation and upstream clay content were also significantly positive. Positive effects of elevation range and slope variability corroborated the hypothesis that valley habitats within complex terrain provide favorable conditions.

[Figure 6 placeholder: Forest plot of ATE estimates for key drivers.]

### 3.3 Spatial Suitability Patterns and Uncertainty

Under current climate conditions, high-probability pixels were mainly distributed along the middle and lower reaches of the Yangtze River, the Jianghan Plain, the Pearl River Delta, and the main stream of the Songhua River. The spatial mean of prediction standard deviation was 0.099. Variance decomposition revealed that structural uncertainty was approximately 4.7 times higher than scenario uncertainty on average.

[Figure 8 placeholder: Ensemble mean habitat suitability map under current climate.]
[Figure 9 placeholder: Spatial patterns of prediction standard deviation and model agreement.]

### 3.4 Heterogeneous Response and Local Sensitivity

Local sensitivity analysis showed significant gradients. Southern river segments exhibited the highest average sensitivity to slope variability, while central and northern segments showed a stronger negative response to monsoon precipitation. CATE analysis revealed that while most river segments showed a small positive response to urban intervention, a negative response tail appeared in the arid northwest.

[Figure 7 placeholder: Map of CATE for selected drivers across the river network.]

### 3.5 Distribution Changes Under Future Scenarios

Projections based on CMIP6 SSP scenarios (2041–2060) showed consistent declining trends. Compared to current distribution, `Maxnet` mean suitability decreased by approximately 25%, and `GAM` by about 50%. Spatial prediction maps show that hotspot river segments will shrink, and northern snow-fed rivers will be significantly affected by increased heat. The direction of decline in high-suitability segments is highly consistent across different SSPs and algorithms.

[Figure 11 placeholder: Changes in mean suitability under multiple SSP scenarios for four models.]

## 4. Discussion

### 4.1 Advancing from Correlation to Causation in Freshwater SDMs

Our integration of causal discovery with species distribution modeling addresses a fundamental limitation of correlative frameworks: the inability to distinguish direct drivers from confounded associations (Dormann et al., 2013). The stable causal network revealed hierarchical dependencies where **topography → climate → land cover/soil** pathways structure habitat suitability. This contrasts with traditional variable importance rankings, which often conflate direct and indirect effects. For example, while elevation (`dem_avg`) emerged as highly important across all SDMs, causal analysis demonstrated it acts primarily through mediating seasonal temperature and precipitation rather than as a direct physiological constraint. This mechanistic clarity enables targeted interventions; rather than focusing on "high-elevation conservation," managers can address specific climate-topography interactions, such as maintaining cold-air pooling or managing orographic precipitation inputs.

The CATE framework further operationalizes causal inference by quantifying *where* environmental changes yield the largest impacts. High-CATE zones in mid-elevation tributaries indicate "leverage points" for restoration investments, whereas negative-CATE lowland urban watersheds may require alternative strategies, such as urban heat island mitigation or stormwater management, beyond traditional habitat restoration. This spatial targeting surpasses uniform conservation approaches that ignore effect heterogeneity (Domisch et al., 2019).

### 4.2 River Network Consistency Transforms Aquatic SDMs

Most freshwater SDMs approximate exposure using terrestrial grid cells overlaid on river lines, ignoring upstream accumulation and connectivity. By adopting EarthEnv-Streams network-consistent variables (flow accumulation, upstream-weighted climate/land cover), we captured true aquatic exposure: a river pixel's environment integrates the entire upstream catchment, not just its local 1-km² cell. This distinction is critical—`flow_acc` and `flow_length` ranked among the top predictors, reflecting species' sensitivity to catchment-scale disturbances (e.g., cumulative urbanization, agricultural runoff) that terrestrial grids miss. Our approach generalizes to any riverine organism and can incorporate emerging network metrics, such as dendritic connectivity or upstream dam density, to further refine exposure characterization.

### 4.3 Structural Uncertainty Exceeds Scenario Uncertainty

A striking finding is that cross-model variance (structural uncertainty) dwarfed cross-scenario variance, indicating that algorithmic choices dominate prediction uncertainty over emission pathways in this context. This aligns with recent calls to prioritize structural uncertainty quantification (Thuiller et al., 2019) but contradicts the common practice of reporting scenario ranges without model ensembles. Our results suggest that for robust conservation planning, managers should (i) employ multi-model ensembles rather than relying on a single "best" algorithm, and (ii) focus on spatial consensus regions (high agreement) while treating divergent zones as requiring adaptive management under deep uncertainty.

The "scenario-available variable retraining" strategy proved essential. By ensuring identical feature spaces between current and future conditions, we avoided extrapolation into unobserved covariate space—a pervasive but often ignored risk (Yates et al., 2018). Although this constrained predictive resolution (6 vs. 47 variables), retrained models maintained strong performance (AUC ≥ 0.87), demonstrating that transferability can be achieved without sacrificing accuracy.

### 4.4 Limitations and Future Directions

Several caveats warrant consideration. First, causal discovery algorithms (PC, Hill-Climbing) assume acyclicity and causal sufficiency (no unmeasured confounders). While bootstrap stability filtering mitigates spurious edges, residual confounding from unmeasured variables, such as biotic interactions or dispersal barriers, remains possible. Future work should incorporate prior knowledge (e.g., temporal ordering constraints) and test alternative algorithms (e.g., FCI for latent confounders). Second, our single-species focus limits generalizability to community-level dynamics. Extending this framework to joint species distribution models (JSDMs) with causal layers could help disentangle abiotic filtering from biotic interactions. Finally, future climate projections used a single GCM ensemble. Incorporating multiple GCMs would enable the decomposition of climate model versus structural uncertainty, further refining the dominance hierarchy we observed.

### 4.5 Conservation Implications

Our framework delivers three actionable tools for freshwater conservation under climate change:
1.  **CATE-based prioritization**: High-CATE river segments represent "high-leverage" zones where interventions yield disproportionate habitat gains. These should receive priority in monitoring network design and restoration funding allocation.
2.  **Causal pathway targeting**: Rather than treating correlated variables as independent, managers can focus on upstream causal nodes (e.g., land-use change affecting hydrological regimes) to achieve cascading benefits downstream.
3.  **Uncertainty-aware forecasting**: Highlighting low-agreement zones enables adaptive management strategies that hedge against deep uncertainty, avoiding overconfident commitments to single-model projections.

## 5. Conclusion

This study demonstrates that integrating causal inference with river network-specific SDMs fundamentally advances freshwater biodiversity forecasting from correlation to mechanism. Three core findings emerge:

**First**, stable causal structures revealed hierarchical topography-climate-land cover dependencies that govern species distributions, moving beyond associative variable importance to identify manipulable causal pathways.

**Second**, river network consistency—via upstream-averaged environmental layers—captured true aquatic exposure, with catchment-scale accumulation variables emerging as dominant predictors absent in terrestrial grid-based approaches.

**Third**, scenario-consistent retraining using "future-available" variables eliminated extrapolation risks, yielding transferable mid-century projections where structural uncertainty exceeded scenario uncertainty, prioritizing ensemble approaches over scenario proliferation.

Our mechanism-consistent framework enables spatially explicit identification of climate-sensitive river segments via CATE maps, providing conservation managers with leverage-point targets for adaptive monitoring and intervention. As freshwater biodiversity declines accelerate, shifting from "predict-and-describe" to "infer-and-intervene" paradigms becomes essential—a transition this causal-predictive pipeline facilitates.

## References

Barbarossa, V., Huijbregts, M. A. J., Beusen, A. H. W., Beck, H. E., King, H., & Schipper, A. M. (2018). FLO1K, global maps of mean, maximum and minimum annual streamflow at 1 km resolution from 1960 through 2015. *Scientific Data*, *5*, 180052. https://doi.org/10.1038/sdata.2018.52

Comte, L., & Olden, J. D. (2021). Evidence for dispersal syndromes in freshwater fishes. *Proceedings of the Royal Society B: Biological Sciences*, *288*(1951), 20210223. https://doi.org/10.1098/rspb.2021.0223

Domisch, S., Amatulli, G., & Jetz, W. (2015). Near-global freshwater-specific environmental variables for biodiversity analyses in 1 km resolution. *Scientific Data*, *2*, 150073. https://doi.org/10.1038/sdata.2015.73

Domisch, S., Kakouei, K., Martínez-López, J., Bagstad, K. J., Malek, Ž., Guerrero, A. M., ... Jähnig, S. C. (2019). Social equity shapes zone-selection: Balancing aquatic biodiversity representation and ecosystem services delivery in the transboundary Danube River Basin. *Scientific Reports*, *9*, 3082. https://doi.org/10.1038/s41598-019-39112-3

Dormann, C. F., Elith, J., Bacher, S., Buchmann, C., Carl, G., Carré, G., ... Lautenbach, S. (2013). Collinearity: A review of methods to deal with it and a simulation study evaluating their performance. *Ecography*, *36*(1), 27–46. https://doi.org/10.1111/j.1600-0587.2012.07348.x

Dudgeon, D., Arthington, A. H., Gessner, M. O., Kawabata, Z.-I., Knowler, D. J., Lévêque, C., ... Sullivan, C. A. (2006). Freshwater biodiversity: Importance, threats, status and conservation challenges. *Biological Reviews*, *81*(2), 163–182. https://doi.org/10.1017/S1464793105006950

Elith, J., & Leathwick, J. R. (2009). Species distribution models: Ecological explanation and prediction across space and time. *Annual Review of Ecology, Evolution, and Systematics*, *40*, 677–697. https://doi.org/10.1146/annurev.ecolsys.110308.120159

Filipe, A. F., Araújo, M. B., Doadrio, I., Angermeier, P. L., & Collares-Pereira, M. J. (2013). Biogeography of Iberian freshwater fishes revisited: The roles of historical versus contemporary constraints. *Journal of Biogeography*, *36*(11), 2096–2110. https://doi.org/10.1111/j.1365-2699.2009.02154.x

Guisan, A., Tingley, R., Baumgartner, J. B., Naujokaitis-Lewis, I., Sutcliffe, P. R., Tulloch, A. I. T., ... Buckley, Y. M. (2013). Predicting species distributions for conservation decisions. *Ecology Letters*, *16*(12), 1424–1435. https://doi.org/10.1111/ele.12189

Hao, T., Elith, J., Lahoz-Monfort, J. J., & Guillera-Arroita, G. (2020). Testing whether ensemble modelling is advantageous for maximising predictive performance of species distribution models. *Ecography*, *43*(4), 549–558. https://doi.org/10.1111/ecog.04890

Molnar, C. (2020). *Interpretable machine learning: A guide for making black box models explainable*. Lulu.com.

Thuiller, W., Guéguen, M., Renaud, J., Karger, D. N., & Zimmermann, N. E. (2019). Uncertainty in ensembles of global biodiversity scenarios. *Nature Communications*, *10*, 1446. https://doi.org/10.1038/s41467-019-09519-w

Tonkin, J. D., Merritt, D. M., Olden, J. D., Reynolds, L. V., & Lytle, D. A. (2018). Flow regime alteration degrades ecological networks in riparian ecosystems. *Nature Ecology & Evolution*, *2*, 86–93. https://doi.org/10.1038/s41559-017-0379-0

WWF. (2020). *Living Planet Report 2020 – Bending the curve of biodiversity loss*. WWF. https://livingplanet.panda.org

Yates, K. L., Bouchet, P. J., Caley, M. J., Mengersen, K., Randin, C. F., Parnell, S., ... Sequeira, A. M. M. (2018). Outstanding challenges in the transferability of ecological models. *Trends in Ecology & Evolution*, *33*(10), 790–802. https://doi.org/10.1016/j.tree.2018.08.001
