---
title: ""
author: ""
date: ""
fontsize: 12pt
output:
  bookdown::gitbook:
    number_sections: yes
    split_by: chapter  
    config:
      toc:
        collapse: section
        scroll_highlight: yes
        before: null
        after: null
      toolbar:
        position: fixed
      edit: null
      download: ["pdf"]
      search: yes
      fontsettings:
        theme: white
        family: sans
        size: 2
      sharing:
        facebook: yes
        twitter: yes
        google: no
        linkedin: no
        weibo: no
        instapper: no
        vk: no
        all: ['facebook', 'google', 'twitter', 'linkedin', 'weibo', 'instapaper']
bibliography: /home/ghislain/Documents/Bibliography/biblio.bib
biblio-style: bib/jae.bst
link-citations: yes
csl: bib/journal-of-applied-ecology.csl
---


\linenumbers
\newpage

<!--chapter:end:index.Rmd-->

# Materials and Methods

## Study-areas

We defined 119 study-areas representing 92 countries (Table \@ref(tab:samp-size)) and covering the integrality of the tropical moist forest in the world, at the exception of some islands (eg. Sao Tome and Principe or Wallis-and-Futuna). Each country was identified by one unique three-letter code following the ISO 3166-1 standard (eg. MDG for Madagascar or GUF for French Guiana). Most of the countries corresponded to one unique study-area, with three exceptions: Brazil, India, and Australia. Brazil, because of its large size, was divided into 26 study-areas corresponding to the 26 administrative states (the state of Goias including the Federal District). For India, which is also a large country, the tropical moist forest is located in three distinct regions far from each other. We thus considered three independent study-areas for India: the Western Ghats, North-East India (including the West Bengal), and the union territory of the Andaman and Nicobar Islands. For Australia, we only considered the Queensland state as a study-area. Data sampling and spatial deforestation modelling were performed independently for each study-area. Study-area borders were obtained from version 3.6 of the Global Administrative Areas database (<https://gadm.org>). We used level-0 data for study-areas corresponding to countries and level-1 data for study-areas corresponding to states or regions. We grouped the study-areas in three continents (Fig. \@ref(fig:study-areas)): America (64 study-areas for 39 countries), Africa (32 study-areas for 32 countries), and Asia (23 study-areas for 21 countries).

## Historical forest cover change maps

For each study-area, we derived historical forest cover change maps on two periods of time: 1$^{\text{st}}$ January 2000 -- 1$^{\text{st}}$ January 2010, and 1$^{\text{st}}$ January 2010 -- 1$^{\text{st}}$ January 2020 from the forest cover change annual product by @Vancutsem2020. The annual product by @Vancutsem2020 classifies Landsat image pixels at 30 m resolution in 16 categories for each year (on the 31$^{\text{st}}$ of December) between 1982 and 2019 and allows identifying tropical moist forest pixels at each date (Table \@ref(tab:cat-annual-product)). This classification is based on an expert model analyzing time-series data at the pixel level extracted from the full Landsat satellite image archive on the period 1982--2019. The expert model was built on Google Earth Engine [@Gorelick2017]. For our forest definition, we only considered natural old-growth tropical moist forests, disregarding plantations and regrowths. We included degraded forests (not yet deforested) in our forest definition. As a consequence, we considered pixels in the following categories: 1, 2, 3, 4, 5, 13, or 14 in the annual product, to be natural old-growth tropical moist forest pixels (simply abbreviated "forest" in this manuscript). Because several decades are usually necessary to reach the state of old-growth forest, we assumed every pixel classified as "forest" at a given date between 1999 and 2019 to be also classified as "forest" in the previous years of that period of time. We thus obtained three forest cover maps for the dates 1$^{\text{st}}$ January 2000, 1$^{\text{st}}$ January 2010, and 1st January 2020. We combined these three maps to obtain high-resolution forest cover change maps in the periods 2000--2010--2020 at 30 m resolution in the humid tropics (Fig. \@ref(fig:fcc-maps)). We used Google Earth Engine to process the annual product by @Vancutsem2020 and derive the historical forest cover change map for each country. An interactive forest cover change map for the humid tropics is available at <https://forestatrisk.cirad.fr/tropics>.

We did not consider potential forest regrowth in our forest definition for three main reasons. First, throughout the humid tropics, forest regeneration involves much smaller areas than deforestation [@Vancutsem2020, see also 2000-2012 tree cover gain in @Hansen2013]. Second, there is little evidence of natural forest regeneration in the long term in the tropics [@Grouzis2001]. This can be explained by several ecological processes following deforestation such as soil erosion [@Grinand2017] and reduced seed bank due to fire-induced deforestation and soil loss [@Grouzis2001]. Moreover, in areas where forest regeneration is ecologically possible, young forest regrowths are more easily re-burnt for agriculture and pasture [@Vieilledent2020]. Third, young secondary forests generally provide more limited ecosystem services compared to old-growth natural forests in terms of biodiversity [@Gibson2011] and carbon storage [@Blanc2009].

## Spatial explanatory variables

To explain the observed deforestation on the period 2010--2020, we considered a set of spatial explanatory variables describing: topography (altitude and slope), accessibility (distances to nearest road, town, and river), forest landscape (distance to forest edge), deforestation history (distance to past deforestation), and land conservation status (presence of a protected area). This set of variables were selected based on an _a priori_ knowledge of the deforestation process [@Gorenflo2011; @Vieilledent2013; @Geist2002; @Brown1994]. For example, the risk of deforestation is supposed to decrease with the distance to road and forest edge (lower accessibility), to increase at lower elevation and slope (higher probability to find arable lands), and to decrease in protected areas (higher level of protection). Characteristics of each explanatory variable are summarized in Table \@ref(tab:variables).

Elevation (in m) and slope (in degree) at 90 m resolution were obtained from the SRTM Digital Elevation Database v4.1 (<http://srtm.csi.cgiar.org/>). Distances (in m) to nearest road, town and river at 150 m resolution were computed from the road, town and river network which were obtained from the OpenStreetMap (OSM) project (<https://www.openstreetmap.org/>). OSM country data were downloaded from two web-sites: Geofabric (<http://download.geofabrik.de/>) and OpenStreetMap.fr (<https://download.openstreetmap.fr/extracts/>) depending on the availability of the data for each country. To obtain the road network in each country, we considered the "motorway", "trunk", "primary", "secondary" and "tertiary" categories for the "highway" key in OSM. To obtain the network of populated places in each country (that we simply call "towns" in the present study), we considered the "city", "town" and "village" categories for the "place" key in OSM. To obtain the river network, we considered the "river" and "canal" categories for the "waterway" key in OSM. For a more detailed description of each category, see the OSM wiki page (<https://wiki.openstreetmap.org/wiki/Tags>). OSM data have been downloaded during the period January--March 2020 for all countries. Distance to forest edge was computed at 30 m resolution from the forest cover map in 2010. Distance to past deforestation in 2010 was computed at 30 m resolution from the 2000--2010 forest cover change map. To minimize border effect for the computation of distance to forest edge and distance to past deforestation, a buffer of 10 km around each study-area extent was considered. Data on protected areas were obtained from the World Database on Protected Areas (<https://www.protectedplanet.net>, @WDPA2020) using the `pywdpa` Python package (<https://pypi.org/project/pywdpa/>). WDPA data have been downloaded during the period January--March 2020 for all countries. All protected areas defined by at least one polygon were considered in the analysis, but protected areas defined by a point were not taken into account. Data included protected areas of all IUCN categories (from Ia to VI) and of all types defined at the national level (e.g. National Parks, Reserves), even if the type and IUCN category were not reported. Polygons representing protected areas were rasterized at 30 m resolution.

In total, we obtained 8 spatial explanatory variables to model the spatial probability of deforestation (Table \@ref(tab:variables)).

## Data sampling for spatial modelling of the deforestation

With the spatial model of deforestation, our aim was to estimate the effects of a set of variables in determining the location of the deforestation (or "allocation" census @Pontius2011) and compute the relative probability of deforestation for each forest pixel. With the spatial model, our objective was not to estimate the intensity of the deforestation (or "quantity" census @Pontius2011), that could be expressed in %/year or in ha/year for example. A balanced sampling between deforested and non-deforested pixels is preferable in this case [@Dezecache2017; @Vieilledent2013]. Because deforestation events are rare ($\approx$ 1 %/yr), a non-stratified random sampling would lead to very few observations of deforestation events, rendering difficult a good estimation of the effects of the explanatory variables. Stratified balanced sampling provided unbiased estimates of the model's parameters, except for the model's intercept (estimated average deforestation). Having a biased model intercept (which has the same value for all forest pixels) is not a problem as we are interested in estimating a _relative_ probability of deforestation for all forest pixels.

As a consequence, we performed a stratified balanced sampling between (i) forest pixels in 2010 which have been deforested on the period 2010--2020 ("deforested" pixels), and (ii) forest pixels in 2010 which have not been deforested on that period of time and which represent the remaining forest in 2020 ("non-deforested" pixels). Forest pixels in each category were sampled randomly (Fig. \@ref(fig:sampling)). To maximize the representativity of the data, the total number of forest pixels sampled in each study-area for the year 2010 was chosen proportionally to the area of forest in 2010 in that study-area (2000 points for 1 Mha of forest), with the condition that this number had to be between 20,000 (to be representative of the deforestation process) and 100,000 (to limit computation time). When, for a specific study-area, the total number of pixels in one of the two categories (deforested vs. non-deforested pixels) was $\leq$ 10,000, all the pixels of that category were included in the sample. This could happen for study-areas with low moist forest cover such as small islands (eg. Antigua and Barbuda). For each sampled pixel, we retrieved information regarding the 8 computed explanatory variables at their original spatial resolution. When the information was not complete for a given pixel (eg. elevation and slope data missing for a forest pixel located close to the sea border), the observation was removed from the data-set. Missing information affected a minority of pixels. The global data-set included a total of 3,186,698 observations (1,601,810 of non-deforested pixels and 1,584,888 of deforested pixels, corresponding to an area of 144,163 ha and 142,647 ha, respectively).

## Spatial deforestation model

Using observations of forest cover change in the period 2010--2020, we modelled the spatial probability of deforestation as a function of the $n$ explanatory variables using a logistic regression. We considered the random variable $y_i$ which takes value 1 if the forest pixel $i$ was deforested in the period 2010--2020 and 0 if it was not. We assumed that $y_i$ follows a Bernoulli distribution of parameter $\theta_i$ (Eq. \@ref(eq:icar)). In our model, $\theta_i$ represents the spatial relative probability of deforestation for pixel $i$. We assumed that $\theta_i$ is linked, through a logit function, to a linear combination of the explanatory variables $X_i \beta$, where $X_i$ is the vector of explanatory variables for pixel $i$, and $\beta$ is the vector of effects $[\beta_1, \ldots, \beta_n]$ associated to the $n$ variables. All the continuous explanatory variables were normalized before fitting the model. The model includes an intercept $\alpha$. To account for the residual spatial variation in the deforestation process, we included an additional random effect $\rho_{j(i)}$ for each spatial cell $j$ of a 10 $\times$ 10 km grid covering each study-area (Fig. \@ref(fig:grid)). This grid resolution was chosen in order to have a reasonable balance between a good representation of the spatial variability of the deforestation process and a limited number of parameters to estimate. A sampled forest pixel $i$ was associated to one cell $j$ and one random effect $\rho_{j(i)}$. We assumed that random effects were spatially autocorrelated through an intrinsic conditional autoregressive (iCAR) model [@Besag1991; @Banerjee2014]. This model is denoted "icar" in subsequent sections and results. In an iCAR model, the random effect $\rho_j$ associated to cell $j$ depends on the values of the random effects $\rho_{j^{\prime}}$ associated to neighbouring cells $j^{\prime}$. In our case, the neighbouring cells are connected to the target cell $j$ through a common border or corner (cells defined by the "king move" in chess, see Fig. \@ref(fig:grid)). The variance of the spatial random effects $\rho_j$ was denoted $V_{\rho}$. The number of neighbouring cells for cell $j$, which might vary, was denoted $n_j$. Spatial random effects $\rho_j$ account for unmeasured or unmeasurable variables [@Clark2005] that explain a part of the residual spatial variation in the deforestation process that is not explained by the fixed spatial explanatory variables ($X_i$).

\begin{equation}
\begin{split}
  y_i \sim \mathcal{B}ernoulli(\theta_i)\\
  \text{logit}(\theta_i) = \alpha + X_i \beta + \rho_{j(i)}\\
  \rho_{j(i)} \sim \mathcal{N}ormal(\sum_{j^{\prime}} \rho_{j^{\prime}} / n_j,V_{\rho} / n_j)
\end{split}
(\#eq:icar)
\end{equation}

## Variable selection

Variable selection was performed using a backward elimination procedure. All the spatial explanatory variables in our data-set should decrease the deforestation risk (having a negative effect on the probability of deforestation). For example, the probability of deforestation should decrease with the distance to the forest edge and should also decrease inside a protected area. Our variable selection procedure was thus not based on statistical significance but on background knowledge regarding the deforestation process and on the interpretability of the variable effects [@Heinze2018]. For each study-area, we started to fit a model with the full set of explanatory variables (8 variables). At each step of the procedure, we removed the variables having a positive effect on the probability of deforestation. This was done in order to avoid unrealistic predictions of the spatial probability of deforestation at the scale of the study-area. For example, it is not realistic to observe, at a country scale, a decrease of the deforestation with the distance to forest edge (higher deforestation risk in the core of the forest compared with forest edge). This might happen in a particular context for a very specific region but is very unlikely at a country scale. In most of the cases, when we found a positive effect for a given explanatory variable, it was non-significant (95% credible interval including zero). On the contrary, when we found a negative but non-significant effect for a given explanatory variable, we kept this variable in the model. This effect, albeit non-significant, was interpretable given our background knowledge of the deforestation process and was relatively lower than the effects associated with the other variables.

## Parameter inference

Parameter inference was done in a hierarchical Bayesian framework. Non-informative priors were used for all parameters: $\alpha \sim \mathcal{N}ormal(\text{mean}=0,\text{var}=10^6)$, $\beta \sim \mathcal{N}ormal(\text{mean}=0,\text{var}=10^6)$, and $V_{\rho} \sim 1/\mathcal{G}amma(\text{shape}=0.05,\text{rate}=0.0005)$. During the variable selection procedure, we run a Markov Chain Monte Carlo (MCMC) of 2000 iterations, discarding the first 1000 iterations (burn-in phase). For the final model, we repeated the parameter inference using a longer MCMC of 10,000 iterations. We discarded the first 5000 iterations (burn-in phase), and we thinned the chain each 5 iterations (to reduce autocorrelation between samples). We obtained 1000 estimates for each parameter. MCMC convergence was visually checked looking at MCMC traces and parameter posterior distributions. Function `model_binomial_iCAR()` from the `forestatrisk` Python package was used for parameter inference. This function calls an adaptive Metropolis-within-Gibbs algorithm [@Rosenthal2011] written in C for maximum computation speed.

## Model comparison

### Alternative models

We compared the performance of the "icar" model at predicting deforestation with three other models: a null model (denoted "null"), a simple generalized linear model ("glm"), and a random forest model ("rf"). The "null" model assumes that all the slope parameters have value zero (no effect of explanatory variables), and that all the spatial random effects have also value zero (no residual regional variability in the deforestation process). For the "null" model, the probability of deforestation is only determined by a mean intercept common to every forest pixel ($\text{logit}(\theta_i) = \alpha$). The simple "glm" is a logistic regression which does not include spatial random effects (no residual regional variability in the deforestation process). For the "glm" model, the probability of deforestation is only determined by the mean intercept $\alpha$ and the parameters $\beta$ associated to each explanatory variables ($\text{logit}(\theta_i) = \alpha + X_i \beta$). Using simple "glm" model is a commonly proposed approach for spatial modelling of deforestation [@Verburg2002; @Mas2007; @Rosa2014; @Ludeke1990; @Soares-Filho2002; @Mertens1997; @Soares-Filho2001; @Eastman2017]. The random forest model [@Breiman2001] is a machine learning approach using an ensemble of random classification trees (where both observations and features are chosen at random to build the classification trees) to predict the deforestation probability for a forest pixel. Random forest has been intensively used for species distribution modelling [@Thuiller2009] and is now also commonly used for spatial modelling of deforestation [@Grinand2020; @Zanella2017; @Santos2019]. The "glm" and "rf" models were fitted using functions `LinearRegression` and `RandomForestClassifier` respectively, both available in the `scikit-learn` Python package [@Pedregosa2011]. We used the same set of selected explanatory variables for the "glm" and "rf" models as those used for the final "icar" model. For the "rf" model, we set the number of random classification trees to 500.

### Percentage of deviance explained

We computed the deviance $\mathcal{D}$ of the four models ("icar", "null", "glm", and "rf") with the formula $\mathcal{D}=-2 \log \mathcal{L}$, $\mathcal{L}$ being the likelihood of the model, i.e. the probability of observing the data given the model and estimated parameters. The deviance is a measure of error and a model with a lower deviance fits better the data. We also considered the deviance of the "full" model (also called the saturated model) which has as many parameters as there are observations. When $y_i=0$, the deforestation probability predicted by the full model is 0. When $y_i=1$ the deforestation probability predicted by the full model is 1. The deviance of the full model is then equal to 0. Considering that the "null" model explained 0% of the deviance and the "full" model explains 100% of the deviance, we then computed the percentage of deviance explained by each of the three other models: "icar", "glm", and "rf".

### Cross-validation procedure

To compare the performance of the "icar", "glm", and "rf" models at predicting correctly the relative probability of deforestation on independent observations, we also performed a five-fold cross-validation procedure. We used 70% of the observations for the model training and 30% of the observations for the model validation. We used the fitted models to predict the deforestation probability of all the forest pixels of the validation data-set. To transform the deforestation probabilities into binary values, we identified the probability threshold respecting the percentage of deforested pixels in the validation data-set (eg. the mode of the predicted probabilities for a percentage of 50% of deforested pixels). Consequently, the predicted number of deforested pixels was equal to the observed number of deforested pixels in the validation data-set. This implies that there was no "quantity disagreement" (sensus @Pontius2008) for any of the three models in the cross-validation procedure. Through this cross-validation, we only compared the ability of the models to correctly identify the pixels to be deforested, given a particular mean deforestation rate. This corresponds to estimating the "allocation disagreement" (sensus @Pontius2008) for each of the three models. Using model predictions and observations in the validation data-set, we computed several accuracy indices: the Area Under the ROC Curve (AUC), the Figure of Merit (FOM), the Overall Accuracy (OA), the Expected Accuracy (EA), the Kappa of Cohen (K), the Specificity (Spe), the Sensitivity(Sen), and the True Skill Statistics (TSS). A detailed description of these indices can be found in @Pontius2008 (for the FOM) and @Liu2011 (for all the other indices). Formulas used to compute these indices are presented in Table \@ref(tab:confusion-matrix) and Table \@ref(tab:accuracy-indices).

### Model selection

For all the study-areas, we found that the percentage of deviance explained for the "rf" model was much higher than for the "icar" and "glm" models (Table \@ref(tab:accuracy)), suggesting that the "rf" model was fitting much better the data than the two other models. Nevertheless, when looking at the results of the cross-validation procedure, the "rf" model had a lower accuracy in comparison with the "icar" model (Table \@ref(tab:accuracy)). These two results show clearly that the "rf" model was overfitting the data and was less performant at predicting the probability of deforestation at new sites than the "icar" model. This can be explained by the fact that the strong nonlinear relationship between explanatory variables and the spatial probability of deforestation, which is estimated by the "rf" model based on the training data-set, does not represent the true relationship at the landscape scale. This limitation associated with machine-learning techniques such as random forest has already been raised in previous scientific articles dealing with large scale mapping of ecological variables [@Ploton2020]. On the contrary, the "icar" model had the highest accuracy indices of the three models suggesting a better predictive performance than the "glm" and "rf" models. Also, the "icar" model increased the explained deviance from 38.8 to 52.2% in average compared with the "glm" model. This shows that environmental explanatory variables alone explain a relative small part of the spatial deforestation process and that including spatial random effects to account for unexplained residual spatial variability strongly improves model fit (+13.4% of deviance explained in average) and model predictive performance (+6.9% for the TSS for example). Same results were obtained when comparing accuracy indices for the three statistical models per continent (Table \@ref(tab:accuracy-cont)). We thus selected the "icar" model for predicting the spatial probability of deforestation for all the study-areas.

## Computing the spatial probability of deforestation for the year 2020

Before computing the predictions of the deforestation probability, the spatial random effects at 10 km were interpolated at 1 km using a bicubic interpolation method. This was done in order to obtain spatial random effects at a resolution closer to the original forest raster resolution of 30 m, and to smooth the deforestation probability predictions spatially.

Distance to forest edge in 2020 was recomputed at 30 m resolution from the forest cover map in 2020. Distance to past deforestation in 2020 was computed at 30 m resolution from the 2010--2020 forest cover change map. All other explanatory variables (protected areas, distance to nearest road, town and river, elevation, and slope) were supposed unchanged between years 2010 and 2020. Using rasters of explanatory variables at their original resolution, interpolated spatial random effects at 1 km resolution, and the fitted "icar" model for each study-area, we computed the spatial probability of deforestation at 30 m resolution for the year 2020 for each study-area.

Deforestation probabilities (float values in the interval $[0, 1]$) were rescaled and transformed as integer values on the interval $[\![1, 65535]\!]$. This allowed us to record the large rasters of probabilities as UInt16 type (using zero as no-data value) and save space on disk. We then obtained a map of the relative probability of deforestation for the year 2020 at 30 m resolution. An interactive global map of the spatial probability of deforestation is available at <https://forestatrisk.cirad.fr/tropics>.

## Forecasting forest cover change on the period 2020--2100 {#forecast}

For each study-area, we computed the observed mean annual deforested area $d$ (in ha/yr) on the recent period 2010--2020 from the forest cover maps at these two dates (Tables \@ref(tab:fcc-hist) and \@ref(tab:fcc-hist-reg)). To forecast the forest cover change at a particular date $y$ in the future, we computed the estimated total deforestation $D_y$ (in ha) between years 2020 and $y$ assuming a "business-as-usual" scenario. The "business-as-usual" scenario makes the assumption of an absence of change in the deforestation intensity in the future (no increase in the deforestation intensity that could be attributable to a future increase in the demand of agricultural commodities for example, nor decrease in the deforestation intensity that could be attributable to new conservation policies or increase in agricultural yields for example). The "business-as-usual" scenario also makes the assumption that the spatial deforestation process will remain the same in the future (assuming a constant effect of the spatial explanatory variables in the future) and that areas with a higher relative probability of deforestation will remain the same in the future. To compute the total deforestation $D_y$ (in ha) between years 2020 and $y$, we projected the mean annual deforestation $d$ on the time interval $y-2020$: $D_y=d \times (y-2020)$.

For Brazil, which was divided into 26 study-areas, the mean annual deforested area $d$ was supposed constant at the country level, not at the study-area level. This is because we assumed that deforestation inside Brazil should spread between study-areas and should not stop at the study-area administrative borders. As a consequence, when all the forest of a specific study-area in Brazil was deforested, the corresponding residual deforested area for that study-area was redistributed to the other study-areas of Brazil still having forest. This ensured that the annual deforested area was constant for Brazil, but imply an increase of the annual deforested area with time for some study areas. We assumed that this contagious deforestation between study-areas was only valuable for connected study areas inside a specific country (here, Brazil). On the contrary, we assumed that it was much less likely to observe contagious deforestation between two neighbouring countries. First, because of the presence of less permeable international borders, and second, because the socio-economic factors driving the intensity of deforestation between two neighbouring countries can be significantly different (see contrasting historical deforestation intensity for Haiti--Dominican Republic, Congo--DRC, French Guiana--Suriname, or Indonesia--Papua New Guinea for example). We assumed no contagious deforestation for the three study-areas in India which are not connected (Fig. \@ref(fig:study-areas)).

The map of the relative spatial probability of deforestation in 2020 has a resolution of 30 m equivalent to an area of $r_{\text{ha}}=0.09$ ha. Using this map, we computed a probability threshold $p_y$ in the interval $[\![1, 65535]\!]$ identifying the $n_y$ forest pixels in 2020 with the highest probability of deforestation so that $n_y \times r_{\text{ha}} = D_y + \epsilon$. Because deforestation probabilities have finite values in the interval $[\![1, 65535]\!]$, some forest pixels might have the same deforestation probability and it might not be possible to identify $p_y$ such that $\epsilon=0$. We thus selected the threshold $p_y$ minimizing $\epsilon$. Because $D_y$ represents the total deforested area for several years and that few pixels had the same probability of deforestation, we always obtained negligible $\epsilon$ compared to $D_y$ ($\epsilon << D_y$). We considered those $n_y$ forest pixels in 2020 as deforested between years 2020 and $y$, and we derived the corresponding forest cover change map for the period 2020--$y$.

We projected the future forest cover for several years in the period 2030--2100 using a base time-interval of 10 years (Tables \@ref(tab:fcc-proj) and \@ref(tab:fcc-proj-reg)). We also projected the future forest cover for years 2055 and 2085 as these two years are often used as pivot years in studies on future climate change. They correspond to 30-year climate averages on the periods 2040--2070 and 2070--2100, and to mid-term and long-term climate projections, respectively, for the 21st century [@IPCC2014]. We then computed the percentage of forest cover loss in 2100 in comparison with the forest cover in 2000 for each study-area and continent (Tables \@ref(tab:fcc-proj) and \@ref(tab:fcc-proj-reg)). We also computed the year at which all the forest will have disappeared for each study area (Tables \@ref(tab:fcc-proj)). At the continental level, it makes less sense to compute the year at which all the forest will have disappeared, as some countries might conserve forest for a very long time, even though they account for a very small proportion of the total forest area at the continental scale. Instead, we computed the estimated year at which 75% of the forest cover in 2000 will have disappeared (Table \@ref(tab:fcc-proj-reg) and Fig. \@ref(fig:perc-loss)).

## Carbon emissions associated to deforestation

We estimated the carbon emissions associated with historical deforestation (2010--2020) and projected deforestation (2030--2100). To do so, we used the pantropical 1 km resolution aboveground dry biomass (AGB, in Mg/ha) map (Fig. \@ref(fig:agb)) by @Avitabile2016. This map is a combination of two pantropical aboveground biomass maps by @Saatchi2011 and @Baccini2012. The fusion map is representative of the aboveground biomass for the years 2000--2010. While the RMSE (root mean square error) of the fused map by @Avitabile2016 is still substantial (87--98 Mg/ha), the fused map achieved a lower RMSE (a decrease of 5--74%) and bias (a decrease of 90--153%) than the two input maps for all continents. We used the Intergovernmental Panel on Climate Change (IPCC) default carbon fraction of 0.47 [@McGroddy2004] to convert aboveground dry biomass to carbon stocks. We assumed no change of the forest carbon stocks in the future while computing carbon emissions associated with projected deforestation. We estimated average annual carbon emissions for ten-year periods from 2010 to 2100 (Tab. XXX). Under a "business-as-usual" scenario of deforestation (no change in the annual deforested area, in ha/yr, in the future), the change in mean annual carbon emissions in the future is only attributable to the spatial variation of the forest carbon stocks and to the location of the future deforestation.

## Software used: the `forestatrisk` Python module

We developed a specific module called `forestatrisk` to model and forecast tropical deforestation spatially using the Python programming language [@Python2020]. The module is available either on GitHub at <https://github.com/ghislainv/forestatrisk> or PyPI (The Python Package Index) at <https://pypi.org/project/forestatrisk>. It can be easily installed using `pip` (the package installer for Python) in any Python virtual environment created with either `virtualenv` or `conda`. The `forestatrisk` module includes functions (i) to build a data-set from deforestation observations (function `.sample`), (ii) to estimate the parameters of several spatial deforestation models (functions `.model*`, including function `.model_binomial_iCAR` for the "icar" model), (iii) to assess the model performance (function `.cross_validation`), (iv) to derive predictive maps of the probability of deforestation (functions `.predict_raster*`), and (v) to forecast future forest cover under a given intensity of deforestation (function `.deforest`). Using functions from the `forestatrisk` Python package makes computation fast and efficient (with low memory usage) by treating large raster data by blocks. Numerical computations on blocks of data are performed with the NumPy (<https://numpy.org>) Python module whose core is mostly made of optimized and compiled C code which runs fast [@Harris2020].

We also developed another smaller Python module called `pywdpa` that allows downloading the shapefile of the protected areas for each country using the API of the World Database on Protected Areas (<https://api.protectedplanet.net>). The module is also available on GitHub at <https://github.com/ghislainv/pywdpa> or PyPI at <https://pypi.org/project/pywdpa>, and can also be easily installed using `pip`.

Computations for each country were run in parallel on the computing cluster of the Montpellier Bioinformatics Biodiversity (MBB) platform (https://mbb.univ-montp2.fr) provided by LabEx CeMEB (http://www.labex-cemeb.org/). Google Earth Engine [@Gorelick2017] was used to process the annual product by @Vancutsem2020 and derive the historical forest cover change map for each country.

While the raw results of our study (which implied intensive computations on large raster data) were obtained using the Python programming language, the summarized results (tables and figures) presented in this article (main text and supplementary materials) were obtained using the R software [@R2020].

## Repeating the analysis with the tree cover loss product by @Hansen2013

We performed the exact same analysis replacing historical forest cover change maps derived from the annual product of @Vancutsem2020 by forest cover change maps derived from the tree cover products of @Hansen2013. We assumed that the extent of the natural old-growth tropical moist forest cover in 2000 could be approximated by using a tree cover >70% at this date [@Vieilledent2018; @Aleman2017a]. Then, we used the tree cover loss provided by @Hansen2013 on the periods 2001--2009 and 2010--2019 to derive pantropical forest cover change maps on the periods 2000--2010--2020 at 30 m resolution. All the results we obtained were robust to the change of the historical forest cover change map. Results obtained with the forest cover change maps derived from the tree cover products of @Hansen2013 are available on the CIRAD Dataverse repository (see below).

## Reproducibility of the results

All the data and code used for this study have been made permanently and publicly available on the CIRAD Dataverse repository so that the results are entirely reproducible.

- Input data: XXX
- Code: XXX
- Output data: XXX

Note that the code used for this study is also available on GitHub at <https://github.com/ghislainv/forestatrisk-tropics>.

\nolinenumbers
\newpage

<!--chapter:end:01-Supplementary_Materials.Rmd-->

# Figures

<!--------------------------------------------->
<!-- Study-areas -->
<!--------------------------------------------->

## Study-areas in the three continents

(ref:cap-study-areas) **Study-areas in the three continents: America, Africa, and Asia**. America included 64 study-areas (39 countries), Africa included 32 study areas (32 countries), and Asia included 23 study-areas (21 countries). Each country was identified by one unique three-letter code following the ISO 3166-1 standard (eg. MDG for Madagascar or GUF for French Guiana). In America, Brazil was divided in 26 study areas corresponding to the 26 Brazilian states. Each Brazilian state was defined by one unique two-letter code (eg. AM for Amazonas). For India, three study areas were considered: the Whestern Ghats (WG), the North-East India (NE), and the Andaman and Nicobar Islands (AN). For Australia, we only considered the Queensland (QLD) state as a study-area. In the three figures, each study-area is identified by one unique code and a set of polygons with the same color. The horizontal lines on each figure indicate the position of the Equator (plain line) and the two tropics (Cancer at the North and Capricorn at the South, dashed lines).


\begin{center}\includegraphics[width=\textwidth]{figures/study_areas_America} \end{center}


\begin{center}\includegraphics[width=\textwidth]{figures/study_areas_Africa} \end{center}

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/study_areas_Asia} 

}

\caption{(ref:cap-study-areas)}(\#fig:study-areas)
\end{figure}

<!--------------------------------------------->
<!-- Historical forest cover change map -->
<!--------------------------------------------->

## Historical forest cover change map

(ref:cap-fcc-maps) **Historical forest cover change map**. Forest cover change map on the period 2000--2010--2020 for the Democratic Republic of the Congo in central Africa. \textcolor{orange}{orange}: 2000--2010 deforestation, \textcolor{red}{red}: 2010--2020 deforestation, \textcolor{darkgreen}{green}: forest cover in 2020. Forest cover change map was derived from the forest cover change annual product by @Vancutsem2020. Original resolution of the forest cover change map is 30 m. The inset at the bottom left shows a zoom of the map for an area at the North-East of the country which is close to the city of Beni and the Virunga national park. An interactive global forest cover change map is available at <https://forestatrisk.cirad.fr/tropics>.

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/fcc123} 

}

\caption{(ref:cap-fcc-maps)}(\#fig:fcc-maps)
\end{figure}

<!--------------------------------------------->
<!-- Spatial explanatory variables -->
<!--------------------------------------------->

## Spatial explanatory variables used for spatial modelling of deforestation

(ref:cap-var) **Spatial explanatory variables**. Spatial explanatory variables for the Democratic Republic of the Congo in central Africa. Elevation (in m) and slope (in degree) at 90 m resolution were obtained from the SRTM Digital Elevation Database v4.1 (<http://srtm.csi.cgiar.org/>). Distances (in m) to nearest road, town and river at 150 m resolution were computed from the road, town and river network obtained from OpenStreetMap (OSM) (<https://www.openstreetmap.org/>). Roads include "motorway", "trunk", "primary", "secondary" and "tertiary" roads from OSM. Towns include "city", "town" and "village"  categories from OSM. Rivers include "river" and "canal" categories from OSM. Protected areas were obtained from the World Database on Protected Areas (<https://www.protectedplanet.net>, @WDPA2020). Data included protected areas of all IUCN categories (from Ia to VI) and of all types defined at the national level (e.g. National Parks, Reserves). Two additional spatial explanatory variables (distance to forest edge and distance to past deforestation) were obtained from the historical forest cover change map (Fig. \@ref(fig:fcc-maps)).

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/var} 

}

\caption{(ref:cap-var)}(\#fig:var)
\end{figure}

(ref:cap-corr-var) **Correlation between explanatory variables**. To compute the correlations, we used a representative data-set at the global scale where the number of observations for each study-area was proportional to its forest cover in 2010. The total number of observations changed from 3,186,698 to 813,796. We computed the Pearson's correlation matrix for the seven continuous explanatory variables used to model the deforestation: elevation ("elev"), slope ("slope"), distance to nearest road, town, and river ("droad", "dtown", and "driver", respectively), distance to forest edge ("dedge"), and distance to past deforestation ("ddefor"). For the protected areas ("pa"), which is a categorical variable for which a Pearson's correlation coefficient cannot be computed, we reported the slope coefficient of simple logistic regressions where the probability of presence of a protected area was a function of one intercept and one of the continuous variable (which was normalized).

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/corr-var} 

}

\caption{(ref:cap-corr-var)}(\#fig:corr-var)
\end{figure}

(ref:cap-data-pa) **Pantropical data-set on protected areas**. Protected areas were downloaded from the World Database on Protected Areas (<https://www.protectedplanet.net>, @WDPA2020) using the `pywdpa` Python package during the period January-March 2020. Data included protected areas of all IUCN categories (from Ia to VI) and of all types defined at the national level (e.g. National Parks, Reserves).

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/pa} 

}

\caption{(ref:cap-data-pa)}(\#fig:data-pa)
\end{figure}

(ref:cap-data-roads) **Pantropical road network**. The road network was obtained from OpenStreetMap (OSM) (<https://www.openstreetmap.org/>). Roads included "motorway", "trunk", "primary", "secondary" and "tertiary" roads from OSM.

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/roads} 

}

\caption{(ref:cap-data-roads)}(\#fig:data-roads)
\end{figure}

<!--------------------------------------------->
<!-- Data sampling -->
<!--------------------------------------------->

## Data sampling

(ref:cap-sampling) **Data sampling for spatial modelling of the deforestation**. Map on the left corresponds to the top left inset in Fig. \@ref(fig:fcc-maps) representing a zoom of the forest cover change map on the period 2000--2010--2020 for an area at the North-East of the Democratic Republic of the Congo. Map on the right presents an inner zoom showing the delimitation of the 30 m forest pixels with three sample points. We used a stratified balanced sampling between (i) forest pixels in 2010 which have been deforested on the period 2010--2020 ("deforested" pixels in \textcolor{red}{red}), and (ii) forest pixels in 2010 which have not been deforested on that period of time and which represent the remaining forest in 2020 ("non-deforested" pixels in \textcolor{darkgreen}{green}). Forest pixels in each category were sampled randomly.

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/sample} 

}

\caption{(ref:cap-sampling)}(\#fig:sampling)
\end{figure}

<!--------------------------------------------->
<!-- Grid for spatial random effects -->
<!--------------------------------------------->

## Grid for spatial random effects

(ref:cap-grid) **Grid used to compute the spatial random effects**. _Main figure_: $10 \times 10$ km grid covering the Democratic Republic of the Congo (DRC). The grid over DRC includes 45,154 $10 \times 10$ km cells (214 cells on the x axis by 211 cells on the y axis). The background map shows the historical forest cover change on the periods 2000--2010--2020 (see Fig. \@ref(fig:fcc-maps)). _Top inset_: Zoom for an area at the North-East of the country (black square) showing specific grid cells. One grid cell can include several sample points (see Fig. \@ref(fig:sampling)). _Bottom inset_: One random effect $\rho_j$ is estimated for each grid cell $j$. Spatial autocorrelation is taken into account through an intrinsic CAR process: the value of the random effect for one cell depends on the values of the random effects $\rho_{j^{\prime}}$ for the eight neighbouring cells $j^{\prime}$ (see Eq. \@ref(eq:icar)).

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/grid} 

}

\caption{(ref:cap-grid)}(\#fig:grid)
\end{figure}

<!--------------------------------------------->
<!-- Spatial random effects -->
<!--------------------------------------------->

## Estimated spatial random effects

(ref:cap-rho) **Estimated spatial random effects**. _Left_: Estimated spatial random effects at 10 km resolution for the Democratic Republic of the Congo (DRC). _Right_: Interpolated spatial random effects at 1 km resolution. A bicubic interpolation method was used. _Bottom_: Zoom for an area at the North-East of the country (black square) which is close to the city of Beni and the Virunga national park. Due to the structure of the intrinsic CAR model (see Eq. \@ref(eq:icar)), spatial random effects are also estimated for cells without sampled points. This includes cells for which there was no forest cover in the period 2000--2010--2020, and also cells outside the country's borders.

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/rho} 

}

\caption{(ref:cap-rho)}(\#fig:rho)
\end{figure}

<!--------------------------------------------->
<!-- Spatial probability of deforestation -->
<!--------------------------------------------->

## Relative spatial probability of deforestation

(ref:cap-prob) **Predicted relative spatial probability of deforestation**. _Main figure_: Map of the spatial probability of deforestation computed for each forest pixel in 2020 for the Democratic Republic of the Congo. On the map, we clearly see the effect of the distance to nearest town, road, and river, and the effect of the distance to forest edge on the spatial probability of deforestation. Also, we clearly see the importance of the spatial random effects in structuring the spatial variability of the deforestation probability. For example, the area at the North of the zoom (black square) shows very high deforestation probabilities (in black). This area is politically unstable and is home to a large number of militias who survive at the expense of the forest. _Inset_: Zoom of the map for an area at the North-East of the country which is close to the city of Beni and the Virunga national park.

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/prob} 

}

\caption{(ref:cap-prob)}(\#fig:prob)
\end{figure}

<!--------------------------------------------->
<!-- Projected forest cover change -->
<!--------------------------------------------->

## Projected forest cover change

(ref:cap-proj) **Projected forest cover change**. _Main figures_: Maps of the projected forest cover change (left: 2020--2050, right: 2020--2100) for the Democratic Republic of the Congo (DRC). \textcolor{red}{red}: projected deforestation, \textcolor{darkgreen}{green}: remaining forest cover. Besides the loss of forest cover, we show a progressive fragmentation of the forest in the future, with an increasing number of forest patches of smaller size in DRC. _Insets_: Zoom of the map for an area at the North-East of the country (black square) which is close to the city of Beni and the Virunga national park. Interactive global maps of the projected forest cover change for years 2050 and 2100 are available at <https://forestatrisk.cirad.fr/tropics>.

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/fcc2050_2100} 

}

\caption{(ref:cap-proj)}(\#fig:proj)
\end{figure}

<!--------------------------------------------->
<!-- Percentage of forest cover loss -->
<!--------------------------------------------->

## Percentage of forest cover loss

(ref:cap-perc-loss) **Projected percentage of forest cover loss per continent**. Points represent the observed percentage of forest cover loss (in comparison with the year 2000) for years 2000 (0%), 2010, and 2020, for the three continents: America, Africa, and Asia. Lines represent the projected percentage of forest cover loss (in comparison with the year 2000) from year 2020 to 2400 per continent. For the deforestation projections, we assumed no diffusion of the deforestation between countries (see section \@ref(forecast) in Methods). When large countries with high annual deforested areas (Brazil for America, DRC for Africa, and Indonesia for Asia) have nor more forest (in 2242, 2160, and 2122, respectively, Table \@ref(tab:fcc-proj)), deforestation at the continent scale is rapidly decreasing. The horizontal black line indicates a loss of 75% of the forest cover in comparison with the year 2000. Under a "business-as-usual" scenario, this should happen in 2099, 2145, and 2201 for Asia, Africa, and America, respectively (see Table \@ref(tab:fcc-proj-reg)).

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/perc_loss_cont} 

}

\caption{(ref:cap-perc-loss)}(\#fig:perc-loss)
\end{figure}

<!--------------------------------------------->
<!-- Aboveground biomass map -->
<!--------------------------------------------->

## Aboveground biomass map

(ref:cap-agb) **Aboveground biomass map**. _Main figure_: Aboveground biomass (AGB in Mg/ha) map for the Democratic Republic of Congo (DRC) at 1 km resolution produced by @Avitabile2016. This map is a combination of two pantropical aboveground biomass maps by @Saatchi2011 and @Baccini2012. While the RMSE of the fused map by @Avitabile2016 is still substantial (87--98 Mg/ha), the fused map achieved a lower RMSE (a decrease of 5--74%) and bias (a decrease of 90--153%) than the two input maps for all continents. The fusion map is representative of the aboveground biomass for the years 2000--2010. _Inset_: Zoom of the map for an area at the North-East of the country which is close to the city of Beni and the Virunga national park.

\begin{figure}[H]

{\centering \includegraphics[width=\textwidth]{figures/AGB} 

}

\caption{(ref:cap-agb)}(\#fig:agb)
\end{figure}

\newpage

<!--chapter:end:02-Figures-S.Rmd-->

# Tables



<!--------------------------------------------------->
<!-- Categories of the forest cover annual product -->
<!--------------------------------------------------->

## Categories of the forest cover annual product

(ref:cap-cat-annual-product) **Categories of the forest cover annual product by @Vancutsem2020**. The forest cover annual product classifies Landsat image pixels in 16 categories for each year (on the 31$^{\text{st}}$ of December) between 1982 and 2019 and allows identifying moist tropical forest pixels at each date. Water data comes from @Pekel2016.\vspace{0.5cm}

\begin{table}[H]

\caption{(\#tab:cat-annual-product)(ref:cap-cat-annual-product)}
\centering
\begin{tabular}[t]{>{\raggedright\arraybackslash}p{2cm}>{\raggedleft\arraybackslash}p{10cm}}
\toprule
Class & Definition\\
\midrule
\cellcolor{gray!6}{1} & \cellcolor{gray!6}{Tropical moist forest (TMF including bamboo-dominated forest and mangroves)}\\
2 & TMF converted later in a tree plantation\\
\cellcolor{gray!6}{3} & \cellcolor{gray!6}{NEW degradation}\\
4 & Ongoing degradation (disturbances still detected)\\
\cellcolor{gray!6}{5} & \cellcolor{gray!6}{Degraded forest (former degradation, no disturbances detected anymore)}\\
6 & NEW deforestation (may follow degradation)\\
\cellcolor{gray!6}{7} & \cellcolor{gray!6}{Ongoing deforestation (disturbances still detected)}\\
8 & NEW Regrowth\\
\cellcolor{gray!6}{9} & \cellcolor{gray!6}{Regrowthing}\\
10 & Other land cover (not water)\\
\cellcolor{gray!6}{11} & \cellcolor{gray!6}{Permanent Water (Pekel et al. 2016)}\\
12 & Seasonal Water (Pekel et al. 2016)\\
\cellcolor{gray!6}{13} & \cellcolor{gray!6}{Init period without valid data - Init class = TMF}\\
14 & Init period with min 1 valid obs - Init class = TMF\\
\cellcolor{gray!6}{15} & \cellcolor{gray!6}{Nodata - Init class = other LC}\\
16 & Init period without valid data - Init class = Plantation\\
\bottomrule
\end{tabular}
\end{table}

<!--------------->
<!-- Variables -->
<!--------------->

## Variables

(ref:cap-variables) **Set of explicative variables used to model the spatial probability of deforestation**. A total of height variables were tested. They indicate topography, forest accessibility, forest landscape, deforestation history, and conservation status.\vspace{0.5cm}

\begin{table}[H]

\caption{(\#tab:variables)(ref:cap-variables)}
\centering
\fontsize{11}{13}\selectfont
\begin{tabular}[t]{>{\raggedright\arraybackslash}p{2.5cm}>{\raggedright\arraybackslash}p{2.5cm}>{\raggedright\arraybackslash}p{2.5cm}>{\raggedleft\arraybackslash}p{1cm}>{\raggedleft\arraybackslash}p{2cm}>{\raggedleft\arraybackslash}p{2cm}}
\toprule
Product & Source & Variable derived & Unit & Resolution (m) & Date\\
\midrule
\cellcolor{gray!6}{Forest maps (2000-2010-2020)} & \cellcolor{gray!6}{Vancutsem et al. 2020} & \cellcolor{gray!6}{distance to forest edge} & \cellcolor{gray!6}{m} & \cellcolor{gray!6}{30} & \cellcolor{gray!6}{--}\\
 &  & distance to past deforestation & m & 30 & --\\
\cellcolor{gray!6}{Digital Elevation Model} & \cellcolor{gray!6}{SRTM v4.1 CSI-CGIAR} & \cellcolor{gray!6}{elevation} & \cellcolor{gray!6}{m} & \cellcolor{gray!6}{90} & \cellcolor{gray!6}{--}\\
 &  & slope & degree & 90 & --\\
\cellcolor{gray!6}{Highways} & \cellcolor{gray!6}{OSM-Geofabrik} & \cellcolor{gray!6}{distance to road} & \cellcolor{gray!6}{m} & \cellcolor{gray!6}{150} & \cellcolor{gray!6}{Jan-Mar 2020}\\
Places &  & distance to town & m & 150 & Jan-Mar 2020\\
\cellcolor{gray!6}{Waterways} & \cellcolor{gray!6}{} & \cellcolor{gray!6}{distance to river} & \cellcolor{gray!6}{m} & \cellcolor{gray!6}{150} & \cellcolor{gray!6}{Jan-Mar 2020}\\
Protected areas & WDPA & presence of protected area & -- & 30 & Jan-Mar 2020\\
\bottomrule
\end{tabular}
\end{table}

<!----------------->
<!-- Sample size -->
<!----------------->

## Sample size

(ref:cap-samp-size) **Number of observations used for the spatial model of deforestation for each study area**. The table includes the number of non-deforested (nfor) and deforested (ndef) pixels per study area. These numbers include the forest pixels with full information regarding the explanatory variables. The corresponding number of hectares is also provided (nfHa and ndHa, respectively).\vspace{0.5cm}

\begingroup\fontsize{11}{13}\selectfont

\begin{longtable}[t]{llrrrr}
\caption{(\#tab:samp-size)(ref:cap-samp-size)}\\
\toprule
Country -- Study-area & Code & nfor & ndef & nfHa & ndHa\\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
Country -- Study-area & Code & nfor & ndef & nfHa & ndHa\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{America}}\\
\cellcolor{gray!6}{\hspace{1em}Antigua and B.} & \cellcolor{gray!6}{ATG} & \cellcolor{gray!6}{9,840} & \cellcolor{gray!6}{6,336} & \cellcolor{gray!6}{886} & \cellcolor{gray!6}{570}\\
\hspace{1em}Bahamas & BHS & 9,923 & 9,946 & 893 & 895\\
\cellcolor{gray!6}{\hspace{1em}Barbados} & \cellcolor{gray!6}{BRB} & \cellcolor{gray!6}{9,958} & \cellcolor{gray!6}{8,251} & \cellcolor{gray!6}{896} & \cellcolor{gray!6}{743}\\
\hspace{1em}Belize & BLZ & 9,995 & 9,999 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Bolivia} & \cellcolor{gray!6}{BOL} & \cellcolor{gray!6}{30,476} & \cellcolor{gray!6}{30,476} & \cellcolor{gray!6}{2,743} & \cellcolor{gray!6}{2,743}\\
\hspace{1em}Brazil – Acre & AC & 13,164 & 13,164 & 1,185 & 1,185\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Alagoas} & \cellcolor{gray!6}{AL} & \cellcolor{gray!6}{9,996} & \cellcolor{gray!6}{9,998} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Brazil – Amapa & AP & 11,466 & 11,469 & 1,032 & 1,032\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Amazonas} & \cellcolor{gray!6}{AM} & \cellcolor{gray!6}{50,000} & \cellcolor{gray!6}{50,000} & \cellcolor{gray!6}{4,500} & \cellcolor{gray!6}{4,500}\\
\hspace{1em}Brazil – Bahia & BA & 9,986 & 9,998 & 899 & 900\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Ceara} & \cellcolor{gray!6}{CE} & \cellcolor{gray!6}{9,996} & \cellcolor{gray!6}{9,999} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Brazil – Espirito Santo & ES & 9,989 & 9,999 & 899 & 900\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Goias} & \cellcolor{gray!6}{GO} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Brazil – Maranhao & MA & 9,987 & 9,997 & 899 & 900\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Mato Grosso} & \cellcolor{gray!6}{MT} & \cellcolor{gray!6}{31,678} & \cellcolor{gray!6}{31,678} & \cellcolor{gray!6}{2,851} & \cellcolor{gray!6}{2,851}\\
\hspace{1em}Brazil – Mato Grosso do Sul & MS & 10,000 & 10,000 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Minas Gerais} & \cellcolor{gray!6}{MG} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Brazil – Para & PA & 49,999 & 49,999 & 4,500 & 4,500\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Paraiba} & \cellcolor{gray!6}{PB} & \cellcolor{gray!6}{9,973} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{898} & \cellcolor{gray!6}{900}\\
\hspace{1em}Brazil – Parana & PR & 9,996 & 10,000 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Pernambouco} & \cellcolor{gray!6}{PE} & \cellcolor{gray!6}{9,964} & \cellcolor{gray!6}{9,999} & \cellcolor{gray!6}{897} & \cellcolor{gray!6}{900}\\
\hspace{1em}Brazil – Piaui & PI & 10,000 & 10,000 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio de Janeiro} & \cellcolor{gray!6}{RJ} & \cellcolor{gray!6}{9,994} & \cellcolor{gray!6}{9,991} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{899}\\
\hspace{1em}Brazil – Rio Grande do Norte & RN & 9,944 & 9,992 & 895 & 899\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio Grande do Sul} & \cellcolor{gray!6}{RS} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{9,999} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Brazil – Rondonia & RO & 12,964 & 12,964 & 1,167 & 1,167\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Roraima} & \cellcolor{gray!6}{RR} & \cellcolor{gray!6}{15,548} & \cellcolor{gray!6}{15,548} & \cellcolor{gray!6}{1,399} & \cellcolor{gray!6}{1,399}\\
\hspace{1em}Brazil – Santa Catarina & SC & 9,999 & 10,000 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Sao Paulo} & \cellcolor{gray!6}{SP} & \cellcolor{gray!6}{9,997} & \cellcolor{gray!6}{9,997} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Brazil – Sergipe & SE & 9,970 & 9,999 & 897 & 900\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Tocantins} & \cellcolor{gray!6}{TO} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Colombia & COL & 49,995 & 49,999 & 4,500 & 4,500\\
\cellcolor{gray!6}{\hspace{1em}Costa Rica} & \cellcolor{gray!6}{CRI} & \cellcolor{gray!6}{9,997} & \cellcolor{gray!6}{9,991} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{899}\\
\hspace{1em}Cuba & CUB & 9,955 & 9,958 & 896 & 896\\
\cellcolor{gray!6}{\hspace{1em}Dominica} & \cellcolor{gray!6}{DMA} & \cellcolor{gray!6}{9,991} & \cellcolor{gray!6}{9,878} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{889}\\
\hspace{1em}Dominican Rep. & DOM & 9,992 & 9,997 & 899 & 900\\
\cellcolor{gray!6}{\hspace{1em}Ecuador} & \cellcolor{gray!6}{ECU} & \cellcolor{gray!6}{14,395} & \cellcolor{gray!6}{14,395} & \cellcolor{gray!6}{1,296} & \cellcolor{gray!6}{1,296}\\
\hspace{1em}El Salvador & SLV & 9,965 & 9,982 & 897 & 898\\
\cellcolor{gray!6}{\hspace{1em}French Guiana} & \cellcolor{gray!6}{GUF} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Grenada & GRD & 9,974 & 9,928 & 898 & 894\\
\cellcolor{gray!6}{\hspace{1em}Guadeloupe} & \cellcolor{gray!6}{GLP} & \cellcolor{gray!6}{9,979} & \cellcolor{gray!6}{9,923} & \cellcolor{gray!6}{898} & \cellcolor{gray!6}{893}\\
\hspace{1em}Guatemala & GTM & 9,999 & 9,999 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Guyana} & \cellcolor{gray!6}{GUY} & \cellcolor{gray!6}{18,511} & \cellcolor{gray!6}{18,511} & \cellcolor{gray!6}{1,666} & \cellcolor{gray!6}{1,666}\\
\hspace{1em}Haiti & HTI & 9,952 & 9,970 & 896 & 897\\
\cellcolor{gray!6}{\hspace{1em}Honduras} & \cellcolor{gray!6}{HND} & \cellcolor{gray!6}{9,997} & \cellcolor{gray!6}{9,999} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Jamaica & JAM & 9,988 & 9,996 & 899 & 900\\
\cellcolor{gray!6}{\hspace{1em}Martinique} & \cellcolor{gray!6}{MTQ} & \cellcolor{gray!6}{9,972} & \cellcolor{gray!6}{9,973} & \cellcolor{gray!6}{897} & \cellcolor{gray!6}{898}\\
\hspace{1em}Mexico & MEX & 9,994 & 9,997 & 899 & 900\\
\cellcolor{gray!6}{\hspace{1em}Montserrat} & \cellcolor{gray!6}{MSR} & \cellcolor{gray!6}{9,982} & \cellcolor{gray!6}{1,258} & \cellcolor{gray!6}{898} & \cellcolor{gray!6}{113}\\
\hspace{1em}Nicaragua & NIC & 9,999 & 9,999 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Panama} & \cellcolor{gray!6}{PAN} & \cellcolor{gray!6}{9,996} & \cellcolor{gray!6}{9,994} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{899}\\
\hspace{1em}Paraguay & PRY & 10,000 & 10,000 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Peru} & \cellcolor{gray!6}{PER} & \cellcolor{gray!6}{50,000} & \cellcolor{gray!6}{50,000} & \cellcolor{gray!6}{4,500} & \cellcolor{gray!6}{4,500}\\
\hspace{1em}Puerto Rico & PRI & 9,993 & 9,984 & 899 & 899\\
\cellcolor{gray!6}{\hspace{1em}Saint Kitts and N.} & \cellcolor{gray!6}{KNA} & \cellcolor{gray!6}{9,988} & \cellcolor{gray!6}{4,614} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{415}\\
\hspace{1em}Saint Lucia & LCA & 9,994 & 9,967 & 899 & 897\\
\cellcolor{gray!6}{\hspace{1em}Saint Martin} & \cellcolor{gray!6}{MAF} & \cellcolor{gray!6}{3,156} & \cellcolor{gray!6}{3,808} & \cellcolor{gray!6}{284} & \cellcolor{gray!6}{343}\\
\hspace{1em}Saint Vincent & VCT & 9,981 & 9,792 & 898 & 881\\
\cellcolor{gray!6}{\hspace{1em}Sint Maarten} & \cellcolor{gray!6}{SXM} & \cellcolor{gray!6}{1,305} & \cellcolor{gray!6}{1,396} & \cellcolor{gray!6}{117} & \cellcolor{gray!6}{126}\\
\hspace{1em}Suriname & SUR & 13,699 & 13,698 & 1,233 & 1,233\\
\cellcolor{gray!6}{\hspace{1em}Trinidad and Tobago} & \cellcolor{gray!6}{TTO} & \cellcolor{gray!6}{9,993} & \cellcolor{gray!6}{9,984} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{899}\\
\hspace{1em}Venezuela & VEN & 41,916 & 41,911 & 3,772 & 3,772\\
\cellcolor{gray!6}{\hspace{1em}Virgin Isl. UK} & \cellcolor{gray!6}{VGB} & \cellcolor{gray!6}{9,869} & \cellcolor{gray!6}{9,851} & \cellcolor{gray!6}{888} & \cellcolor{gray!6}{887}\\
\hspace{1em}Virgin Isl. US & VIR & 9,945 & 9,874 & 895 & 889\\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{Africa}}\\
\cellcolor{gray!6}{\hspace{1em}Angola} & \cellcolor{gray!6}{AGO} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Benin & BEN & 9,981 & 9,995 & 898 & 900\\
\cellcolor{gray!6}{\hspace{1em}Burundi} & \cellcolor{gray!6}{BDI} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{9,999} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Cameroon & CMR & 22,983 & 22,987 & 2,068 & 2,069\\
\cellcolor{gray!6}{\hspace{1em}CAR} & \cellcolor{gray!6}{CAF} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Comoros & COM & 9,990 & 9,984 & 899 & 899\\
\cellcolor{gray!6}{\hspace{1em}Congo} & \cellcolor{gray!6}{COG} & \cellcolor{gray!6}{23,412} & \cellcolor{gray!6}{23,411} & \cellcolor{gray!6}{2,107} & \cellcolor{gray!6}{2,107}\\
\hspace{1em}DRC & COD & 50,000 & 50,000 & 4,500 & 4,500\\
\cellcolor{gray!6}{\hspace{1em}Eq. Guinea} & \cellcolor{gray!6}{GNQ} & \cellcolor{gray!6}{9,997} & \cellcolor{gray!6}{9,988} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{899}\\
\hspace{1em}Ethiopia & ETH & 10,000 & 10,000 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Gabon} & \cellcolor{gray!6}{GAB} & \cellcolor{gray!6}{23,986} & \cellcolor{gray!6}{23,966} & \cellcolor{gray!6}{2,159} & \cellcolor{gray!6}{2,157}\\
\hspace{1em}Gambia & GMB & 9,988 & 9,999 & 899 & 900\\
\cellcolor{gray!6}{\hspace{1em}Ghana} & \cellcolor{gray!6}{GHA} & \cellcolor{gray!6}{9,999} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Guinea & GIN & 9,978 & 9,998 & 898 & 900\\
\cellcolor{gray!6}{\hspace{1em}Guinea Bissau} & \cellcolor{gray!6}{GNB} & \cellcolor{gray!6}{9,882} & \cellcolor{gray!6}{9,977} & \cellcolor{gray!6}{889} & \cellcolor{gray!6}{898}\\
\hspace{1em}Ivory Coast & CIV & 10,000 & 10,000 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Kenya} & \cellcolor{gray!6}{KEN} & \cellcolor{gray!6}{9,986} & \cellcolor{gray!6}{9,998} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{900}\\
\hspace{1em}Liberia & LBR & 9,999 & 9,997 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Madagascar} & \cellcolor{gray!6}{MDG} & \cellcolor{gray!6}{9,994} & \cellcolor{gray!6}{9,998} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{900}\\
\hspace{1em}Malawi & MWI & 10,000 & 10,000 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Mauritius} & \cellcolor{gray!6}{MUS} & \cellcolor{gray!6}{9,968} & \cellcolor{gray!6}{9,969} & \cellcolor{gray!6}{897} & \cellcolor{gray!6}{897}\\
\hspace{1em}Mayotte & MYT & 9,963 & 9,983 & 897 & 898\\
\cellcolor{gray!6}{\hspace{1em}Nigeria} & \cellcolor{gray!6}{NGA} & \cellcolor{gray!6}{9,984} & \cellcolor{gray!6}{9,995} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{900}\\
\hspace{1em}Reunion & REU & 9,996 & 9,995 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Rwanda} & \cellcolor{gray!6}{RWA} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Senegal & SEN & 9,892 & 9,976 & 890 & 898\\
\cellcolor{gray!6}{\hspace{1em}Sierra Leone} & \cellcolor{gray!6}{SLE} & \cellcolor{gray!6}{9,997} & \cellcolor{gray!6}{9,999} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}South Sudan & SSD & 5,737 & 7,613 & 516 & 685\\
\cellcolor{gray!6}{\hspace{1em}Tanzania} & \cellcolor{gray!6}{TZA} & \cellcolor{gray!6}{9,983} & \cellcolor{gray!6}{9,970} & \cellcolor{gray!6}{898} & \cellcolor{gray!6}{897}\\
\hspace{1em}Togo & TGO & 10,000 & 9,999 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}Uganda} & \cellcolor{gray!6}{UGA} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Zambia & ZMB & 10,000 & 10,000 & 900 & 900\\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{Asia}}\\
\cellcolor{gray!6}{\hspace{1em}Australia – Queensland} & \cellcolor{gray!6}{QLD} & \cellcolor{gray!6}{9,979} & \cellcolor{gray!6}{9,984} & \cellcolor{gray!6}{898} & \cellcolor{gray!6}{899}\\
\hspace{1em}Bangladesh & BGD & 9,962 & 9,993 & 897 & 899\\
\cellcolor{gray!6}{\hspace{1em}Bhutan} & \cellcolor{gray!6}{BTN} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Brunei & BRN & 9,989 & 9,997 & 899 & 900\\
\cellcolor{gray!6}{\hspace{1em}Cambodia} & \cellcolor{gray!6}{KHM} & \cellcolor{gray!6}{9,993} & \cellcolor{gray!6}{9,999} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{900}\\
\hspace{1em}Fiji & FJI & 9,983 & 9,957 & 898 & 896\\
\cellcolor{gray!6}{\hspace{1em}India – Andaman and N.} & \cellcolor{gray!6}{AN} & \cellcolor{gray!6}{9,950} & \cellcolor{gray!6}{9,908} & \cellcolor{gray!6}{896} & \cellcolor{gray!6}{892}\\
\hspace{1em}India – North-East & NE & 9,997 & 10,000 & 900 & 900\\
\cellcolor{gray!6}{\hspace{1em}India – West. Ghats} & \cellcolor{gray!6}{WG} & \cellcolor{gray!6}{9,989} & \cellcolor{gray!6}{9,990} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{899}\\
\hspace{1em}Indonesia & IDN & 49,967 & 49,984 & 4,497 & 4,499\\
\cellcolor{gray!6}{\hspace{1em}Laos} & \cellcolor{gray!6}{LAO} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\hspace{1em}Malaysia & MYS & 20,199 & 20,209 & 1,818 & 1,819\\
\cellcolor{gray!6}{\hspace{1em}Myanmar} & \cellcolor{gray!6}{MMR} & \cellcolor{gray!6}{16,076} & \cellcolor{gray!6}{16,078} & \cellcolor{gray!6}{1,447} & \cellcolor{gray!6}{1,447}\\
\hspace{1em}New Caledonia & NCL & 9,961 & 9,932 & 896 & 894\\
\cellcolor{gray!6}{\hspace{1em}Papua New Guinea} & \cellcolor{gray!6}{PNG} & \cellcolor{gray!6}{39,738} & \cellcolor{gray!6}{39,691} & \cellcolor{gray!6}{3,576} & \cellcolor{gray!6}{3,572}\\
\hspace{1em}Philippines & PHL & 13,251 & 13,251 & 1,193 & 1,193\\
\cellcolor{gray!6}{\hspace{1em}Singapore} & \cellcolor{gray!6}{SGP} & \cellcolor{gray!6}{9,904} & \cellcolor{gray!6}{9,961} & \cellcolor{gray!6}{891} & \cellcolor{gray!6}{896}\\
\hspace{1em}Solomon Isl. & SLB & 9,939 & 9,853 & 895 & 887\\
\cellcolor{gray!6}{\hspace{1em}Sri Lanka} & \cellcolor{gray!6}{LKA} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{9,993} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{899}\\
\hspace{1em}Thailand & THA & 9,990 & 9,997 & 899 & 900\\
\cellcolor{gray!6}{\hspace{1em}Timor-Leste} & \cellcolor{gray!6}{TLS} & \cellcolor{gray!6}{9,994} & \cellcolor{gray!6}{9,966} & \cellcolor{gray!6}{899} & \cellcolor{gray!6}{897}\\
\hspace{1em}Vanuatu & VUT & 9,977 & 9,925 & 898 & 893\\
\cellcolor{gray!6}{\hspace{1em}Vietnam} & \cellcolor{gray!6}{VNM} & \cellcolor{gray!6}{9,998} & \cellcolor{gray!6}{10,000} & \cellcolor{gray!6}{900} & \cellcolor{gray!6}{900}\\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{All continents}}\\
\hspace{1em}TOTAL &  & 1,601,805 & 1,584,888 & 144,163 & 142,647\\*
\end{longtable}
\endgroup{}

<!------------------------------------------------>
<!-- Mathematical formulas for accuracy indices -->
<!------------------------------------------------>

## Mathematical formulas for accuracy indices

(ref:cap-confusion-matrix) **Confusion matrix used to compute accuracy indices**. A confusion matrix can be computed to compare model predictions with observations.\vspace{0.5cm}

\begin{table}[H]

\caption{(\#tab:confusion-matrix)(ref:cap-confusion-matrix)}
\centering
\begin{tabular}[t]{lllll}
\toprule
 &  & Observations &  & Total\\
\midrule
\cellcolor{gray!6}{} & \cellcolor{gray!6}{} & \cellcolor{gray!6}{0 (non-deforested)} & \cellcolor{gray!6}{1 (deforested)} & \cellcolor{gray!6}{}\\
Predictions & 0 & $n_{00}$ & $n_{01}$ & $n_{0+}$\\
\cellcolor{gray!6}{} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{$n_{10}$} & \cellcolor{gray!6}{$n_{11}$} & \cellcolor{gray!6}{$n_{1+}$}\\
Total &  & $n_{+0}$ & $n_{+1}$ & $n$\\
\bottomrule
\end{tabular}
\end{table}


(ref:cap-accuracy-indices) **Formulas used to compute accuracy indices**. Several accuracy indices can be computed from the confusion matrix to estimate and compare models' predictive skill. We followed the definitions of @Pontius2008 for the FOM and @Liu2011 for the other indices. Note that the AUC relies on the predicted probabilities for observations 0 (non-deforested) and 1 (deforested), not on the confusion matrix.\vspace{0.5cm}

\begin{table}[H]

\caption{(\#tab:accuracy-indices)(ref:cap-accuracy-indices)}
\centering
\begin{tabular}[t]{ll}
\toprule
Index & Formula\\
\midrule
\cellcolor{gray!6}{Overall Accuracy} & \cellcolor{gray!6}{$\text{OA} = (n_{11}+n_{00})/n$}\\
Expected Accuracy & $\text{EA} = (n_{1+} n_{+1}+n_{0+} n_{+0})/n^2$\\
\cellcolor{gray!6}{Figure Of Merit} & \cellcolor{gray!6}{$\text{FOM} = n_{11}/(n_{11}+n_{10}+n_{01})$}\\
Sensitivity & $\text{Sen} = n_{11}/(n_{11}+n_{01})$\\
\cellcolor{gray!6}{Specificity} & \cellcolor{gray!6}{$\text{Spe} = n_{00}/(n_{00}+n_{10})$}\\
True Skill Statistics & $\text{TSS} = \text{Sen}+\text{Spe}-1$\\
\cellcolor{gray!6}{Cohen's Kappa} & \cellcolor{gray!6}{$\text{K} = (\text{OA}-\text{EA})/(1-\text{EA})$}\\
Area Under ROC Curve & $\text{AUC} = 1/(n_{+1} n_{+0}) \sum_{i=1}^{n_{+0}} \sum_{j=1}^{n_{+1}} \phi(\delta_i,\theta_j)$\\
\cellcolor{gray!6}{} & \cellcolor{gray!6}{where $\phi(\delta_i,\theta_j)$ equals 1 if $\theta_j>\delta_i$, 1/2 if $\theta_j=\delta_i$, and 0 otherwise}\\
 & $\delta_i$ and $\theta_j$ are the predicted probabilities for $Y_i=0$ and $Y_j=1$\\
\bottomrule
\end{tabular}
\end{table}

<!---------------------->
<!-- Accuracy indices -->
<!---------------------->

## Accuracy indices

(ref:cap-accuracy) **Accuracy indices' weighted averages for the three statistical models**. Accuracy indices were averaged across study areas using forest cover areas in 2010 as weights. D: percentage of deviance explained, AUC Area Under ROC Curve, OA: overall accuracy, FOM: Figure Of Merit, TSS: True Skill Statistics. Averaged accuracy indices were computed for the three statistical models: "glm", "icar", and "rf" model.\vspace{0.5cm}

\begin{table}[H]

\caption{(\#tab:accuracy)(ref:cap-accuracy)}
\centering
\begin{tabular}[t]{lrrrrr}
\toprule
Model & D & AUC & OA & FOM & TSS\\
\midrule
glm & 28.3 & 83.1 & 75.7 & 61.2 & 51.7\\
\cellcolor{gray!6}{icar} & \cellcolor{gray!6}{38.6} & \cellcolor{gray!6}{86.8} & \cellcolor{gray!6}{79.1} & \cellcolor{gray!6}{65.5} & \cellcolor{gray!6}{58.1}\\
rf & 75.9 & 84.0 & 76.7 & 62.4 & 53.2\\
\bottomrule
\end{tabular}
\end{table}


(ref:cap-accuracy-cont) **Accuracy indices' weighted averages for the three statistical models by continent**. Accuracy indices were averaged across study areas for each continent using forest cover areas in 2010 as weights. D: percentage of deviance explained, AUC Area Under ROC Curve, OA: overall accuracy, FOM: Figure Of Merit, TSS: True Skill Statistics. Averaged accuracy indices were computed for the three statistical models: "glm", "icar", and "rf" model.\vspace{0.5cm}

\begin{table}[H]

\caption{(\#tab:accuracy-cont)(ref:cap-accuracy-cont)}
\centering
\begin{tabular}[t]{>{}lrrrrrl}
\toprule
Continent & Model & D & AUC & OA & FOM & TSS\\
\midrule
 & glm & 30.6 & 83.9 & 76.6 & 62.5 & 53.6\\

\cellcolor{gray!6}{} & \cellcolor{gray!6}{icar} & \cellcolor{gray!6}{40.2} & \cellcolor{gray!6}{87.1} & \cellcolor{gray!6}{79.8} & \cellcolor{gray!6}{66.6} & \cellcolor{gray!6}{59.5}\\

\multirow{-3}{*}{\raggedright\arraybackslash \textbf{America}} & rf & 81.9 & 85.4 & 78.1 & 64.6 & 56.2\\
\cmidrule{1-7}
 & glm & 25.9 & 81.4 & 73.9 & 58.9 & 48.1\\

\cellcolor{gray!6}{} & \cellcolor{gray!6}{icar} & \cellcolor{gray!6}{35.8} & \cellcolor{gray!6}{86.0} & \cellcolor{gray!6}{77.8} & \cellcolor{gray!6}{63.7} & \cellcolor{gray!6}{55.5}\\

\multirow{-3}{*}{\raggedright\arraybackslash \textbf{Africa}} & rf & 65.4 & 81.2 & 74.3 & 58.9 & 48.0\\
\cmidrule{1-7}
 & glm & 27.6 & 84.3 & 77.2 & 63.1 & 54.4\\

\cellcolor{gray!6}{} & \cellcolor{gray!6}{icar} & \cellcolor{gray!6}{40.5} & \cellcolor{gray!6}{88.0} & \cellcolor{gray!6}{80.2} & \cellcolor{gray!6}{66.6} & \cellcolor{gray!6}{59.8}\\

\multirow{-3}{*}{\raggedright\arraybackslash \textbf{Asia}} & rf & 83.4 & 86.3 & 78.5 & 64.6 & 57.0\\
\bottomrule
\end{tabular}
\end{table}

<!------------------------------------>
<!-- Historical forest cover change -->
<!------------------------------------>

## Historical forest cover change

(ref:cap-fcc-hist) **Historical forest cover change for each study-area**. Forest cover areas are given in thousand hectares (Kha) for the years 2000, 2010 and 2020 ("fc2000", "fc2010", and "fc2020", respectively). The mean annual deforested area $d$ for the period 2010--2020 is given in hectare per year (ha/yr). For comparing the deforestation intensity between study-areas, the corresponding mean annual deforestation rate $p$ on the period 2010--2020 is also provided in percent per year (%/yr), with one decimal precision.\vspace{0.5cm}

\begingroup\fontsize{11}{13}\selectfont

\begin{longtable}[t]{lrrrrr}
\caption{(\#tab:fcc-hist)(ref:cap-fcc-hist)}\\
\toprule
Country -- Study-area & fc2000 & fc2010 & fc2020 & $d$ (ha/yr) & $p$ (\%/yr)\\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
Country -- Study-area & fc2000 & fc2010 & fc2020 & $d$ (ha/yr) & $p$ (\%/yr)\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{America}}\\
\cellcolor{gray!6}{\hspace{1em}Antigua and B.} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{58} & \cellcolor{gray!6}{1.6}\\
\hspace{1em}Bahamas & 195 & 148 & 125 & 2,256 & 1.6\\
\cellcolor{gray!6}{\hspace{1em}Barbados} & \cellcolor{gray!6}{5} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{1.9}\\
\hspace{1em}Belize & 1,564 & 1,462 & 1,318 & 14,410 & 1.0\\
\cellcolor{gray!6}{\hspace{1em}Bolivia} & \cellcolor{gray!6}{34,936} & \cellcolor{gray!6}{32,612} & \cellcolor{gray!6}{30,476} & \cellcolor{gray!6}{213,540} & \cellcolor{gray!6}{0.7}\\
\hspace{1em}Brazil – Acre & 14,278 & 13,646 & 13,164 & 48,170 & 0.4\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Alagoas} & \cellcolor{gray!6}{119} & \cellcolor{gray!6}{103} & \cellcolor{gray!6}{91} & \cellcolor{gray!6}{1,169} & \cellcolor{gray!6}{1.2}\\
\hspace{1em}Brazil – Amapa & 11,765 & 11,602 & 11,469 & 13,244 & 0.1\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Amazonas} & \cellcolor{gray!6}{149,433} & \cellcolor{gray!6}{148,106} & \cellcolor{gray!6}{146,494} & \cellcolor{gray!6}{161,139} & \cellcolor{gray!6}{0.1}\\
\hspace{1em}Brazil – Bahia & 2,779 & 2,272 & 2,042 & 22,928 & 1.1\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Ceara} & \cellcolor{gray!6}{62} & \cellcolor{gray!6}{52} & \cellcolor{gray!6}{41} & \cellcolor{gray!6}{1,075} & \cellcolor{gray!6}{2.3}\\
\hspace{1em}Brazil – Espirito Santo & 581 & 484 & 435 & 4,831 & 1.0\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Goias} & \cellcolor{gray!6}{722} & \cellcolor{gray!6}{532} & \cellcolor{gray!6}{377} & \cellcolor{gray!6}{15,560} & \cellcolor{gray!6}{3.4}\\
\hspace{1em}Brazil – Maranhao & 5,846 & 4,156 & 3,281 & 87,533 & 2.3\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Mato Grosso} & \cellcolor{gray!6}{42,257} & \cellcolor{gray!6}{34,846} & \cellcolor{gray!6}{31,678} & \cellcolor{gray!6}{316,838} & \cellcolor{gray!6}{0.9}\\
\hspace{1em}Brazil – Mato Grosso do Sul & 1,008 & 852 & 750 & 10,202 & 1.3\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Minas Gerais} & \cellcolor{gray!6}{2,143} & \cellcolor{gray!6}{1,477} & \cellcolor{gray!6}{1,055} & \cellcolor{gray!6}{42,183} & \cellcolor{gray!6}{3.3}\\
\hspace{1em}Brazil – Para & 101,374 & 93,017 & 88,500 & 451,686 & 0.5\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Paraiba} & \cellcolor{gray!6}{48} & \cellcolor{gray!6}{42} & \cellcolor{gray!6}{39} & \cellcolor{gray!6}{342} & \cellcolor{gray!6}{0.8}\\
\hspace{1em}Brazil – Parana & 3,980 & 3,296 & 2,995 & 30,140 & 1.0\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Pernambouco} & \cellcolor{gray!6}{147} & \cellcolor{gray!6}{124} & \cellcolor{gray!6}{112} & \cellcolor{gray!6}{1,191} & \cellcolor{gray!6}{1.0}\\
\hspace{1em}Brazil – Piaui & 112 & 81 & 55 & 2,634 & 3.9\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio de Janeiro} & \cellcolor{gray!6}{981} & \cellcolor{gray!6}{875} & \cellcolor{gray!6}{809} & \cellcolor{gray!6}{6,630} & \cellcolor{gray!6}{0.8}\\
\hspace{1em}Brazil – Rio Grande do Norte & 33 & 26 & 23 & 339 & 1.4\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio Grande do Sul} & \cellcolor{gray!6}{3,398} & \cellcolor{gray!6}{2,926} & \cellcolor{gray!6}{2,709} & \cellcolor{gray!6}{21,681} & \cellcolor{gray!6}{0.8}\\
\hspace{1em}Brazil – Rondonia & 17,155 & 14,343 & 12,964 & 137,959 & 1.0\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Roraima} & \cellcolor{gray!6}{16,775} & \cellcolor{gray!6}{16,277} & \cellcolor{gray!6}{15,548} & \cellcolor{gray!6}{72,854} & \cellcolor{gray!6}{0.5}\\
\hspace{1em}Brazil – Santa Catarina & 3,785 & 3,171 & 2,956 & 21,442 & 0.7\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Sao Paulo} & \cellcolor{gray!6}{3,663} & \cellcolor{gray!6}{3,336} & \cellcolor{gray!6}{3,158} & \cellcolor{gray!6}{17,783} & \cellcolor{gray!6}{0.5}\\
\hspace{1em}Brazil – Sergipe & 79 & 65 & 57 & 774 & 1.3\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Tocantins} & \cellcolor{gray!6}{1,820} & \cellcolor{gray!6}{1,407} & \cellcolor{gray!6}{958} & \cellcolor{gray!6}{44,871} & \cellcolor{gray!6}{3.8}\\
\hspace{1em}Colombia & 70,622 & 67,322 & 64,483 & 283,945 & 0.4\\
\cellcolor{gray!6}{\hspace{1em}Costa Rica} & \cellcolor{gray!6}{2,507} & \cellcolor{gray!6}{2,359} & \cellcolor{gray!6}{2,181} & \cellcolor{gray!6}{17,772} & \cellcolor{gray!6}{0.8}\\
\hspace{1em}Cuba & 1,807 & 1,538 & 1,387 & 15,088 & 1.0\\
\cellcolor{gray!6}{\hspace{1em}Dominica} & \cellcolor{gray!6}{76} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{73} & \cellcolor{gray!6}{274} & \cellcolor{gray!6}{0.4}\\
\hspace{1em}Dominican Rep. & 1,415 & 1,137 & 968 & 16,879 & 1.6\\
\cellcolor{gray!6}{\hspace{1em}Ecuador} & \cellcolor{gray!6}{15,487} & \cellcolor{gray!6}{14,968} & \cellcolor{gray!6}{14,395} & \cellcolor{gray!6}{57,270} & \cellcolor{gray!6}{0.4}\\
\hspace{1em}El Salvador & 140 & 117 & 102 & 1,499 & 1.4\\
\cellcolor{gray!6}{\hspace{1em}French Guiana} & \cellcolor{gray!6}{8,160} & \cellcolor{gray!6}{8,132} & \cellcolor{gray!6}{8,105} & \cellcolor{gray!6}{2,708} & \cellcolor{gray!6}{0.0}\\
\hspace{1em}Grenada & 27 & 23 & 20 & 319 & 1.5\\
\cellcolor{gray!6}{\hspace{1em}Guadeloupe} & \cellcolor{gray!6}{87} & \cellcolor{gray!6}{84} & \cellcolor{gray!6}{79} & \cellcolor{gray!6}{498} & \cellcolor{gray!6}{0.6}\\
\hspace{1em}Guatemala & 3,801 & 2,959 & 2,476 & 48,333 & 1.8\\
\cellcolor{gray!6}{\hspace{1em}Guyana} & \cellcolor{gray!6}{18,782} & \cellcolor{gray!6}{18,655} & \cellcolor{gray!6}{18,511} & \cellcolor{gray!6}{14,360} & \cellcolor{gray!6}{0.1}\\
\hspace{1em}Haiti & 286 & 187 & 136 & 5,080 & 3.1\\
\cellcolor{gray!6}{\hspace{1em}Honduras} & \cellcolor{gray!6}{3,682} & \cellcolor{gray!6}{3,236} & \cellcolor{gray!6}{2,757} & \cellcolor{gray!6}{47,884} & \cellcolor{gray!6}{1.6}\\
\hspace{1em}Jamaica & 534 & 473 & 438 & 3,487 & 0.8\\
\cellcolor{gray!6}{\hspace{1em}Martinique} & \cellcolor{gray!6}{78} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{70} & \cellcolor{gray!6}{492} & \cellcolor{gray!6}{0.7}\\
\hspace{1em}Mexico & 10,399 & 8,340 & 6,957 & 138,359 & 1.8\\
\cellcolor{gray!6}{\hspace{1em}Montserrat} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{12} & \cellcolor{gray!6}{0.3}\\
\hspace{1em}Nicaragua & 5,244 & 4,549 & 3,562 & 98,723 & 2.4\\
\cellcolor{gray!6}{\hspace{1em}Panama} & \cellcolor{gray!6}{4,545} & \cellcolor{gray!6}{4,311} & \cellcolor{gray!6}{4,090} & \cellcolor{gray!6}{22,120} & \cellcolor{gray!6}{0.5}\\
\hspace{1em}Paraguay & 2,864 & 1,769 & 1,302 & 46,731 & 3.0\\
\cellcolor{gray!6}{\hspace{1em}Peru} & \cellcolor{gray!6}{74,971} & \cellcolor{gray!6}{73,476} & \cellcolor{gray!6}{72,136} & \cellcolor{gray!6}{133,992} & \cellcolor{gray!6}{0.2}\\
\hspace{1em}Puerto Rico & 472 & 406 & 354 & 5,193 & 1.4\\
\cellcolor{gray!6}{\hspace{1em}Saint Kitts and N.} & \cellcolor{gray!6}{10} & \cellcolor{gray!6}{10} & \cellcolor{gray!6}{9} & \cellcolor{gray!6}{42} & \cellcolor{gray!6}{0.4}\\
\hspace{1em}Saint Lucia & 50 & 50 & 47 & 310 & 0.6\\
\cellcolor{gray!6}{\hspace{1em}Saint Martin} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{34} & \cellcolor{gray!6}{7.6}\\
\hspace{1em}Saint Vincent & 30 & 30 & 29 & 114 & 0.4\\
\cellcolor{gray!6}{\hspace{1em}Sint Maarten} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{13} & \cellcolor{gray!6}{7.0}\\
\hspace{1em}Suriname & 13,902 & 13,814 & 13,699 & 11,480 & 0.1\\
\cellcolor{gray!6}{\hspace{1em}Trinidad and Tobago} & \cellcolor{gray!6}{361} & \cellcolor{gray!6}{345} & \cellcolor{gray!6}{317} & \cellcolor{gray!6}{2,818} & \cellcolor{gray!6}{0.8}\\
\hspace{1em}Venezuela & 45,344 & 43,476 & 41,917 & 155,907 & 0.4\\
\cellcolor{gray!6}{\hspace{1em}Virgin Isl. UK} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{2} & \cellcolor{gray!6}{114} & \cellcolor{gray!6}{4.3}\\
\hspace{1em}Virgin Isl. US & 10 & 9 & 7 & 194 & 2.4\\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{Africa}}\\
\cellcolor{gray!6}{\hspace{1em}Angola} & \cellcolor{gray!6}{7,322} & \cellcolor{gray!6}{6,236} & \cellcolor{gray!6}{5,365} & \cellcolor{gray!6}{87,108} & \cellcolor{gray!6}{1.5}\\
\hspace{1em}Benin & 76 & 48 & 33 & 1,508 & 3.7\\
\cellcolor{gray!6}{\hspace{1em}Burundi} & \cellcolor{gray!6}{107} & \cellcolor{gray!6}{65} & \cellcolor{gray!6}{55} & \cellcolor{gray!6}{1,014} & \cellcolor{gray!6}{1.7}\\
\hspace{1em}Cameroon & 24,001 & 23,688 & 22,989 & 69,890 & 0.3\\
\cellcolor{gray!6}{\hspace{1em}CAR} & \cellcolor{gray!6}{9,885} & \cellcolor{gray!6}{9,416} & \cellcolor{gray!6}{8,885} & \cellcolor{gray!6}{53,044} & \cellcolor{gray!6}{0.6}\\
\hspace{1em}Comoros & 90 & 90 & 86 & 490 & 0.6\\
\cellcolor{gray!6}{\hspace{1em}Congo} & \cellcolor{gray!6}{24,126} & \cellcolor{gray!6}{24,006} & \cellcolor{gray!6}{23,412} & \cellcolor{gray!6}{59,470} & \cellcolor{gray!6}{0.3}\\
\hspace{1em}DRC & 131,998 & 126,505 & 118,112 & 839,343 & 0.7\\
\cellcolor{gray!6}{\hspace{1em}Eq. Guinea} & \cellcolor{gray!6}{2,650} & \cellcolor{gray!6}{2,645} & \cellcolor{gray!6}{2,613} & \cellcolor{gray!6}{3,280} & \cellcolor{gray!6}{0.1}\\
\hspace{1em}Ethiopia & 3,870 & 2,907 & 2,168 & 73,914 & 2.9\\
\cellcolor{gray!6}{\hspace{1em}Gabon} & \cellcolor{gray!6}{24,143} & \cellcolor{gray!6}{24,124} & \cellcolor{gray!6}{23,989} & \cellcolor{gray!6}{13,422} & \cellcolor{gray!6}{0.1}\\
\hspace{1em}Gambia & 55 & 48 & 44 & 445 & 1.0\\
\cellcolor{gray!6}{\hspace{1em}Ghana} & \cellcolor{gray!6}{4,987} & \cellcolor{gray!6}{4,517} & \cellcolor{gray!6}{3,363} & \cellcolor{gray!6}{115,463} & \cellcolor{gray!6}{2.9}\\
\hspace{1em}Guinea & 1,961 & 1,269 & 861 & 40,811 & 3.8\\
\cellcolor{gray!6}{\hspace{1em}Guinea Bissau} & \cellcolor{gray!6}{427} & \cellcolor{gray!6}{349} & \cellcolor{gray!6}{303} & \cellcolor{gray!6}{4,658} & \cellcolor{gray!6}{1.4}\\
\hspace{1em}Ivory Coast & 7,800 & 6,439 & 3,914 & 252,474 & 4.9\\
\cellcolor{gray!6}{\hspace{1em}Kenya} & \cellcolor{gray!6}{1,207} & \cellcolor{gray!6}{902} & \cellcolor{gray!6}{750} & \cellcolor{gray!6}{15,265} & \cellcolor{gray!6}{1.8}\\
\hspace{1em}Liberia & 9,023 & 8,818 & 7,902 & 91,554 & 1.1\\
\cellcolor{gray!6}{\hspace{1em}Madagascar} & \cellcolor{gray!6}{7,956} & \cellcolor{gray!6}{6,242} & \cellcolor{gray!6}{5,002} & \cellcolor{gray!6}{123,984} & \cellcolor{gray!6}{2.2}\\
\hspace{1em}Malawi & 120 & 74 & 40 & 3,397 & 6.0\\
\cellcolor{gray!6}{\hspace{1em}Mauritius} & \cellcolor{gray!6}{56} & \cellcolor{gray!6}{54} & \cellcolor{gray!6}{49} & \cellcolor{gray!6}{440} & \cellcolor{gray!6}{0.8}\\
\hspace{1em}Mayotte & 19 & 18 & 13 & 528 & 3.4\\
\cellcolor{gray!6}{\hspace{1em}Nigeria} & \cellcolor{gray!6}{7,862} & \cellcolor{gray!6}{7,323} & \cellcolor{gray!6}{6,327} & \cellcolor{gray!6}{99,580} & \cellcolor{gray!6}{1.5}\\
\hspace{1em}Reunion & 164 & 164 & 155 & 874 & 0.5\\
\cellcolor{gray!6}{\hspace{1em}Rwanda} & \cellcolor{gray!6}{289} & \cellcolor{gray!6}{199} & \cellcolor{gray!6}{160} & \cellcolor{gray!6}{3,932} & \cellcolor{gray!6}{2.2}\\
\hspace{1em}Senegal & 147 & 139 & 130 & 828 & 0.6\\
\cellcolor{gray!6}{\hspace{1em}Sierra Leone} & \cellcolor{gray!6}{3,535} & \cellcolor{gray!6}{2,413} & \cellcolor{gray!6}{1,408} & \cellcolor{gray!6}{100,542} & \cellcolor{gray!6}{5.2}\\
\hspace{1em}South Sudan & 271 & 209 & 165 & 4,340 & 2.3\\
\cellcolor{gray!6}{\hspace{1em}Tanzania} & \cellcolor{gray!6}{1,464} & \cellcolor{gray!6}{1,225} & \cellcolor{gray!6}{1,105} & \cellcolor{gray!6}{11,984} & \cellcolor{gray!6}{1.0}\\
\hspace{1em}Togo & 163 & 106 & 66 & 3,976 & 4.6\\
\cellcolor{gray!6}{\hspace{1em}Uganda} & \cellcolor{gray!6}{1,887} & \cellcolor{gray!6}{1,106} & \cellcolor{gray!6}{745} & \cellcolor{gray!6}{36,029} & \cellcolor{gray!6}{3.9}\\
\hspace{1em}Zambia & 191 & 121 & 83 & 3,785 & 3.7\\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{Asia}}\\
\cellcolor{gray!6}{\hspace{1em}Australia – Queensland} & \cellcolor{gray!6}{2,268} & \cellcolor{gray!6}{2,069} & \cellcolor{gray!6}{1,938} & \cellcolor{gray!6}{13,170} & \cellcolor{gray!6}{0.7}\\
\hspace{1em}Bangladesh & 1,146 & 973 & 916 & 5,652 & 0.6\\
\cellcolor{gray!6}{\hspace{1em}Bhutan} & \cellcolor{gray!6}{2,554} & \cellcolor{gray!6}{2,386} & \cellcolor{gray!6}{2,270} & \cellcolor{gray!6}{11,535} & \cellcolor{gray!6}{0.5}\\
\hspace{1em}Brunei & 518 & 505 & 498 & 784 & 0.2\\
\cellcolor{gray!6}{\hspace{1em}Cambodia} & \cellcolor{gray!6}{5,075} & \cellcolor{gray!6}{4,079} & \cellcolor{gray!6}{2,931} & \cellcolor{gray!6}{114,795} & \cellcolor{gray!6}{3.3}\\
\hspace{1em}Fiji & 1,108 & 1,062 & 1,018 & 4,378 & 0.4\\
\cellcolor{gray!6}{\hspace{1em}India – Andaman and N.} & \cellcolor{gray!6}{636} & \cellcolor{gray!6}{614} & \cellcolor{gray!6}{596} & \cellcolor{gray!6}{1,834} & \cellcolor{gray!6}{0.3}\\
\hspace{1em}India – North-East & 9,043 & 7,537 & 6,899 & 63,805 & 0.9\\
\cellcolor{gray!6}{\hspace{1em}India – West. Ghats} & \cellcolor{gray!6}{3,300} & \cellcolor{gray!6}{2,846} & \cellcolor{gray!6}{2,305} & \cellcolor{gray!6}{54,132} & \cellcolor{gray!6}{2.1}\\
\hspace{1em}Indonesia & 141,692 & 128,119 & 116,712 & 1,140,663 & 0.9\\
\cellcolor{gray!6}{\hspace{1em}Laos} & \cellcolor{gray!6}{13,295} & \cellcolor{gray!6}{10,970} & \cellcolor{gray!6}{9,188} & \cellcolor{gray!6}{178,200} & \cellcolor{gray!6}{1.8}\\
\hspace{1em}Malaysia & 26,062 & 22,594 & 20,213 & 238,167 & 1.1\\
\cellcolor{gray!6}{\hspace{1em}Myanmar} & \cellcolor{gray!6}{21,969} & \cellcolor{gray!6}{18,387} & \cellcolor{gray!6}{16,080} & \cellcolor{gray!6}{230,647} & \cellcolor{gray!6}{1.3}\\
\hspace{1em}New Caledonia & 1,042 & 1,016 & 982 & 3,459 & 0.3\\
\cellcolor{gray!6}{\hspace{1em}Papua New Guinea} & \cellcolor{gray!6}{41,153} & \cellcolor{gray!6}{40,433} & \cellcolor{gray!6}{39,758} & \cellcolor{gray!6}{67,522} & \cellcolor{gray!6}{0.2}\\
\hspace{1em}Philippines & 15,608 & 14,423 & 13,274 & 114,923 & 0.8\\
\cellcolor{gray!6}{\hspace{1em}Singapore} & \cellcolor{gray!6}{17} & \cellcolor{gray!6}{15} & \cellcolor{gray!6}{14} & \cellcolor{gray!6}{172} & \cellcolor{gray!6}{1.2}\\
\hspace{1em}Solomon Isl. & 2,830 & 2,827 & 2,806 & 2,112 & 0.1\\
\cellcolor{gray!6}{\hspace{1em}Sri Lanka} & \cellcolor{gray!6}{2,151} & \cellcolor{gray!6}{1,784} & \cellcolor{gray!6}{1,604} & \cellcolor{gray!6}{18,007} & \cellcolor{gray!6}{1.1}\\
\hspace{1em}Thailand & 7,767 & 6,804 & 6,164 & 64,040 & 1.0\\
\cellcolor{gray!6}{\hspace{1em}Timor-Leste} & \cellcolor{gray!6}{137} & \cellcolor{gray!6}{92} & \cellcolor{gray!6}{77} & \cellcolor{gray!6}{1,500} & \cellcolor{gray!6}{1.8}\\
\hspace{1em}Vanuatu & 1,256 & 1,256 & 1,247 & 910 & 0.1\\
\cellcolor{gray!6}{\hspace{1em}Vietnam} & \cellcolor{gray!6}{12,172} & \cellcolor{gray!6}{9,739} & \cellcolor{gray!6}{8,322} & \cellcolor{gray!6}{141,721} & \cellcolor{gray!6}{1.6}\\*
\end{longtable}
\endgroup{}

(ref:cap-fcc-hist-reg) **Historical forest cover change per region and continent**. Areas of forest cover are given in thousand hectares (Kha) for the years 2000, 2010 and 2020 ("fc2000", "fc2010", and "fc2020", respectively). The mean annual deforested area $d$ for the period 2010--2020 is given in hectare per year (ha/yr). For comparing the deforestation intensity between study-areas, the corresponding mean annual deforestation rate $p$ on the period 2010--2020 is also provided in percent per year (%/yr), with one decimal precision. Estimates for America include Brazil, and estimates for Asia include India. Around 7.5 Mha (75,000 km$^2$, about the size of Scotland or South Carolina) of natural old-growth moist tropical forest have been disappearing each year in the period 2010--2020.\vspace{0.5cm}

\begin{table}[H]

\caption{(\#tab:fcc-hist-reg)(ref:cap-fcc-hist-reg)}
\centering
\fontsize{11}{13}\selectfont
\begin{tabular}[t]{lrrrrr}
\toprule
Region & fc2000 & fc2010 & fc2020 & $d$ (ha/yr) & $p$ (\%/yr)\\
\midrule
\cellcolor{gray!6}{India} & \cellcolor{gray!6}{12,978} & \cellcolor{gray!6}{10,998} & \cellcolor{gray!6}{9,800} & \cellcolor{gray!6}{119,771} & \cellcolor{gray!6}{1.1}\\
Brazil & 384,341 & 357,113 & 341,761 & 1,535,198 & 0.4\\
\cellcolor{gray!6}{America} & \cellcolor{gray!6}{706,748} & \cellcolor{gray!6}{663,278} & \cellcolor{gray!6}{634,302} & \cellcolor{gray!6}{2,897,581} & \cellcolor{gray!6}{0.4}\\
Africa & 277,854 & 261,466 & 240,292 & 2,117,372 & 0.8\\
\cellcolor{gray!6}{Asia} & \cellcolor{gray!6}{312,797} & \cellcolor{gray!6}{280,532} & \cellcolor{gray!6}{255,810} & \cellcolor{gray!6}{2,472,128} & \cellcolor{gray!6}{0.9}\\
All continents & 1,297,398 & 1,205,276 & 1,130,405 & 7,487,081 & 0.6\\
\bottomrule
\end{tabular}
\end{table}

\newpage

<!------------------------------>
<!-- Forest cover projections -->
<!------------------------------>

## Forest cover projections

(ref:cap-fcc-proj) **Forest cover projections for each study-area**. Projected areas of forest cover are given in thousand hectares (Kha) for four years in the future (2040, 2060, 2080, and 2100). Projections were made using the forest cover in 2020 and the mean annual deforested area on the period 2010--2020 ("fc2000" and $d$ respectively in Table \@ref(tab:fcc-hist)), assuming a "business-as-usual" scenario of deforestation. Column "loss21" indicates the projected percentage of forest cover loss during the 21$^\text{st}$ century (2100 vs. 2000). Column "yrdis" indicates the estimated year at which all the forest of the study-area will have disappeared.\vspace{0.5cm}

\begingroup\fontsize{11}{13}\selectfont

\begin{longtable}[t]{lrrrrrr}
\caption{(\#tab:fcc-proj)(ref:cap-fcc-proj)}\\
\toprule
Country -- Study-area & fc2040 & fc2060 & fc2080 & fc2100 & loss21 (\%) & yrdis\\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
Country -- Study-area & fc2040 & fc2060 & fc2080 & fc2100 & loss21 (\%) & yrdis\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{7}{l}{\textbf{America}}\\
\cellcolor{gray!6}{\hspace{1em}Antigua and B.} & \cellcolor{gray!6}{2} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2078}\\
\hspace{1em}Bahamas & 80 & 35 & 0 & 0 & 100 & 2075\\
\cellcolor{gray!6}{\hspace{1em}Barbados} & \cellcolor{gray!6}{2} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2068}\\
\hspace{1em}Belize & 1,030 & 742 & 453 & 165 & 89 & 2111\\
\cellcolor{gray!6}{\hspace{1em}Bolivia} & \cellcolor{gray!6}{26,205} & \cellcolor{gray!6}{21,935} & \cellcolor{gray!6}{17,664} & \cellcolor{gray!6}{13,393} & \cellcolor{gray!6}{62} & \cellcolor{gray!6}{2162}\\
\hspace{1em}Brazil – Acre & 12,201 & 11,116 & 9,851 & 8,443 & 41 & 2148\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Alagoas} & \cellcolor{gray!6}{68} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2053}\\
\hspace{1em}Brazil – Amapa & 11,204 & 10,818 & 10,251 & 9,542 & 19 & 2157\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Amazonas} & \cellcolor{gray!6}{143,271} & \cellcolor{gray!6}{139,927} & \cellcolor{gray!6}{136,403} & \cellcolor{gray!6}{132,736} & \cellcolor{gray!6}{11} & \cellcolor{gray!6}{2242}\\
\hspace{1em}Brazil – Bahia & 1,584 & 1,004 & 244 & 0 & 100 & 2086\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Ceara} & \cellcolor{gray!6}{20} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2047}\\
\hspace{1em}Brazil – Espirito Santo & 339 & 121 & 0 & 0 & 100 & 2067\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Goias} & \cellcolor{gray!6}{65} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2044}\\
\hspace{1em}Brazil – Maranhao & 1,530 & 0 & 0 & 0 & 100 & 2057\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Mato Grosso} & \cellcolor{gray!6}{25,341} & \cellcolor{gray!6}{18,883} & \cellcolor{gray!6}{12,244} & \cellcolor{gray!6}{5,463} & \cellcolor{gray!6}{87} & \cellcolor{gray!6}{2115}\\
\hspace{1em}Brazil – Mato Grosso do Sul & 546 & 221 & 0 & 0 & 100 & 2070\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Minas Gerais} & \cellcolor{gray!6}{212} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2045}\\
\hspace{1em}Brazil – Para & 79,466 & 70,311 & 60,976 & 51,498 & 49 & 2176\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Paraiba} & \cellcolor{gray!6}{32} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2049}\\
\hspace{1em}Brazil – Parana & 2,392 & 1,668 & 764 & 0 & 100 & 2096\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Pernambouco} & \cellcolor{gray!6}{88} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2056}\\
\hspace{1em}Brazil – Piaui & 2 & 0 & 0 & 0 & 100 & 2041\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio de Janeiro} & \cellcolor{gray!6}{677} & \cellcolor{gray!6}{423} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2080}\\
\hspace{1em}Brazil – Rio Grande do Norte & 16 & 0 & 0 & 0 & 100 & 2047\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio Grande do Sul} & \cellcolor{gray!6}{2,276} & \cellcolor{gray!6}{1,721} & \cellcolor{gray!6}{985} & \cellcolor{gray!6}{108} & \cellcolor{gray!6}{97} & \cellcolor{gray!6}{2103}\\
\hspace{1em}Brazil – Rondonia & 10,205 & 7,324 & 4,263 & 1,060 & 94 & 2107\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Roraima} & \cellcolor{gray!6}{14,091} & \cellcolor{gray!6}{12,513} & \cellcolor{gray!6}{10,754} & \cellcolor{gray!6}{8,853} & \cellcolor{gray!6}{47} & \cellcolor{gray!6}{2146}\\
\hspace{1em}Brazil – Santa Catarina & 2,527 & 1,977 & 1,247 & 374 & 90 & 2107\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Sao Paulo} & \cellcolor{gray!6}{2,803} & \cellcolor{gray!6}{2,326} & \cellcolor{gray!6}{1,668} & \cellcolor{gray!6}{869} & \cellcolor{gray!6}{76} & \cellcolor{gray!6}{2114}\\
\hspace{1em}Brazil – Sergipe & 42 & 0 & 0 & 0 & 100 & 2050\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Tocantins} & \cellcolor{gray!6}{61} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2042}\\
\hspace{1em}Colombia & 58,804 & 53,125 & 47,446 & 41,767 & 41 & 2247\\
\cellcolor{gray!6}{\hspace{1em}Costa Rica} & \cellcolor{gray!6}{1,826} & \cellcolor{gray!6}{1,470} & \cellcolor{gray!6}{1,115} & \cellcolor{gray!6}{759} & \cellcolor{gray!6}{70} & \cellcolor{gray!6}{2142}\\
\hspace{1em}Cuba & 1,086 & 784 & 482 & 180 & 90 & 2111\\
\cellcolor{gray!6}{\hspace{1em}Dominica} & \cellcolor{gray!6}{67} & \cellcolor{gray!6}{62} & \cellcolor{gray!6}{56} & \cellcolor{gray!6}{51} & \cellcolor{gray!6}{33} & \cellcolor{gray!6}{2284}\\
\hspace{1em}Dominican Rep. & 631 & 293 & 0 & 0 & 100 & 2077\\
\cellcolor{gray!6}{\hspace{1em}Ecuador} & \cellcolor{gray!6}{13,250} & \cellcolor{gray!6}{12,104} & \cellcolor{gray!6}{10,959} & \cellcolor{gray!6}{9,814} & \cellcolor{gray!6}{37} & \cellcolor{gray!6}{2271}\\
\hspace{1em}El Salvador & 72 & 42 & 12 & 0 & 100 & 2088\\
\cellcolor{gray!6}{\hspace{1em}French Guiana} & \cellcolor{gray!6}{8,051} & \cellcolor{gray!6}{7,997} & \cellcolor{gray!6}{7,942} & \cellcolor{gray!6}{7,888} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{5012}\\
\hspace{1em}Grenada & 14 & 7 & 1 & 0 & 100 & 2082\\
\cellcolor{gray!6}{\hspace{1em}Guadeloupe} & \cellcolor{gray!6}{69} & \cellcolor{gray!6}{59} & \cellcolor{gray!6}{49} & \cellcolor{gray!6}{39} & \cellcolor{gray!6}{55} & \cellcolor{gray!6}{2178}\\
\hspace{1em}Guatemala & 1,509 & 542 & 0 & 0 & 100 & 2071\\
\cellcolor{gray!6}{\hspace{1em}Guyana} & \cellcolor{gray!6}{18,224} & \cellcolor{gray!6}{17,937} & \cellcolor{gray!6}{17,650} & \cellcolor{gray!6}{17,362} & \cellcolor{gray!6}{8} & \cellcolor{gray!6}{3309}\\
\hspace{1em}Haiti & 35 & 0 & 0 & 0 & 100 & 2046\\
\cellcolor{gray!6}{\hspace{1em}Honduras} & \cellcolor{gray!6}{1,800} & \cellcolor{gray!6}{842} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2077}\\
\hspace{1em}Jamaica & 368 & 298 & 229 & 159 & 70 & 2145\\
\cellcolor{gray!6}{\hspace{1em}Martinique} & \cellcolor{gray!6}{61} & \cellcolor{gray!6}{51} & \cellcolor{gray!6}{41} & \cellcolor{gray!6}{31} & \cellcolor{gray!6}{60} & \cellcolor{gray!6}{2163}\\
\hspace{1em}Mexico & 4,189 & 1,422 & 0 & 0 & 100 & 2070\\
\cellcolor{gray!6}{\hspace{1em}Montserrat} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{30} & \cellcolor{gray!6}{2318}\\
\hspace{1em}Nicaragua & 1,587 & 0 & 0 & 0 & 100 & 2056\\
\cellcolor{gray!6}{\hspace{1em}Panama} & \cellcolor{gray!6}{3,647} & \cellcolor{gray!6}{3,205} & \cellcolor{gray!6}{2,763} & \cellcolor{gray!6}{2,320} & \cellcolor{gray!6}{49} & \cellcolor{gray!6}{2204}\\
\hspace{1em}Paraguay & 367 & 0 & 0 & 0 & 100 & 2047\\
\cellcolor{gray!6}{\hspace{1em}Peru} & \cellcolor{gray!6}{69,457} & \cellcolor{gray!6}{66,777} & \cellcolor{gray!6}{64,097} & \cellcolor{gray!6}{61,417} & \cellcolor{gray!6}{18} & \cellcolor{gray!6}{2558}\\
\hspace{1em}Puerto Rico & 250 & 147 & 43 & 0 & 100 & 2088\\
\cellcolor{gray!6}{\hspace{1em}Saint Kitts and N.} & \cellcolor{gray!6}{9} & \cellcolor{gray!6}{8} & \cellcolor{gray!6}{7} & \cellcolor{gray!6}{6} & \cellcolor{gray!6}{40} & \cellcolor{gray!6}{2244}\\
\hspace{1em}Saint Lucia & 40 & 34 & 28 & 22 & 57 & 2170\\
\cellcolor{gray!6}{\hspace{1em}Saint Martin} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2028}\\
\hspace{1em}Saint Vincent & 27 & 24 & 22 & 20 & 35 & 2274\\
\cellcolor{gray!6}{\hspace{1em}Sint Maarten} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2029}\\
\hspace{1em}Suriname & 13,470 & 13,240 & 13,010 & 12,781 & 8 & 3213\\
\cellcolor{gray!6}{\hspace{1em}Trinidad and Tobago} & \cellcolor{gray!6}{261} & \cellcolor{gray!6}{204} & \cellcolor{gray!6}{148} & \cellcolor{gray!6}{92} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{2132}\\
\hspace{1em}Venezuela & 38,799 & 35,681 & 32,562 & 29,444 & 35 & 2288\\
\cellcolor{gray!6}{\hspace{1em}Virgin Isl. UK} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2037}\\
\hspace{1em}Virgin Isl. US & 3 & 0 & 0 & 0 & 100 & 2056\\
\addlinespace[0.3em]
\multicolumn{7}{l}{\textbf{Africa}}\\
\cellcolor{gray!6}{\hspace{1em}Angola} & \cellcolor{gray!6}{3,623} & \cellcolor{gray!6}{1,880} & \cellcolor{gray!6}{138} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2081}\\
\hspace{1em}Benin & 3 & 0 & 0 & 0 & 100 & 2041\\
\cellcolor{gray!6}{\hspace{1em}Burundi} & \cellcolor{gray!6}{35} & \cellcolor{gray!6}{15} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2074}\\
\hspace{1em}Cameroon & 21,591 & 20,194 & 18,796 & 17,398 & 28 & 2348\\
\cellcolor{gray!6}{\hspace{1em}CAR} & \cellcolor{gray!6}{7,824} & \cellcolor{gray!6}{6,763} & \cellcolor{gray!6}{5,703} & \cellcolor{gray!6}{4,642} & \cellcolor{gray!6}{53} & \cellcolor{gray!6}{2187}\\
\hspace{1em}Comoros & 76 & 66 & 56 & 46 & 49 & 2194\\
\cellcolor{gray!6}{\hspace{1em}Congo} & \cellcolor{gray!6}{22,222} & \cellcolor{gray!6}{21,033} & \cellcolor{gray!6}{19,844} & \cellcolor{gray!6}{18,654} & \cellcolor{gray!6}{23} & \cellcolor{gray!6}{2413}\\
\hspace{1em}DRC & 101,325 & 84,538 & 67,751 & 50,965 & 61 & 2160\\
\cellcolor{gray!6}{\hspace{1em}Eq. Guinea} & \cellcolor{gray!6}{2,547} & \cellcolor{gray!6}{2,481} & \cellcolor{gray!6}{2,416} & \cellcolor{gray!6}{2,350} & \cellcolor{gray!6}{11} & \cellcolor{gray!6}{2816}\\
\hspace{1em}Ethiopia & 690 & 0 & 0 & 0 & 100 & 2049\\
\cellcolor{gray!6}{\hspace{1em}Gabon} & \cellcolor{gray!6}{23,721} & \cellcolor{gray!6}{23,452} & \cellcolor{gray!6}{23,184} & \cellcolor{gray!6}{22,916} & \cellcolor{gray!6}{5} & \cellcolor{gray!6}{3807}\\
\hspace{1em}Gambia & 35 & 26 & 17 & 8 & 85 & 2118\\
\cellcolor{gray!6}{\hspace{1em}Ghana} & \cellcolor{gray!6}{1,053} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2049}\\
\hspace{1em}Guinea & 45 & 0 & 0 & 0 & 100 & 2041\\
\cellcolor{gray!6}{\hspace{1em}Guinea Bissau} & \cellcolor{gray!6}{210} & \cellcolor{gray!6}{117} & \cellcolor{gray!6}{23} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2085}\\
\hspace{1em}Ivory Coast & 0 & 0 & 0 & 0 & 100 & 2035\\
\cellcolor{gray!6}{\hspace{1em}Kenya} & \cellcolor{gray!6}{444} & \cellcolor{gray!6}{139} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2069}\\
\hspace{1em}Liberia & 6,071 & 4,240 & 2,409 & 578 & 94 & 2106\\
\cellcolor{gray!6}{\hspace{1em}Madagascar} & \cellcolor{gray!6}{2,522} & \cellcolor{gray!6}{42} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2060}\\
\hspace{1em}Malawi & 0 & 0 & 0 & 0 & 100 & 2031\\
\cellcolor{gray!6}{\hspace{1em}Mauritius} & \cellcolor{gray!6}{41} & \cellcolor{gray!6}{32} & \cellcolor{gray!6}{23} & \cellcolor{gray!6}{14} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{2132}\\
\hspace{1em}Mayotte & 2 & 0 & 0 & 0 & 100 & 2044\\
\cellcolor{gray!6}{\hspace{1em}Nigeria} & \cellcolor{gray!6}{4,335} & \cellcolor{gray!6}{2,344} & \cellcolor{gray!6}{352} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2083}\\
\hspace{1em}Reunion & 138 & 120 & 103 & 86 & 48 & 2197\\
\cellcolor{gray!6}{\hspace{1em}Rwanda} & \cellcolor{gray!6}{81} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2060}\\
\hspace{1em}Senegal & 114 & 97 & 81 & 64 & 56 & 2177\\
\cellcolor{gray!6}{\hspace{1em}Sierra Leone} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2034}\\
\hspace{1em}South Sudan & 78 & 0 & 0 & 0 & 100 & 2058\\
\cellcolor{gray!6}{\hspace{1em}Tanzania} & \cellcolor{gray!6}{866} & \cellcolor{gray!6}{626} & \cellcolor{gray!6}{386} & \cellcolor{gray!6}{146} & \cellcolor{gray!6}{90} & \cellcolor{gray!6}{2112}\\
\hspace{1em}Togo & 0 & 0 & 0 & 0 & 100 & 2036\\
\cellcolor{gray!6}{\hspace{1em}Uganda} & \cellcolor{gray!6}{25} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2040}\\
\hspace{1em}Zambia & 8 & 0 & 0 & 0 & 100 & 2042\\
\addlinespace[0.3em]
\multicolumn{7}{l}{\textbf{Asia}}\\
\cellcolor{gray!6}{\hspace{1em}Australia – Queensland} & \cellcolor{gray!6}{1,674} & \cellcolor{gray!6}{1,411} & \cellcolor{gray!6}{1,147} & \cellcolor{gray!6}{884} & \cellcolor{gray!6}{61} & \cellcolor{gray!6}{2167}\\
\hspace{1em}Bangladesh & 803 & 690 & 577 & 464 & 59 & 2182\\
\cellcolor{gray!6}{\hspace{1em}Bhutan} & \cellcolor{gray!6}{2,040} & \cellcolor{gray!6}{1,809} & \cellcolor{gray!6}{1,578} & \cellcolor{gray!6}{1,348} & \cellcolor{gray!6}{47} & \cellcolor{gray!6}{2216}\\
\hspace{1em}Brunei & 482 & 466 & 450 & 435 & 16 & 2654\\
\cellcolor{gray!6}{\hspace{1em}Cambodia} & \cellcolor{gray!6}{635} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2045}\\
\hspace{1em}Fiji & 931 & 843 & 756 & 668 & 40 & 2252\\
\cellcolor{gray!6}{\hspace{1em}India – Andaman and N.} & \cellcolor{gray!6}{559} & \cellcolor{gray!6}{523} & \cellcolor{gray!6}{486} & \cellcolor{gray!6}{449} & \cellcolor{gray!6}{29} & \cellcolor{gray!6}{2345}\\
\hspace{1em}India – North-East & 5,623 & 4,347 & 3,071 & 1,795 & 80 & 2128\\
\cellcolor{gray!6}{\hspace{1em}India – West. Ghats} & \cellcolor{gray!6}{1,222} & \cellcolor{gray!6}{139} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2062}\\
\hspace{1em}Indonesia & 93,899 & 71,086 & 48,273 & 25,459 & 82 & 2122\\
\cellcolor{gray!6}{\hspace{1em}Laos} & \cellcolor{gray!6}{5,624} & \cellcolor{gray!6}{2,060} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2071}\\
\hspace{1em}Malaysia & 15,449 & 10,686 & 5,923 & 1,159 & 96 & 2104\\
\cellcolor{gray!6}{\hspace{1em}Myanmar} & \cellcolor{gray!6}{11,467} & \cellcolor{gray!6}{6,854} & \cellcolor{gray!6}{2,241} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2089}\\
\hspace{1em}New Caledonia & 912 & 843 & 774 & 705 & 32 & 2303\\
\cellcolor{gray!6}{\hspace{1em}Papua New Guinea} & \cellcolor{gray!6}{38,408} & \cellcolor{gray!6}{37,057} & \cellcolor{gray!6}{35,707} & \cellcolor{gray!6}{34,356} & \cellcolor{gray!6}{17} & \cellcolor{gray!6}{2608}\\
\hspace{1em}Philippines & 10,975 & 8,677 & 6,378 & 4,080 & 74 & 2135\\
\cellcolor{gray!6}{\hspace{1em}Singapore} & \cellcolor{gray!6}{10} & \cellcolor{gray!6}{7} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2098}\\
\hspace{1em}Solomon Isl. & 2,763 & 2,721 & 2,679 & 2,637 & 7 & 3348\\
\cellcolor{gray!6}{\hspace{1em}Sri Lanka} & \cellcolor{gray!6}{1,244} & \cellcolor{gray!6}{884} & \cellcolor{gray!6}{523} & \cellcolor{gray!6}{163} & \cellcolor{gray!6}{92} & \cellcolor{gray!6}{2109}\\
\hspace{1em}Thailand & 4,883 & 3,602 & 2,322 & 1,041 & 87 & 2116\\
\cellcolor{gray!6}{\hspace{1em}Timor-Leste} & \cellcolor{gray!6}{47} & \cellcolor{gray!6}{17} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2071}\\
\hspace{1em}Vanuatu & 1,229 & 1,210 & 1,192 & 1,174 & 7 & 3390\\
\cellcolor{gray!6}{\hspace{1em}Vietnam} & \cellcolor{gray!6}{5,487} & \cellcolor{gray!6}{2,653} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{100} & \cellcolor{gray!6}{2078}\\*
\end{longtable}
\endgroup{}

(ref:cap-fcc-proj-reg) **Forest cover projections per region and continent**. Projected areas of forest cover are given in thousand hectares (Kha) for four dates in the future (2040, 2060, 2080, and 2100). Projections were made using the forest cover in 2020 and the mean annual deforested area on the period 2010--2020 ("fc2000" and $d$ respectively in Table \@ref(tab:fcc-hist-reg)), assuming a "business-as-usual" scenario of deforestation. Column "loss21" indicates the projected percentage of forest cover loss during the 21$^\text{st}$ century (2100 vs. 2000). At the continental level, it makes less sense to compute the year at which all the forest will have disappeared, as some countries might conserve forest for a very long time, even though they account for a very small proportion of the total forest area at the continental scale. Instead, we computed the estimated year at which 75% of the forest cover in 2000 will have disappeared ("yr75dis").\vspace{0.5cm}

\begin{table}[H]

\caption{(\#tab:fcc-proj-reg)(ref:cap-fcc-proj-reg)}
\centering
\fontsize{11}{13}\selectfont
\begin{tabular}[t]{lrrrrrr}
\toprule
Region & fc2040 & fc2060 & fc2080 & fc2100 & loss21 (\%) & yr75dis\\
\midrule
\cellcolor{gray!6}{India} & \cellcolor{gray!6}{7,405} & \cellcolor{gray!6}{5,009} & \cellcolor{gray!6}{3,557} & \cellcolor{gray!6}{2,244} & \cellcolor{gray!6}{83} & \cellcolor{gray!6}{2064}\\
Brazil & 311,057 & 280,353 & 249,649 & 218,945 & 43 & 2187\\
\cellcolor{gray!6}{America} & \cellcolor{gray!6}{576,351} & \cellcolor{gray!6}{519,425} & \cellcolor{gray!6}{466,432} & \cellcolor{gray!6}{416,659} & \cellcolor{gray!6}{41} & \cellcolor{gray!6}{2201}\\
Africa & 199,725 & 168,209 & 141,282 & 117,867 & 58 & 2145\\
\cellcolor{gray!6}{Asia} & \cellcolor{gray!6}{206,368} & \cellcolor{gray!6}{158,586} & \cellcolor{gray!6}{114,082} & \cellcolor{gray!6}{76,818} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{2099}\\
All continents & 982,444 & 846,220 & 721,796 & 611,344 & 53 & 2173\\
\bottomrule
\end{tabular}
\end{table}

<!------------------------------------------------->
<!-- Parameter estimates and variable importance -->
<!------------------------------------------------->

## Parameter estimates and variable importance

(ref:cap-par) **Parameter estimates for each study area**. For each study area, we computed the posterior mean of each parameter ("int": intercept, "pa": protected area effect, "elev", "slope", "ddefor","dedge", "driver", "droad", "dtown": slope parameters associated to elevation, slope, distance to past deforestation, distance to forest edge, distance to nearest river, distance to nearest road, and distance to nearest town, respectively, "Vrho": variance of the spatial random effects). Continuous explanatory variables were normalized (mean=0 and standard-deviation=1), allowing us to estimate relative variable importance in determining the spatial probability of deforestation from parameter values. Comparison can be done within and between study-areas.\vspace{0.5cm}

\begingroup\fontsize{10}{12}\selectfont

\begin{longtable}[t]{lrrrrrrrrrr}
\caption{(\#tab:par)(ref:cap-par)}\\
\toprule
Study-area & int & pa & elev & slope & ddefor & dedge & driver & droad & dtown & Vrho\\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
Study-area & int & pa & elev & slope & ddefor & dedge & driver & droad & dtown & Vrho\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{11}{l}{\textbf{America}}\\
\cellcolor{gray!6}{\hspace{1em}ATG} & \cellcolor{gray!6}{0.599} & \cellcolor{gray!6}{-0.259} & \cellcolor{gray!6}{-0.247} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.719} & \cellcolor{gray!6}{-1.090} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{10.000}\\
\hspace{1em}BHS & -1.010 & -- & -- & -0.050 & -1.890 & -0.616 & -- & -- & -0.083 & 6.930\\
\cellcolor{gray!6}{\hspace{1em}BRB} & \cellcolor{gray!6}{-0.477} & \cellcolor{gray!6}{-0.277} & \cellcolor{gray!6}{-0.484} & \cellcolor{gray!6}{-0.067} & \cellcolor{gray!6}{-0.308} & \cellcolor{gray!6}{-1.560} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{0.727}\\
\hspace{1em}BLZ & -1.300 & -0.828 & -0.020 & -0.323 & -1.560 & -1.710 & -- & -0.169 & -0.301 & 6.640\\
\cellcolor{gray!6}{\hspace{1em}BOL} & \cellcolor{gray!6}{-0.717} & \cellcolor{gray!6}{-0.204} & \cellcolor{gray!6}{-0.303} & \cellcolor{gray!6}{-0.122} & \cellcolor{gray!6}{-0.684} & \cellcolor{gray!6}{-3.430} & \cellcolor{gray!6}{-0.047} & \cellcolor{gray!6}{-0.333} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.080}\\
\hspace{1em}AC & -4.170 & -0.723 & -- & -0.006 & -5.460 & -5.170 & -0.004 & -- & -0.308 & 3.930\\
\cellcolor{gray!6}{\hspace{1em}AL} & \cellcolor{gray!6}{0.718} & \cellcolor{gray!6}{-0.248} & \cellcolor{gray!6}{-0.427} & \cellcolor{gray!6}{-0.010} & \cellcolor{gray!6}{-0.480} & \cellcolor{gray!6}{-1.650} & \cellcolor{gray!6}{-0.114} & \cellcolor{gray!6}{-0.116} & \cellcolor{gray!6}{-0.207} & \cellcolor{gray!6}{3.510}\\
\hspace{1em}AP & -5.800 & -0.157 & -0.638 & -- & -2.430 & -9.690 & -- & -0.293 & -0.377 & 6.290\\
\cellcolor{gray!6}{\hspace{1em}AM} & \cellcolor{gray!6}{-4.000} & \cellcolor{gray!6}{-0.558} & \cellcolor{gray!6}{-0.076} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-2.710} & \cellcolor{gray!6}{-4.650} & \cellcolor{gray!6}{-0.196} & \cellcolor{gray!6}{-0.721} & \cellcolor{gray!6}{-0.472} & \cellcolor{gray!6}{12.400}\\
\hspace{1em}BA & 1.320 & -0.171 & -0.381 & -0.174 & -0.645 & -1.240 & -- & -0.050 & -0.019 & 3.710\\
\cellcolor{gray!6}{\hspace{1em}CE} & \cellcolor{gray!6}{1.100} & \cellcolor{gray!6}{-0.654} & \cellcolor{gray!6}{-1.010} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.058} & \cellcolor{gray!6}{-1.680} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{12.600}\\
\hspace{1em}ES & -1.230 & -0.056 & -0.973 & -0.007 & -0.638 & -0.996 & -0.012 & -0.068 & -0.008 & 2.780\\
\cellcolor{gray!6}{\hspace{1em}GO} & \cellcolor{gray!6}{0.016} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.571} & \cellcolor{gray!6}{-0.529} & \cellcolor{gray!6}{-0.112} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.023} & \cellcolor{gray!6}{3.920}\\
\hspace{1em}MA & -2.070 & -0.056 & -- & -0.099 & -3.690 & -4.630 & -- & -- & -0.133 & 3.270\\
\cellcolor{gray!6}{\hspace{1em}MT} & \cellcolor{gray!6}{-0.438} & \cellcolor{gray!6}{-0.477} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.995} & \cellcolor{gray!6}{-1.140} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.252} & \cellcolor{gray!6}{-0.291} & \cellcolor{gray!6}{9.680}\\
\hspace{1em}MS & -0.249 & -- & -0.229 & -- & -0.294 & -1.280 & -- & -- & -- & 3.430\\
\cellcolor{gray!6}{\hspace{1em}MG} & \cellcolor{gray!6}{0.228} & \cellcolor{gray!6}{-0.113} & \cellcolor{gray!6}{-0.444} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.504} & \cellcolor{gray!6}{-0.571} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.069} & \cellcolor{gray!6}{-0.032} & \cellcolor{gray!6}{3.730}\\
\hspace{1em}PA & -1.470 & -0.992 & -- & -0.045 & -1.520 & -2.530 & -- & -0.544 & -0.446 & 7.140\\
\cellcolor{gray!6}{\hspace{1em}PB} & \cellcolor{gray!6}{2.360} & \cellcolor{gray!6}{-0.453} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.707} & \cellcolor{gray!6}{-1.310} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.174} & \cellcolor{gray!6}{4.060}\\
\hspace{1em}PR & -0.857 & -0.206 & -- & -0.239 & -0.858 & -2.570 & -0.028 & -- & -- & 2.910\\
\cellcolor{gray!6}{\hspace{1em}PE} & \cellcolor{gray!6}{3.710} & \cellcolor{gray!6}{-0.405} & \cellcolor{gray!6}{-0.513} & \cellcolor{gray!6}{-0.093} & \cellcolor{gray!6}{-0.596} & \cellcolor{gray!6}{-2.140} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.148} & \cellcolor{gray!6}{-0.076} & \cellcolor{gray!6}{3.350}\\
\hspace{1em}PI & 0.033 & -0.144 & -- & -- & -0.448 & -0.789 & -0.052 & -0.041 & -0.002 & 3.460\\
\cellcolor{gray!6}{\hspace{1em}RJ} & \cellcolor{gray!6}{-0.460} & \cellcolor{gray!6}{-0.205} & \cellcolor{gray!6}{-0.177} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.070} & \cellcolor{gray!6}{-3.010} & \cellcolor{gray!6}{-0.056} & \cellcolor{gray!6}{-0.100} & \cellcolor{gray!6}{-0.010} & \cellcolor{gray!6}{2.400}\\
\hspace{1em}RN & 1.120 & -0.301 & -- & -- & -1.110 & -0.657 & -- & -0.094 & -- & 17.100\\
\cellcolor{gray!6}{\hspace{1em}RS} & \cellcolor{gray!6}{0.072} & \cellcolor{gray!6}{-0.143} & \cellcolor{gray!6}{-0.075} & \cellcolor{gray!6}{-0.525} & \cellcolor{gray!6}{-0.551} & \cellcolor{gray!6}{-1.500} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{1.730}\\
\hspace{1em}RO & -0.882 & -1.340 & -- & -0.086 & -0.662 & -1.840 & -- & -0.240 & -0.310 & 7.120\\
\cellcolor{gray!6}{\hspace{1em}RR} & \cellcolor{gray!6}{-1.610} & \cellcolor{gray!6}{-0.670} & \cellcolor{gray!6}{-0.445} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.550} & \cellcolor{gray!6}{-2.800} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.989} & \cellcolor{gray!6}{7.840}\\
\hspace{1em}SC & -0.690 & -0.232 & -- & -0.479 & -0.688 & -1.450 & -- & -- & -- & 2.100\\
\cellcolor{gray!6}{\hspace{1em}SP} & \cellcolor{gray!6}{-0.700} & \cellcolor{gray!6}{-0.185} & \cellcolor{gray!6}{-0.177} & \cellcolor{gray!6}{-0.232} & \cellcolor{gray!6}{-1.190} & \cellcolor{gray!6}{-2.760} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.095} & \cellcolor{gray!6}{-0.015} & \cellcolor{gray!6}{3.960}\\
\hspace{1em}SE & 0.682 & -- & -- & -- & -0.797 & -1.070 & -0.045 & -- & -- & 3.840\\
\cellcolor{gray!6}{\hspace{1em}TO} & \cellcolor{gray!6}{-0.309} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.394} & \cellcolor{gray!6}{-0.204} & \cellcolor{gray!6}{-0.139} & \cellcolor{gray!6}{0.010} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{4.630}\\
\hspace{1em}COL & -3.480 & -0.535 & -0.687 & -0.310 & -4.210 & -3.720 & -- & -0.817 & -0.286 & 5.320\\
\cellcolor{gray!6}{\hspace{1em}CRI} & \cellcolor{gray!6}{-2.990} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.031} & \cellcolor{gray!6}{-0.370} & \cellcolor{gray!6}{-2.330} & \cellcolor{gray!6}{-6.030} & \cellcolor{gray!6}{-0.031} & \cellcolor{gray!6}{-0.060} & \cellcolor{gray!6}{-0.086} & \cellcolor{gray!6}{2.090}\\
\hspace{1em}CUB & -0.308 & -0.077 & -0.131 & -0.064 & -0.681 & -0.930 & -- & -0.147 & -- & 5.700\\
\cellcolor{gray!6}{\hspace{1em}DMA} & \cellcolor{gray!6}{-1.040} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.258} & \cellcolor{gray!6}{-0.221} & \cellcolor{gray!6}{-0.809} & \cellcolor{gray!6}{-2.360} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.228} & \cellcolor{gray!6}{2.160}\\
\hspace{1em}DOM & -0.211 & -0.357 & -0.392 & -0.029 & -1.110 & -1.410 & -0.021 & -- & -- & 2.670\\
\cellcolor{gray!6}{\hspace{1em}ECU} & \cellcolor{gray!6}{-4.290} & \cellcolor{gray!6}{-0.867} & \cellcolor{gray!6}{-0.124} & \cellcolor{gray!6}{-0.359} & \cellcolor{gray!6}{-1.110} & \cellcolor{gray!6}{-11.500} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.455} & \cellcolor{gray!6}{-0.159} & \cellcolor{gray!6}{2.470}\\
\hspace{1em}SLV & -0.077 & -0.426 & -0.322 & -0.172 & -0.766 & -1.570 & -0.012 & -0.068 & -0.111 & 5.860\\
\cellcolor{gray!6}{\hspace{1em}GUF} & \cellcolor{gray!6}{-1.940} & \cellcolor{gray!6}{-1.070} & \cellcolor{gray!6}{-0.635} & \cellcolor{gray!6}{-0.025} & \cellcolor{gray!6}{-1.430} & \cellcolor{gray!6}{-3.160} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.894} & \cellcolor{gray!6}{-0.730} & \cellcolor{gray!6}{21.400}\\
\hspace{1em}GRD & -0.599 & -0.287 & -0.991 & -- & -1.610 & -1.900 & -- & -0.232 & -- & 5.780\\
\cellcolor{gray!6}{\hspace{1em}GLP} & \cellcolor{gray!6}{-3.980} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.520} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-2.120} & \cellcolor{gray!6}{-7.090} & \cellcolor{gray!6}{-0.012} & \cellcolor{gray!6}{-0.039} & \cellcolor{gray!6}{-0.092} & \cellcolor{gray!6}{3.110}\\
\hspace{1em}GTM & -0.312 & -0.172 & -0.377 & -0.234 & -0.796 & -1.160 & -0.150 & -0.156 & -0.002 & 2.700\\
\cellcolor{gray!6}{\hspace{1em}GUY} & \cellcolor{gray!6}{-1.750} & \cellcolor{gray!6}{-0.475} & \cellcolor{gray!6}{-0.553} & \cellcolor{gray!6}{-0.181} & \cellcolor{gray!6}{-0.736} & \cellcolor{gray!6}{-3.920} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.221} & \cellcolor{gray!6}{-0.391} & \cellcolor{gray!6}{9.980}\\
\hspace{1em}HTI & -0.087 & -0.152 & -- & -0.131 & -0.437 & -0.916 & -0.089 & -0.124 & -0.028 & 2.840\\
\cellcolor{gray!6}{\hspace{1em}HND} & \cellcolor{gray!6}{-1.050} & \cellcolor{gray!6}{-0.403} & \cellcolor{gray!6}{-0.776} & \cellcolor{gray!6}{-0.286} & \cellcolor{gray!6}{-0.529} & \cellcolor{gray!6}{-1.350} & \cellcolor{gray!6}{-0.001} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{3.580}\\
\hspace{1em}JAM & -1.320 & -0.043 & -0.323 & -0.090 & -1.950 & -3.470 & 0.001 & -0.067 & -0.119 & 2.570\\
\cellcolor{gray!6}{\hspace{1em}MTQ} & \cellcolor{gray!6}{-1.750} & \cellcolor{gray!6}{-0.124} & \cellcolor{gray!6}{-0.392} & \cellcolor{gray!6}{-0.213} & \cellcolor{gray!6}{-1.270} & \cellcolor{gray!6}{-4.940} & \cellcolor{gray!6}{-0.011} & \cellcolor{gray!6}{-0.020} & \cellcolor{gray!6}{-0.006} & \cellcolor{gray!6}{1.850}\\
\hspace{1em}MEX & -0.773 & -0.328 & -0.185 & -0.257 & -0.746 & -1.600 & -- & -0.056 & -- & 3.600\\
\cellcolor{gray!6}{\hspace{1em}MSR} & \cellcolor{gray!6}{-6.510} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.243} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.470} & \cellcolor{gray!6}{-5.270} & \cellcolor{gray!6}{-0.329} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{7.120}\\
\hspace{1em}NIC & -1.140 & -- & -0.299 & -0.092 & -0.584 & -1.440 & -- & -0.178 & -- & 3.450\\
\cellcolor{gray!6}{\hspace{1em}PAN} & \cellcolor{gray!6}{-2.230} & \cellcolor{gray!6}{-0.344} & \cellcolor{gray!6}{-0.147} & \cellcolor{gray!6}{-0.230} & \cellcolor{gray!6}{-2.310} & \cellcolor{gray!6}{-3.580} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.521} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{3.170}\\
\hspace{1em}PRY & -0.497 & -0.167 & -0.117 & -0.182 & -0.673 & -0.394 & -- & -0.058 & -- & 3.460\\
\cellcolor{gray!6}{\hspace{1em}PER} & \cellcolor{gray!6}{-2.360} & \cellcolor{gray!6}{-0.684} & \cellcolor{gray!6}{-0.549} & \cellcolor{gray!6}{-0.338} & \cellcolor{gray!6}{-2.940} & \cellcolor{gray!6}{-4.880} & \cellcolor{gray!6}{-0.008} & \cellcolor{gray!6}{-0.531} & \cellcolor{gray!6}{-0.457} & \cellcolor{gray!6}{4.980}\\
\hspace{1em}PRI & -0.034 & -- & -0.264 & -0.090 & -0.650 & -0.741 & -- & -- & -0.061 & 2.470\\
\cellcolor{gray!6}{\hspace{1em}KNA} & \cellcolor{gray!6}{-1.330} & \cellcolor{gray!6}{-0.137} & \cellcolor{gray!6}{-0.924} & \cellcolor{gray!6}{-0.405} & \cellcolor{gray!6}{-0.032} & \cellcolor{gray!6}{-1.220} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.780}\\
\hspace{1em}LCA & -1.350 & -- & -0.278 & -0.152 & -0.766 & -3.290 & -- & -0.044 & -0.077 & 3.900\\
\cellcolor{gray!6}{\hspace{1em}MAF} & \cellcolor{gray!6}{-0.505} & \cellcolor{gray!6}{-1.180} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.334} & \cellcolor{gray!6}{-0.249} & \cellcolor{gray!6}{-0.609} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{7.700}\\
\hspace{1em}VCT & -0.339 & -- & -0.806 & -0.140 & -1.790 & -1.170 & -- & -0.098 & -- & 2.320\\
\cellcolor{gray!6}{\hspace{1em}SXM} & \cellcolor{gray!6}{-0.032} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.533} & \cellcolor{gray!6}{-0.629} & \cellcolor{gray!6}{-0.461} & \cellcolor{gray!6}{-0.523} & \cellcolor{gray!6}{-0.186} & \cellcolor{gray!6}{-0.290} & \cellcolor{gray!6}{30.400}\\
\hspace{1em}SUR & -1.380 & -0.171 & -0.584 & -- & -0.716 & -1.230 & -- & -0.553 & -0.510 & 14.200\\
\cellcolor{gray!6}{\hspace{1em}TTO} & \cellcolor{gray!6}{-2.840} & \cellcolor{gray!6}{-0.292} & \cellcolor{gray!6}{-0.083} & \cellcolor{gray!6}{-0.032} & \cellcolor{gray!6}{-2.420} & \cellcolor{gray!6}{-5.560} & \cellcolor{gray!6}{0.011} & \cellcolor{gray!6}{-0.074} & \cellcolor{gray!6}{-0.163} & \cellcolor{gray!6}{1.790}\\
\hspace{1em}VEN & -3.160 & -0.120 & -0.378 & -0.109 & -1.660 & -7.350 & -0.090 & -0.444 & -0.356 & 5.860\\
\cellcolor{gray!6}{\hspace{1em}VGB} & \cellcolor{gray!6}{0.511} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.386} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.050} & \cellcolor{gray!6}{-0.937} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{12.900}\\
\hspace{1em}VIR & 0.008 & -0.037 & -- & -0.127 & -0.162 & -0.647 & -- & -0.037 & -- & 6.460\\
\addlinespace[0.3em]
\multicolumn{11}{l}{\textbf{Africa}}\\
\cellcolor{gray!6}{\hspace{1em}AGO} & \cellcolor{gray!6}{-0.125} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.080} & \cellcolor{gray!6}{-0.062} & \cellcolor{gray!6}{-3.080} & \cellcolor{gray!6}{-3.050} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.231} & \cellcolor{gray!6}{-0.133} & \cellcolor{gray!6}{3.780}\\
\hspace{1em}BEN & 1.960 & -0.218 & -- & -0.158 & -0.474 & -1.170 & 0.009 & -- & -- & 6.730\\
\cellcolor{gray!6}{\hspace{1em}BDI} & \cellcolor{gray!6}{-1.010} & \cellcolor{gray!6}{-1.860} & \cellcolor{gray!6}{-0.624} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.660} & \cellcolor{gray!6}{-4.010} & \cellcolor{gray!6}{-0.143} & \cellcolor{gray!6}{-0.078} & \cellcolor{gray!6}{-0.035} & \cellcolor{gray!6}{5.060}\\
\hspace{1em}CMR & -0.583 & -0.933 & -0.036 & -0.172 & -1.310 & -3.260 & -0.047 & -0.449 & -0.232 & 5.740\\
\cellcolor{gray!6}{\hspace{1em}CAF} & \cellcolor{gray!6}{-0.998} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.021} & \cellcolor{gray!6}{-4.830} & \cellcolor{gray!6}{-3.850} & \cellcolor{gray!6}{-0.122} & \cellcolor{gray!6}{-0.086} & \cellcolor{gray!6}{-0.273} & \cellcolor{gray!6}{4.270}\\
\hspace{1em}COM & -2.580 & -0.038 & -- & -0.293 & -- & -8.680 & -- & -- & -- & 6.220\\
\cellcolor{gray!6}{\hspace{1em}COG} & \cellcolor{gray!6}{-2.370} & \cellcolor{gray!6}{-0.112} & \cellcolor{gray!6}{-0.076} & \cellcolor{gray!6}{-0.022} & \cellcolor{gray!6}{-0.877} & \cellcolor{gray!6}{-5.930} & \cellcolor{gray!6}{-0.016} & \cellcolor{gray!6}{-0.708} & \cellcolor{gray!6}{-0.325} & \cellcolor{gray!6}{9.550}\\
\hspace{1em}COD & -3.950 & -0.261 & -- & -- & -5.030 & -6.070 & -- & -0.368 & -0.206 & 3.740\\
\cellcolor{gray!6}{\hspace{1em}GNQ} & \cellcolor{gray!6}{-0.980} & \cellcolor{gray!6}{-0.368} & \cellcolor{gray!6}{-0.412} & \cellcolor{gray!6}{-0.382} & \cellcolor{gray!6}{-0.156} & \cellcolor{gray!6}{-2.050} & \cellcolor{gray!6}{-0.117} & \cellcolor{gray!6}{-1.070} & \cellcolor{gray!6}{-0.122} & \cellcolor{gray!6}{10.600}\\
\hspace{1em}ETH & -0.775 & -0.241 & -0.160 & -0.215 & -0.919 & -2.170 & -- & -0.056 & -0.157 & 3.280\\
\cellcolor{gray!6}{\hspace{1em}GAB} & \cellcolor{gray!6}{-2.480} & \cellcolor{gray!6}{-0.206} & \cellcolor{gray!6}{-0.532} & \cellcolor{gray!6}{-0.246} & \cellcolor{gray!6}{-1.030} & \cellcolor{gray!6}{-4.750} & \cellcolor{gray!6}{-0.024} & \cellcolor{gray!6}{-0.720} & \cellcolor{gray!6}{-0.128} & \cellcolor{gray!6}{14.800}\\
\hspace{1em}GMB & 0.626 & -- & -- & -0.110 & -0.517 & -0.713 & -- & -0.341 & -- & 3.980\\
\cellcolor{gray!6}{\hspace{1em}GHA} & \cellcolor{gray!6}{0.298} & \cellcolor{gray!6}{-0.174} & \cellcolor{gray!6}{-0.177} & \cellcolor{gray!6}{-0.131} & \cellcolor{gray!6}{-0.636} & \cellcolor{gray!6}{-1.980} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.093} & \cellcolor{gray!6}{-0.064} & \cellcolor{gray!6}{2.110}\\
\hspace{1em}GIN & -0.935 & -0.272 & -0.027 & -- & -4.170 & -2.580 & -- & -0.073 & -- & 2.140\\
\cellcolor{gray!6}{\hspace{1em}GNB} & \cellcolor{gray!6}{0.665} & \cellcolor{gray!6}{-0.502} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.240} & \cellcolor{gray!6}{-0.966} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.039} & \cellcolor{gray!6}{-0.100} & \cellcolor{gray!6}{5.200}\\
\hspace{1em}CIV & -0.120 & -- & -- & -0.143 & -1.590 & -1.150 & -- & -- & -0.006 & 2.370\\
\cellcolor{gray!6}{\hspace{1em}KEN} & \cellcolor{gray!6}{-0.950} & \cellcolor{gray!6}{-0.073} & \cellcolor{gray!6}{-0.611} & \cellcolor{gray!6}{-0.033} & \cellcolor{gray!6}{-0.680} & \cellcolor{gray!6}{-2.180} & \cellcolor{gray!6}{-0.126} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.145} & \cellcolor{gray!6}{6.550}\\
\hspace{1em}LBR & -1.370 & -0.347 & -- & -0.170 & -1.330 & -2.870 & -- & -0.279 & -0.273 & 1.700\\
\cellcolor{gray!6}{\hspace{1em}MDG} & \cellcolor{gray!6}{-1.450} & \cellcolor{gray!6}{-0.393} & \cellcolor{gray!6}{-0.439} & \cellcolor{gray!6}{-0.155} & \cellcolor{gray!6}{-2.240} & \cellcolor{gray!6}{-1.380} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.033} & \cellcolor{gray!6}{3.320}\\
\hspace{1em}MWI & -0.020 & -0.428 & -0.265 & -0.128 & -0.698 & -0.240 & -0.114 & -0.540 & -0.107 & 13.800\\
\cellcolor{gray!6}{\hspace{1em}MUS} & \cellcolor{gray!6}{-1.490} & \cellcolor{gray!6}{-0.486} & \cellcolor{gray!6}{-0.154} & \cellcolor{gray!6}{-0.128} & \cellcolor{gray!6}{-0.688} & \cellcolor{gray!6}{-3.850} & \cellcolor{gray!6}{0.000} & \cellcolor{gray!6}{-0.083} & \cellcolor{gray!6}{0.001} & \cellcolor{gray!6}{0.870}\\
\hspace{1em}MYT & -0.399 & -1.320 & -- & -- & -0.738 & -1.760 & -0.100 & -- & -0.027 & 1.060\\
\cellcolor{gray!6}{\hspace{1em}NGA} & \cellcolor{gray!6}{0.563} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.056} & \cellcolor{gray!6}{-1.580} & \cellcolor{gray!6}{-1.130} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.255} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.010}\\
\hspace{1em}REU & -1.010 & -0.391 & -- & -0.301 & -0.207 & -3.820 & -- & -0.070 & -- & 2.230\\
\cellcolor{gray!6}{\hspace{1em}RWA} & \cellcolor{gray!6}{-3.080} & \cellcolor{gray!6}{-0.808} & \cellcolor{gray!6}{-0.524} & \cellcolor{gray!6}{-0.139} & \cellcolor{gray!6}{-1.790} & \cellcolor{gray!6}{-6.260} & \cellcolor{gray!6}{-0.024} & \cellcolor{gray!6}{-0.278} & \cellcolor{gray!6}{-0.032} & \cellcolor{gray!6}{3.460}\\
\hspace{1em}SEN & 0.698 & -0.080 & -- & -0.048 & -0.626 & -0.474 & -- & -- & -- & 8.300\\
\cellcolor{gray!6}{\hspace{1em}SLE} & \cellcolor{gray!6}{-0.577} & \cellcolor{gray!6}{-0.231} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.750} & \cellcolor{gray!6}{-1.260} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.116} & \cellcolor{gray!6}{-0.015} & \cellcolor{gray!6}{1.120}\\
\hspace{1em}SSD & 1.150 & -- & -0.283 & -- & -0.284 & -0.533 & -- & -0.066 & -0.125 & 3.790\\
\cellcolor{gray!6}{\hspace{1em}TZA} & \cellcolor{gray!6}{-0.274} & \cellcolor{gray!6}{-0.553} & \cellcolor{gray!6}{-0.253} & \cellcolor{gray!6}{-0.046} & \cellcolor{gray!6}{-0.977} & \cellcolor{gray!6}{-3.090} & \cellcolor{gray!6}{-0.082} & \cellcolor{gray!6}{-0.248} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.320}\\
\hspace{1em}TGO & 1.300 & -0.751 & -0.123 & -0.117 & -0.406 & -0.530 & -- & -0.173 & -- & 2.970\\
\cellcolor{gray!6}{\hspace{1em}UGA} & \cellcolor{gray!6}{-1.650} & \cellcolor{gray!6}{-0.952} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-3.010} & \cellcolor{gray!6}{-2.610} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.079} & \cellcolor{gray!6}{-0.109} & \cellcolor{gray!6}{4.830}\\
\hspace{1em}ZMB & -0.422 & -0.196 & -- & -0.051 & -1.160 & -0.499 & -- & -0.008 & -0.102 & 8.340\\
\addlinespace[0.3em]
\multicolumn{11}{l}{\textbf{Asia}}\\
\cellcolor{gray!6}{\hspace{1em}QLD} & \cellcolor{gray!6}{-1.060} & \cellcolor{gray!6}{-0.428} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.088} & \cellcolor{gray!6}{-0.794} & \cellcolor{gray!6}{-2.600} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.025} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.390}\\
\hspace{1em}BGD & -0.323 & -0.511 & -0.100 & -0.143 & -0.625 & -1.330 & -- & -0.185 & -0.024 & 4.130\\
\cellcolor{gray!6}{\hspace{1em}BTN} & \cellcolor{gray!6}{-1.220} & \cellcolor{gray!6}{-0.009} & \cellcolor{gray!6}{-0.033} & \cellcolor{gray!6}{-0.009} & \cellcolor{gray!6}{-6.340} & \cellcolor{gray!6}{-1.580} & \cellcolor{gray!6}{-0.024} & \cellcolor{gray!6}{-0.066} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{0.821}\\
\hspace{1em}BRN & -0.491 & -1.100 & -0.863 & -0.353 & -2.100 & -0.404 & -- & -0.744 & -0.256 & 16.000\\
\cellcolor{gray!6}{\hspace{1em}KHM} & \cellcolor{gray!6}{-0.413} & \cellcolor{gray!6}{-0.255} & \cellcolor{gray!6}{-1.190} & \cellcolor{gray!6}{-0.217} & \cellcolor{gray!6}{-0.962} & \cellcolor{gray!6}{-0.423} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.419} & \cellcolor{gray!6}{-0.273} & \cellcolor{gray!6}{9.430}\\
\hspace{1em}FJI & -2.030 & -0.381 & -0.454 & -0.125 & -- & -5.280 & -- & -0.199 & -0.052 & 5.340\\
\cellcolor{gray!6}{\hspace{1em}AN} & \cellcolor{gray!6}{-6.390} & \cellcolor{gray!6}{0.029} & \cellcolor{gray!6}{-0.080} & \cellcolor{gray!6}{-0.175} & \cellcolor{gray!6}{-0.096} & \cellcolor{gray!6}{-19.100} & \cellcolor{gray!6}{0.025} & \cellcolor{gray!6}{-0.102} & \cellcolor{gray!6}{0.029} & \cellcolor{gray!6}{4.330}\\
\hspace{1em}NE & -0.939 & -0.163 & -0.321 & -0.219 & -4.860 & -1.690 & -0.030 & -0.159 & -- & 1.730\\
\cellcolor{gray!6}{\hspace{1em}WG} & \cellcolor{gray!6}{-0.041} & \cellcolor{gray!6}{-0.305} & \cellcolor{gray!6}{-0.473} & \cellcolor{gray!6}{-0.036} & \cellcolor{gray!6}{-0.630} & \cellcolor{gray!6}{-2.050} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.163} & \cellcolor{gray!6}{-0.046} & \cellcolor{gray!6}{2.340}\\
\hspace{1em}IDN & -1.290 & -0.738 & -0.378 & -0.541 & -1.160 & -2.090 & -- & -0.389 & -0.115 & 7.260\\
\cellcolor{gray!6}{\hspace{1em}LAO} & \cellcolor{gray!6}{-0.366} & \cellcolor{gray!6}{-0.427} & \cellcolor{gray!6}{-0.376} & \cellcolor{gray!6}{-0.315} & \cellcolor{gray!6}{-1.250} & \cellcolor{gray!6}{-0.860} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.118} & \cellcolor{gray!6}{-0.208} & \cellcolor{gray!6}{2.180}\\
\hspace{1em}MYS & -0.733 & -1.080 & -0.507 & -0.459 & -0.665 & -1.960 & -- & -0.227 & -0.039 & 6.830\\
\cellcolor{gray!6}{\hspace{1em}MMR} & \cellcolor{gray!6}{-0.351} & \cellcolor{gray!6}{-0.210} & \cellcolor{gray!6}{-0.404} & \cellcolor{gray!6}{-0.183} & \cellcolor{gray!6}{-1.020} & \cellcolor{gray!6}{-1.740} & \cellcolor{gray!6}{-0.094} & \cellcolor{gray!6}{-0.166} & \cellcolor{gray!6}{-0.039} & \cellcolor{gray!6}{2.550}\\
\hspace{1em}NCL & -1.580 & -- & -0.374 & -0.080 & -0.415 & -7.730 & -- & -0.122 & -0.035 & 3.740\\
\cellcolor{gray!6}{\hspace{1em}PNG} & \cellcolor{gray!6}{-1.540} & \cellcolor{gray!6}{-0.288} & \cellcolor{gray!6}{-0.582} & \cellcolor{gray!6}{-0.411} & \cellcolor{gray!6}{-0.158} & \cellcolor{gray!6}{-3.830} & \cellcolor{gray!6}{-0.079} & \cellcolor{gray!6}{-0.254} & \cellcolor{gray!6}{-0.144} & \cellcolor{gray!6}{6.340}\\
\hspace{1em}PHL & -1.410 & -0.180 & -0.184 & -0.410 & -1.520 & -2.700 & -- & -0.156 & -0.029 & 2.220\\
\cellcolor{gray!6}{\hspace{1em}SGP} & \cellcolor{gray!6}{0.346} & \cellcolor{gray!6}{-1.230} & \cellcolor{gray!6}{-0.166} & \cellcolor{gray!6}{-0.145} & \cellcolor{gray!6}{-1.140} & \cellcolor{gray!6}{-1.260} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.882} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{10.900}\\
\hspace{1em}SLB & -0.275 & -0.191 & -0.553 & -0.477 & -0.390 & -1.210 & -0.375 & -0.485 & -- & 5.680\\
\cellcolor{gray!6}{\hspace{1em}LKA} & \cellcolor{gray!6}{-0.090} & \cellcolor{gray!6}{-0.307} & \cellcolor{gray!6}{-0.088} & \cellcolor{gray!6}{-0.249} & \cellcolor{gray!6}{-0.828} & \cellcolor{gray!6}{-1.440} & \cellcolor{gray!6}{-0.028} & \cellcolor{gray!6}{-0.108} & \cellcolor{gray!6}{-0.075} & \cellcolor{gray!6}{1.820}\\
\hspace{1em}THA & -1.390 & -0.202 & -0.761 & -0.004 & -1.720 & -3.550 & -0.026 & -0.054 & -- & 1.900\\
\cellcolor{gray!6}{\hspace{1em}TLS} & \cellcolor{gray!6}{-0.474} & \cellcolor{gray!6}{-0.003} & \cellcolor{gray!6}{-0.283} & \cellcolor{gray!6}{-0.047} & \cellcolor{gray!6}{-0.673} & \cellcolor{gray!6}{-1.590} & \cellcolor{gray!6}{-0.076} & \cellcolor{gray!6}{-0.022} & \cellcolor{gray!6}{-0.077} & \cellcolor{gray!6}{1.930}\\
\hspace{1em}VUT & -2.790 & -0.480 & -0.042 & -0.493 & -0.597 & -2.830 & -0.184 & -- & -- & 18.000\\
\cellcolor{gray!6}{\hspace{1em}VNM} & \cellcolor{gray!6}{-0.790} & \cellcolor{gray!6}{-0.525} & \cellcolor{gray!6}{-0.291} & \cellcolor{gray!6}{-0.299} & \cellcolor{gray!6}{-1.630} & \cellcolor{gray!6}{-1.550} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.161} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{3.270}\\*
\end{longtable}
\endgroup{}

(ref:cap-par-region) **Parameter estimate weighted means per region**. We used the forest cover in 2010 to compute the parameter estimate weighted mean per region. Continuous explanatory variables were normalized (mean=0 and standard-deviation=1), allowing us to estimate relative variable importance in determining the spatial probability of deforestation from parameter values. Comparison can be done within and between regions. At the pantropical scale (when considering all continents together), explanatory variables can be classified in the following decreasing order of importance: distance to forest edge, distance to past deforestation, presence of a protected area, distance to nearest road, distance to nearest town, altitude, slope, and distance to nearest river. \vspace{0.5cm}

\begingroup\fontsize{10}{12}\selectfont

\begin{longtable}[t]{lrrrrrrrrrr}
\caption{(\#tab:par-region)(ref:cap-par-region)}\\
\toprule
Region & int & pa & elev & slope & ddefor & dedge & driver & droad & dtown & Vrho\\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
Region & int & pa & elev & slope & ddefor & dedge & driver & droad & dtown & Vrho\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\cellcolor{gray!6}{India} & \cellcolor{gray!6}{-1.011} & \cellcolor{gray!6}{-0.189} & \cellcolor{gray!6}{-0.347} & \cellcolor{gray!6}{-0.169} & \cellcolor{gray!6}{-3.499} & \cellcolor{gray!6}{-2.756} & \cellcolor{gray!6}{-0.019} & \cellcolor{gray!6}{-0.157} & \cellcolor{gray!6}{-0.010} & \cellcolor{gray!6}{2.033}\\
Brazil & -2.578 & -0.663 & -0.082 & -0.031 & -2.088 & -3.567 & -0.083 & -0.486 & -0.424 & 9.168\\
\cellcolor{gray!6}{America} & \cellcolor{gray!6}{-2.502} & \cellcolor{gray!6}{-0.567} & \cellcolor{gray!6}{-0.266} & \cellcolor{gray!6}{-0.123} & \cellcolor{gray!6}{-2.154} & \cellcolor{gray!6}{-4.028} & \cellcolor{gray!6}{-0.054} & \cellcolor{gray!6}{-0.495} & \cellcolor{gray!6}{-0.367} & \cellcolor{gray!6}{7.734}\\
Africa & -2.554 & -0.284 & -0.085 & -0.067 & -3.241 & -4.767 & -0.014 & -0.392 & -0.191 & 5.453\\
\cellcolor{gray!6}{Asia} & \cellcolor{gray!6}{-1.136} & \cellcolor{gray!6}{-0.554} & \cellcolor{gray!6}{-0.419} & \cellcolor{gray!6}{-0.420} & \cellcolor{gray!6}{-1.129} & \cellcolor{gray!6}{-2.326} & \cellcolor{gray!6}{-0.024} & \cellcolor{gray!6}{-0.285} & \cellcolor{gray!6}{-0.094} & \cellcolor{gray!6}{5.797}\\
All continents & -2.195 & -0.503 & -0.262 & -0.180 & -2.151 & -3.792 & -0.039 & -0.424 & -0.265 & 6.788\\*
\end{longtable}
\endgroup{}

\newpage

<!--------------------------------->
<!-- Back-transformed parameters -->
<!--------------------------------->

## Back-transformed parameters

(ref:cap-bt-par) **Back-transformed parameters for each study area**. We back-transformed the parameters using the mean and standard-deviation of each continuous variable for each study-area. Doing so, we can use Eq. \@ref(eq:icar) to compute the change in the probability of deforestation associated to a particular change in the explanatory variables, in their original units. To use this table of parameters, distances and elevation must be expressed in kilometers (Km), and slope must be expressed in hecto-degrees ($10^2$°). Note that the intercept is affected by the back-transformation but that the effect associated to protected areas ("pa") and the variance of the spatial random effects ("Vrho") are left unchanged.\vspace{0.5cm}

\begingroup\fontsize{10}{12}\selectfont

\begin{longtable}[t]{lrrrrrrrrrr}
\caption{(\#tab:bt-par)(ref:cap-bt-par)}\\
\toprule
\multicolumn{1}{l}{Study-area} & \multicolumn{1}{r}{int} & \multicolumn{1}{r}{pa} & \multicolumn{1}{r}{elev} & \multicolumn{1}{r}{slope} & \multicolumn{1}{r}{ddefor} & \multicolumn{1}{r}{dedge} & \multicolumn{1}{r}{driver} & \multicolumn{1}{r}{droad} & \multicolumn{1}{r}{dtown} & \multicolumn{1}{r}{Vrho} \\
 &  &  & (Km) & ($10^2$°) & (Km) & (Km) & (Km) & (Km) & (Km) & \\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
\multicolumn{1}{l}{Study-area} & \multicolumn{1}{r}{int} & \multicolumn{1}{r}{pa} & \multicolumn{1}{r}{elev} & \multicolumn{1}{r}{slope} & \multicolumn{1}{r}{ddefor} & \multicolumn{1}{r}{dedge} & \multicolumn{1}{r}{driver} & \multicolumn{1}{r}{droad} & \multicolumn{1}{r}{dtown} & \multicolumn{1}{r}{Vrho} \\
 &  &  & (Km) & ($10^2$°) & (Km) & (Km) & (Km) & (Km) & (Km) & \\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{11}{l}{\textbf{America}}\\
\cellcolor{gray!6}{\hspace{1em}ATG} & \cellcolor{gray!6}{2.678} & \cellcolor{gray!6}{-0.259} & \cellcolor{gray!6}{-2.607} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-3.604} & \cellcolor{gray!6}{-14.700} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{10.000}\\
\hspace{1em}BHS & 0.539 & -- & -- & -5.067 & -6.109 & -6.728 & -- & -- & -0.003 & 6.930\\
\cellcolor{gray!6}{\hspace{1em}BRB} & \cellcolor{gray!6}{2.449} & \cellcolor{gray!6}{-0.277} & \cellcolor{gray!6}{-5.647} & \cellcolor{gray!6}{-1.681} & \cellcolor{gray!6}{-2.168} & \cellcolor{gray!6}{-24.767} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{0.727}\\
\hspace{1em}BLZ & 1.372 & -0.828 & -0.093 & -5.875 & -1.290 & -2.356 & -- & -0.027 & -0.033 & 6.640\\
\cellcolor{gray!6}{\hspace{1em}BOL} & \cellcolor{gray!6}{2.306} & \cellcolor{gray!6}{-0.204} & \cellcolor{gray!6}{-0.593} & \cellcolor{gray!6}{-1.629} & \cellcolor{gray!6}{-0.608} & \cellcolor{gray!6}{-3.887} & \cellcolor{gray!6}{-0.005} & \cellcolor{gray!6}{-0.017} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.080}\\
\hspace{1em}AC & 2.534 & -0.723 & -- & -0.373 & -2.089 & -3.516 & -0.001 & -- & -0.012 & 3.930\\
\cellcolor{gray!6}{\hspace{1em}AL} & \cellcolor{gray!6}{3.973} & \cellcolor{gray!6}{-0.248} & \cellcolor{gray!6}{-2.488} & \cellcolor{gray!6}{-0.168} & \cellcolor{gray!6}{-3.648} & \cellcolor{gray!6}{-25.258} & \cellcolor{gray!6}{-0.039} & \cellcolor{gray!6}{-0.047} & \cellcolor{gray!6}{-0.057} & \cellcolor{gray!6}{3.510}\\
\hspace{1em}AP & 2.753 & -0.157 & -6.392 & -- & -0.654 & -3.667 & -- & -0.007 & -0.007 & 6.290\\
\cellcolor{gray!6}{\hspace{1em}AM} & \cellcolor{gray!6}{2.435} & \cellcolor{gray!6}{-0.558} & \cellcolor{gray!6}{-0.980} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.640} & \cellcolor{gray!6}{-1.312} & \cellcolor{gray!6}{-0.023} & \cellcolor{gray!6}{-0.010} & \cellcolor{gray!6}{-0.009} & \cellcolor{gray!6}{12.400}\\
\hspace{1em}BA & 3.310 & -0.171 & -1.466 & -2.797 & -2.939 & -8.872 & -- & -0.017 & -0.002 & 3.710\\
\cellcolor{gray!6}{\hspace{1em}CE} & \cellcolor{gray!6}{4.207} & \cellcolor{gray!6}{-0.654} & \cellcolor{gray!6}{-2.708} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.191} & \cellcolor{gray!6}{-28.850} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{12.600}\\
\hspace{1em}ES & 0.929 & -0.056 & -2.378 & -0.074 & -2.769 & -7.403 & -0.004 & -0.047 & -0.002 & 2.780\\
\cellcolor{gray!6}{\hspace{1em}GO} & \cellcolor{gray!6}{1.167} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-2.364} & \cellcolor{gray!6}{-11.867} & \cellcolor{gray!6}{-0.018} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.002} & \cellcolor{gray!6}{3.920}\\
\hspace{1em}MA & 1.336 & -0.056 & -- & -2.749 & -3.139 & -4.937 & -- & -- & -0.014 & 3.270\\
\cellcolor{gray!6}{\hspace{1em}MT} & \cellcolor{gray!6}{1.416} & \cellcolor{gray!6}{-0.477} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.788} & \cellcolor{gray!6}{-1.027} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.012} & \cellcolor{gray!6}{-0.012} & \cellcolor{gray!6}{9.680}\\
\hspace{1em}MS & 1.456 & -- & -1.422 & -- & -0.773 & -18.365 & -- & -- & -- & 3.430\\
\cellcolor{gray!6}{\hspace{1em}MG} & \cellcolor{gray!6}{2.541} & \cellcolor{gray!6}{-0.113} & \cellcolor{gray!6}{-1.413} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-3.584} & \cellcolor{gray!6}{-9.224} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.021} & \cellcolor{gray!6}{-0.006} & \cellcolor{gray!6}{3.730}\\
\hspace{1em}PA & 1.773 & -0.992 & -- & -1.104 & -0.834 & -1.614 & -- & -0.012 & -0.007 & 7.140\\
\cellcolor{gray!6}{\hspace{1em}PB} & \cellcolor{gray!6}{4.422} & \cellcolor{gray!6}{-0.453} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-3.672} & \cellcolor{gray!6}{-11.859} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.061} & \cellcolor{gray!6}{4.060}\\
\hspace{1em}PR & 1.264 & -0.206 & -- & -4.249 & -2.332 & -9.678 & -0.007 & -- & -- & 2.910\\
\cellcolor{gray!6}{\hspace{1em}PE} & \cellcolor{gray!6}{7.310} & \cellcolor{gray!6}{-0.405} & \cellcolor{gray!6}{-2.041} & \cellcolor{gray!6}{-1.616} & \cellcolor{gray!6}{-4.883} & \cellcolor{gray!6}{-31.274} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.071} & \cellcolor{gray!6}{-0.032} & \cellcolor{gray!6}{3.350}\\
\hspace{1em}PI & 1.538 & -0.144 & -- & -- & -2.239 & -23.456 & -0.006 & -0.008 & 0.000 & 3.460\\
\cellcolor{gray!6}{\hspace{1em}RJ} & \cellcolor{gray!6}{2.740} & \cellcolor{gray!6}{-0.205} & \cellcolor{gray!6}{-0.432} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-4.532} & \cellcolor{gray!6}{-17.406} & \cellcolor{gray!6}{-0.019} & \cellcolor{gray!6}{-0.038} & \cellcolor{gray!6}{-0.002} & \cellcolor{gray!6}{2.400}\\
\hspace{1em}RN & 2.766 & -0.301 & -- & -- & -6.202 & -8.393 & -- & -0.058 & -- & 17.100\\
\cellcolor{gray!6}{\hspace{1em}RS} & \cellcolor{gray!6}{2.443} & \cellcolor{gray!6}{-0.143} & \cellcolor{gray!6}{-0.289} & \cellcolor{gray!6}{-7.504} & \cellcolor{gray!6}{-3.560} & \cellcolor{gray!6}{-21.518} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{1.730}\\
\hspace{1em}RO & 1.480 & -1.340 & -- & -3.024 & -0.498 & -1.686 & -- & -0.013 & -0.016 & 7.120\\
\cellcolor{gray!6}{\hspace{1em}RR} & \cellcolor{gray!6}{2.461} & \cellcolor{gray!6}{-0.670} & \cellcolor{gray!6}{-2.015} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.632} & \cellcolor{gray!6}{-1.362} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.021} & \cellcolor{gray!6}{7.840}\\
\hspace{1em}SC & 1.467 & -0.232 & -- & -7.218 & -3.370 & -10.375 & -- & -- & -- & 2.100\\
\cellcolor{gray!6}{\hspace{1em}SP} & \cellcolor{gray!6}{2.162} & \cellcolor{gray!6}{-0.185} & \cellcolor{gray!6}{-0.510} & \cellcolor{gray!6}{-3.381} & \cellcolor{gray!6}{-2.229} & \cellcolor{gray!6}{-7.245} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.024} & \cellcolor{gray!6}{-0.002} & \cellcolor{gray!6}{3.960}\\
\hspace{1em}SE & 2.556 & -- & -- & -- & -6.064 & -14.803 & -0.017 & -- & -- & 3.840\\
\cellcolor{gray!6}{\hspace{1em}TO} & \cellcolor{gray!6}{0.272} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.324} & \cellcolor{gray!6}{-1.193} & \cellcolor{gray!6}{-0.015} & \cellcolor{gray!6}{0.001} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{4.630}\\
\hspace{1em}COL & 2.185 & -0.535 & -0.915 & -3.401 & -1.121 & -1.217 & -- & -0.011 & -0.016 & 5.320\\
\cellcolor{gray!6}{\hspace{1em}CRI} & \cellcolor{gray!6}{1.848} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.042} & \cellcolor{gray!6}{-4.093} & \cellcolor{gray!6}{-2.599} & \cellcolor{gray!6}{-8.827} & \cellcolor{gray!6}{-0.017} & \cellcolor{gray!6}{-0.009} & \cellcolor{gray!6}{-0.012} & \cellcolor{gray!6}{2.090}\\
\hspace{1em}CUB & 1.030 & -0.077 & -0.499 & -0.801 & -2.248 & -5.478 & -- & -0.015 & -- & 5.700\\
\cellcolor{gray!6}{\hspace{1em}DMA} & \cellcolor{gray!6}{2.164} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.146} & \cellcolor{gray!6}{-2.607} & \cellcolor{gray!6}{-0.787} & \cellcolor{gray!6}{-4.124} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.130} & \cellcolor{gray!6}{2.160}\\
\hspace{1em}DOM & 1.846 & -0.357 & -0.714 & -0.351 & -4.776 & -9.541 & -0.005 & -- & -- & 2.670\\
\cellcolor{gray!6}{\hspace{1em}ECU} & \cellcolor{gray!6}{2.617} & \cellcolor{gray!6}{-0.867} & \cellcolor{gray!6}{-0.134} & \cellcolor{gray!6}{-3.742} & \cellcolor{gray!6}{-0.448} & \cellcolor{gray!6}{-5.973} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.022} & \cellcolor{gray!6}{-0.025} & \cellcolor{gray!6}{2.470}\\
\hspace{1em}SLV & 2.659 & -0.426 & -0.609 & -2.064 & -3.650 & -18.572 & -0.004 & -0.047 & -0.041 & 5.860\\
\cellcolor{gray!6}{\hspace{1em}GUF} & \cellcolor{gray!6}{3.797} & \cellcolor{gray!6}{-1.070} & \cellcolor{gray!6}{-8.284} & \cellcolor{gray!6}{-0.718} & \cellcolor{gray!6}{-0.418} & \cellcolor{gray!6}{-1.431} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.011} & \cellcolor{gray!6}{-0.024} & \cellcolor{gray!6}{21.400}\\
\hspace{1em}GRD & 3.447 & -0.287 & -6.263 & -- & -3.760 & -11.809 & -- & -0.231 & -- & 5.780\\
\cellcolor{gray!6}{\hspace{1em}GLP} & \cellcolor{gray!6}{2.080} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-2.164} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.696} & \cellcolor{gray!6}{-16.219} & \cellcolor{gray!6}{-0.006} & \cellcolor{gray!6}{-0.026} & \cellcolor{gray!6}{-0.040} & \cellcolor{gray!6}{3.110}\\
\hspace{1em}GTM & 1.756 & -0.172 & -0.568 & -2.559 & -2.151 & -3.959 & -0.018 & -0.015 & 0.000 & 2.700\\
\cellcolor{gray!6}{\hspace{1em}GUY} & \cellcolor{gray!6}{2.666} & \cellcolor{gray!6}{-0.475} & \cellcolor{gray!6}{-2.425} & \cellcolor{gray!6}{-3.873} & \cellcolor{gray!6}{-0.239} & \cellcolor{gray!6}{-1.779} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.005} & \cellcolor{gray!6}{-0.018} & \cellcolor{gray!6}{9.980}\\
\hspace{1em}HTI & 1.682 & -0.152 & -- & -1.315 & -5.259 & -18.411 & -0.052 & -0.037 & -0.012 & 2.840\\
\cellcolor{gray!6}{\hspace{1em}HND} & \cellcolor{gray!6}{1.063} & \cellcolor{gray!6}{-0.403} & \cellcolor{gray!6}{-1.493} & \cellcolor{gray!6}{-3.362} & \cellcolor{gray!6}{-0.394} & \cellcolor{gray!6}{-1.254} & \cellcolor{gray!6}{0.000} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{3.580}\\
\hspace{1em}JAM & 2.095 & -0.043 & -1.209 & -1.189 & -4.395 & -11.087 & 0.000 & -0.032 & -0.049 & 2.570\\
\cellcolor{gray!6}{\hspace{1em}MTQ} & \cellcolor{gray!6}{2.591} & \cellcolor{gray!6}{-0.124} & \cellcolor{gray!6}{-2.206} & \cellcolor{gray!6}{-3.127} & \cellcolor{gray!6}{-1.739} & \cellcolor{gray!6}{-20.515} & \cellcolor{gray!6}{-0.010} & \cellcolor{gray!6}{-0.023} & \cellcolor{gray!6}{-0.004} & \cellcolor{gray!6}{1.850}\\
\hspace{1em}MEX & 0.853 & -0.328 & -0.218 & -2.562 & -2.096 & -5.676 & -- & -0.007 & -- & 3.600\\
\cellcolor{gray!6}{\hspace{1em}MSR} & \cellcolor{gray!6}{2.005} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.525} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-2.812} & \cellcolor{gray!6}{-32.264} & \cellcolor{gray!6}{-0.274} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{7.120}\\
\hspace{1em}NIC & 0.245 & -- & -1.280 & -1.680 & -0.270 & -1.158 & -- & -0.011 & -- & 3.450\\
\cellcolor{gray!6}{\hspace{1em}PAN} & \cellcolor{gray!6}{1.897} & \cellcolor{gray!6}{-0.344} & \cellcolor{gray!6}{-0.364} & \cellcolor{gray!6}{-3.031} & \cellcolor{gray!6}{-1.330} & \cellcolor{gray!6}{-3.680} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.030} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{3.170}\\
\hspace{1em}PRY & 0.952 & -0.167 & -1.348 & -6.938 & -2.918 & -3.210 & -- & -0.010 & -- & 3.460\\
\cellcolor{gray!6}{\hspace{1em}PER} & \cellcolor{gray!6}{3.557} & \cellcolor{gray!6}{-0.684} & \cellcolor{gray!6}{-0.796} & \cellcolor{gray!6}{-3.810} & \cellcolor{gray!6}{-0.847} & \cellcolor{gray!6}{-1.734} & \cellcolor{gray!6}{-0.001} & \cellcolor{gray!6}{-0.012} & \cellcolor{gray!6}{-0.024} & \cellcolor{gray!6}{4.980}\\
\hspace{1em}PRI & 1.481 & -- & -1.108 & -1.321 & -2.556 & -4.709 & -- & -- & -0.022 & 2.470\\
\cellcolor{gray!6}{\hspace{1em}KNA} & \cellcolor{gray!6}{1.927} & \cellcolor{gray!6}{-0.137} & \cellcolor{gray!6}{-4.480} & \cellcolor{gray!6}{-5.042} & \cellcolor{gray!6}{-0.048} & \cellcolor{gray!6}{-3.941} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.780}\\
\hspace{1em}LCA & 2.185 & -- & -2.269 & -2.127 & -0.953 & -9.779 & -- & -0.033 & -0.064 & 3.900\\
\cellcolor{gray!6}{\hspace{1em}MAF} & \cellcolor{gray!6}{1.384} & \cellcolor{gray!6}{-1.180} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-4.541} & \cellcolor{gray!6}{-2.964} & \cellcolor{gray!6}{-12.450} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{7.700}\\
\hspace{1em}VCT & 3.086 & -- & -3.736 & -1.947 & -1.691 & -4.476 & -- & -0.060 & -- & 2.320\\
\cellcolor{gray!6}{\hspace{1em}SXM} & \cellcolor{gray!6}{6.198} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-8.470} & \cellcolor{gray!6}{-6.827} & \cellcolor{gray!6}{-12.379} & \cellcolor{gray!6}{-1.219} & \cellcolor{gray!6}{-0.808} & \cellcolor{gray!6}{-0.660} & \cellcolor{gray!6}{30.400}\\
\hspace{1em}SUR & 1.737 & -0.171 & -5.316 & -- & -0.355 & -0.720 & -- & -0.008 & -0.018 & 14.200\\
\cellcolor{gray!6}{\hspace{1em}TTO} & \cellcolor{gray!6}{2.060} & \cellcolor{gray!6}{-0.292} & \cellcolor{gray!6}{-0.761} & \cellcolor{gray!6}{-0.507} & \cellcolor{gray!6}{-2.880} & \cellcolor{gray!6}{-8.056} & \cellcolor{gray!6}{0.003} & \cellcolor{gray!6}{-0.036} & \cellcolor{gray!6}{-0.062} & \cellcolor{gray!6}{1.790}\\
\hspace{1em}VEN & 3.050 & -0.120 & -0.823 & -1.359 & -0.737 & -4.136 & -0.010 & -0.005 & -0.014 & 5.860\\
\cellcolor{gray!6}{\hspace{1em}VGB} & \cellcolor{gray!6}{3.121} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-3.312} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-7.725} & \cellcolor{gray!6}{-13.468} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{12.900}\\
\hspace{1em}VIR & 0.888 & -0.037 & -- & -1.927 & -0.613 & -5.077 & -- & -0.067 & -- & 6.460\\
\addlinespace[0.3em]
\multicolumn{11}{l}{\textbf{Africa}}\\
\cellcolor{gray!6}{\hspace{1em}AGO} & \cellcolor{gray!6}{2.429} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.207} & \cellcolor{gray!6}{-1.247} & \cellcolor{gray!6}{-1.490} & \cellcolor{gray!6}{-9.294} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.020} & \cellcolor{gray!6}{-0.008} & \cellcolor{gray!6}{3.780}\\
\hspace{1em}BEN & 3.855 & -0.218 & -- & -13.778 & -1.267 & -29.995 & 0.001 & -- & -- & 6.730\\
\cellcolor{gray!6}{\hspace{1em}BDI} & \cellcolor{gray!6}{5.004} & \cellcolor{gray!6}{-1.860} & \cellcolor{gray!6}{-1.228} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-5.560} & \cellcolor{gray!6}{-15.233} & \cellcolor{gray!6}{-0.043} & \cellcolor{gray!6}{-0.047} & \cellcolor{gray!6}{-0.017} & \cellcolor{gray!6}{5.060}\\
\hspace{1em}CMR & 2.795 & -0.933 & -0.122 & -3.696 & -0.480 & -1.869 & -0.005 & -0.046 & -0.034 & 5.740\\
\cellcolor{gray!6}{\hspace{1em}CAF} & \cellcolor{gray!6}{2.968} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.038} & \cellcolor{gray!6}{-3.281} & \cellcolor{gray!6}{-3.153} & \cellcolor{gray!6}{-0.011} & \cellcolor{gray!6}{-0.005} & \cellcolor{gray!6}{-0.020} & \cellcolor{gray!6}{4.270}\\
\hspace{1em}COM & 2.147 & -0.038 & -- & -3.728 & -- & -32.087 & -- & -- & -- & 6.220\\
\cellcolor{gray!6}{\hspace{1em}COG} & \cellcolor{gray!6}{2.560} & \cellcolor{gray!6}{-0.112} & \cellcolor{gray!6}{-0.626} & \cellcolor{gray!6}{-0.607} & \cellcolor{gray!6}{-0.174} & \cellcolor{gray!6}{-1.555} & \cellcolor{gray!6}{-0.001} & \cellcolor{gray!6}{-0.029} & \cellcolor{gray!6}{-0.018} & \cellcolor{gray!6}{9.550}\\
\hspace{1em}COD & 2.244 & -0.261 & -- & -- & -1.606 & -2.387 & -- & -0.023 & -0.012 & 3.740\\
\cellcolor{gray!6}{\hspace{1em}GNQ} & \cellcolor{gray!6}{2.574} & \cellcolor{gray!6}{-0.368} & \cellcolor{gray!6}{-1.407} & \cellcolor{gray!6}{-7.729} & \cellcolor{gray!6}{-0.046} & \cellcolor{gray!6}{-1.728} & \cellcolor{gray!6}{-0.016} & \cellcolor{gray!6}{-0.284} & \cellcolor{gray!6}{-0.014} & \cellcolor{gray!6}{10.600}\\
\hspace{1em}ETH & 2.227 & -0.241 & -0.274 & -2.810 & -2.028 & -6.497 & -- & -0.008 & -0.016 & 3.280\\
\cellcolor{gray!6}{\hspace{1em}GAB} & \cellcolor{gray!6}{2.894} & \cellcolor{gray!6}{-0.206} & \cellcolor{gray!6}{-2.339} & \cellcolor{gray!6}{-6.273} & \cellcolor{gray!6}{-0.177} & \cellcolor{gray!6}{-1.807} & \cellcolor{gray!6}{-0.003} & \cellcolor{gray!6}{-0.034} & \cellcolor{gray!6}{-0.008} & \cellcolor{gray!6}{14.800}\\
\hspace{1em}GMB & 2.754 & -- & -- & -12.034 & -4.609 & -14.050 & -- & -0.145 & -- & 3.980\\
\cellcolor{gray!6}{\hspace{1em}GHA} & \cellcolor{gray!6}{2.187} & \cellcolor{gray!6}{-0.174} & \cellcolor{gray!6}{-1.712} & \cellcolor{gray!6}{-3.660} & \cellcolor{gray!6}{-1.006} & \cellcolor{gray!6}{-3.892} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.037} & \cellcolor{gray!6}{-0.019} & \cellcolor{gray!6}{2.110}\\
\hspace{1em}GIN & 1.664 & -0.272 & -0.108 & -- & -6.581 & -5.747 & -- & -0.015 & -- & 2.140\\
\cellcolor{gray!6}{\hspace{1em}GNB} & \cellcolor{gray!6}{2.807} & \cellcolor{gray!6}{-0.502} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-6.527} & \cellcolor{gray!6}{-11.743} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.004} & \cellcolor{gray!6}{-0.018} & \cellcolor{gray!6}{5.200}\\
\hspace{1em}CIV & 0.937 & -- & -- & -4.434 & -0.983 & -0.902 & -- & -- & -0.002 & 2.370\\
\cellcolor{gray!6}{\hspace{1em}KEN} & \cellcolor{gray!6}{2.859} & \cellcolor{gray!6}{-0.073} & \cellcolor{gray!6}{-0.763} & \cellcolor{gray!6}{-0.459} & \cellcolor{gray!6}{-1.341} & \cellcolor{gray!6}{-7.598} & \cellcolor{gray!6}{-0.019} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.029} & \cellcolor{gray!6}{6.550}\\
\hspace{1em}LBR & 1.665 & -0.347 & -- & -5.328 & -0.857 & -3.130 & -- & -0.034 & -0.035 & 1.700\\
\cellcolor{gray!6}{\hspace{1em}MDG} & \cellcolor{gray!6}{1.390} & \cellcolor{gray!6}{-0.393} & \cellcolor{gray!6}{-0.971} & \cellcolor{gray!6}{-2.186} & \cellcolor{gray!6}{-3.595} & \cellcolor{gray!6}{-2.599} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.005} & \cellcolor{gray!6}{3.320}\\
\hspace{1em}MWI & 3.021 & -0.428 & -0.680 & -2.018 & -4.069 & -3.727 & -0.023 & -0.112 & -0.027 & 13.800\\
\cellcolor{gray!6}{\hspace{1em}MUS} & \cellcolor{gray!6}{2.095} & \cellcolor{gray!6}{-0.486} & \cellcolor{gray!6}{-0.824} & \cellcolor{gray!6}{-1.682} & \cellcolor{gray!6}{-2.184} & \cellcolor{gray!6}{-18.004} & \cellcolor{gray!6}{0.000} & \cellcolor{gray!6}{-0.091} & \cellcolor{gray!6}{0.001} & \cellcolor{gray!6}{0.870}\\
\hspace{1em}MYT & 2.137 & -1.320 & -- & -- & -3.170 & -20.554 & -0.025 & -- & -0.032 & 1.060\\
\cellcolor{gray!6}{\hspace{1em}NGA} & \cellcolor{gray!6}{2.205} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.310} & \cellcolor{gray!6}{-1.711} & \cellcolor{gray!6}{-2.168} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.027} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.010}\\
\hspace{1em}REU & 1.972 & -0.391 & -- & -2.883 & -0.103 & -14.670 & -- & -0.041 & -- & 2.230\\
\cellcolor{gray!6}{\hspace{1em}RWA} & \cellcolor{gray!6}{5.019} & \cellcolor{gray!6}{-0.808} & \cellcolor{gray!6}{-1.142} & \cellcolor{gray!6}{-1.719} & \cellcolor{gray!6}{-3.851} & \cellcolor{gray!6}{-17.521} & \cellcolor{gray!6}{-0.013} & \cellcolor{gray!6}{-0.070} & \cellcolor{gray!6}{-0.004} & \cellcolor{gray!6}{3.460}\\
\hspace{1em}SEN & 1.841 & -0.080 & -- & -8.912 & -1.856 & -8.943 & -- & -- & -- & 8.300\\
\cellcolor{gray!6}{\hspace{1em}SLE} & \cellcolor{gray!6}{1.068} & \cellcolor{gray!6}{-0.231} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-5.263} & \cellcolor{gray!6}{-4.807} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.023} & \cellcolor{gray!6}{-0.004} & \cellcolor{gray!6}{1.120}\\
\hspace{1em}SSD & 2.945 & -- & -1.557 & -- & -1.250 & -5.962 & -- & -0.004 & -0.012 & 3.790\\
\cellcolor{gray!6}{\hspace{1em}TZA} & \cellcolor{gray!6}{2.919} & \cellcolor{gray!6}{-0.553} & \cellcolor{gray!6}{-0.338} & \cellcolor{gray!6}{-0.545} & \cellcolor{gray!6}{-2.066} & \cellcolor{gray!6}{-13.659} & \cellcolor{gray!6}{-0.008} & \cellcolor{gray!6}{-0.039} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.320}\\
\hspace{1em}TGO & 3.174 & -0.751 & -0.639 & -1.643 & -5.875 & -18.033 & -- & -0.041 & -- & 2.970\\
\cellcolor{gray!6}{\hspace{1em}UGA} & \cellcolor{gray!6}{1.718} & \cellcolor{gray!6}{-0.952} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-3.655} & \cellcolor{gray!6}{-4.427} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.022} & \cellcolor{gray!6}{-0.033} & \cellcolor{gray!6}{4.830}\\
\hspace{1em}ZMB & 0.765 & -0.196 & -- & -2.309 & -1.720 & -15.636 & -- & -0.001 & -0.005 & 8.340\\
\addlinespace[0.3em]
\multicolumn{11}{l}{\textbf{Asia}}\\
\cellcolor{gray!6}{\hspace{1em}QLD} & \cellcolor{gray!6}{1.130} & \cellcolor{gray!6}{-0.428} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-1.132} & \cellcolor{gray!6}{-1.500} & \cellcolor{gray!6}{-5.703} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.002} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{5.390}\\
\hspace{1em}BGD & 1.378 & -0.511 & -0.943 & -2.891 & -1.654 & -9.399 & -- & -0.014 & -0.002 & 4.130\\
\cellcolor{gray!6}{\hspace{1em}BTN} & \cellcolor{gray!6}{1.832} & \cellcolor{gray!6}{-0.009} & \cellcolor{gray!6}{-0.029} & \cellcolor{gray!6}{-0.101} & \cellcolor{gray!6}{-5.078} & \cellcolor{gray!6}{-13.930} & \cellcolor{gray!6}{-0.008} & \cellcolor{gray!6}{-0.008} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{0.821}\\
\hspace{1em}BRN & 2.819 & -1.100 & -6.585 & -7.202 & -1.620 & -0.333 & -- & -0.106 & -0.021 & 16.000\\
\cellcolor{gray!6}{\hspace{1em}KHM} & \cellcolor{gray!6}{3.099} & \cellcolor{gray!6}{-0.255} & \cellcolor{gray!6}{-4.695} & \cellcolor{gray!6}{-3.603} & \cellcolor{gray!6}{-1.656} & \cellcolor{gray!6}{-0.905} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.043} & \cellcolor{gray!6}{-0.031} & \cellcolor{gray!6}{9.430}\\
\hspace{1em}FJI & 1.718 & -0.381 & -2.130 & -1.898 & -- & -6.085 & -- & -0.007 & -0.002 & 5.340\\
\cellcolor{gray!6}{\hspace{1em}AN} & \cellcolor{gray!6}{1.938} & \cellcolor{gray!6}{0.029} & \cellcolor{gray!6}{-1.344} & \cellcolor{gray!6}{-3.915} & \cellcolor{gray!6}{-0.026} & \cellcolor{gray!6}{-17.004} & \cellcolor{gray!6}{0.002} & \cellcolor{gray!6}{-0.009} & \cellcolor{gray!6}{0.002} & \cellcolor{gray!6}{4.330}\\
\hspace{1em}NE & 2.148 & -0.163 & -0.324 & -2.074 & -4.556 & -8.213 & -0.008 & -0.021 & -- & 1.730\\
\cellcolor{gray!6}{\hspace{1em}WG} & \cellcolor{gray!6}{2.669} & \cellcolor{gray!6}{-0.305} & \cellcolor{gray!6}{-1.007} & \cellcolor{gray!6}{-0.423} & \cellcolor{gray!6}{-3.975} & \cellcolor{gray!6}{-19.336} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.071} & \cellcolor{gray!6}{-0.014} & \cellcolor{gray!6}{2.340}\\
\hspace{1em}IDN & 1.507 & -0.738 & -0.792 & -6.830 & -0.604 & -1.836 & -- & -0.022 & -0.016 & 7.260\\
\cellcolor{gray!6}{\hspace{1em}LAO} & \cellcolor{gray!6}{2.608} & \cellcolor{gray!6}{-0.427} & \cellcolor{gray!6}{-1.020} & \cellcolor{gray!6}{-3.801} & \cellcolor{gray!6}{-4.575} & \cellcolor{gray!6}{-4.047} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.015} & \cellcolor{gray!6}{-0.034} & \cellcolor{gray!6}{2.180}\\
\hspace{1em}MYS & 1.727 & -1.080 & -1.518 & -6.167 & -0.777 & -3.046 & -- & -0.021 & -0.003 & 6.830\\
\cellcolor{gray!6}{\hspace{1em}MMR} & \cellcolor{gray!6}{2.306} & \cellcolor{gray!6}{-0.210} & \cellcolor{gray!6}{-0.622} & \cellcolor{gray!6}{-2.001} & \cellcolor{gray!6}{-3.163} & \cellcolor{gray!6}{-6.421} & \cellcolor{gray!6}{-0.018} & \cellcolor{gray!6}{-0.011} & \cellcolor{gray!6}{-0.003} & \cellcolor{gray!6}{2.550}\\
\hspace{1em}NCL & 3.180 & -- & -1.483 & -0.868 & -0.385 & -23.014 & -- & -0.019 & -0.004 & 3.740\\
\cellcolor{gray!6}{\hspace{1em}PNG} & \cellcolor{gray!6}{1.912} & \cellcolor{gray!6}{-0.288} & \cellcolor{gray!6}{-0.745} & \cellcolor{gray!6}{-4.497} & \cellcolor{gray!6}{-0.052} & \cellcolor{gray!6}{-2.906} & \cellcolor{gray!6}{-0.010} & \cellcolor{gray!6}{-0.006} & \cellcolor{gray!6}{-0.006} & \cellcolor{gray!6}{6.340}\\
\hspace{1em}PHL & 1.627 & -0.180 & -0.463 & -4.875 & -2.494 & -6.668 & -- & -0.036 & -0.006 & 2.220\\
\cellcolor{gray!6}{\hspace{1em}SGP} & \cellcolor{gray!6}{3.579} & \cellcolor{gray!6}{-1.230} & \cellcolor{gray!6}{-9.055} & \cellcolor{gray!6}{-6.720} & \cellcolor{gray!6}{-5.960} & \cellcolor{gray!6}{-18.039} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.894} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{10.900}\\
\hspace{1em}SLB & 2.104 & -0.191 & -2.584 & -6.107 & -0.077 & -0.988 & -0.015 & -0.006 & -- & 5.680\\
\cellcolor{gray!6}{\hspace{1em}LKA} & \cellcolor{gray!6}{1.889} & \cellcolor{gray!6}{-0.307} & \cellcolor{gray!6}{-0.276} & \cellcolor{gray!6}{-3.465} & \cellcolor{gray!6}{-4.531} & \cellcolor{gray!6}{-10.719} & \cellcolor{gray!6}{-0.009} & \cellcolor{gray!6}{-0.032} & \cellcolor{gray!6}{-0.018} & \cellcolor{gray!6}{1.820}\\
\hspace{1em}THA & 2.528 & -0.202 & -2.052 & -0.051 & -4.471 & -12.968 & -0.004 & -0.007 & -- & 1.900\\
\cellcolor{gray!6}{\hspace{1em}TLS} & \cellcolor{gray!6}{2.037} & \cellcolor{gray!6}{-0.003} & \cellcolor{gray!6}{-0.609} & \cellcolor{gray!6}{-0.494} & \cellcolor{gray!6}{-5.618} & \cellcolor{gray!6}{-18.042} & \cellcolor{gray!6}{-0.022} & \cellcolor{gray!6}{-0.008} & \cellcolor{gray!6}{-0.018} & \cellcolor{gray!6}{1.930}\\
\hspace{1em}VUT & 0.432 & -0.480 & -0.115 & -6.206 & -0.055 & -3.348 & -0.014 & -- & -- & 18.000\\
\cellcolor{gray!6}{\hspace{1em}VNM} & \cellcolor{gray!6}{2.197} & \cellcolor{gray!6}{-0.525} & \cellcolor{gray!6}{-0.658} & \cellcolor{gray!6}{-3.420} & \cellcolor{gray!6}{-5.369} & \cellcolor{gray!6}{-6.549} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{-0.042} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{3.270}\\*
\end{longtable}
\endgroup{}

(ref:cap-bt-par-reg) **Back-transformed parameters per region**. We back-transformed the parameters using the mean and standard-deviation of each continuous variable for each study-area. We then used the forest cover in 2010 to compute the back-transformed parameter estimate weighted mean per region. Doing so, we can use Eq. \@ref(eq:icar) to compute the change in the probability of deforestation associated with a particular change in the explanatory variables, in their original units. To use this table of parameters, distances and elevation must be expressed in kilometers (Km), and slope must be expressed in hecto-degrees ($10^2$°).\vspace{0.5cm}

\begingroup\fontsize{10}{12}\selectfont

\begin{longtable}[t]{lrrrrrrrrrr}
\caption{(\#tab:bt-par-region)(ref:cap-bt-par-reg)}\\
\toprule
\multicolumn{1}{l}{Region} & \multicolumn{1}{r}{int} & \multicolumn{1}{r}{pa} & \multicolumn{1}{r}{elev} & \multicolumn{1}{r}{slope} & \multicolumn{1}{r}{ddefor} & \multicolumn{1}{r}{dedge} & \multicolumn{1}{r}{driver} & \multicolumn{1}{r}{droad} & \multicolumn{1}{r}{dtown} & \multicolumn{1}{r}{Vrho} \\
 &  &  & (Km) & ($10^2$°) & (Km) & (Km) & (Km) & (Km) & (Km) & \\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
\multicolumn{1}{l}{Region} & \multicolumn{1}{r}{int} & \multicolumn{1}{r}{pa} & \multicolumn{1}{r}{elev} & \multicolumn{1}{r}{slope} & \multicolumn{1}{r}{ddefor} & \multicolumn{1}{r}{dedge} & \multicolumn{1}{r}{driver} & \multicolumn{1}{r}{droad} & \multicolumn{1}{r}{dtown} & \multicolumn{1}{r}{Vrho} \\
 &  &  & (Km) & ($10^2$°) & (Km) & (Km) & (Km) & (Km) & (Km) & \\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\cellcolor{gray!6}{India} & \cellcolor{gray!6}{2.271} & \cellcolor{gray!6}{-0.189} & \cellcolor{gray!6}{-0.558} & \cellcolor{gray!6}{-1.750} & \cellcolor{gray!6}{-4.152} & \cellcolor{gray!6}{-11.582} & \cellcolor{gray!6}{-0.005} & \cellcolor{gray!6}{-0.033} & \cellcolor{gray!6}{-0.003} & \cellcolor{gray!6}{2.033}\\
Brazil & 2.100 & -0.663 & -0.738 & -0.670 & -0.912 & -2.177 & -0.010 & -0.010 & -0.009 & 9.168\\
\cellcolor{gray!6}{America} & \cellcolor{gray!6}{2.338} & \cellcolor{gray!6}{-0.567} & \cellcolor{gray!6}{-0.978} & \cellcolor{gray!6}{-1.644} & \cellcolor{gray!6}{-0.898} & \cellcolor{gray!6}{-2.400} & \cellcolor{gray!6}{-0.006} & \cellcolor{gray!6}{-0.010} & \cellcolor{gray!6}{-0.012} & \cellcolor{gray!6}{7.734}\\
Africa & 2.336 & -0.284 & -0.367 & -1.608 & -1.365 & -2.649 & -0.001 & -0.028 & -0.015 & 5.453\\
\cellcolor{gray!6}{Asia} & \cellcolor{gray!6}{1.799} & \cellcolor{gray!6}{-0.554} & \cellcolor{gray!6}{-0.916} & \cellcolor{gray!6}{-5.190} & \cellcolor{gray!6}{-1.438} & \cellcolor{gray!6}{-3.834} & \cellcolor{gray!6}{-0.003} & \cellcolor{gray!6}{-0.020} & \cellcolor{gray!6}{-0.011} & \cellcolor{gray!6}{5.797}\\
All continents & 2.212 & -0.503 & -0.831 & -2.461 & -1.125 & -2.788 & -0.005 & -0.017 & -0.012 & 6.788\\*
\end{longtable}
\endgroup{}

\newpage

<!------------------------------------------------>
<!-- Effect of protected areas on deforestation -->
<!------------------------------------------------>

## Effect of protected areas on deforestation



(ref:cap-pa) **Effect of protected areas on deforestation**. We show here the estimated effect of the presence of a protected area on the probability of deforestation for each study-area. We computed the mean ("Mean"), the standard-deviation ("Sd"), and the bayesian 95% credible interval ("CI 95%") of the estimated parameter. Column "signif" indicates (with a star) that the estimated effect was negative and significantly different from zero (zero not included in the credible interval). Out of the 119 study-areas, 68 showed a significant negative effect (57% of the countries). These 68 study-areas accounted for 88% of the moist tropical forest in 2010 ("fc2010" in Kha).\vspace{0.5cm}

\begingroup\fontsize{11}{13}\selectfont

\begin{longtable}[t]{lrrrrc}
\caption{(\#tab:pa)(ref:cap-pa)}\\
\toprule
Country -- Study-area & fc2010 & Mean & Sd & CI 95\% & signif\\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
Country -- Study-area & fc2010 & Mean & Sd & CI 95\% & signif\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{America}}\\
\cellcolor{gray!6}{\hspace{1em}Antigua and B.} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{-0.259} & \cellcolor{gray!6}{0.213} & \cellcolor{gray!6}{(-0.669,  0.165)} & \cellcolor{gray!6}{}\\
\hspace{1em}Bahamas & 148 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Barbados} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{-0.277} & \cellcolor{gray!6}{0.133} & \cellcolor{gray!6}{(-0.532, -0.011)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Belize & 1,462 & -0.828 & 0.098 & (-1.020, -0.632) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Bolivia} & \cellcolor{gray!6}{32,612} & \cellcolor{gray!6}{-0.204} & \cellcolor{gray!6}{0.054} & \cellcolor{gray!6}{(-0.321, -0.094)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Acre & 13,646 & -0.723 & 0.109 & (-0.937, -0.521) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Alagoas} & \cellcolor{gray!6}{103} & \cellcolor{gray!6}{-0.248} & \cellcolor{gray!6}{0.131} & \cellcolor{gray!6}{(-0.526,  0.014)} & \cellcolor{gray!6}{}\\
\hspace{1em}Brazil – Amapa & 11,602 & -0.157 & 0.147 & (-0.439,  0.124) & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Amazonas} & \cellcolor{gray!6}{148,106} & \cellcolor{gray!6}{-0.558} & \cellcolor{gray!6}{0.066} & \cellcolor{gray!6}{(-0.690, -0.431)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Bahia & 2,272 & -0.171 & 0.071 & (-0.307, -0.037) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Ceara} & \cellcolor{gray!6}{52} & \cellcolor{gray!6}{-0.654} & \cellcolor{gray!6}{0.117} & \cellcolor{gray!6}{(-0.881, -0.427)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Espirito Santo & 484 & -0.056 & 0.089 & (-0.230,  0.115) & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Goias} & \cellcolor{gray!6}{532} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Brazil – Maranhao & 4,156 & -0.056 & 0.090 & (-0.238,  0.115) & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Mato Grosso} & \cellcolor{gray!6}{34,846} & \cellcolor{gray!6}{-0.477} & \cellcolor{gray!6}{0.073} & \cellcolor{gray!6}{(-0.615, -0.339)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Mato Grosso do Sul & 852 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Minas Gerais} & \cellcolor{gray!6}{1,477} & \cellcolor{gray!6}{-0.113} & \cellcolor{gray!6}{0.093} & \cellcolor{gray!6}{(-0.296,  0.070)} & \cellcolor{gray!6}{}\\
\hspace{1em}Brazil – Para & 93,017 & -0.992 & 0.063 & (-1.120, -0.868) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Paraiba} & \cellcolor{gray!6}{42} & \cellcolor{gray!6}{-0.453} & \cellcolor{gray!6}{0.112} & \cellcolor{gray!6}{(-0.684, -0.240)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Parana & 3,296 & -0.206 & 0.083 & (-0.363, -0.044) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Pernambouco} & \cellcolor{gray!6}{124} & \cellcolor{gray!6}{-0.405} & \cellcolor{gray!6}{0.073} & \cellcolor{gray!6}{(-0.554, -0.260)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Piaui & 81 & -0.144 & 0.212 & (-0.582,  0.259) & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio de Janeiro} & \cellcolor{gray!6}{875} & \cellcolor{gray!6}{-0.205} & \cellcolor{gray!6}{0.067} & \cellcolor{gray!6}{(-0.331, -0.071)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Rio Grande do Norte & 26 & -0.301 & 0.099 & (-0.508, -0.117) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio Grande do Sul} & \cellcolor{gray!6}{2,926} & \cellcolor{gray!6}{-0.143} & \cellcolor{gray!6}{0.124} & \cellcolor{gray!6}{(-0.386,  0.098)} & \cellcolor{gray!6}{}\\
\hspace{1em}Brazil – Rondonia & 14,343 & -1.340 & 0.081 & (-1.500, -1.190) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Roraima} & \cellcolor{gray!6}{16,277} & \cellcolor{gray!6}{-0.670} & \cellcolor{gray!6}{0.103} & \cellcolor{gray!6}{(-0.898, -0.465)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Santa Catarina & 3,171 & -0.232 & 0.118 & (-0.473, -0.012) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Sao Paulo} & \cellcolor{gray!6}{3,336} & \cellcolor{gray!6}{-0.185} & \cellcolor{gray!6}{0.077} & \cellcolor{gray!6}{(-0.336, -0.035)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Sergipe & 65 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Tocantins} & \cellcolor{gray!6}{1,407} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Colombia & 67,322 & -0.535 & 0.053 & (-0.632, -0.433) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Costa Rica} & \cellcolor{gray!6}{2,359} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Cuba & 1,538 & -0.077 & 0.072 & (-0.229,  0.058) & \\
\cellcolor{gray!6}{\hspace{1em}Dominica} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Dominican Rep. & 1,137 & -0.357 & 0.074 & (-0.502, -0.203) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Ecuador} & \cellcolor{gray!6}{14,968} & \cellcolor{gray!6}{-0.867} & \cellcolor{gray!6}{0.098} & \cellcolor{gray!6}{(-1.050, -0.675)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}El Salvador & 117 & -0.426 & 0.116 & (-0.643, -0.186) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}French Guiana} & \cellcolor{gray!6}{8,132} & \cellcolor{gray!6}{-1.070} & \cellcolor{gray!6}{0.217} & \cellcolor{gray!6}{(-1.470, -0.639)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Grenada & 23 & -0.287 & 0.195 & (-0.647,  0.144) & \\
\cellcolor{gray!6}{\hspace{1em}Guadeloupe} & \cellcolor{gray!6}{84} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Guatemala & 2,959 & -0.172 & 0.086 & (-0.340, -0.003) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Guyana} & \cellcolor{gray!6}{18,655} & \cellcolor{gray!6}{-0.475} & \cellcolor{gray!6}{0.225} & \cellcolor{gray!6}{(-0.909, -0.047)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Haiti & 187 & -0.152 & 0.093 & (-0.325,  0.024) & \\
\cellcolor{gray!6}{\hspace{1em}Honduras} & \cellcolor{gray!6}{3,236} & \cellcolor{gray!6}{-0.403} & \cellcolor{gray!6}{0.068} & \cellcolor{gray!6}{(-0.539, -0.273)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Jamaica & 473 & -0.043 & 0.103 & (-0.252,  0.155) & \\
\cellcolor{gray!6}{\hspace{1em}Martinique} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{-0.124} & \cellcolor{gray!6}{0.223} & \cellcolor{gray!6}{(-0.517,  0.316)} & \cellcolor{gray!6}{}\\
\hspace{1em}Mexico & 8,340 & -0.328 & 0.074 & (-0.477, -0.189) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Montserrat} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Nicaragua & 4,549 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Panama} & \cellcolor{gray!6}{4,311} & \cellcolor{gray!6}{-0.344} & \cellcolor{gray!6}{0.078} & \cellcolor{gray!6}{(-0.491, -0.192)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Paraguay & 1,769 & -0.167 & 0.096 & (-0.355,  0.019) & \\
\cellcolor{gray!6}{\hspace{1em}Peru} & \cellcolor{gray!6}{73,476} & \cellcolor{gray!6}{-0.684} & \cellcolor{gray!6}{0.064} & \cellcolor{gray!6}{(-0.812, -0.563)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Puerto Rico & 406 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Saint Kitts and N.} & \cellcolor{gray!6}{10} & \cellcolor{gray!6}{-0.137} & \cellcolor{gray!6}{0.276} & \cellcolor{gray!6}{(-0.699,  0.415)} & \cellcolor{gray!6}{}\\
\hspace{1em}Saint Lucia & 50 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Saint Martin} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{-1.180} & \cellcolor{gray!6}{0.505} & \cellcolor{gray!6}{(-2.190, -0.148)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Saint Vincent & 30 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Sint Maarten} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Suriname & 13,814 & -0.171 & 0.143 & (-0.455,  0.098) & \\
\cellcolor{gray!6}{\hspace{1em}Trinidad and Tobago} & \cellcolor{gray!6}{345} & \cellcolor{gray!6}{-0.292} & \cellcolor{gray!6}{0.077} & \cellcolor{gray!6}{(-0.435, -0.142)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Venezuela & 43,476 & -0.120 & 0.059 & (-0.237, -0.006) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Virgin Isl. UK} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Virgin Isl. US & 9 & -0.037 & 0.256 & (-0.520,  0.474) & \\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{Africa}}\\
\cellcolor{gray!6}{\hspace{1em}Angola} & \cellcolor{gray!6}{6,236} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Benin & 48 & -0.218 & 0.199 & (-0.595,  0.174) & \\
\cellcolor{gray!6}{\hspace{1em}Burundi} & \cellcolor{gray!6}{65} & \cellcolor{gray!6}{-1.860} & \cellcolor{gray!6}{0.171} & \cellcolor{gray!6}{(-2.200, -1.540)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Cameroon & 23,688 & -0.933 & 0.122 & (-1.150, -0.697) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}CAR} & \cellcolor{gray!6}{9,416} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Comoros & 90 & -0.038 & 0.180 & (-0.396,  0.315) & \\
\cellcolor{gray!6}{\hspace{1em}Congo} & \cellcolor{gray!6}{24,006} & \cellcolor{gray!6}{-0.112} & \cellcolor{gray!6}{0.101} & \cellcolor{gray!6}{(-0.302,  0.100)} & \cellcolor{gray!6}{}\\
\hspace{1em}DRC & 126,505 & -0.261 & 0.078 & (-0.409, -0.096) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Eq. Guinea} & \cellcolor{gray!6}{2,645} & \cellcolor{gray!6}{-0.368} & \cellcolor{gray!6}{0.127} & \cellcolor{gray!6}{(-0.616, -0.112)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Ethiopia & 2,907 & -0.241 & 0.062 & (-0.354, -0.114) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Gabon} & \cellcolor{gray!6}{24,124} & \cellcolor{gray!6}{-0.206} & \cellcolor{gray!6}{0.143} & \cellcolor{gray!6}{(-0.486,  0.097)} & \cellcolor{gray!6}{}\\
\hspace{1em}Gambia & 48 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Ghana} & \cellcolor{gray!6}{4,517} & \cellcolor{gray!6}{-0.174} & \cellcolor{gray!6}{0.062} & \cellcolor{gray!6}{(-0.298, -0.050)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Guinea & 1,269 & -0.272 & 0.097 & (-0.469, -0.083) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Guinea Bissau} & \cellcolor{gray!6}{349} & \cellcolor{gray!6}{-0.502} & \cellcolor{gray!6}{0.166} & \cellcolor{gray!6}{(-0.813, -0.212)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Ivory Coast & 6,439 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Kenya} & \cellcolor{gray!6}{902} & \cellcolor{gray!6}{-0.073} & \cellcolor{gray!6}{0.066} & \cellcolor{gray!6}{(-0.196,  0.049)} & \cellcolor{gray!6}{}\\
\hspace{1em}Liberia & 8,818 & -0.347 & 0.122 & (-0.588, -0.106) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Madagascar} & \cellcolor{gray!6}{6,242} & \cellcolor{gray!6}{-0.393} & \cellcolor{gray!6}{0.095} & \cellcolor{gray!6}{(-0.576, -0.210)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Malawi & 74 & -0.428 & 0.138 & (-0.695, -0.149) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Mauritius} & \cellcolor{gray!6}{54} & \cellcolor{gray!6}{-0.486} & \cellcolor{gray!6}{0.105} & \cellcolor{gray!6}{(-0.676, -0.269)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Mayotte & 18 & -1.320 & 0.132 & (-1.570, -1.060) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Nigeria} & \cellcolor{gray!6}{7,323} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Reunion & 164 & -0.391 & 0.108 & (-0.624, -0.193) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Rwanda} & \cellcolor{gray!6}{199} & \cellcolor{gray!6}{-0.808} & \cellcolor{gray!6}{0.198} & \cellcolor{gray!6}{(-1.210, -0.450)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Senegal & 139 & -0.080 & 0.136 & (-0.346,  0.184) & \\
\cellcolor{gray!6}{\hspace{1em}Sierra Leone} & \cellcolor{gray!6}{2,413} & \cellcolor{gray!6}{-0.231} & \cellcolor{gray!6}{0.109} & \cellcolor{gray!6}{(-0.442, -0.014)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}South Sudan & 209 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Tanzania} & \cellcolor{gray!6}{1,225} & \cellcolor{gray!6}{-0.553} & \cellcolor{gray!6}{0.072} & \cellcolor{gray!6}{(-0.695, -0.414)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Togo & 106 & -0.751 & 0.097 & (-0.934, -0.548) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Uganda} & \cellcolor{gray!6}{1,106} & \cellcolor{gray!6}{-0.952} & \cellcolor{gray!6}{0.069} & \cellcolor{gray!6}{(-1.090, -0.813)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Zambia & 121 & -0.196 & 0.070 & (-0.319, -0.049) & $\star$\\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{Asia}}\\
\cellcolor{gray!6}{\hspace{1em}Australia – Queensland} & \cellcolor{gray!6}{2,069} & \cellcolor{gray!6}{-0.428} & \cellcolor{gray!6}{0.113} & \cellcolor{gray!6}{(-0.643, -0.202)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Bangladesh & 973 & -0.511 & 0.162 & (-0.798, -0.188) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Bhutan} & \cellcolor{gray!6}{2,386} & \cellcolor{gray!6}{-0.009} & \cellcolor{gray!6}{0.056} & \cellcolor{gray!6}{(-0.122,  0.098)} & \cellcolor{gray!6}{}\\
\hspace{1em}Brunei & 505 & -1.100 & 0.228 & (-1.550, -0.650) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Cambodia} & \cellcolor{gray!6}{4,079} & \cellcolor{gray!6}{-0.255} & \cellcolor{gray!6}{0.125} & \cellcolor{gray!6}{(-0.490,  0.006)} & \cellcolor{gray!6}{}\\
\hspace{1em}Fiji & 1,062 & -0.381 & 0.139 & (-0.640, -0.080) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}India – Andaman and N.} & \cellcolor{gray!6}{614} & \cellcolor{gray!6}{0.029} & \cellcolor{gray!6}{0.183} & \cellcolor{gray!6}{(-0.321,  0.394)} & \cellcolor{gray!6}{}\\
\hspace{1em}India – North-East & 7,537 & -0.163 & 0.096 & (-0.358,  0.022) & \\
\cellcolor{gray!6}{\hspace{1em}India – West. Ghats} & \cellcolor{gray!6}{2,846} & \cellcolor{gray!6}{-0.305} & \cellcolor{gray!6}{0.088} & \cellcolor{gray!6}{(-0.494, -0.142)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Indonesia & 128,119 & -0.738 & 0.081 & (-0.891, -0.579) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Laos} & \cellcolor{gray!6}{10,970} & \cellcolor{gray!6}{-0.427} & \cellcolor{gray!6}{0.072} & \cellcolor{gray!6}{(-0.558, -0.282)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Malaysia & 22,594 & -1.080 & 0.106 & (-1.290, -0.876) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Myanmar} & \cellcolor{gray!6}{18,387} & \cellcolor{gray!6}{-0.210} & \cellcolor{gray!6}{0.086} & \cellcolor{gray!6}{(-0.376, -0.044)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}New Caledonia & 1,016 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Papua New Guinea} & \cellcolor{gray!6}{40,433} & \cellcolor{gray!6}{-0.288} & \cellcolor{gray!6}{0.114} & \cellcolor{gray!6}{(-0.509, -0.068)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Philippines & 14,423 & -0.180 & 0.077 & (-0.336, -0.033) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Singapore} & \cellcolor{gray!6}{15} & \cellcolor{gray!6}{-1.230} & \cellcolor{gray!6}{0.494} & \cellcolor{gray!6}{(-2.220, -0.319)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Solomon Isl. & 2,827 & -0.191 & 0.319 & (-0.799,  0.438) & \\
\cellcolor{gray!6}{\hspace{1em}Sri Lanka} & \cellcolor{gray!6}{1,784} & \cellcolor{gray!6}{-0.307} & \cellcolor{gray!6}{0.069} & \cellcolor{gray!6}{(-0.435, -0.160)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Thailand & 6,804 & -0.202 & 0.055 & (-0.306, -0.093) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Timor-Leste} & \cellcolor{gray!6}{92} & \cellcolor{gray!6}{-0.003} & \cellcolor{gray!6}{0.075} & \cellcolor{gray!6}{(-0.145,  0.148)} & \cellcolor{gray!6}{}\\
\hspace{1em}Vanuatu & 1,256 & -0.480 & 0.174 & (-0.813, -0.148) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Vietnam} & \cellcolor{gray!6}{9,739} & \cellcolor{gray!6}{-0.525} & \cellcolor{gray!6}{0.086} & \cellcolor{gray!6}{(-0.680, -0.364)} & \cellcolor{gray!6}{$\star$}\\*
\end{longtable}
\endgroup{}

\newpage

<!----------------------------------------------------->
<!-- Effect of the distance to road on deforestation -->
<!----------------------------------------------------->

## Effect of the distance to road on deforestation



(ref:cap-road) **Effect of the distance to road on deforestation**. We show here the estimated effect of the distance to the nearest road on the probability of deforestation for each study-area. We computed the mean ("Mean"), the standard-deviation ("Sd"), and the bayesian 95% credible interval ("CI 95%") of the estimated parameter. Column "signif" indicates (with a star) that the estimated effect was negative and significantly different from zero (zero not included in the credible interval). Out of the 119 study-areas, 61 showed a significant negative effect (51% of the countries). These 61 study-areas accounted for 87% of the moist tropical forest in 2010 ("fc2010" in Kha).\vspace{0.5cm}

\begingroup\fontsize{11}{13}\selectfont

\begin{longtable}[t]{lrrrrc}
\caption{(\#tab:road)(ref:cap-road)}\\
\toprule
Country -- Study-area & fc2010 & Mean & Sd & CI 95\% & signif\\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
Country -- Study-area & fc2010 & Mean & Sd & CI 95\% & signif\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{America}}\\
\cellcolor{gray!6}{\hspace{1em}Antigua and B.} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Bahamas & 148 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Barbados} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Belize & 1,462 & -0.169 & 0.061 & (-0.299, -0.051) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Bolivia} & \cellcolor{gray!6}{32,612} & \cellcolor{gray!6}{-0.333} & \cellcolor{gray!6}{0.044} & \cellcolor{gray!6}{(-0.421, -0.252)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Acre & 13,646 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Alagoas} & \cellcolor{gray!6}{103} & \cellcolor{gray!6}{-0.116} & \cellcolor{gray!6}{0.031} & \cellcolor{gray!6}{(-0.176, -0.055)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Amapa & 11,602 & -0.293 & 0.181 & (-0.635,  0.012) & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Amazonas} & \cellcolor{gray!6}{148,106} & \cellcolor{gray!6}{-0.721} & \cellcolor{gray!6}{0.119} & \cellcolor{gray!6}{(-0.917, -0.493)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Bahia & 2,272 & -0.050 & 0.025 & (-0.101, -8e-04) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Ceara} & \cellcolor{gray!6}{52} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Brazil – Espirito Santo & 484 & -0.068 & 0.026 & (-0.118, -0.017) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Goias} & \cellcolor{gray!6}{532} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Brazil – Maranhao & 4,156 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Mato Grosso} & \cellcolor{gray!6}{34,846} & \cellcolor{gray!6}{-0.252} & \cellcolor{gray!6}{0.058} & \cellcolor{gray!6}{(-0.369, -0.145)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Mato Grosso do Sul & 852 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Minas Gerais} & \cellcolor{gray!6}{1,477} & \cellcolor{gray!6}{-0.069} & \cellcolor{gray!6}{0.030} & \cellcolor{gray!6}{(-0.127, -0.015)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Para & 93,017 & -0.544 & 0.105 & (-0.740, -0.309) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Paraiba} & \cellcolor{gray!6}{42} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Brazil – Parana & 3,296 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Pernambouco} & \cellcolor{gray!6}{124} & \cellcolor{gray!6}{-0.148} & \cellcolor{gray!6}{0.029} & \cellcolor{gray!6}{(-0.205, -0.091)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Piaui & 81 & -0.041 & 0.040 & (-0.123,  0.031) & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio de Janeiro} & \cellcolor{gray!6}{875} & \cellcolor{gray!6}{-0.100} & \cellcolor{gray!6}{0.037} & \cellcolor{gray!6}{(-0.176, -0.030)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Rio Grande do Norte & 26 & -0.094 & 0.029 & (-0.153, -0.034) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio Grande do Sul} & \cellcolor{gray!6}{2,926} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Brazil – Rondonia & 14,343 & -0.240 & 0.094 & (-0.422, -0.060) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Roraima} & \cellcolor{gray!6}{16,277} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Brazil – Santa Catarina & 3,171 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Sao Paulo} & \cellcolor{gray!6}{3,336} & \cellcolor{gray!6}{-0.095} & \cellcolor{gray!6}{0.039} & \cellcolor{gray!6}{(-0.170, -0.018)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Brazil – Sergipe & 65 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Brazil – Tocantins} & \cellcolor{gray!6}{1,407} & \cellcolor{gray!6}{0.010} & \cellcolor{gray!6}{0.046} & \cellcolor{gray!6}{(-0.075,  0.095)} & \cellcolor{gray!6}{}\\
\hspace{1em}Colombia & 67,322 & -0.817 & 0.118 & (-0.973, -0.506) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Costa Rica} & \cellcolor{gray!6}{2,359} & \cellcolor{gray!6}{-0.060} & \cellcolor{gray!6}{0.058} & \cellcolor{gray!6}{(-0.177,  0.050)} & \cellcolor{gray!6}{}\\
\hspace{1em}Cuba & 1,538 & -0.147 & 0.061 & (-0.279, -0.045) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Dominica} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Dominican Rep. & 1,137 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Ecuador} & \cellcolor{gray!6}{14,968} & \cellcolor{gray!6}{-0.455} & \cellcolor{gray!6}{0.095} & \cellcolor{gray!6}{(-0.639, -0.253)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}El Salvador & 117 & -0.068 & 0.031 & (-0.130, -0.008) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}French Guiana} & \cellcolor{gray!6}{8,132} & \cellcolor{gray!6}{-0.894} & \cellcolor{gray!6}{0.574} & \cellcolor{gray!6}{(-1.960,  0.031)} & \cellcolor{gray!6}{}\\
\hspace{1em}Grenada & 23 & -0.232 & 0.088 & (-0.457, -0.096) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Guadeloupe} & \cellcolor{gray!6}{84} & \cellcolor{gray!6}{-0.039} & \cellcolor{gray!6}{0.066} & \cellcolor{gray!6}{(-0.165,  0.089)} & \cellcolor{gray!6}{}\\
\hspace{1em}Guatemala & 2,959 & -0.156 & 0.079 & (-0.304, -0.020) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Guyana} & \cellcolor{gray!6}{18,655} & \cellcolor{gray!6}{-0.221} & \cellcolor{gray!6}{0.181} & \cellcolor{gray!6}{(-0.492,  0.147)} & \cellcolor{gray!6}{}\\
\hspace{1em}Haiti & 187 & -0.124 & 0.030 & (-0.186, -0.069) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Honduras} & \cellcolor{gray!6}{3,236} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Jamaica & 473 & -0.067 & 0.032 & (-0.134, -0.005) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Martinique} & \cellcolor{gray!6}{75} & \cellcolor{gray!6}{-0.020} & \cellcolor{gray!6}{0.040} & \cellcolor{gray!6}{(-0.099,  0.058)} & \cellcolor{gray!6}{}\\
\hspace{1em}Mexico & 8,340 & -0.056 & 0.036 & (-0.129,  0.018) & \\
\cellcolor{gray!6}{\hspace{1em}Montserrat} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Nicaragua & 4,549 & -0.178 & 0.064 & (-0.298, -0.050) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Panama} & \cellcolor{gray!6}{4,311} & \cellcolor{gray!6}{-0.521} & \cellcolor{gray!6}{0.065} & \cellcolor{gray!6}{(-0.649, -0.395)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Paraguay & 1,769 & -0.058 & 0.028 & (-0.114, -0.003) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Peru} & \cellcolor{gray!6}{73,476} & \cellcolor{gray!6}{-0.531} & \cellcolor{gray!6}{0.096} & \cellcolor{gray!6}{(-0.692, -0.364)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Puerto Rico & 406 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Saint Kitts and N.} & \cellcolor{gray!6}{10} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Saint Lucia & 50 & -0.044 & 0.062 & (-0.168,  0.084) & \\
\cellcolor{gray!6}{\hspace{1em}Saint Martin} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Saint Vincent & 30 & -0.098 & 0.062 & (-0.222,  0.025) & \\
\cellcolor{gray!6}{\hspace{1em}Sint Maarten} & \cellcolor{gray!6}{0} & \cellcolor{gray!6}{-0.186} & \cellcolor{gray!6}{0.086} & \cellcolor{gray!6}{(-0.390, -0.048)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Suriname & 13,814 & -0.553 & 0.467 & (-1.190,  0.295) & \\
\cellcolor{gray!6}{\hspace{1em}Trinidad and Tobago} & \cellcolor{gray!6}{345} & \cellcolor{gray!6}{-0.074} & \cellcolor{gray!6}{0.040} & \cellcolor{gray!6}{(-0.150,  0.008)} & \cellcolor{gray!6}{}\\
\hspace{1em}Venezuela & 43,476 & -0.444 & 0.132 & (-0.734, -0.241) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Virgin Isl. UK} & \cellcolor{gray!6}{3} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Virgin Isl. US & 9 & -0.037 & 0.032 & (-0.099,  0.032) & \\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{Africa}}\\
\cellcolor{gray!6}{\hspace{1em}Angola} & \cellcolor{gray!6}{6,236} & \cellcolor{gray!6}{-0.231} & \cellcolor{gray!6}{0.040} & \cellcolor{gray!6}{(-0.309, -0.150)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Benin & 48 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Burundi} & \cellcolor{gray!6}{65} & \cellcolor{gray!6}{-0.078} & \cellcolor{gray!6}{0.041} & \cellcolor{gray!6}{(-0.158,  6e-04)} & \cellcolor{gray!6}{}\\
\hspace{1em}Cameroon & 23,688 & -0.449 & 0.055 & (-0.564, -0.352) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}CAR} & \cellcolor{gray!6}{9,416} & \cellcolor{gray!6}{-0.086} & \cellcolor{gray!6}{0.065} & \cellcolor{gray!6}{(-0.201,  0.039)} & \cellcolor{gray!6}{}\\
\hspace{1em}Comoros & 90 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Congo} & \cellcolor{gray!6}{24,006} & \cellcolor{gray!6}{-0.708} & \cellcolor{gray!6}{0.100} & \cellcolor{gray!6}{(-0.892, -0.510)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}DRC & 126,505 & -0.368 & 0.025 & (-0.416, -0.319) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Eq. Guinea} & \cellcolor{gray!6}{2,645} & \cellcolor{gray!6}{-1.070} & \cellcolor{gray!6}{0.062} & \cellcolor{gray!6}{(-1.190, -0.951)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Ethiopia & 2,907 & -0.056 & 0.040 & (-0.131,  0.025) & \\
\cellcolor{gray!6}{\hspace{1em}Gabon} & \cellcolor{gray!6}{24,124} & \cellcolor{gray!6}{-0.720} & \cellcolor{gray!6}{0.091} & \cellcolor{gray!6}{(-0.872, -0.540)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Gambia & 48 & -0.341 & 0.028 & (-0.394, -0.284) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Ghana} & \cellcolor{gray!6}{4,517} & \cellcolor{gray!6}{-0.093} & \cellcolor{gray!6}{0.025} & \cellcolor{gray!6}{(-0.143, -0.043)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Guinea & 1,269 & -0.073 & 0.033 & (-0.137, -0.012) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Guinea Bissau} & \cellcolor{gray!6}{349} & \cellcolor{gray!6}{-0.039} & \cellcolor{gray!6}{0.079} & \cellcolor{gray!6}{(-0.223,  0.107)} & \cellcolor{gray!6}{}\\
\hspace{1em}Ivory Coast & 6,439 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Kenya} & \cellcolor{gray!6}{902} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Liberia & 8,818 & -0.279 & 0.047 & (-0.375, -0.196) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Madagascar} & \cellcolor{gray!6}{6,242} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{}\\
\hspace{1em}Malawi & 74 & -0.540 & 0.061 & (-0.663, -0.421) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Mauritius} & \cellcolor{gray!6}{54} & \cellcolor{gray!6}{-0.083} & \cellcolor{gray!6}{0.023} & \cellcolor{gray!6}{(-0.130, -0.036)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Mayotte & 18 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Nigeria} & \cellcolor{gray!6}{7,323} & \cellcolor{gray!6}{-0.255} & \cellcolor{gray!6}{0.060} & \cellcolor{gray!6}{(-0.373, -0.132)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Reunion & 164 & -0.070 & 0.034 & (-0.141, -0.003) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Rwanda} & \cellcolor{gray!6}{199} & \cellcolor{gray!6}{-0.278} & \cellcolor{gray!6}{0.086} & \cellcolor{gray!6}{(-0.458, -0.115)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Senegal & 139 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Sierra Leone} & \cellcolor{gray!6}{2,413} & \cellcolor{gray!6}{-0.116} & \cellcolor{gray!6}{0.028} & \cellcolor{gray!6}{(-0.170, -0.061)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}South Sudan & 209 & -0.066 & 0.084 & (-0.233,  0.106) & \\
\cellcolor{gray!6}{\hspace{1em}Tanzania} & \cellcolor{gray!6}{1,225} & \cellcolor{gray!6}{-0.248} & \cellcolor{gray!6}{0.040} & \cellcolor{gray!6}{(-0.329, -0.170)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Togo & 106 & -0.173 & 0.041 & (-0.251, -0.086) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Uganda} & \cellcolor{gray!6}{1,106} & \cellcolor{gray!6}{-0.079} & \cellcolor{gray!6}{0.045} & \cellcolor{gray!6}{(-0.170,  0.007)} & \cellcolor{gray!6}{}\\
\hspace{1em}Zambia & 121 & -0.008 & 0.056 & (-0.128,  0.095) & \\
\addlinespace[0.3em]
\multicolumn{6}{l}{\textbf{Asia}}\\
\cellcolor{gray!6}{\hspace{1em}Australia – Queensland} & \cellcolor{gray!6}{2,069} & \cellcolor{gray!6}{-0.025} & \cellcolor{gray!6}{0.058} & \cellcolor{gray!6}{(-0.130,  0.100)} & \cellcolor{gray!6}{}\\
\hspace{1em}Bangladesh & 973 & -0.185 & 0.101 & (-0.388,  0.014) & \\
\cellcolor{gray!6}{\hspace{1em}Bhutan} & \cellcolor{gray!6}{2,386} & \cellcolor{gray!6}{-0.066} & \cellcolor{gray!6}{0.036} & \cellcolor{gray!6}{(-0.138,  0.004)} & \cellcolor{gray!6}{}\\
\hspace{1em}Brunei & 505 & -0.744 & 0.208 & (-1.210, -0.388) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Cambodia} & \cellcolor{gray!6}{4,079} & \cellcolor{gray!6}{-0.419} & \cellcolor{gray!6}{0.076} & \cellcolor{gray!6}{(-0.582, -0.285)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Fiji & 1,062 & -0.199 & 0.202 & (-0.552,  0.243) & \\
\cellcolor{gray!6}{\hspace{1em}India – Andaman and N.} & \cellcolor{gray!6}{614} & \cellcolor{gray!6}{-0.102} & \cellcolor{gray!6}{0.113} & \cellcolor{gray!6}{(-0.336,  0.113)} & \cellcolor{gray!6}{}\\
\hspace{1em}India – North-East & 7,537 & -0.159 & 0.036 & (-0.224, -0.088) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}India – West. Ghats} & \cellcolor{gray!6}{2,846} & \cellcolor{gray!6}{-0.163} & \cellcolor{gray!6}{0.038} & \cellcolor{gray!6}{(-0.237, -0.090)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Indonesia & 128,119 & -0.389 & 0.032 & (-0.449, -0.320) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Laos} & \cellcolor{gray!6}{10,970} & \cellcolor{gray!6}{-0.118} & \cellcolor{gray!6}{0.032} & \cellcolor{gray!6}{(-0.184, -0.056)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Malaysia & 22,594 & -0.227 & 0.045 & (-0.326, -0.144) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Myanmar} & \cellcolor{gray!6}{18,387} & \cellcolor{gray!6}{-0.166} & \cellcolor{gray!6}{0.034} & \cellcolor{gray!6}{(-0.233, -0.106)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}New Caledonia & 1,016 & -0.122 & 0.049 & (-0.232, -0.029) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Papua New Guinea} & \cellcolor{gray!6}{40,433} & \cellcolor{gray!6}{-0.254} & \cellcolor{gray!6}{0.074} & \cellcolor{gray!6}{(-0.393, -0.102)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Philippines & 14,423 & -0.156 & 0.031 & (-0.222, -0.100) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Singapore} & \cellcolor{gray!6}{15} & \cellcolor{gray!6}{-0.882} & \cellcolor{gray!6}{0.080} & \cellcolor{gray!6}{(-1.060, -0.745)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Solomon Isl. & 2,827 & -0.485 & 0.126 & (-0.719, -0.255) & $\star$\\
\cellcolor{gray!6}{\hspace{1em}Sri Lanka} & \cellcolor{gray!6}{1,784} & \cellcolor{gray!6}{-0.108} & \cellcolor{gray!6}{0.046} & \cellcolor{gray!6}{(-0.195, -0.017)} & \cellcolor{gray!6}{$\star$}\\
\hspace{1em}Thailand & 6,804 & -0.054 & 0.038 & (-0.129,  0.023) & \\
\cellcolor{gray!6}{\hspace{1em}Timor-Leste} & \cellcolor{gray!6}{92} & \cellcolor{gray!6}{-0.022} & \cellcolor{gray!6}{0.031} & \cellcolor{gray!6}{(-0.088,  0.039)} & \cellcolor{gray!6}{}\\
\hspace{1em}Vanuatu & 1,256 & -- & -- & -- & \\
\cellcolor{gray!6}{\hspace{1em}Vietnam} & \cellcolor{gray!6}{9,739} & \cellcolor{gray!6}{-0.161} & \cellcolor{gray!6}{0.028} & \cellcolor{gray!6}{(-0.213, -0.104)} & \cellcolor{gray!6}{$\star$}\\*
\end{longtable}
\endgroup{}

\newpage

<!------------------------------------------------->
<!-- Carbon emissions associated with deforestation -->
<!------------------------------------------------->

## Cumulative carbon emissions associated to deforestation

(ref:cap-c-em) **Cumulative carbon emissions associated to future deforestation**. We computed the cumulative carbon emissions associated with future deforestation from 2020 for each study-area (C in Gg, $10^9$ g). To do so, we used our maps of projected forest cover change together with the aboveground biomass map by @Avitabile2016. The aboveground biomass map was not covering some islands in the tropics (Reunion island or Fiji for example). \vspace{0.5cm}

\begingroup\fontsize{11}{13}\selectfont

\begin{longtable}[t]{lrrrr}
\caption{(\#tab:c-em)(ref:cap-c-em)}\\
\toprule
Country -- Study-area & C2040 & C2060 & C2080 & C2100\\
\midrule
\endfirsthead
\caption[]{\textit{(continued)}}\\
\toprule
Country -- Study-area & C2040 & C2060 & C2080 & C2100\\
\midrule
\endhead

\endfoot
\bottomrule
\endlastfoot
\addlinespace[0.3em]
\multicolumn{5}{l}{\textbf{America}}\\
\cellcolor{gray!6}{\hspace{1em}Antigua and B.} & \cellcolor{gray!6}{12} & \cellcolor{gray!6}{37} & \cellcolor{gray!6}{61} & \cellcolor{gray!6}{61}\\
\hspace{1em}Bahamas & 955 & 2,013 & 2,686 & 2,686\\
\cellcolor{gray!6}{\hspace{1em}Barbados} & \cellcolor{gray!6}{15} & \cellcolor{gray!6}{38} & \cellcolor{gray!6}{51} & \cellcolor{gray!6}{51}\\
\hspace{1em}Belize & 11,614 & 24,624 & 40,394 & 59,273\\
\cellcolor{gray!6}{\hspace{1em}Bolivia} & \cellcolor{gray!6}{396,729} & \cellcolor{gray!6}{809,521} & \cellcolor{gray!6}{1,253,118} & \cellcolor{gray!6}{1,717,256}\\
\hspace{1em}Brazil – Acre & 99,775 & 228,000 & 386,518 & 561,680\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Alagoas} & \cellcolor{gray!6}{419} & \cellcolor{gray!6}{2,390} & \cellcolor{gray!6}{2,390} & \cellcolor{gray!6}{2,390}\\
\hspace{1em}Brazil – Amapa & 12,825 & 40,047 & 92,561 & 182,303\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Amazonas} & \cellcolor{gray!6}{324,816} & \cellcolor{gray!6}{646,720} & \cellcolor{gray!6}{990,323} & \cellcolor{gray!6}{1,359,616}\\
\hspace{1em}Brazil – Bahia & 24,118 & 62,243 & 120,978 & 143,914\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Ceara} & \cellcolor{gray!6}{1,092} & \cellcolor{gray!6}{2,033} & \cellcolor{gray!6}{2,033} & \cellcolor{gray!6}{2,033}\\
\hspace{1em}Brazil – Espirito Santo & 3,966 & 14,979 & 23,955 & 23,955\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Goias} & \cellcolor{gray!6}{5,071} & \cellcolor{gray!6}{5,990} & \cellcolor{gray!6}{5,990} & \cellcolor{gray!6}{5,990}\\
\hspace{1em}Brazil – Maranhao & 110,426 & 259,741 & 259,741 & 259,741\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Mato Grosso} & \cellcolor{gray!6}{475,985} & \cellcolor{gray!6}{994,408} & \cellcolor{gray!6}{1,602,327} & \cellcolor{gray!6}{2,314,341}\\
\hspace{1em}Brazil – Mato Grosso do Sul & 5,624 & 15,123 & 22,553 & 22,553\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Minas Gerais} & \cellcolor{gray!6}{34,483} & \cellcolor{gray!6}{49,066} & \cellcolor{gray!6}{49,066} & \cellcolor{gray!6}{49,066}\\
\hspace{1em}Brazil – Para & 944,271 & 1,916,871 & 2,981,477 & 4,130,806\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Paraiba} & \cellcolor{gray!6}{137} & \cellcolor{gray!6}{925} & \cellcolor{gray!6}{925} & \cellcolor{gray!6}{925}\\
\hspace{1em}Brazil – Parana & 34,247 & 73,118 & 130,637 & 199,490\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Pernambouco} & \cellcolor{gray!6}{527} & \cellcolor{gray!6}{3,158} & \cellcolor{gray!6}{3,158} & \cellcolor{gray!6}{3,158}\\
\hspace{1em}Brazil – Piaui & 2,232 & 2,291 & 2,291 & 2,291\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio de Janeiro} & \cellcolor{gray!6}{5,472} & \cellcolor{gray!6}{22,554} & \cellcolor{gray!6}{61,872} & \cellcolor{gray!6}{61,872}\\
\hspace{1em}Brazil – Rio Grande do Norte & 85 & 410 & 410 & 410\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Rio Grande do Sul} & \cellcolor{gray!6}{12,126} & \cellcolor{gray!6}{32,926} & \cellcolor{gray!6}{69,307} & \cellcolor{gray!6}{127,534}\\
\hspace{1em}Brazil – Rondonia & 250,464 & 526,196 & 871,950 & 1,259,579\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Roraima} & \cellcolor{gray!6}{133,850} & \cellcolor{gray!6}{278,847} & \cellcolor{gray!6}{450,149} & \cellcolor{gray!6}{649,834}\\
\hspace{1em}Brazil – Santa Catarina & 26,710 & 59,984 & 108,273 & 175,923\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Sao Paulo} & \cellcolor{gray!6}{12,299} & \cellcolor{gray!6}{42,212} & \cellcolor{gray!6}{94,584} & \cellcolor{gray!6}{168,393}\\
\hspace{1em}Brazil – Sergipe & 357 & 1,552 & 1,552 & 1,552\\
\cellcolor{gray!6}{\hspace{1em}Brazil – Tocantins} & \cellcolor{gray!6}{31,271} & \cellcolor{gray!6}{34,696} & \cellcolor{gray!6}{34,696} & \cellcolor{gray!6}{34,696}\\
\hspace{1em}Colombia & 425,150 & 884,687 & 1,365,656 & 1,921,767\\
\cellcolor{gray!6}{\hspace{1em}Costa Rica} & \cellcolor{gray!6}{11,705} & \cellcolor{gray!6}{27,354} & \cellcolor{gray!6}{45,887} & \cellcolor{gray!6}{68,208}\\
\hspace{1em}Cuba & 7,759 & 16,913 & 26,852 & 37,744\\
\cellcolor{gray!6}{\hspace{1em}Dominica} & \cellcolor{gray!6}{157} & \cellcolor{gray!6}{362} & \cellcolor{gray!6}{586} & \cellcolor{gray!6}{831}\\
\hspace{1em}Dominican Rep. & 9,554 & 21,400 & 35,478 & 35,478\\
\cellcolor{gray!6}{\hspace{1em}Ecuador} & \cellcolor{gray!6}{69,637} & \cellcolor{gray!6}{156,589} & \cellcolor{gray!6}{253,094} & \cellcolor{gray!6}{356,106}\\
\hspace{1em}El Salvador & 1,041 & 2,094 & 2,963 & 3,315\\
\cellcolor{gray!6}{\hspace{1em}French Guiana} & \cellcolor{gray!6}{7,148} & \cellcolor{gray!6}{13,210} & \cellcolor{gray!6}{19,838} & \cellcolor{gray!6}{26,804}\\
\hspace{1em}Grenada & 134 & 308 & 497 & 536\\
\cellcolor{gray!6}{\hspace{1em}Guadeloupe} & \cellcolor{gray!6}{202} & \cellcolor{gray!6}{382} & \cellcolor{gray!6}{588} & \cellcolor{gray!6}{844}\\
\hspace{1em}Guatemala & 41,054 & 93,454 & 128,709 & 128,709\\
\cellcolor{gray!6}{\hspace{1em}Guyana} & \cellcolor{gray!6}{40,770} & \cellcolor{gray!6}{77,898} & \cellcolor{gray!6}{116,111} & \cellcolor{gray!6}{153,850}\\
\hspace{1em}Haiti & 2,278 & 3,071 & 3,071 & 3,071\\
\cellcolor{gray!6}{\hspace{1em}Honduras} & \cellcolor{gray!6}{44,057} & \cellcolor{gray!6}{88,430} & \cellcolor{gray!6}{138,387} & \cellcolor{gray!6}{138,387}\\
\hspace{1em}Jamaica & 2,210 & 4,749 & 7,483 & 10,632\\
\cellcolor{gray!6}{\hspace{1em}Martinique} & \cellcolor{gray!6}{189} & \cellcolor{gray!6}{408} & \cellcolor{gray!6}{645} & \cellcolor{gray!6}{900}\\
\hspace{1em}Mexico & 108,573 & 226,238 & 295,430 & 295,430\\
\cellcolor{gray!6}{\hspace{1em}Montserrat} & \cellcolor{gray!6}{4} & \cellcolor{gray!6}{8} & \cellcolor{gray!6}{12} & \cellcolor{gray!6}{16}\\
\hspace{1em}Nicaragua & 89,423 & 181,712 & 181,712 & 181,712\\
\cellcolor{gray!6}{\hspace{1em}Panama} & \cellcolor{gray!6}{15,984} & \cellcolor{gray!6}{34,300} & \cellcolor{gray!6}{55,548} & \cellcolor{gray!6}{79,297}\\
\hspace{1em}Paraguay & 45,755 & 62,617 & 62,617 & 62,617\\
\cellcolor{gray!6}{\hspace{1em}Peru} & \cellcolor{gray!6}{291,512} & \cellcolor{gray!6}{597,020} & \cellcolor{gray!6}{902,950} & \cellcolor{gray!6}{1,206,302}\\
\hspace{1em}Puerto Rico & 2,908 & 7,075 & 12,338 & 14,899\\
\cellcolor{gray!6}{\hspace{1em}Saint Kitts and N.} & \cellcolor{gray!6}{15} & \cellcolor{gray!6}{35} & \cellcolor{gray!6}{55} & \cellcolor{gray!6}{78}\\
\hspace{1em}Saint Lucia & 133 & 275 & 436 & 615\\
\cellcolor{gray!6}{\hspace{1em}Saint Martin} & \cellcolor{gray!6}{5} & \cellcolor{gray!6}{5} & \cellcolor{gray!6}{5} & \cellcolor{gray!6}{5}\\
\hspace{1em}Saint Vincent & 24 & 56 & 97 & 144\\
\cellcolor{gray!6}{\hspace{1em}Sint Maarten} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{1} & \cellcolor{gray!6}{1}\\
\hspace{1em}Suriname & 31,792 & 60,582 & 89,994 & 119,111\\
\cellcolor{gray!6}{\hspace{1em}Trinidad and Tobago} & \cellcolor{gray!6}{1,441} & \cellcolor{gray!6}{3,192} & \cellcolor{gray!6}{5,171} & \cellcolor{gray!6}{7,455}\\
\hspace{1em}Venezuela & 214,597 & 479,878 & 781,171 & 1,109,571\\
\cellcolor{gray!6}{\hspace{1em}Virgin Isl. UK} & \cellcolor{gray!6}{65} & \cellcolor{gray!6}{65} & \cellcolor{gray!6}{65} & \cellcolor{gray!6}{65}\\
\hspace{1em}Virgin Isl. US & 135 & 228 & 228 & 228\\
\addlinespace[0.3em]
\multicolumn{5}{l}{\textbf{Africa}}\\
\cellcolor{gray!6}{\hspace{1em}Angola} & \cellcolor{gray!6}{97,155} & \cellcolor{gray!6}{202,787} & \cellcolor{gray!6}{344,003} & \cellcolor{gray!6}{368,734}\\
\hspace{1em}Benin & 153 & 171 & 171 & 171\\
\cellcolor{gray!6}{\hspace{1em}Burundi} & \cellcolor{gray!6}{925} & \cellcolor{gray!6}{2,797} & \cellcolor{gray!6}{4,779} & \cellcolor{gray!6}{4,779}\\
\hspace{1em}Cameroon & 103,883 & 244,503 & 408,857 & 589,883\\
\cellcolor{gray!6}{\hspace{1em}CAR} & \cellcolor{gray!6}{57,822} & \cellcolor{gray!6}{117,651} & \cellcolor{gray!6}{187,150} & \cellcolor{gray!6}{272,626}\\
\hspace{1em}Comoros & 308 & 725 & 1,231 & 1,826\\
\cellcolor{gray!6}{\hspace{1em}Congo} & \cellcolor{gray!6}{83,458} & \cellcolor{gray!6}{184,765} & \cellcolor{gray!6}{309,011} & \cellcolor{gray!6}{453,102}\\
\hspace{1em}DRC & 1,455,501 & 3,300,764 & 5,510,506 & 8,090,246\\
\cellcolor{gray!6}{\hspace{1em}Eq. Guinea} & \cellcolor{gray!6}{7,686} & \cellcolor{gray!6}{16,250} & \cellcolor{gray!6}{25,450} & \cellcolor{gray!6}{34,939}\\
\hspace{1em}Ethiopia & 118,854 & 199,649 & 199,649 & 199,649\\
\cellcolor{gray!6}{\hspace{1em}Gabon} & \cellcolor{gray!6}{23,948} & \cellcolor{gray!6}{48,024} & \cellcolor{gray!6}{74,963} & \cellcolor{gray!6}{103,796}\\
\hspace{1em}Gambia & 25 & 54 & 85 & 119\\
\cellcolor{gray!6}{\hspace{1em}Ghana} & \cellcolor{gray!6}{77,157} & \cellcolor{gray!6}{186,999} & \cellcolor{gray!6}{186,999} & \cellcolor{gray!6}{186,999}\\
\hspace{1em}Guinea & 47,192 & 54,141 & 54,141 & 54,141\\
\cellcolor{gray!6}{\hspace{1em}Guinea Bissau} & \cellcolor{gray!6}{1,377} & \cellcolor{gray!6}{2,394} & \cellcolor{gray!6}{3,479} & \cellcolor{gray!6}{3,750}\\
\hspace{1em}Ivory Coast & 233,120 & 233,120 & 233,120 & 233,120\\
\cellcolor{gray!6}{\hspace{1em}Kenya} & \cellcolor{gray!6}{19,633} & \cellcolor{gray!6}{42,648} & \cellcolor{gray!6}{55,935} & \cellcolor{gray!6}{55,935}\\
\hspace{1em}Liberia & 118,813 & 319,490 & 596,951 & 915,819\\
\cellcolor{gray!6}{\hspace{1em}Madagascar} & \cellcolor{gray!6}{237,253} & \cellcolor{gray!6}{546,918} & \cellcolor{gray!6}{546,918} & \cellcolor{gray!6}{546,918}\\
\hspace{1em}Malawi & 2,502 & 2,502 & 2,502 & 2,502\\
\cellcolor{gray!6}{\hspace{1em}Mauritius} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--} & \cellcolor{gray!6}{--}\\
\hspace{1em}Mayotte & 568 & 723 & 723 & 723\\
\cellcolor{gray!6}{\hspace{1em}Nigeria} & \cellcolor{gray!6}{62,674} & \cellcolor{gray!6}{148,764} & \cellcolor{gray!6}{274,599} & \cellcolor{gray!6}{316,487}\\
\hspace{1em}Reunion & -- & -- & -- & --\\
\cellcolor{gray!6}{\hspace{1em}Rwanda} & \cellcolor{gray!6}{4,266} & \cellcolor{gray!6}{13,997} & \cellcolor{gray!6}{13,997} & \cellcolor{gray!6}{13,997}\\
\hspace{1em}Senegal & 100 & 191 & 284 & 379\\
\cellcolor{gray!6}{\hspace{1em}Sierra Leone} & \cellcolor{gray!6}{83,910} & \cellcolor{gray!6}{83,910} & \cellcolor{gray!6}{83,910} & \cellcolor{gray!6}{83,910}\\
\hspace{1em}South Sudan & 2,006 & 2,392 & 2,392 & 2,392\\
\cellcolor{gray!6}{\hspace{1em}Tanzania} & \cellcolor{gray!6}{11,230} & \cellcolor{gray!6}{25,324} & \cellcolor{gray!6}{42,633} & \cellcolor{gray!6}{63,506}\\
\hspace{1em}Togo & 1,325 & 1,325 & 1,325 & 1,325\\
\cellcolor{gray!6}{\hspace{1em}Uganda} & \cellcolor{gray!6}{58,257} & \cellcolor{gray!6}{65,092} & \cellcolor{gray!6}{65,092} & \cellcolor{gray!6}{65,092}\\
\hspace{1em}Zambia & 3,856 & 4,124 & 4,124 & 4,124\\
\addlinespace[0.3em]
\multicolumn{5}{l}{\textbf{Asia}}\\
\cellcolor{gray!6}{\hspace{1em}Australia – Queensland} & \cellcolor{gray!6}{38,639} & \cellcolor{gray!6}{73,531} & \cellcolor{gray!6}{109,650} & \cellcolor{gray!6}{148,682}\\
\hspace{1em}Bangladesh & 5,795 & 9,804 & 14,401 & 19,507\\
\cellcolor{gray!6}{\hspace{1em}Bhutan} & \cellcolor{gray!6}{29,039} & \cellcolor{gray!6}{59,618} & \cellcolor{gray!6}{90,537} & \cellcolor{gray!6}{121,877}\\
\hspace{1em}Brunei & 1,170 & 2,517 & 3,883 & 5,289\\
\cellcolor{gray!6}{\hspace{1em}Cambodia} & \cellcolor{gray!6}{240,635} & \cellcolor{gray!6}{322,672} & \cellcolor{gray!6}{322,672} & \cellcolor{gray!6}{322,672}\\
\hspace{1em}Fiji & -- & -- & -- & --\\
\cellcolor{gray!6}{\hspace{1em}India – Andaman and N.} & \cellcolor{gray!6}{2,333} & \cellcolor{gray!6}{4,683} & \cellcolor{gray!6}{7,018} & \cellcolor{gray!6}{9,422}\\
\hspace{1em}India – North-East & 129,518 & 278,660 & 442,761 & 622,153\\
\cellcolor{gray!6}{\hspace{1em}India – West. Ghats} & \cellcolor{gray!6}{60,673} & \cellcolor{gray!6}{152,780} & \cellcolor{gray!6}{166,639} & \cellcolor{gray!6}{166,639}\\
\hspace{1em}Indonesia & 1,819,929 & 4,520,972 & 7,879,538 & 11,499,558\\
\cellcolor{gray!6}{\hspace{1em}Laos} & \cellcolor{gray!6}{345,510} & \cellcolor{gray!6}{748,769} & \cellcolor{gray!6}{1,035,521} & \cellcolor{gray!6}{1,035,521}\\
\hspace{1em}Malaysia & 451,545 & 1,119,286 & 1,919,819 & 2,766,275\\
\cellcolor{gray!6}{\hspace{1em}Myanmar} & \cellcolor{gray!6}{504,967} & \cellcolor{gray!6}{1,110,067} & \cellcolor{gray!6}{1,798,752} & \cellcolor{gray!6}{2,167,216}\\
\hspace{1em}New Caledonia & -- & -- & -- & --\\
\cellcolor{gray!6}{\hspace{1em}Papua New Guinea} & \cellcolor{gray!6}{136,568} & \cellcolor{gray!6}{264,397} & \cellcolor{gray!6}{399,319} & \cellcolor{gray!6}{542,116}\\
\hspace{1em}Philippines & 156,258 & 331,913 & 524,414 & 755,230\\
\cellcolor{gray!6}{\hspace{1em}Singapore} & \cellcolor{gray!6}{50} & \cellcolor{gray!6}{143} & \cellcolor{gray!6}{258} & \cellcolor{gray!6}{405}\\
\hspace{1em}Solomon Isl. & -- & -- & -- & --\\
\cellcolor{gray!6}{\hspace{1em}Sri Lanka} & \cellcolor{gray!6}{12,956} & \cellcolor{gray!6}{31,479} & \cellcolor{gray!6}{52,847} & \cellcolor{gray!6}{77,075}\\
\hspace{1em}Thailand & 100,290 & 220,984 & 353,560 & 498,116\\
\cellcolor{gray!6}{\hspace{1em}Timor-Leste} & \cellcolor{gray!6}{1,959} & \cellcolor{gray!6}{3,999} & \cellcolor{gray!6}{5,172} & \cellcolor{gray!6}{5,172}\\
\hspace{1em}Vanuatu & -- & -- & -- & --\\
\cellcolor{gray!6}{\hspace{1em}Vietnam} & \cellcolor{gray!6}{205,809} & \cellcolor{gray!6}{450,873} & \cellcolor{gray!6}{737,542} & \cellcolor{gray!6}{737,542}\\*
\end{longtable}
\endgroup{}

\newpage

<!--chapter:end:03-Tables-S.Rmd-->

# References {-}

<!--chapter:end:04-References.Rmd-->

