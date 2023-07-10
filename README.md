SL2P-CCRS

Simplified Level 2 Processor - CCRS


The Simplified Level 2 Processor – CCRS (SL2P-CCRS) produces estimates of canopy biophysical variables (able 1) given inputs of either a top-of-atmosphere (TOA) or top-of-canopy (TOC) bi-directional reflectance spectrum together with the illumination, view and relative azimuth angles. Separate non-linear regression models are used to estimate the expected value and the expected root mean square error (RMSE) of each output. The regression estimators are optimized for multi-spectral (i.e. <10 bands with >10nm bandwidth) reflectance inputs but can be applied to arbitrary spectra as long as a radiative transfer model with sufficient accuracy to simulate such spectra is included in the processor.  This processor includes separate predictors for North American Broadleaf Forests and Needleleaf Forests in addition to a predictor for other vegetated land cover defined using the SL2P algorithm (Weiss and Baret, 2016).


Table 1.  Mapped variables.
{|class="wikitable"
|-
! Variable
| Definition
|Units
|-
| Surface_Reflectance
| Surface reflectance (at local overpass time)
| Ratio 
|-
| ALBEDO
| Black sky albedo (at local overpass time)
| Ratio
|-
| CCC
| Canopy chlorophyll content
| g/m<sup>2</sup>
|-
| CWC
| Canopy water content 
| g/m<sup>2</sup>
|-
| FAPAR
| fraction of absorbed photosynthetically active radiation (at local overpass time)
| 0-1
|-
| FCOVER
| fraction canopy cover (vertical)
| 0-1
|-
| LAI
| leaf area index (half total green foliage area per unit horizontal ground area)
| m<sup>2</sup> foliage/m<sup>2</sup> ground
|-
| LAND COVER
| Land cover based on a specified legend.
| Legend categories
|} 

Contact: richard.fernandes@canada.ca


