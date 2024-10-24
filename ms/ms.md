---
title: Wall-to-wall mapping of peat depth from Lidar terrain and airborne radiometrics in Norwegian landscapes
journal: "soil"
author:
  - given_name: Julien
    surname: Vollering
    affiliation: 1
    email: julien.vollering@hvl.no
    corresponding: true
  - given_name: Naomi
    surname: Gatis
    affiliation: 2
  - given_name: Mette
    surname: Kusk Gillespie
    affiliation: 1
  - given_name: Karl-Kristian
    surname: Muggerud
    affiliation: 1
  - given_name: Sigurd Daniel
    surname: Nerhus
    affiliation: 1
  - given_name: Knut
    surname: Rydgren
    affiliation: 1
  - given_name: Mikko
    surname: Sparf
    affiliation: 1
affiliation:
  - code: 1
    address: Department of Civil Engineering and Environmental Sciences, Western Norway University of Applied Sciences, Norway
  - code: 2
    address: Department of Geography, University of Exeter, United Kingdom
abstract: |
  The abstract goes here.
  It can also be on _multiple lines_.
bibliography: ms.bib
running:
  title: Peat depth from terrain and radiometrics
  author: Vollering et al.
# This section is mandatory even if you declare that no competing interests are present.
competinginterests: |
  *Competing interests.* The authors declare that they have no conflict of interest.
# See https://publications.copernicus.org/for_authors/licence_and_copyright.html, normally used for transferring the copyright, if needed. 
# Note: additional copyright statements for affiliated software or data need to be placed in the data availability section. 
### The following commands are for the statements about the availability of data sets and/or software code corresponding to the manuscript.
### It is strongly recommended to make use of these sections in case data sets and/or software code have been part of your research the article is based on.
### Note: unless stated otherwise, software and data affiliated with the manuscript are assumed to be published under the same licence as the article (currently Creative Commons 4.0)
availability:
  #code: |
  #  use this to add a statement when having only software code available
  #data: |
  #  use this to add a statement when having only data sets available
  codedata: |
    *Code and data availability.* Use this to add a statement when having data sets and software code available
authorcontribution: |
  *Author contributions.*
  JV: Conceptualization, Investigation, Data curation, Formal analysis, Writing – original draft.
  NG: Conceptualization, Methodology, Writing - review & editing.
  MKG: Investigation, Writing - review & editing.
  KKM: Investigation, Data curation, Writing - review & editing.
  SDN: Investigation, Writing - review & editing.
  KR: Conceptualization, Investigation, Writing - review & editing.
  MS: Investigation, Data curation, Writing - review & editing.
disclaimer: |
  *Disclaimer.* The authors declare that the results, discussions, and interpretations presented in this study are solely their own. The views expressed herein do not necessarily reflect those of their respective institutions or funding agencies.
acknowledgements: |
  *Acknowledgements.* We thank the Norwegian Public Roads Administration for sharing data from ground-penetrating radar surveys. We also thank Vikas Baranwal from the Geological Survey of Norway for helping us access the radiometric data from Skrim. This work contains data under the following licenses: (1) Creative Commons Attribution 4.0 International, © Kartverket, (2) *Norge digitalt* license, Norwegian Institute of Bioeconomy Research (NIBIO), © Geovekst, and (3) the Norwegian License for Public Data (NLOD), made available by the Geological Survey of Norway (NGU).
appendix: |
  \section{For submission}
  "Appendices: all material required to understand the essential aspects of the paper
  such as experimental methods, data, and interpretation 
  should preferably be included in the main text. 
  Additional figures, tables, as well as technical and theoretical developments 
  which are not critical to support the conclusion of the paper,
  but which provide extra detail and/or support useful for experts in the field
  and whose inclusion in the main text would disrupt the flow of descriptions or demonstrations 
  may be presented as appendices.
  These should be labelled with capital letters: Appendix A, Appendix B etc.
  Equations, figures and tables should be numbered as (A1), Fig. B5 or Table C6, respectively.
  Please keep in mind that appendices are part of the manuscript
  whereas supplements (see below) are published along with the manuscript."
  
  \section{Figures and tables in appendices}
  Please also sort the appendix figures and appendix tables into the respective appendix sections.
  They will be correctly named automatically.
  
  \section{Copernicus from Rmarkdown}
  **Please note:** Per [their guidelines](https://publications.copernicus.org/for_authors/manuscript_preparation.html),
  Copernicus does not support additional \LaTeX{} packages 
  or new \LaTeX{} commands than those defined in their `.cls` file. 
  This means that you cannot add any extra dependencies 
  and a warning will be thrown if so.
  **Important**: Always double-check with the official manuscript preparation guidelines 
  at [https://publications.copernicus.org/for_authors/manuscript_preparation.html](https://publications.copernicus.org/for_authors/manuscript_preparation.html),
  especially the sections "Technical instructions for LaTeX" and "Manuscript composition". 
  Please contact Daniel Nüst, `daniel.nuest@uni-muenster.de`, with any problems.
output:
  bookdown::pdf_book:
    keep_md: true
    base_format: rticles::copernicus_article # for using bookdown features like \@ref()
---

# Introduction

Introduction text goes here.
Read Gatis et al. [-@gatisMappingUplandPeat2019] and related work [@minasnyDigitalMappingPeatlands2019].

# Materials and methods

## Sites

We assessed how well we could predict peat depth at two sites with conspicuously different physical geography: Skrimfjella in eastern Norway and Ørskogfjellet in western Norway (Fig. \@ref(fig:sites)c).
These sites were chosen because they were covered by radiometric data from airborne surveys, relatively little built-up area, and road access.  

\begin{figure}
\includegraphics[height=0.9\textheight]{figures/sites-patchwork} \caption{Study areas at Ørskogfjellet (a) and Skrimfjella (b) within southern Norway (c). Land cover shown here is from the AR50 national land resource database and has simplified geometry with respect to the AR5 database used in the study.}(\#fig:sites)
\end{figure}

At Skrimfjella we delineated a study area of \unit{34\,km^{2}} based on radiometric coverage and accessibility (Fig. \@ref(fig:sites)b).
The study area has a diverse bedrock, with 32 % alkali feltspat granite, 26 % mergelstein, 10 % granite, and eight other types with > 1 % coverage (NGU, 1:250 000 dataset).
The landscape within our delineation is classified as *inland hills and mountains* [@simensenDiversityDistributionLandscape2021].
It is almost without human infrastructure, dominated by forest, and borders on a large nature reserve.
The study area has a mean elevation of \unit{438\,m} above sea level (range 223--711, IQR 351--509), and its mean slope at \unit{10\,m} resolution is 10.8° (IQR 4.6--15.1°).
In Norway's AR5 national land capability dataset [@ahlstromAR5Klassifikasjonssystem2019], \unit{1.5\,km^{2}} (4.5 %) of the study area is classified as mire --- defined as areas with mire vegetation and at least \unit{30\,cm} of peat depth.

At Ørskogfjellet we defined a study area of \unit{124\,km^{2}} which basically followed the footprint of the radiometric survey (Fig. \@ref(fig:sites)a).
According to the Geological Survey of Norway, bedrock in the area is 84 % granitic gneiss,  11 % granite, and 5 % aluminium silicate gneiss (NGU, 1:250 000 dataset). 
This study area comprises a wide range of major landscape types: *coastal plains*, *coastal fjord*, *inland valleys*, as well as *inland hills and mountains* [@simensenDiversityDistributionLandscape2021].
It is mostly forested, but also contains considerable farmland and open upland, and has several large lakes.
Its mean elevation is \unit{211\,m} above sea level (range 0--807, IQR 73--310), and its mean slope at \unit{10\,m} resolution is 13.0° (IQR 4.7--18.3°).
The AR5 dataset counts \unit{15.3\,km^{2}} (12.4 %) of the study area as mire.

## Peat depth measurements

At both study sites, our measurements of peat depth were made for the purpose of training a Random Forest model of peat depth, and we designed our sampling with this in mind [@brusSamplingDigitalSoil2019]. 
Broadly, we aimed for a sample that was representative of the predictor space defined by the most important predictors of peat depth [@wadouxSamplingDesignOptimization2019; @maComparisonConditionedLatin2020].
A sample that preserves the properties of the multivariate distribution of predictor and outcome variables is most likely to maintain any complex, non-linear relationships that exist in the population while avoiding spurious ones [@brusSamplingDigitalSoil2019].
We chose for our sampling and modelling a spatial resolution of \unit{10\,m}.
We considered this a reasonable compromise between digital terrain model (DTM) resolution (\unit{1\,m}) and small mires on the one hand, and airborne radiometric resolution (\unit{50\,m}) on the other.

### Skrimfjella

We measured peat depth in selected locations (\unit{10\,m} raster cells) at Skrimfjella.
The locations were chosen only from areas delineated as mire in the AR5 national land capability dataset.
Within this mire area, we stratified our sample across values of elevation, slope, and potassium ground concentration [from processed airborne gamma ray spectrometry, @baranwalHelicopterborneMagneticElectromagnetic2013].
Specifically, we used the *eSample* function in the *iSDM* R package (v.1.0) to chose an environmentally systematic sample.
This function defines the environmental space as a two-dimensional convex hull around the ordinated data, then creates a regular grid across that space, and lastly finds for each grid cell the datum that is nearest [@hattabUnifiedFrameworkModel2017].
Elevation was extracted from the \unit{10\,m} national DTM, slope calculated in degrees, and potassium ground concentration downscaled with bilinear resampling.
We set a target sample size of 100, excluded the top and bottom percentile from the convex hull, and with these parameters *eSample* returned 105 raster cells.

In addition to the peat depth locations, we had another arm of our sampling design for measuring peatland occurrence, as binary variable.
We wanted to measure peatland occurrence outside of mapped mire areas because the AR5 dataset is known to underestimate peatland coverage [especially in forests, @brynLandCoverNorway2018], and because airborne radiometrics may help identify unmapped peatland [@gatisMappingUplandPeat2019; @olearyDigitalSoilMapping2022].
The occurrence locations were sampled from the part of the study area that (1) was mapped as something other than mire in the AR5 database and (2) had a slope < 20°. 
We performed environmentally systematic sampling of this population with the same procedure as for the depth locations, and *eSample* returned 106 raster cells.

Field work at Skrimfjella was conducted in August 2020.
We navigated to the centers of the raster cells in the depth and occurrence samples by handheld GPS, checking that positional error was below \unit{3\,m}.
For each depth sample location, we measured peat depth three times (at the vertices of a triangle with \unit{2\,m} sides) to get a more representative value for the \unit{10\,m} raster cell, and to dampen the effect of outlying measurements [@parryEvaluatingApproachesEstimating2014].
We used a metal probe pushed downward until resistance indicated the base of the peat column.
Probe locations were adjusted up to \unit{20\,cm} if the base of the peat column seemed to be blocked by an obvious artifact.
For each occurrence sample location, we recorded the presence or absence of peatland --- primarily by digging and examining the top \unit{20\,cm} of soil (where this was possible).
We judged whether the soil was a peat soil based on its density, texture, and color.
Occasionally, when the soil itself was difficult to judge, we made our determination also based on the presence or absence of mire vegetation.
Although peat soil is strictly defined by organic content (which we did not analyse), we believe our protocol produced reasonable determinations of presence or absence that would generally satisfy most of the varying definitions of peatland [@minasnyMappingMonitoringPeatland2023].

Besides the depth and occurrence measurements described above, we also measured peat depth in three subjectively-chosen, individual mires, using ground-penetrating radar (GPR).
We used the Malå ProEx GPR system (Guideline Geo AB, Sweden) with its \unit{500\,MHz} shielded antenna mounted in a plastic sledge, and its control unit connected to a GNSS receiver.
At each of the three mires we recorded GPR traces along walking transects that covered the extent of the mire, mostly in traversing, zigzag patterns with between \unit{5\,m} and \unit{20\,m} spacing at their vertices.
Along the GPR transects we also probed peat depth at marked trace locations, to be able to calibrate the GPR wave speed velocity.
We processed the GPR data with Reflex2DQuick software (v.3.0, Sandmeier Scientific Software, Germany), applying a time-zero correction, a dewow filter, and a gain filter based on observed energy decay.
Then we picked strong reflectors in the radargrams that we interpreted as the base of the peat column.
We used picks at marked trace locations to calibrate wave speed velocity; we pooled calibration points across the three mires and fitted a linear regression of depth on one-way travel time with the intercept fixed at zero.
In total we had 46 calibration points along \unit{3.5\,km} of GPR transects.
Finally, we used the calculated wave velocity (\unit{0.0387\,m\,ns^{-1}}, $R^2 = 0.874$) to convert the travel times of all picks to calibrated peat depths.

### Ørskogfjellet

At Ørskogfjellet we also measured peat depth in a sample of \unit{10\,m} raster cells, selected from the part of the study area classified in the AR5 dataset as mire.
Before selecting locations, we determined a minimal sample size that would adequately capture the terrain and radiometric properties of the entire mire area.
Specifically, we aimed to identify the size at which adding locations produced diminishing decreases in divergence between sample and population distributions --- i.e. the elbow point in a curve of similarity between sample and population [@maloneMethodsImproveUtility2019].
This approach has been found to identify sample sizes that correspond with diminishing returns in predictive model performance on external evaluation data [@sauretteDivergenceMetricsDetermining2023].
We defined a sequence of sample sizes (50--500) and for ten replicate samples at each size [drawn by conditioned latin hypercube sampling, @minasnyConditionedLatinHypercube2006; @roudierClhsPackageConditioned2011], we calculated the mean Kullback-Leibler divergence between sample and population distributions [@maloneMethodsImproveUtility2019; @sauretteDivergenceMetricsDetermining2023].
The variables in the divergence calculation were terrain slope and four radiometrics: potassium, thorium, uranium, and total count.
Next, we fitted a asymptotic regression of mean divergence on sample size, and identified the sample size at which the curve reached 95 % of the fitted asymptote.
Through this procedure we found that we could adequately capture the population distribution with a sample of 160 locations.

To choose 160 locations, we performed feature space coverage sampling.
This approach has been found to produce higher accuracy in Random Forest models than conditioned latin hypercube sampling [@wadouxSamplingDesignOptimization2019; @maComparisonConditionedLatin2020].
Feature space coverage sampling aims to disperse samples as uniformly as possible in multidimensional predictor space, and is implemented by choosing locations that are closest to cluster centers in a k-means clustering of the standardized predictor space [@brusSamplingDigitalSoil2019].
Feature space coverage sampling works best when all dimensions are important predictors of the outcome [@wadouxSamplingDesignOptimization2019], and we used the same five predictors that we used to choose sample size: terrain slope and four radiometrics.
The radiometric variables were downscaled to \unit{10\,m} resolution with cubic B-spline resampling in QGIS [v.3, @QGISsoftware].
We adjusted the feature space coverage sampling to ensure that locations were accessible within time constraints, and assessed how this changed our sample from an ideal feature space coverage sample.
Adjusting for accessibility is justified because the smaller sample size that would result if accessibility were ignored can degrade model accuracy as much or more as deviations from ideal sampling designs [@wadouxSamplingDesignOptimization2019; @maComparisonConditionedLatin2020].
To adjust, we first restricted the sampling population to mire areas that were within an arbitrary cost distance of publicly-accessible roads.
Cost distance was calculated using GRASS's *r.walk* function, with friction costs defined by AR5 land classes [@GRASSv8-2].
After creating a feature space coverage sample with this restriction, we also inspected a map of the sample and substituted 16 inaccessible locations with accessible locations from the same or a nearby cluster.
Our two accessibility adjustments increased the distance in standardized predictor space between sample locations and cluster centers by 78 % (with respect to the ideal sample), but distance in our sample was still only 46 % of the mean distance to cluster centers --- i.e., accessibility did not force locations far from cluster centers relative to the size of the clusters.

Field work at Ørskogfjellet was conducted in August 2023.
We navigated to the centers of the raster cells in the sample using real time kinematic differential GNSS (Topcon Positioning Systems, USA), to ensure sub-meter positional accuracy.
At each location we measured peat depth three times by manual probing, with probe locations spaced approximately \unit{2.5\,m} apart.
In areas with dense sampling locations, we also measured peat depth with GPR along snaking transects passing through the centers of the sampling cells (seven transects, \unit{6.2\,km} total length).
We used the same GPR system as at Skrimfjella, but with a \unit{100\,MHz} Malå rough terrain antenna (Guideline Geo AB, Sweden) at some transects.
To navigate the GPR transects, we placed flags at the cell centers of the sample locations, and used a handheld GNSS receiver to guide the GPR operator.
At sampling locations crossed by a GPR transect, we arranged the manual probe positions along the transect (for better calibration of the GPR wave speed velocity), while other locations were probed in a triangular pattern around the cell center like at Skrimfjella.

We processed the GPR data with Reflexw software (v.8.5, Sandmeier Scientific Software, Germany), applying a dewow filter, time-zero correction, bandpass filter, gain filter, and a dynamic correction that accounts for the non-vertical wave path between offset transmitter and receiver antennae.
The last correction is important for the rough terrain antenna, which has an antenna separation (\unit{2.2\,m}) --- comparable to typical peat depths.
As with the data from Skrimfjella, we picked the base of the peat column from strong reflectors in the radargram, and calibrated wave velocity with manual probe measurements in a linear regression.
The points in the regression were created by joining to each probe measurement the travel time of the nearest pick, but only if these were within \unit{2\,m} of each other.
In total we had 78 calibration points along \unit{7.8\,km} of interpretable GPR traces (transect length exceeded because of extra GPR data).
Finally, we used the calculated wave velocity (\unit{0.0427\,m\,ns^{-1}}, $R^2 = 0.946$) to convert the travel times of all picks to calibrated peat depths.

We also used two sets of existing depth measurements from Ørskogfjellet.
The first set was provided by the Norwegian Public Roads Administration, who commissioned GPR surveys of particular peatland areas in 2020 and 2021.
The surveys were conducted with a dual channel system (\unit{70\,MHz} and \unit{300\,MHz}; ImpulseRadar AB, Sweden), connected to GNSS with CPOS correction.
We used interpreted and calibrated traces from these surveys, and discarded some data where multiple depths were interpreted for the same locations.
This summed to \unit{7.4\,km} of interpreted traces.
The second set of existing depth data we extracted from a paper map made by the Norwegian Soil and Mire Company in 1984.
This map presents 44 borehole depths (in decimeters) across a \unit{9\,ha} peatland area.
We georeferenced the map and digitized the borehole locations and depths.

## Peat depth predictors

We created the same suite of 25 quantitative peat depth predictors for both sites (Table \@ref(tab:preds)).
All quantitative predictors were derived either from an airborne radiometric survey or from a DTM.
From the radiometric surveys we simply used the four variables produced by the surveyors (Geological Survey of Norway): ground concentration of Potassium, Thorium, Uranium, as well as total count.
From the DTMs we calculated several land surface parameters, ranging from simple terrain indices to more complex geomorphometric and hydrological variables [@maxwellLandsurfaceParametersSpatial2022].
\begin{table}

\caption{(\#tab:preds)Quantitative predictors of peat depth.}
\centering
\begin{tabular}[t]{l|l}
\hline
name & description\\
\hline
radK & Potassium ground concentration\\
\hline
radTh & Thorium ground concentration\\
\hline
radU & Uranium ground concentration\\
\hline
radTC & Total count of gamma radiation\\
\hline
elevation & Mean elevation\\
\hline
slope1m & Mean of 1 m slope\\
\hline
TPI1m & Mean of 1 m topographic position index\\
\hline
TRI1m & Mean of 1 m terrain ruggedness index\\
\hline
roughness1m & Mean of 1 m roughness\\
\hline
slope10m & 10 m slope\\
\hline
TPI10m & 10 m topographic position index\\
\hline
TRI10m & 10 m terrain ruggedness index\\
\hline
roughness10m & 10 m roughness\\
\hline
MRVBF & Multi-resolution valley bottom flatness\\
\hline
TWI5m & Mean of 5 m topographic wetness index\\
\hline
TWI10m & 10 m topographic wetness index\\
\hline
TWI20m & Bilinear interpolation of 20 m topographic wetness index\\
\hline
TWI50m & Bilinear interpolation of 50 m topographic wetness index\\
\hline
DTW2500 & Depth-to-water index, flow initiation area of 0.25 ha\\
\hline
DTW5000 & Depth-to-water index, flow initiation area of 0.5 ha\\
\hline
DTW10000 & Depth-to-water index, flow initiation area of 1 ha\\
\hline
DTW20000 & Depth-to-water index, flow initiation area of 2 ha\\
\hline
DTW40000 & Depth-to-water index, flow initiation area of 4 ha\\
\hline
DTW80000 & Depth-to-water index, flow initiation area of 8 ha\\
\hline
DTW160000 & Depth-to-water index, flow initiation area of 16 ha\\
\hline
\end{tabular}
\end{table}

The radiometric survey covering Skrimfjella was conducted in 2008--2011. 
The survey was flown at an average altitude of \unit{75\,m} and average speed of \unit{108\,km\,h^{-1}}, with flight lines spaced \unit{200\,m} apart.
Spectrometer count rates were calibrated annually to known concentrations of Potassium, Thorium, and Uranium in mobile pads.
The Geological Survey of Norway processed data from the spectrometer following standard procedures outlined by the International Atomic Energy Association, and the processing included: correction for aircraft and cosmic background radiation, correction for radon in the air, window stripping of the gamma ray spectrum, correction for flying height, conversion of count rates to ground concentrations, and finally gridding to \unit{50\,m} resolution with micro-leveling.
Further details about the survey and data processing are provided in Baranwal et al. [-@baranwalHelicopterborneMagneticElectromagnetic2013].
We downscaled the processed data to \unit{10\,m} resolution by cubic spline resampling, using the *terra* package in R. 

A very similar radiometric survey covering Ørskogfjellet was conducted in December 2014 and January 2015. 
This survey was flown at an average altitude of \unit{80\,m} and average speed of \unit{88\,km\,h^{-1}}, with flight lines also spaced \unit{200\,m} apart.
Spectrometer count rates were calibrated in 2013 to known concentrations of Potassium, Thorium, and Uranium in mobile pads.
The spectrometer data were processed following the same procedure as for the survey at Skrimfjella, except that a convolution filter was added to smooth the gridded data.
Further details about the survey and data processing are provided in Ofstad [-@ofstadHelicopterborneMagneticRadiometric2015].
The \unit{10\,m} resolution predictors from this survey were identical to the layers used in the sampling design (resampled from \unit{50\,m} resolution with cubic B-splines).

For terrain-derived variables, we obtained \unit{1\,m} resolution DTMs from the Norwegian Mapping Authority.
The DTM for Skrimfjella was produced from airborne laser scanning surveys in 2015 and 2022, with laser point density of \unit{5\,pts\,m^{-2}}.
For Ørskogfjellet, the DTM was produced from a 2015 survey with \unit{2\,pts\,m^{-2}}.
Where necessary, DTMs were resampled to the coordinate reference system of the radiometric data.

We used the *terra* R package to calculate from the DTMs: slope, topographic position index (difference from mean of eight neighbors), terrain ruggedness index (mean of absolute differences from eight neighbors), and roughness (range in the nine-cell neighborhood).
These variables were derived at two scales to produce different predictors; we either calculated from \unit{1\,m} DTM resolution and then aggregated, or aggregated to \unit{10\,m} DTM resolution and then calculated the indices.
This kind of multiscale feature engineering of land surface parameters has been found to improve machine learning predictions of soil properties [@millerImpactMultiscalePredictor2015; @dornikOptimalScalingPredictors2022; @newmanAssessingSpatiallyHeterogeneous2023].
We chose \unit{1\,m} and \unit{10\,m} resolutions because we know that peat depth tends to vary at fine scales in Norway [@maxwellLandsurfaceParametersSpatial2022].
We also calculated the the Multi-Resolution Valley Bottom Flatness index, which indicates the degree of valley bottom flatness at a given location via a multiscale algorithm [@gallantMultiresolutionIndexValley2003].
We calculated this index in SAGA GIS [v.9.3.2, @conradSystemAutomatedGeoscientific2015] with default parameters.

Next, we calculated the Topographic Wetness Index [@quinnPredictionHillslopeFlow1991].
This index is notoriously scale-dependent and often matches real hydrological conditions best when calculated from moderate-to-coarse resolution DTMs [@agrenEvaluatingDigitalTerrain2014; @riihimakiTopographicWetnessIndex2021], so we calculated it from \unit{5\,m}, \unit{10\,m}, \unit{20\,m}, and \unit{50\,m} DTM resolution.
The calculations were performed with Whitebox software [@lindsayWhiteboxGATCase2016], accessed through the *whitebox* R package [v2.4, @wuWhiteboxWhiteboxToolsFrontend2022].
We filled depressions in the DTM with the algorithm in Wang & Liu [-@wangEfficientMethodIdentifying2006], and used the deterministic infinity flow accumulation algorithm [@tarbotonNewMethodDetermination1997].

The last terrain-based predictor we included was the depth-to-water index [@murphyMappingWetlandsComparison2007].
This index approximates a location's vertical height above the surface water feature that it is likely to drain towards.
It is calculated as the minimum cumulative slope (scaled by cell size) to a surface water feature [eq. 5 in @murphyTopographicModellingSoil2009].
We calculated unitless slope from the \unit{1\,m} DTM using the Whitebox software.
Also using Whitebox, we defined surface water features from the DTM by filling depressions and then calculating flow accumulation to define catchment areas for each cell [@schonauerSpatiotemporalPredictionSoil2021; @schonauerRcodeCalculatingDepthwater2021].
This catchments area layer was then thresholded at seven different levels (*flow initiation areas* from \unit{0.25\,ha} to \unit{16\,ha}) to estimate surface water features under moisture scenarios varying from wet to dry [@murphyModellingMappingTopographic2011; @agrenEvaluatingDigitalTerrain2014; @schonauerSpatiotemporalPredictionSoil2021].
In addition, all surface water features mapped in the AR5 dataset were also transferred to the raster layer.
For each of the seven surface water layers, we derived the depth-to-water index using the *Distance Accumulation* tool in ArcGIS Pro (v.3.1, ESRI, USA), which has an efficient algorithm to find the cumulative distance over a cost surface to the least-cost source.

In addition to the quantitative predictors described above, we prepared one categorical predictor --- peat depth class --- from a historical national map dataset called *DMK* [@ahlstromAR5Klassifikasjonssystem2019].
The DMK peat depth classes are: < \unit{1\,m} (*shallow*), > \unit{1\,m} (*deep*), and *unknown*.
This dataset stems from field measurements made in 1964--2007 as part of a wider land cover mapping in Norway.
The mapping was primarily for identifying agricultural and silvicultural resources, so it covers productive areas below the tree line and has greater coverage of peat depth class in peatlands judged to be arable or afforestable [@bjordalMarkslagsklassifikasjonOkonomiskKartverk2007; @ahlstromAR5Klassifikasjonssystem2019].
Mappers generally assigned peat depth classes to polygons of at least \unit{0.5\,ha}, although delineating polygons down to \unit{0.2\,ha} was allowed if peat depth showed a "particularly marked difference" [@bjordalMarkslagsklassifikasjonOkonomiskKartverk2007].
We rasterized the peat depth class attribute to our \unit{10\,m} grid.

## Predictive models of peat depth 

We used Random Forests (RF) to predict peat depth at both sites.
RF is a tree-based ensemble machine learning algorithm that builds many decision trees on bootstrapped samples of the training data, randomly subsets predictors in the trees, and averages the predictions of the trees [@breimanRandomForests2001].
We chose RF because it can handle complex interactions between predictors, is robust to overfitting, and generally shows higher performance in digital soil mapping applications than other algorithms [@beguinPredictingSoilProperties2017; @nussbaumEvaluationDigitalSoil2018; @lamichhaneDigitalSoilMapping2019].
It is suited for use on relatively small training data sets and its predictions can be interrogated to learn about predictor importance [@khaledianSelectingAppropriateMachine2020].
Evaluating variable importance in a maximally-predictive model aligns with the aim of this study.

RF by itself is not a spatial model, and it will only predict spatial structure in the outcome to the degree that the structure is captured by predictors.
We considered using regression kriging --- a hybrid between non-spatial and spatial techniques that would be achieved by adding to the RF predictions a geostatistically-interpolated surface of RF residuals [@henglGenericFrameworkSpatial2004].
The spatial component in regression kriging often improves map accuracy compared to a non-spatial model [@beguinPredictingSoilProperties2017; @lamichhaneDigitalSoilMapping2019; @mollaMachineLearningGeostatistical2023], but it can do so only if the spatial autocorrelation range in the non-spatial residuals is large compared to distances between samples and prediction locations [@henglGenericFrameworkSpatial2004; @szaboMappingSoilHydraulic2019; @takoutsingComparingPredictionPerformance2022].
If the outcome varies at fine scales and the samples are clustered in small parts of the study area, a spatial component will hardly improve overall map accuracy.
We used semivariograms to assess the spatial structure in the residuals of the RF predictions, and found that (non-spatial) RF rather regression kriging was justified at both sites. 

For each site we followed the same modelling workflow.
First we trained a RF to predict peat depth based on our 25 quantitative predictors, to simulate applications where radiometric and terrain data are available, but no maps of peat depth.
Next we trained a RF with all 26 predictors (DMK peat depth class included), to fully leverage the relevant available data and check the added value of the national DMK dataset.
Finally, we fit a simple linear model to only DMK depth class, to provide a fair comparison between the accuracy of the RF models and the existing national map of peat depth calibrated on the same data.
We implemented models in the *tidymodels* framework in R [@kuhnTidymodelsCollectionPackages2020], with the *ranger* R package for RFs [v.0.16, @wrightRangerFastImplementation2017].
RFs were fit with 1000 trees, minimum node size of 5, and the number of predictors randomly sampled at each split was the square root of the total number of predictors (*ranger* default).
We did not tune these hyperparameters because RFs are relatively insensitive to tuning [@probstHyperparametersTuningStrategies2019], and because it would require nested spatial cross-validation to prevent data leakage [@schratzHyperparameterTuningPerformance2019].

The performance of digital soil mapping must be evaluated with reference to a specific purpose [i.e., map vs. model validation, interpolation vs. extrapolation, @robertsCrossvalidationStrategiesData2017; @milaNearestNeighbourDistance2022], and here we aimed to evaluate maps covering the study areas.
In the absence of additional field work to collect a design-based independent validation set, we used a spatial cross-validation scheme to evaluate model performance [@wadouxSpatialCrossvalidationNot2021; @meyerMachineLearningbasedGlobal2022].
Specifically, we used k-Means Nearest Neighbor Distance Matching (kNNDM), which creates cross-validation folds that mimic the spatial prediction task that is defined as the goal [@linnenbrinkKNNDMCVKfold2024].
In particular, kNNDM looks for the spatial assignment of training data to folds that minimizes the difference between two distributions: nearest neighbor distances between training and test locations in the cross-validation, and nearest neighbor distances between training and prediction locations for the model.
That way, the spatial separation between folds is similar to the separation between training and prediction locations --- which increases the quality of the map accuracy estimate [@linnenbrinkKNNDMCVKfold2024].
For spatially clustered training data, this approach strikes a balance between the risk of optimistic metrics from random cross-validation and the risk of pessimistic metrics from other forms of spatial cross-validation [@wadouxSpatialCrossvalidationNot2021].
We implemented the kNNDM with the *CAST* R package [v.1.0.2, @meyerCASTPackageTraining2024], setting prediction locations to all AR5 mire cells in the study area, and choosing a number of folds (5--20) that produced the best match between the two NND distributions.
From the cross-validation we quantified *root mean squared error* (accuracy), *R^2^* (correlation), and *Lin's concordance correlation coefficient* (accuracy + correlation).

We quantified global variable importance for the RFs with the 25 quantitative predictors.
Global variable importance measures the influence of a given predictor on the output of the model, aggregated across all locations.
We calculated variable importance with the *vip* R package (v.0.4.1), by three different methods: *FIRM*, *permutation*, and *Shapley* [@greenwellVariableImportancePlots2020].
*Permutation* values were aggregated from ten iterations with RMSE as the performance measure.

# Results

Include a 12cm width figure of Nikolaus Copernicus from [Wikipedia](https://en.wikipedia.org/wiki/File:Nikolaus_Kopernikus.jpg) with caption using R Markdown (Fig. \@ref(fig:portrait)).

\begin{figure}
\includegraphics[width=8.3cm]{Nikolaus_Kopernikus} \caption{one column figure}(\#fig:portrait)
\end{figure}

## Tables

You can add \LaTeX table in an R Markdown document to meet the template requirements (Table \ref{tab:latextable}).


\begin{table}[t]
\caption{TEXT}
\begin{tabular}{l c r}
\tophline

a & b & c \\
\middlehline
1 & 2 & 3 \\

\bottomhline
\end{tabular}
\belowtable{Table Footnotes}
\label{tab:latextable}
\end{table}

Or you can use markdown to create the table with booktabs = FALSE (<https://github.com/rstudio/rticles/issues/558#issuecomment-1907981541>).

See Table \@ref(tab:test).

\begin{table}

\caption{(\#tab:test)My caption}
\centering
\begin{tabular}[t]{l|r|r|r}
\hline
  & mpg & cyl & disp\\
\hline
Mazda RX4 & 21.0 & 6 & 160\\
\hline
Mazda RX4 Wag & 21.0 & 6 & 160\\
\hline
Datsun 710 & 22.8 & 4 & 108\\
\hline
\end{tabular}
\end{table}

# Discussion

Lorem ipsum dolor sit amet, consectetur adipiscing elit.
Sed do eiusmod tempor incididunt ut labore et dolore magna aliqua.
Ut enim ad minim veniam, quis nostrud exercitation ullamco laboris nisi ut aliquip ex ea commodo consequat.
Duis aute irure dolor in reprehenderit in voluptate velit esse cillum dolore eu fugiat nulla pariatur.
Excepteur sint occaecat cupidatat non proident,
sunt in culpa qui officia deserunt mollit anim id est laborum.

# Conclusions

Nulla facilisi.
Maecenas vel nunc nec purus tincidunt congue.
Proin auctor, lectus eu pharetra malesuada,
nisi nunc bibendum nunc,
eget tincidunt nunc nisi id nunc.
Sed euismod, nunc sit amet aliquam tincidunt,
nunc nunc tincidunt nunc,
nec tincidunt nunc nunc nec nunc.
Donec auctor, nunc sit amet aliquam tincidunt,
nunc nunc tincidunt nunc,
nec tincidunt nunc nunc nec nunc.
