# Proposal for Semester Project


<!-- 
Please render a pdf version of this Markdown document with the command below (in your bash terminal) and push this file to Github. Please do not Rename this file (Readme.md has a special meaning on GitHub).

quarto render Readme.md --to pdf
-->

**Patterns & Trends in Environmental Data / Computational Movement
Analysis Geo 880**

| Semester:      | FS25                                     |
|:---------------|:---------------------------------------- |
| **Data:**      | Selected datasets from the Movebank animal tracking database  |
| **Title:**     | How does human activity affect the movement patterns of wild animals?   |
| **Student 1:** | Jannis Bolzern                        |
| **Student 2:** | Elke Michlmayr                        |

## Abstract 
<!-- (50-60 words) -->
In this project, we aim to investigate how human activity influences the movement patterns of wild animals. Using GPS and Argos tracking data from red foxes, bobcats, and coyotes across rural and remote areas in England, Canada, and the US, we will analyze home range sizes, temporal activity shifts, and habitat selection in relation to human footprint and land use data.

## Research Questions
<!-- (50-60 words) -->
1. **Access to Anthropogenic Food Sources:** Do animals in high human-impact areas exhibit smaller home ranges due to easier access to food resources?
2. **Temporal Shifts in Activity Patterns:** Do animals become more nocturnal in high human-impact areas to avoid direct human encounters?
3. **Habitat Selection in Human-Dominated Landscapes:** How do animals select habitats (e.g., forests, agriculture) under varying levels of human influence?

## Results / products
<!-- (50-100 words) -->
<!-- What do you expect, anticipate? -->
The project will produce:

1. A **research paper** in Quarto format, including an abstract, methodology, results, visualizations, and discussion.
2. **R code** for data processing, analysis, and visualization, hosted on GitHub under a Creative Commons license.
3. **Interactive maps** and visualizations of animal movements, home ranges, and habitat selection in relation to human footprint and land use data.

## Data
<!-- (100-150 words) -->
<!-- What data will you use? Will you require additional context data? Where do you get this data from? Do you already have all the data? -->

We will use the following wild animal tracking data published on Movebank under a Creative Commons license by three different research groups (see references):

* **Red fox** data for rural areas (GPS-based, Wiltshire, UK)
* Red fox data for remote areas (Argos-based, from the uninhabitated islands Bylot and Herschel, Canada)
* **Bobcat** and **coyote** data for remote areas with some rural structures (GPS-based, northern Washington, US)

For the human footprint data, we will use the **global 100 meter resolution terrestrial human footprint data (HFP-100)** by Joe Mazzariello et al. The data can be read in Python, R, or any other script that has libraries that can interpret geospatial data (such as folium).

If required, for **land use** in Washington, US, we will rely on the General Land Use Final Dataset published by Washington Spatial Data ([link](https://geo.wa.gov/datasets/a0ddbd4e0e2141b3841a6a42ff5aff46_0/explore?location=48.347066%2C-118.420235%2C9.91)) and for land use in the UK on gov.uk data ([link](https://www.data.gov.uk/dataset/946ce540-de76-441e-bac8-624f30cace8a/land-cover-map-2021-10m-classified-pixels-gb)).

## Analytical concepts
<!-- (100-200 words) -->
<!-- Which analytical concepts will you use? What conceptual movement spaces and respective modelling approaches of trajectories will you be using? What additional spatial analysis methods will you be using? -->
1. **Trajectory Analysis:** Movement paths will be analyzed to identify patterns in speed, direction, and habitat use. Step lengths and turning angles will help infer behavioral states. Step-selection functions (SSFs) will model habitat preferences relative to availability, allowing us to quantify how animals respond to environmental covariates such as human footprint and land use.

2. **Home Range Assessment:** Home range sizes will be calculated using kernel density estimation (KDE) and minimum convex polygons (MCPs). This will provide estimates of the area used by each individual. Home ranges will be compared across regions with varying levels of human activity to test whether access to anthropogenic food sources (e.g., garbage, crops) leads to smaller home ranges due to resource concentration.

3. **Temporal Activity Patterns:** Movement rates will be used to quantify diel activity shifts. We will test whether animals in high human-impact areas exhibit increased nocturnality, potentially as a strategy to avoid direct human encounters. This analysis will reveal how temporal behavior adapts to human presence.

4. **Habitat Selection:** SSFs will quantify selection for human-modified habitats (e.g., agricultural areas, urban edges) relative to natural habitats. By comparing selection patterns across species and regions, we will assess how habitat preferences vary with human influence.


## R concepts
<!-- (50-100 words) -->
<!-- Which R concepts, functions, packages will you mainly use. What additional spatial analysis methods will you be using? -->
We will use the following libraries:

* readr, tidyr, dplyr library for data processing
* ggplot2 and tmap library for visualization
* sf library for spatial data handling
* move2 library for trajectory handling
* folium for reading the human footprint data

We are not yet familiar with using folium.

## Risk analysis
<!-- (100-150 words) -->
<!-- What could be the biggest challenges/problems you might face? What is your plan B? -->
* We have no experience with the human footprint index data and also not the land use data. We have not yet tried to parse it. We'll look into it soon to understand the risk better, and look for alternative datasets if necessary.
* We have a couple of research questions and might not get to all of them. We have ordered them by expected degree of complexity.
* We might run into issues caused by comparing datasets from different sources (such as the Wiltshire and the Canadian studies).
* There might be other unrelated influences on the animals home range choices, nocturnal patterns, and habitat selections that have nothing to do with human activity and are influencing our results. 

## Questions? 
<!-- (100-150 words) -->
<!-- Which questions would you like to discuss at the coaching session? -->

* Are you happy with this plan in general?
* Do you think the level of ambition is appropriate? Are we aiming too low or too high?
* Do you see any obstacles? What are we missing?
* Do you agree with the expected degree of complexity ranking that we chose?
* If it were preferable to focus on fewer research questions, which ones would be the most feasible?
* We don’t have a plan B since we think we can address at least two of the research questions. Is that acceptable? 
* How likely are we going to run into the last two risks mentioned in the risk analysis?

## References
* Porteus TA, Short MJ, Hoodless AN, Reynolds JC. 2024. Movement ecology and minimum density estimates of red foxes in wet grassland habitats used by breeding wading birds. Eur J Wildlife Res. 70:8. https://doi.org/10.1007/s10344-023-01759-y
* Sandra Lai, Chloé Warret Rodrigues, Daniel Gallant, James D Roth, Dominique Berteaux, Red foxes at their northern edge: competition with the Arctic fox and winter movements, Journal of Mammalogy, Volume 103, Issue 3, June 2022, Pages 586–597, https://doi.org/10.1093/jmammal/gyab164
* Laura R. Prugh et al., Fear of large carnivores amplifies human-caused mortality for mesopredators. Science 380,754-758(2023). DOI:10.1126/science.adf2472
* Gassert F., Venter O., Watson J.E.M., Brumby S.P., Mazzariello J.C., Atkinson S.C. and Hyde S., An Operational Approach to Near Real Time Global High Resolution Mapping of the Terrestrial Human Footprint. Front. Remote Sens. 4:1130896 doi: 10.3389/frsen.2023.1130896 (2023)
