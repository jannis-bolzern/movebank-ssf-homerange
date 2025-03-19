# Proposal for Semester Project


<!-- 
Please render a pdf version of this Markdown document with the command below (in your bash terminal) and push this file to Github. Please do not Rename this file (Readme.md has a special meaning on GitHub).

quarto render Readme.md --to pdf
-->

**Patterns & Trends in Environmental Data / Computational Movement
Analysis Geo 880**

| Semester:      | FS25                                     |
|:---------------|:---------------------------------------- |
| **Data:**      | Selected animal tracking dataset from the Movebank database  |
| **Title:**     |                 |
| **Student 1:** | Jannis Bolzern                        |
| **Student 2:** | Elke Michlmayr                        |

## Abstract 
<!-- (50-60 words) -->
In our student project we aim to demonstrate that wild animal movement patterns differ in rural and remote areas based on data from the Movebank animal tracking database. We selected three publicly available datasets with location data from red fox, bobcat, and coyote movements in England, Canada, and the US for comparison and analysis.


## Research Questions
<!-- (50-60 words) -->
Two hypotheses will be tested: (1) animals with access to anthropogenic food sources have shorter home ranges, and (2) animal mortality rates are linked to the human footprint index of their home range area. To prove these, an assessment of home ranges for individual animals will be performed and linked to the human activity found in the respective area.


## Results / products
<!-- (50-100 words) -->
<!-- What do you expect, anticipate? -->
A Quarto document with related data analysis and description will be produced and published. The data analysis part will be hosted on Github under a Creative Commons license.


## Data
<!-- (100-150 words) -->
<!-- What data will you use? Will you require additional context data? Where do you get this data from? Do you already have all the data? -->

Human footprint data:
We will use the global 100 meter resolution terrestrial human footprint data from here. The data can be read in Python, R, or any other script that has libraries that can interpret geospatial data (such as folium). It is described in this publication: Gassert F., Venter O., Watson J.E.M., Brumby S.P., Mazzariello J.C., Atkinson S.C. and Hyde S., An Operational Approach to Near Real Time Global High Resolution Mapping of the Terrestrial Human Footprint. Front. Remote Sens. 4:1130896 doi: 10.3389/frsen.2023.1130896 (2023)

We will use the following animal data:
Red fox data for rural areas (GPS-based, Wiltshire, UK) 
Red fox data for remote areas (Argos-based, from Bylot island and from Herschel island, Canada)
Bobcat and Coyote data for remote areas with some rural structures (GPS-based, northern Washington, US) 
These datasets were created and published on Movebank under a Creative Commons license by three different research groups, and described in the following publications:
Porteus TA, Short MJ, Hoodless AN, Reynolds JC. 2024. Movement ecology and minimum density estimates of red foxes in wet grassland habitats used by breeding wading birds. Eur J Wildlife Res. 70:8. https://doi.org/10.1007/s10344-023-01759-y
Sandra Lai, Chloé Warret Rodrigues, Daniel Gallant, James D Roth, Dominique Berteaux, Red foxes at their northern edge: competition with the Arctic fox and winter movements, Journal of Mammalogy, Volume 103, Issue 3, June 2022, Pages 586–597, https://doi.org/10.1093/jmammal/gyab164
Laura R. Prugh et al., Fear of large carnivores amplifies human-caused mortality for mesopredators. Science 380,754-758(2023). DOI:10.1126/science.adf2472

## Analytical concepts
<!-- (100-200 words) -->
<!-- Which analytical concepts will you use? What conceptual movement spaces and respective modelling approaches of trajectories will you be using? What additional spatial analysis methods will you be using? -->
Trajectory analysis 
Assessment of home ranges

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
Argos tracking and GPS tracking might not be easily comparable. A first analysis reveals that they both provide Latitude and Longitude information so a comparison seems possible in principle.
Movement patterns on (small?) islands might be fundamentally different and make a comparison between remote and rural infeasible.
There could be other reasons for differences in home range sizes that have nothing to do with the patterns we are investigating.
Parsing the global 100 meter resolution terrestrial human footprint data might be tricky.
It’s not clear if the human footprint index for the home range is the relevant factor. It could be that the human footprint index for the animal’s final location is more relevant.
There seem to be some data issues in the mortality data where the paper and data do not seem to have the exact same numbers. It’s not clear if the data is complete.
There could be other reasons for differences in mortality rates that have nothing to do with the patterns we are investigating.

Plan B
Based on first investigations we are confident that the home range comparison based analysis should be doable. We expect the mortality rate analysis to be a much more complex topic and see it as a stretch goal.

## Questions? 
<!-- (100-150 words) -->
<!-- Which questions would you like to discuss at the coaching session? -->

Are you happy with this plan?
Do you think the level of ambition is appropriate?
Do you see any obstacles?

## Bibliography
