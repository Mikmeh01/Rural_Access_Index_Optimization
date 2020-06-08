## Rural_Access_Index_Optimization

This project was used for the [Paper](https://openknowledge.worldbank.org/bitstream/handle/10986/31309/WPS8746.pdf?sequence=5) that aims to quantify how much it will cost for developing countries to increase their Rural Access Index.

## Getting Started

Clone the repo, change the ISO3 country parameter (BLZ for Belize for example) of the RAI_Opti notebook and run the remaining code

## Output

The code exports a shapefile of roads added throughout the process with the corresponding increase of RAI at each iteration and the cost for upgrading the road (making it all-weather road)

## Installation

The requirement file lists all the necessary libraries used in this project

## References

Datesets used:
-  Urban Areas : Global Rural-Urban Mapping Project [GRUMP](https://sedac.ciesin.columbia.edu/data/collection/grump-v1)
- Population: [WorldPop Database](https://www.worldpop.org/)
- Roads : Open Street Map

## Issues

An error might occur at the first block below Country Roads if a new category of roads is added OpenStreetMap contributors. In this case, the user should add this new category inside the "map_roads" function to the functions python file
