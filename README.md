# Preserve ocean gateways and landbridges with skeletons

Paleoclimate models often require the use of low resolution boundary conditions to facilitate the high computational costs of modelling. This can create problems when coarsening elevation maps from high to low resolution. When this is done, land and sea connections can be inadvertently created or lost. For example, the Panama Isthmus is a thin strip of land that connects North and South America. If an elevation map is coarsened by taking the mean or median of the elevations, the resulting map will not have a land bridge between North and South America. 

![skeletonImplementation](https://github.com/JonathonLeonard/ocean-Gateways/blob3c60b4bf007917190fdf9729bed85c55707cfff7/skeletenImplementation.png)

Ocean gateways and landbridges have been shown to have a large influence on climate and oceans, therefore we need to do our best to preserve key ocean and land connections. For present day elevation maps, or paleoelevations of a single age snapshot, this can be done manually. However, for paleoclimate modelling over long time series, this can become burdensome. 

This repository provides a solution to the problem of losing gateway and landbridge information when coarsening elevation grids. It uses skeletonisation to trace contiguous land and ocean connections then enforces those connections on a coarsened grid.

Included in this repository is a workflow that creates skeletons and imposes them for 26 paleoelevation grids from 250 Ma to the present. The outline of the workflow is:

* Loop over the high resolution elevation maps. These are provided as netcdf grids in 'paleoElevation'
* Turn the elevation grid into a land-sea mask
* Identify water bodies with a flood fill algorithm and remove those that are too small (this number can be adjusted in the code)
* Remove small islands with the same method as above (again, this number can be adjusted by setting the minimum). Islands and lakes need to be removed to simplify the resulting skeleton
* Create the skeleton of land and ocean (seaprately) with the scipy function
* Prune and shorten the skeleton using the fil_finder library
* Apply the skeleton over coarsened grids (included here too). If the land skeleton is occupied over an ocean cell, it will become land, and vice versa