###Summary###

The function of this code is to explore co-identification of ID in different organs, and classify according to the number. At the same time, the classification of different IDs can be displayed, such as the glycan type, core type, whether it contains sialic acid or fucose to determine the characteristics of co-identified IDs.


###Depedencies###
R version 4.4.1 is required.
Before using it, packages "ComplexHeatmap", "ggplot2", "tidyverse", "circlize" need to be installed.


###Protocol Steps###
# It does not need to be installed on the desktop, and can be used once R is configured.
# Set the working directory to the folder named "Heatmap" in R 
environment (or RStudio).
# Select all codes and click "run", these results are generated automatically.

# After running the program, you will get: 
1) a complex heatmap graph.svg (IDs are segmented according to the number of identifications in each organ),
2) a .csv of IDs being identified by several organs.
