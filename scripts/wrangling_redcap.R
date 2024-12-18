###################################################################################################
# Load necessary libraries
###################################################################################################
library(conflicted)
library(tidyverse)  # This loads dplyr, ggplot2, and other core tidyverse packages
library(data.table) # For high-performance data manipulation

###################################################################################################
# Set up
###################################################################################################
# Specify preferences for conflicting functions
# I use dplyr functions like these ones. When you load other libraries, sometimes these function names are overridden by some other library's function with the same names.
# So I use this to ensure that the dplyr functions are the ones that are called (otherwise you can just write e.g. dplyr::select).
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("first", "dplyr")
# Remove all variables in current environment
rm(list = ls())
# Set the working directory
setwd("/Users/jamesguevara/sebatlab Dropbox/James Guevara/g2mh")



