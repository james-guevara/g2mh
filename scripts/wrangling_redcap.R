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

###################################################################################################
# Load data
###################################################################################################
RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv = fread("data/RareCNVPhenotypicCol_DATA_2024-12-18_2059.csv") %>%
  filter(redcap_data_access_group != "test")
genomatrix_freeze1_g2mh_04222024.csv = fread("data/genomatrix_freeze1_g2mh_04222024.csv") %>%
  mutate(across(everything(), ~ ifelse(trimws(.) == "", NA, .)))

###################################################################################################
# Inspecting sex discrepancies
###################################################################################################
RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv.sex = RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv %>%
  select(rarecnv_id, dem_sex, dem_sex_fam) %>%
  group_by(rarecnv_id) %>%
  summarize(
    across(where(is.numeric), ~first(na.omit(.))),
    .groups = "drop"
  ) %>%
  mutate(redcapSex = coalesce(dem_sex, dem_sex_fam))

## Testing...
#x = RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv.sex %>% 
#  filter(!is.na(dem_sex) & !is.na(dem_sex_fam)) %>% # 12 obs. (so 12 samples with both columns filled)
#  filter(dem_sex != dem_sex_fam) # 0 obs. (so the values match for samples with both values filled)

genomatrix_freeze1_g2mh_04222024.csv.sex = genomatrix_freeze1_g2mh_04222024.csv %>%
  select(rarecnv_id, GUID, sex_redcap, sex_genetic, sex_aneuploidy_genetic) %>%
  mutate(in_genomatrix = 1)

redcap_genomatrix.sex.concordance = full_join(RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv.sex, genomatrix_freeze1_g2mh_04222024.csv.sex, by = "rarecnv_id") %>%
  mutate(
    oldredcap_newredcap_match = case_when(
      (sex_redcap == "Male" & redcapSex == 1) | (sex_redcap == "Female" & redcapSex == 2) ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(
    geneticSex_newredcap_match = case_when(
      (sex_genetic == "Male" & redcapSex == 1) | (sex_genetic == "Female" & redcapSex == 2) ~ 1,
      TRUE ~ 0
    )
  ) %>%
  select(in_genomatrix, rarecnv_id, GUID, redcapSex, sex_redcap, sex_genetic, oldredcap_newredcap_match, geneticSex_newredcap_match)

# Looking at discordant samples...
x = redcap_genomatrix.sex.concordance %>%
  filter(in_genomatrix == 1) %>%
  filter(oldredcap_newredcap_match == 0)
# 34 observations, all of which have sex_redcap as "Unknown". But the new RedCap sex matches the genetic sex for all of them now. 
x = redcap_genomatrix.sex.concordance %>%
  filter(in_genomatrix == 1) %>%
  filter(geneticSex_newredcap_match == 0)
# 0 observations, so they all match...

# Conclusion:
# The genetic sex matches the sex on RedCap now.


###################################################################################################
# Inspecting affected status discrepancies
###################################################################################################
RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv.affected = RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv %>%
  select(rarecnv_id, affected_unknown, gen_affected, gen_group_class, affected, group_class, group_fish, group_16p_type, group_22q_type, group_22q_type_other, group_22q_type_tbx1, gen_status_gsa_qc, gen_status_wgs_qc) %>%
  group_by(rarecnv_id) %>%
  summarize(
    across(where(is.numeric), ~first(na.omit(.))),
    .groups = "drop"
  )
# RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv.affected is 2347 obs. of 9 variables

## Testing...
## Just checking if any rarecnv_id is duplicated
#RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv.affected %>% filter(duplicated(rarecnv_id)) # No duplicate rarecnv_ids

genomatrix_freeze1_g2mh_04222024.csv.affected = genomatrix_freeze1_g2mh_04222024.csv %>%
  select(rarecnv_id, GUID, group_class_redcap, group_class_genetic) %>%
  mutate(in_genomatrix = 1)

redcap_genomatrix.affected.concordance = full_join(RareCNVPhenotypicCol_DATA_2024_12_18_2059.csv.affected, genomatrix_freeze1_g2mh_04222024.csv.affected, by = "rarecnv_id") %>%
  mutate(
    oldGroupClassRedcap_newGroupClassRedcap_matches = case_when(
      0 |
        (group_class_redcap == "22q11.2_Deletion" & group_class == 1) |  
        (group_class_redcap == "22q11.2_Duplication" & group_class == 2) | 
        (group_class_redcap == "16p11.2_Deletion" & group_class == 3) | 
        (group_class_redcap == "16p11.2_Duplication" & group_class == 4) | 
        (group_class_redcap == "22q11.2 Duplication & 16p11.2 Deletion" & group_class == 5) | 
        (group_class_redcap == "22q11.2 Duplication & 16p11.2 Duplication" & group_class == 6) | 
        (group_class_redcap == "22q11.2 Triplication" & group_class == 11) |
        (is.na(group_class_redcap) & is.na(group_class)) ~ 1,
      TRUE ~ 0
    ),
    geneticGroupClass_newGroupClassRedcap_matches = case_when(
      0 |
        (group_class_genetic == "22q11.2_Deletion" & group_class == 1) |  
        (group_class_genetic == "22q11.2_Duplication" & group_class == 2) | 
        (group_class_genetic == "16p11.2_Deletion" & group_class == 3) | 
        (group_class_genetic == "16p11.2_Duplication" & group_class == 4) | 
        (group_class_genetic == "22q11.2 Duplication & 16p11.2 Deletion" & group_class == 5) | 
        (group_class_genetic == "22q11.2_Duplication_&_16p11.2_Duplication" & group_class == 6) | 
        (group_class_genetic == "22q11.2 Triplication" & group_class == 11) |
        (is.na(group_class_genetic) & is.na(group_class)) ~ 1,
      TRUE ~ 0
    ),
    geneticGroupClass_genGroupClass_matches = case_when(
      0 |
        (group_class_genetic == "22q11.2_Deletion" & gen_group_class == 1) |  
        (group_class_genetic == "22q11.2_Duplication" & gen_group_class == 2) | 
        (group_class_genetic == "16p11.2_Deletion" & gen_group_class == 3) | 
        (group_class_genetic == "16p11.2_Duplication" & gen_group_class == 4) | 
        (group_class_genetic == "22q11.2 Duplication & 16p11.2 Deletion" & gen_group_class == 5) | 
        (group_class_genetic == "22q11.2_Duplication_&_16p11.2_Duplication" & gen_group_class == 6) | 
        (group_class_genetic == "22q11.2 Triplication" & gen_group_class == 11) |
        (is.na(group_class_genetic) & is.na(gen_group_class)) ~ 1,
      TRUE ~ 0
    )
  ) %>%
  mutate(genetic_data_type = case_when
         (
           is.na(gen_status_gsa_qc) & (gen_status_wgs_qc == 1) ~ "WGS",
           is.na(gen_status_wgs_qc) & (gen_status_gsa_qc == 1) ~ "GSA",
           TRUE ~ NA
         )
  ) %>%
  select(in_genomatrix, rarecnv_id, GUID, gen_group_class, group_class, group_class_redcap, group_class_genetic,
         oldGroupClassRedcap_newGroupClassRedcap_matches,
         geneticGroupClass_newGroupClassRedcap_matches,
         geneticGroupClass_genGroupClass_matches,
         genetic_data_type)

    
# Looking at discordant samples...
x = redcap_genomatrix.affected.concordance %>%
  filter(in_genomatrix == 1) %>%
  filter(oldGroupClassRedcap_newGroupClassRedcap_matches == 0)
# 15 observations, though 6 of new RedCap assignments match the genetic group class.
x = redcap_genomatrix.affected.concordance %>%
  filter(in_genomatrix == 1) %>%
  filter(geneticGroupClass_newGroupClassRedcap_matches == 0)
# Looking at just WGS samples (16 of them)
x = redcap_genomatrix.affected.concordance %>%
  filter(in_genomatrix == 1) %>%
  filter(geneticGroupClass_newGroupClassRedcap_matches == 0) %>%
  filter(genetic_data_type == "WGS")
# 20 observations
x = redcap_genomatrix.affected.concordance %>%
  filter(in_genomatrix == 1) %>%
  filter(geneticGroupClass_genGroupClass_matches == 0)
# 2 observations (but group_class_genetic matches the new RedCap group class for both of them)

# Summary:
# I'll want to determine why there are 2 discrepancies between the genomatrix genetic group class and the gen_group_class in RedCap... but I suspect they are simply clerical errors since the genomatrix genetic group class matches the RedCap group class (presumably detected from FISH data).
# The main problem are the 20 discrepancies between RedCap's new group class (from FISH) and the genomatrix genetic group class.
# I've looked at the 4 GSA samples already, but need to write out some information about them and see if there was a sample swap or something like that.
# The rest are WGS samples... so I'll look at the VCFs first, then the IGVs in the region, and then mosdepth coverage output.
# Make some summary slides for these samples with findings. (Both GSA and WGS)

