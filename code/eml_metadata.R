## Goal: create EML metadata for project using the Environmental Data Initiative's EML Assembly Line (https://github.com/EDIorg/EMLassemblyline)

# Put data and code in a folder together to be grabbed by make_eml
# Generate metadata files by editing current ones or generating them (see Github page for tutorial)
# Edit and run this script
# Keywords: https://vocab.lternet.edu/vocab/vocab/index.php


#### set up ####

# clear environment
rm(list=ls())

# load libraries
library(EMLassemblyline)
library(tidyverse)
library(knitr)


#### import templates ####

# list of data files
dlist <- list.files(path = "edi/data",
                    pattern = ".csv")

# create high level text files
template_core_metadata(path = "edi/metadata",
                       license = "CCBY")

# create an attribute file for each data table
template_table_attributes(path = "edi/metadata",
                          data.path = "edi/data",
                          data.table = dlist)

# create a categorical value file for each data table
template_categorical_variables(path = "edi/metadata",
                               data.path = "edi/data")

# look at units
view_unit_dictionary()


#### data file descriptors ####

dlist

# description list
ddlist <- "infection, biomass, and chlorophyll data"

# name list
dnlist <- "experiment data"

# print table
# dtable <- data.frame(data = dlist, description = ddlist)
# kable(dtable)


#### code descriptors ####

# list of code files
clist <- list.files(path = "edi/code",
                    pattern = ".R")

# code descripions
cdlist <- c(
  "code to analyze biomass data and create figure",
  "code to analyze chlorophyll data and create figure",
  "initial data processing",
  "code to analyze infection data and create figures"
)

# name list
cnlist <- c("Biomass analysis",
            "Chlorophyll analysis",
            "Data processing",
            "Infection analysis"
            )

# print table
# ctable <- data.frame(code = clist, desription = cdlist)
# kable(ctable)


#### make eml ####

# copied data and code from the respective folders and put into edi folder

make_eml(path = "metadata",
         data.path = "edi",
         dataset.title = "Long-term nitrogen enrichment mediates the effects of nitrogen supply and co-inoculation on a viral pathogen",
         data.table = dlist,
         data.table.name = dnlist,
         data.table.description = ddlist,
         data.table.quote.character = rep("\"", length(dlist)),
         other.entity = clist,
         other.entity.name = cnlist,
         other.entity.description = cdlist,
         temporal.coverage = c("2014-06-01", "2016-06-26"),
         geographic.description = "St. Paul, MN, USA",
         geographic.coordinates = c(44.98, -93.18, 44.98, -93.18),
         maintenance.description = "completed", 
         user.id = "aekendig",
         user.domain = "EDI",
         package.id = "edi.209.2")


#### check warnings ####

eml <- EML::read_eml("metadata/edi.209.2.xml")
EML::eml_validate(eml)

# https://portal-s.edirepository.org/nis/mapbrowse?scope=edi&identifier=209&revision=1