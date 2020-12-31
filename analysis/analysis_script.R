# This is the R script used to conduct the MixSIAR analysis within the
# manuscript, McCluskey et. al. 2020. Foraging preferences of an apex marine predator
# revealed through stable isotope and stomach content analyses.

library(tidyverse)
library(here)
library(MixSIAR)
library(patchwork)

# read in the raw data stored in the data-raw folder

# read in the consumer data
consumer_data <- readr::read_csv(here('data-raw','MeanConsumerSI.csv'))

# read in the source data and rename the meand/sd columns to match the
# expectations of MixSIAR. Also, we put source (species name) in the first
# column.
source_data <- readr::read_csv(here('data-raw','sourceIsotopeNLatin.csv')) %>%
  dplyr::select(source, everything())

# since the trophic enrichment values are all the same, we just create this
# table from the source data species names.
trophic_enrich <- tibble(
  source = source_data %>% distinct(source) %>% pull(source),
  Meand13C = 1.01,
  SDd13C = 0.37,
  Meand15N = 1.57,
  SDd15N = 0.52
)


# we only want to include species for which there are samples in both the
# estuary and the coastal habitats. So, the easiest thing to do is pull out
# a list of all the species to NOT include (those with only 1 entry).
unique_species <- source_data %>%
  dplyr::group_by(source) %>%
  summarise(count = n()) %>%
  dplyr::filter(count ==1) %>%
  dplyr::pull(source)

# create a temporary file to store our csv
filename <- tempfile()

# write the consumer_data to the temp file and set the other function call
# arguments. Note we are not including dolphin name as an individual random
# effect. Habitat is a fixed effect so that is what we are interested in
# comparing
consumer_data %>%
  readr::write_csv(filename)
iso_names <- c("d13C","d15N")
factors <- c("Habitat")
fac_random <- c(FALSE)
fac_nested <- NULL

# read mix data in for MixSIAR analysis
my_mix <- load_mix_data(filename, iso_names, factors, fac_random, fac_nested,
                        cont_effects = NULL)

# create a temporary file to store our csv
filename <- tempfile()

# before writing out to the csv, remove the species that only have one entry
source_data %>%
  dplyr::filter(!source %in% unique_species) %>%
  readr::write_csv(filename)

# read source data in for MixSIAR analysis
my_source <- load_source_data(filename, source_factors="Habitat", conc_dep = FALSE,
                              data_type = "means", my_mix)

# create a temporary file to store our csv
filename <- tempfile()

# before writing out to the csv, remove the species that only have one entry
trophic_enrich %>%
  dplyr::filter(!source %in% unique_species) %>%
  readr::write_csv(filename)
my_discr <- load_discr_data(filename, my_mix)

# this will create an iso plot of the mix and source data
plot_data(filename = "iso_plot", plot_save_pdf=TRUE, plot_save_png=FALSE,
          mix=my_mix, source = my_source, discr = my_discr)

# create our model text file for JAGS
model_filename <- "MixSIAR_model.txt"   # Name of the JAGS model file

# we want to use multiplicative error Process*Residual which corresponds to
# the preferred default Model 4 in Stock and Semmens paper
resid_err <- TRUE
process_err <- TRUE
write_JAGS_model(model_filename, resid_err, process_err, my_mix, my_source)

# this runs the model and the results are stored in jags.1.
jags.1 <- run_model(run="test", my_mix, my_source, my_discr, model_filename,
                    alpha.prior = 1, resid_err, process_err)

# this creates all of the output plots, summary text, an diagnostic test.
# pdf versions of the plots are stored in the root directory. We source a
# custom modification of `output_JAGS.R` that includes customizations to the
# ggplot code to match our desired aesthetic for the manuscript figure
source(here::here('R','output_JAGS.R'))
output_JAGS(jags.1,my_mix,my_source)

