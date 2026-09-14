library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
library(vegan)
library(GGally)
library(leaflet)
library(pals)
library(stringr)
library(lubridate)

#### 1. load csv files with length and biomasses per DRN and per individual animal ####

cr <- read.csv("Dryver Croatia length and biomass data all campaigns.csv")
fr <- read.csv("Dryver France length and biomass data all campaigns.csv")
hu <- read.csv("Dryver Hungary length and biomass data all campaigns.csv")
fi <- read.csv("Dryver Finland length and biomass data all campaigns.csv")
sp <- read.csv("Dryver Spain length and biomass data all campaigns.csv")
cz <- read.csv("Dryver Czech length and biomass data all campaigns.csv")

data <- rbind(cr,  hu, fi, sp, cz, fr)

data <- data[!is.na(data$country),] # remove empty rows

str(data)

data[data$country == "NA",] 
data[data$stream_type == "NA",] 

data <- data[data$site != "BUK36",] #deleted only done in 1st campaign

write.csv(data, file = "/Dryver length and biomass data per individuals all campaigns 6 DRN's.csv")

d <- read.csv(data, file = "Dryver length and biomass data per individuals all campaigns 6 DRN's.csv")

## Calculates secondary productivity based on Morin & Dumont 1994

#### 2. combine with the abundance datafile ####

d$order
d[d$stream_type == "NA",] # check that all have a stream type
d$length_mm <- as.numeric(d$length_mm, na.rm = TRUE)
str(d)

# calculate means per taxonomic group, campaign and country
# don't include site because we did not measure lengths for every specific site
# as lengths differ between perennial and intermittent streams, separately compute mean lengths for these stream types
# sites receive mean values from means across the DRN for that genus and campaign at that mmoment

env <- read.csv("env_tot.csv") # load environmental data

env <- env %>% rename(temp_c = temperature_C,
                      country = DRN, campaign_number = campaign)
env$campaign_number <- as.numeric(env$campaign_number)
str(env)
str(data)
levels(factor(env$site))
levels(factor(data$site))
levels(factor(d$site))

d <- merge(d, env, by = c("site", "campaign_number", "country"), all.x = TRUE)
str(d)


#### here interpolation of the lengths starts
### not for all campaigns, sites, genera, lengths are known
### this part interpolates lengths based on sites TR or PR
#### assumption: lengths of a specific genus with unknown lengths are the same as 
#### the mean lengths of that genus in other sites of the DRN during the same campaign
### but ONLY lengths of sites that are ALSO TR or also PR, so this is taken separately

d2 <- d %>% group_by(genus, campaign_number, country, stream_type) %>%
  summarize(mean_length_mm = mean(length_mm, na.rm = TRUE), 
            mean_mass_per_ind_in_mg = mean(mass_per_ind_in_mg, na.rm = TRUE),
            temp_c = mean(temp_c) )  %>%
  as.data.frame()

str(d2)

# d2 can be merge into the lengths file, but first load from the abundance file all the
# datapoint that are there. the "d" dataframe contains al actual measurements; the abundance dataframe 
# that will now be loaded contains all the samples. Then next step is to add all the measurements to the sample overview,
# and then interpolate the lengths to the samples not really measured

######################################################

#### 3. load abundance data: number of individuals from the abundance Dryver file

abun <- read_excel("DRYvER minv FULL genus data ver10-3.xlsx", sheet = "MINV genus data")
abun <- abun %>% rename(sample_ID = "Sample ID", 
                        country = "Country", 
                        site = "Site", 
                        order = "Taxagroup",
                        family = "Family",
                        no_ind = "No. of ind.")

str(abun)

# add campaign numbers to the abundance data: load campaign overview file with all campaigns

campaigns <- read_excel("Campaigns overview file.xlsx", sheet = "campaigns")

# check for country information
levels(factor(abun$country))
levels(factor(campaigns$country))
levels(factor(d$country))

# merge abundance data with campaign data
abun <- merge(abun, campaigns, all.x = TRUE)
str(abun)
levels(factor(abun$campaign_number))

abun <- abun %>% rename(genus = "Genus / Higher taxonomic group") %>%
  as.data.frame()
abun <- abun %>% select(no_ind, sample_ID, genus, family, country, site, order, sample_ID, campaign_number)
str(abun)
levels(factor(abun$campaign_number))

abun <- na.omit(abun) #remove missing values
#abun <- abun[abun$campaign_number < 7,] # remove extra campaigns only done in France
levels(factor(abun$campaign_number))

abun$genus_coded <- abun$genus # create the coded version for the data upload

# clean up genus names 
abun$genus<-gsub(" sp.","",as.character(abun$genus))
abun$genus<-gsub(" Lv.","",as.character(abun$genus))
abun$genus<-gsub(" Ad.","",as.character(abun$genus))
abun$genus<-gsub(" Gen.","",as.character(abun$genus))

# check if names of the genera are the same for the abundance data (abun) and the length dataframe (d)
levels(factor(abun$genus)) # genera in the abundance file 
levels(factor(d$genus)) # genera in the lengths database

# load multiplier for the abundance, to convert from "per sample" to "per m2"

mult <- read.csv("Multipliers for 6 DRNs.csv")
str(mult)
mult$country <- factor(mult$country)
plot(multiplier ~ as.factor(country), data = mult)
str(abun)

abun2 <- merge(abun, mult, by = c("site", "campaign_number", "country"), all.x = TRUE)
str(abun2)
# number of individuals per m2
abun2[is.na(abun2$multiplier),]$multiplier <- 1 # replace unknown sampling surfaces with 1 m2

# apply the multiplier to the abundance data
abun2$no_ind_m2 <- abun2$multiplier * abun2$no_ind
plot(no_ind_m2 ~ as.factor(country), data = abun2)

abun <- abun2 # back to the original dataframe
abun[abun$stream_type == "NA",] 

# d is dataset with per campaign, per genus the mean length of the specimens
# creating d ensures that also samples that have not been processed receive a mean length
# we did not measure lengths for all individuals but do have distributions for all samples 
# here we assign a mean to every site / campaign / DRN combination for further use.

# add mean lengths per campaign per genus (d2) to the abundance datafile (abun)

abun$campaign_number <- as.numeric(abun$campaign_number)

# merge abundance and d2 files by campaign number, genus and country
# don't merge per site, because we work with mean lengths per DRN per campaign
# we don't have mean lengths for every site, but measured sites are representatively chosen
# based on spatial distribution over the DRN and variation in densities of genera 

### first link the stream type to abundance dataframe

stream_type <- read.csv("stream_types_final.csv")
stream_type <- stream_type %>% select(site, stream_type)
abun <- merge(abun, stream_type, by = "site", all.x = TRUE)

# first load the actually measured values into the abundance dataframe

d2_per_site <- d %>% group_by(genus, campaign_number, country, site) %>%
  summarize(mean_length_mm = mean(length_mm, na.rm = TRUE), 
            mean_mass_per_ind_in_mg = mean(mass_per_ind_in_mg, na.rm = TRUE),
            temp_c = mean(temp_c) )  %>%
  as.data.frame()


d5 <- merge(abun, d2_per_site, by = c("campaign_number", "genus", "country", "site"), all.x = TRUE)
str(d5)

# add the means originally measured per site

d5$data_type <- ifelse(is.na(d5$mean_length_mm), "Interpolation", "Measurement")

# opknippen datafile in met en zonder lengtes

d6 <- d5[is.na(d5$mean_length_mm),] # d6 is with unknown lengths

#remove the columns that need to be added from d2
d6 <- d6 %>% select(-mean_length_mm, -mean_mass_per_ind_in_mg, -temp_c)

# merge d6 with d2 (which was the dataframe from above in which mean lengths are available based on stream types)
# adding these columns with the interpolated data
d3 <- merge(d6, d2, by = c("campaign_number", "genus", "country", "stream_type"), all.x = TRUE)
str(d3) # including campaign-specific lengths

# combine the d5 with the original data and the d6 with the interpolated data
d7 <- rbind(d5[d5$data_type == "Measurement",], d3)
str(d7)


######################################################################

abun1 <- d7[!is.na(d7$mean_length_mm),] # dataframe without missing lengths
abun2 <- d7[is.na(d7$mean_length_mm),] # dataframe in which lengths are missing
# these lengths are missing because the genera were only measured in other campaigns, not exactly in that camppaign
# take the mean lengths of the genera across all campaigns for each DRN, and add that length to the datafile
# separate between perennial and intermittent streams

d2_means <- d2 %>% group_by(genus, stream_type, campaign_number, country) %>%
  summarize(mean_length_per_genus = mean(mean_length_mm, na.rm=TRUE),
            mean_mass_per_genus = mean(mean_mass_per_ind_in_mg, na.rm = TRUE))

d4 <- merge(abun2, d2_means, by = c("genus", "campaign_number","country", "stream_type"), all.x = TRUE)
str(d4)

# remove the columns with NA's and replace these with estimates based on whole DRN's
abun2 <- d4 %>% select(-mean_length_mm, -mean_mass_per_ind_in_mg) %>% 
  rename(mean_length_mm = "mean_length_per_genus", 
         mean_mass_per_ind_in_mg = "mean_mass_per_genus")

# bind the rows with the original data and those with the filled in missing values
d5 <- rbind(abun1, abun2)

hist(d5$no_ind)
data <- d5[!is.na(d5$country),] # remove empty rows if any

head(data) # final dataset including lengths and biomasses 

missing_values <- data[is.na(data$mean_length_mm),]
str(missing_values) # 317 missing observations on 16000 is acceptable

#### 3. plot size differences among PR and TR streams per DRN ####

str(data) # environmental data already added earlier
# check whether any errors with data

plot(as.numeric(temp_c) ~ as.factor(country), data = d)

# check which countries are there
levels(factor(d$country))
ff <- data[is.na(data$stream_type),]
ff # no missing values

length(ff$campaign_number) # 2250 missing a stream type: not measured for lengths
length(data$campaign_number)

# empty rows were generated in the file, need to delete
data[data$stream_type == "NA",] 

data$stream_type <- as.factor(data$stream_type)

data <- data[data$stream_type != "NA",] # remove empty rows

write.csv(data, file = "Data incl lengths and biomass measured and interpolated.csv")

#### 4. calculate secondary productivity using Morin's formula ####

d <- data
#Morin & Dumont 1994:
"The growth models can be combined with information on 
the biomass size distribution of populations of single species 
or communities to estimate secondary production."

# formula: log10(growth) = a + b*log10(mass) + c*temperature
# load parameters for general model of growth from publication
# load parameters from table 2 in Morin & Dumont-paper per order
a_general <- -2.09
b_general <- -0.27
c_general <- 0.025

a_diptera <- -1.60
b_diptera <- -0.07
c_diptera <- 0.031

a_ephem <- -2.07
b_ephem <- -0.14
c_ephem <- 0.038

a_plecop <- -1.90
b_plecop <- -0.28
c_plecop <- 0

a_trichop <- -2.29
b_trichop <- -0.21
c_trichop <- 0.032

#create datasets for each of these orders for which parameters are available

d <- d[!is.na(d$country),]

d$order <- as.factor(d$order)
rest <- d %>% filter(!order %in% c("Diptera", "Ephemeroptera", "Plecoptera","Trichoptera")) 
dip <- d %>% filter(order %in% c("Diptera")) 
ephem <- d %>% filter(order %in% c("Ephemeroptera")) 
plecop <- d %>% filter(order %in% c("Plecoptera")) 
trichop <- d %>% filter(order %in% c("Trichoptera")) 

# check if no NA's introduced
ff <- d[is.na(d$country),]
ff <- trichop[is.na(trichop$country),]
ff <- rest[is.na(rest$country),]


# appy the Morin & Dumont formula to the data
rest$growth_per_ind_log_in_mg <- a_general + b_general * log10(rest$mean_mass_per_ind_in_mg) + c_general * rest$temp_c
dip$growth_per_ind_log_in_mg <- a_diptera + b_diptera * log10(dip$mean_mass_per_ind_in_mg) + c_diptera * dip$temp_c
ephem$growth_per_ind_log_in_mg <- a_ephem + b_ephem * log10(ephem$mean_mass_per_ind_in_mg) + c_ephem * ephem$temp_c
plecop$growth_per_ind_log_in_mg <- a_plecop + b_plecop * log10(plecop$mean_mass_per_ind_in_mg) + c_plecop * plecop$temp_c
trichop$growth_per_ind_log_in_mg <- a_trichop + b_trichop * log10(trichop$mean_mass_per_ind_in_mg) + c_trichop * trichop$temp_c
str(ephem)
str(rest)
d <- rbind(rest,dip, ephem, plecop, trichop)
str(d)

ff <- d[d$country == "NA",] # no missing values

d$mean_daily_growth_per_ind_in_mg <- 10^d$growth_per_ind_log_in_mg # growth per individual in grams
hist(d$mean_daily_growth_per_ind_in_mg) # in mg
max(d$mean_daily_growth_per_ind_in_mg, na.rm = TRUE)

#### 5. move from individual growth to secondary productivity per surface area ####
# growth in mg dry mass increase per day
# follow Morin formula:

# mass_per_ind # in mg 
# growth_per_ind # in mg per day

#Biomass is the product of mean individual body mass 
# and density (B = M  * N), and the product of biomass 
# and individual growth rate (g) is another simple way 
# to estimate production (P = mean(B) *  g). Benke 2010

# multiply with the abundance data 

# Production = mass*growth* number of individuals per m2
d$sec_prod_mg <- d$no_ind_m2 * d$mean_mass_per_ind_in_mg * d$mean_daily_growth_per_ind_in_mg 
hist(d$sec_prod_mg)
max(d$sec_prod_mg, na.rm = TRUE)
median(d$sec_prod_mg, na.rm = TRUE)
mean(d$sec_prod_mg, na.rm = TRUE)

levels(factor(d$country))

write.csv(d, "Final dataset.csv")

