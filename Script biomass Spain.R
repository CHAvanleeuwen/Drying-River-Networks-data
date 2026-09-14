library(ggplot2)
library(dplyr)
library(reshape2)
library(readxl)
library(vegan)
library(GGally)
library(leaflet)

#### PRODUCES csv file with biomass per each individual measured for its length ####

#### 1. load lengths per taxonomic group per sampling location+time ####

cam <- read_xlsx("Campaigns overview file.xlsx",sheet = "campaigns")

#add biomass per campaign per taxonomic group


#load the xlsx datasheets into lists
#change the lists into dataframes
#concatenate the dataframes into one dataframe called "data"

#Campaign1
cam1 <- lapply(excel_sheets("/length-biomass datafiles Spain/SPAIN_1st campaign_FINAL.xlsx"), function(x)
  read_excel("/length-biomass datafiles Spain/SPAIN_1st campaign_FINAL.xlsx", sheet = x, range = cell_cols("A:Y")))

data <- as.data.frame(cam1[1]) # start dataframe to which rest can be concatenated

for ( i in 2:length(cam1) ) {
  ff <- as.data.frame(cam1[i])
  data <- rbind(data, ff) 
}

cam1 <- data
str(cam1)

#Campaign2
cam2 <- lapply(excel_sheets("/length-biomass datafiles Spain/SPAIN_2nd campaign_FINAL.xlsx"), function(x)
  read_excel("/length-biomass datafiles Spain/SPAIN_2nd campaign_FINAL.xlsx", sheet = x,range = cell_cols("A:Y")))

data <- as.data.frame(cam2[1]) # start dataframe to which rest can be concatenated
str(data)

for ( i in 2:length(cam2) ) {
  ff <- as.data.frame(cam2[i])
  data <- rbind(data, ff)
}

cam2 <- data
str(cam2)

#Campaign3
cam3 <- lapply(excel_sheets("/length-biomass datafiles Spain/Spain_3rd campaign_Final.xlsx"), function(x)
  read_excel("/length-biomass datafiles Spain/Spain_3rd campaign_Final.xlsx", sheet = x,range = cell_cols("A:Y")))

data <- as.data.frame(cam3[1]) # start dataframe to which rest can be concatenated

ff2 <- as.data.frame(cam3[2])

for ( i in 2:length(cam3)) {
  ff <- as.data.frame(cam3[i])
  data <- rbind(data, ff) 
}

cam3 <- data
str(cam3)

#Campaign4
cam4 <- lapply(excel_sheets("/length-biomass datafiles Spain/SPAIN_4th campaign_FINAL.xlsx"), function(x)
  read_excel("/length-biomass datafiles Spain/SPAIN_4th campaign_FINAL.xlsx", sheet = x,range = cell_cols("A:Y")))

data <- as.data.frame(cam4[1]) # start dataframe to which rest can be concatenated

for ( i in 2:length(cam4)) {
  ff <- as.data.frame(cam4[i])
  data <- rbind(data, ff) 
}

cam4 <- data
str(cam4)

#Campaign5
cam5 <- lapply(excel_sheets("/length-biomass datafiles Spain/SPAIN_campaign5_FINAL.xlsx"), function(x)
  read_excel("/length-biomass datafiles Spain/SPAIN_campaign5_FINAL.xlsx", sheet = x,range = cell_cols("A:Y")))

data <- as.data.frame(cam5[1]) # start dataframe to which rest can be concatenated

for ( i in 2:length(cam5)) {
  ff <- as.data.frame(cam5[i])
  data <- rbind(data, ff) 
}

cam5 <- data
str(cam5)

#Campaign6
cam6 <- lapply(excel_sheets("/length-biomass datafiles Spain/Spain_campaign6_FINAL.xlsx"), function(x)
  read_excel("/length-biomass datafiles Spain/Spain_campaign6_FINAL.xlsx", sheet = x,range = cell_cols("A:Y")))

data <- as.data.frame(cam6[1]) # start dataframe to which rest can be concatenated

for ( i in 2:length(cam6)) {
  ff <- as.data.frame(cam6[i])
  data <- rbind(data, ff) 
}

cam6 <- data
str(cam6)

# add campaign numbers to the datafiles before merging
cam1$campaign <- rep(1, length(cam1$Country))
cam2$campaign <- rep(2, length(cam2$Country))
cam3$campaign <- rep(3, length(cam3$Country))
cam4$campaign <- rep(4, length(cam4$Country))
cam5$campaign <- rep(5, length(cam5$Country))
cam6$campaign <- rep(6, length(cam6$Country))


cam <- rbind(cam1, cam2, cam3, cam4, cam5, cam6)
str(cam)
cam <- cam %>% rename(	country = Country, 
                       genus_coded = Genus, 
                       campaign_number = campaign,
                       image_ID = Sample.ID,
                       length_mm = Length..mm.)
cam <- cam %>% select(campaign_number, country, genus_coded, length_mm, image_ID)
str(cam)

# add phylogenetic information to the cam datafile

phylo <- read_excel("Phylo_data_all_taxa.xlsx")

phylo <- unique(phylo)

#create a new variable with the genus name for adding the parameter values later
# keep genus_coded for uploading the data later in IRBAS
phylo$genus_coded <- phylo$genus

phylo$genus<-gsub(" sp.","",as.character(phylo$genus))
phylo$genus<-gsub(" Lv.","",as.character(phylo$genus))
phylo$genus<-gsub(" Ad.","",as.character(phylo$genus))
phylo$genus<-gsub(" Gen.","",as.character(phylo$genus))


data <- merge(cam, phylo, by.x = "genus_coded", all.x = TRUE)


write.csv(data, "Campaign data incl taxonomic information Spain.csv")

#### 2: add taxonomic information for all based on genus ####

data <- read.csv("Campaign data incl taxonomic information Spain.csv")

str(data)


# 2: calculate biomass based on lengths and regression parameters for the taxa from regression file ####
# first by genus level, if not available by family level, if not available by order

para <- read_excel("/Regressions database.xlsx", sheet = "regressions")

para <- para %>% select(order, family, genus, a, b) %>% 
  mutate(genus =genus) %>%
  as.data.frame()
para$order <- as.factor(para$order)
para$family <- as.factor(para$family)
para$genus <- as.factor(para$genus)
levels(para$order)
levels(para$family)
levels(para$genus)

levels(factor(para$order))
# calculate the mean values of a and b per genus, family or order to improve para-dataframe

str(para)
po <- para %>% 	group_by(order) %>%
  summarise(a_ord = mean(a, na.rm = TRUE), b_ord = mean(b,na.rm = TRUE)) %>%
  as.data.frame()
pf <- para %>% 	group_by(family) %>%
  summarise(a_fam = mean(a, na.rm = TRUE), b_fam = mean(b, na.rm = TRUE)) %>%
  as.data.frame() %>%
  na.omit()
pg <- para %>% 	group_by(genus) %>%
  summarise(a_gen = mean(a, na.rm = TRUE), b_gen = mean(b, na.rm = TRUE)) %>%
  as.data.frame() %>%
  na.omit()

# add a and b values at as high possible taxonomic resolution to d1 dataframe

d1 <- merge(data, po, all.x = TRUE, by = "order")
d1 <- left_join(d1, pf, by = "family", all.x = TRUE)
d1 <- left_join(d1, pg, by = "genus")

# for loop that adds a and b values to the dataframe based on availability by gen, fam, order, or "class" level, indicating the level of accuracy in the level_regression column
for (i in 1:length(d1$country) ) {
  if ( !is.na(d1$a_gen[i]) ) {
    d1$a[i] <- d1$a_gen[i]
    d1$b[i] <- d1$b_gen[i]
    d1$level_regression[i] <- "genus"
    
  } else if ( !is.na(d1$a_fam[i])) {
    d1$a[i] <- d1$a_fam[i]
    d1$b[i] <- d1$b_fam[i]
    d1$level_regression[i] <- "family"
    
  } else if ( !is.na(d1$a_ord[i])) {
    d1$a[i] <- d1$a_ord[i]
    d1$b[i] <- d1$b_ord[i]
    d1$level_regression[i] <- "order"
  }
}

str(d1)
d <- d1 %>% select(country, campaign_number, length_mm, image_ID, order, family, genus, genus_coded, a, b, level_regression)
str(d)


levels(factor(d$order)) # all orders represented?

#### 4. calculate invert biomass data ####

#formula 

d$order <- factor(d$order)
d$mass_per_ind_in_mg <- d$a * d$length_mm^d$b # mass in mg dry mass per individual
plot(mass_per_ind_in_mg ~ length_mm, data = d, col = as.numeric(factor(d$order)))

d$site <- sub(".*GEN", "", d$image_ID)
d$site <- substr(d$site, start = 1, stop = 2)
d$site <- gsub(" ", "", paste("GEN",d$site))
d$site

d <- d[d$country != "NA",] # remove missing data

write.csv(d, "Dryver Spain length and biomass data all campaigns.csv")

