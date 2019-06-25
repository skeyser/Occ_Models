#Mangrove 
setwd("C:/Users/Spencer/Box Sync/MSI_Research/Data/CBC Data/CBC_GIS/Extracted_LULC")
library(pacman)
pacman::p_load("reshape2", "lme4", "dplyr", "tidyverse", "cowplot", "readr", "data.table", "tfplot")

##Need to have the mangrove data transforme to sites as rows and veg types as cols
#Reshape should do the trick 


#Load in data
#files <- list.files(pattern = "*.csv")
#tbl <- lapply(files, read_csv) %>% bind_rows()
#tbl <- tbl[complete.cases(tbl), ]


fl.1980 <- read.csv(file = "mangrove_fl_1980.csv.csv", header = TRUE)
fl.1985 <- read.csv(file = "mangrove_fl_1985.csv.csv", header = TRUE)
fl.1990 <- read.csv(file = "mangrove_fl_1990.csv.csv", header = TRUE)
fl.1995 <- read.csv(file = "mangrove_fl_1995.csv.csv", header = TRUE)
fl.2000 <- read.csv(file = "mangrove_fl_2000.csv.csv", header = TRUE)
fl.2005 <- read.csv(file = "mangrove_fl_2005.csv.csv", header = TRUE)
fl.2010 <- read.csv(file = "mangrove_fl_2010.csv.csv", header = TRUE)
fl.2015 <- read.csv(file = "mangrove_fl_2015.csv.csv", header = TRUE)

tx.1980 <- read.csv(file = "mangrove_tx_1980.csv", header = TRUE)
tx.1985 <- read.csv(file = "mangrove_tx_1985.csv", header = TRUE)
tx.1990 <- read.csv(file = "mangrove_tx_1990.csv", header = TRUE)
tx.1995 <- read.csv(file = "mangrove_tx_1995.csv", header = TRUE)
tx.2000 <- read.csv(file = "mangrove_tx_2000.csv", header = TRUE)                    
tx.2005 <- read.csv(file = "mangrove_tx_2005.csv", header = TRUE)
tx.2010 <- read.csv(file = "mangrove_tx_2010.csv", header = TRUE)
tx.2015 <- read.csv(file = "mangrove_tx_2015.csv", header = TRUE)

la.1980 <- read.csv(file = "mangrove_la_1980.csv", header = TRUE)
la.1985 <- read.csv(file = "mangrove_la_1985.csv", header = TRUE)
la.1990 <- read.csv(file = "mangrove_la_1990.csv", header = TRUE)
la.1995 <- read.csv(file = "mangrove_la_1995.csv", header = TRUE)
la.2000 <- read.csv(file = "mangrove_la_2000.csv", header = TRUE)                    
la.2005 <- read.csv(file = "mangrove_la_2005.csv", header = TRUE)
la.2010 <- read.csv(file = "mangrove_la_2010.csv", header = TRUE)
la.2015 <- read.csv(file = "mangrove_la_2015.csv", header = TRUE)


##Re-lable and rearrange to sum up all mangrove classes 
fl.1980 <- fl.1980[, -1]
fl.1980$LABEL <- as.character(fl.1980$LABEL)
class(fl.1980$LABEL)
fl.1980[[2, 1]] <- "Mangrove"
fl.1980[[3, 1]] <- "Mangrove"
fl.1980[[8, 1]] <- "Mangrove"
fl.1980_sum <- aggregate(fl.1980[, 2:39], by = list(LULC = fl.1980$LABEL), FUN = sum)
fl.1980_sum <- as.data.frame(t(fl.1980_sum), stringsAsFactors = FALSE)
fl.1980_sum$year <- "1980"
rownames(fl.1980_sum) <- paste0(rownames(fl.1980_sum), "_", fl.1980_sum$year)
colnames(fl.1980_sum) <- c("UNK","grassland/wetland", "mangrove", "open water", "other veg", "urban", "year")
fl.1980_sum <- fl.1980_sum[-1, ]
fl.1980_sum <- fl.1980_sum[, -1]

fl.1985 <- fl.1985[, -1]
fl.1985$LABEL <- as.character(fl.1985$LABEL)
class(fl.1985$LABEL)
fl.1985[[2, 1]] <- "Mangrove"
fl.1985[[3, 1]] <- "Mangrove"
fl.1985[[8, 1]] <- "Mangrove"
fl.1985_sum <- aggregate(fl.1985[, 2:39], by = list(LULC = fl.1985$LABEL), FUN = sum)
fl.1985_sum <- as.data.frame(t(fl.1985_sum), stringsAsFactors = FALSE)
fl.1985_sum$year <- "1985"
rownames(fl.1985_sum) <- paste0(rownames(fl.1985_sum), "_", fl.1985_sum$year)
colnames(fl.1985_sum) <- c("UNK","grassland/wetland", "mangrove", "open water", "other veg", "urban", "year")
fl.1985_sum <- fl.1985_sum[-1, ]
fl.1985_sum <- fl.1985_sum[, -1]

fl.1990 <- fl.1990[, -1]
fl.1990$LABEL <- as.character(fl.1990$LABEL)
class(fl.1990$LABEL)
fl.1990[[2, 1]] <- "Mangrove"
fl.1990[[3, 1]] <- "Mangrove"
fl.1990[[8, 1]] <- "Mangrove"
fl.1990_sum <- aggregate(fl.1990[, 2:39], by = list(LULC = fl.1990$LABEL), FUN = sum)
fl.1990_sum <- as.data.frame(t(fl.1990_sum), stringsAsFactors = FALSE)
fl.1990_sum$year <- "1990"
rownames(fl.1990_sum) <- paste0(rownames(fl.1990_sum), "_", fl.1990_sum$year)
colnames(fl.1990_sum) <- c("UNK","grassland/wetland", "mangrove", "open water", "other veg", "urban", "year")
fl.1990_sum <- fl.1990_sum[-1, ]
fl.1990_sum <- fl.1990_sum[, -1]

fl.1995 <- fl.1995[, -1]
fl.1995$LABEL <- as.character(fl.1995$LABEL)
class(fl.1995$LABEL)
fl.1995[[2, 1]] <- "Mangrove"
fl.1995[[3, 1]] <- "Mangrove"
fl.1995[[8, 1]] <- "Mangrove"
fl.1995_sum <- aggregate(fl.1995[, 2:39], by = list(LULC = fl.1995$LABEL), FUN = sum)
fl.1995_sum <- as.data.frame(t(fl.1995_sum), stringsAsFactors = FALSE)
fl.1995_sum$year <- "1995"
rownames(fl.1995_sum) <- paste0(rownames(fl.1995_sum), "_", fl.1995_sum$year)
colnames(fl.1995_sum) <- c("UNK","grassland/wetland", "mangrove", "open water", "other veg", "urban", "year")
fl.1995_sum <- fl.1995_sum[-1, ]
fl.1995_sum <- fl.1995_sum[, -1]

fl.2000 <- fl.2000[, -1]
fl.2000$LABEL <- as.character(fl.2000$LABEL)
class(fl.2000$LABEL)
fl.2000[[2, 1]] <- "Mangrove"
fl.2000[[3, 1]] <- "Mangrove"
fl.2000[[8, 1]] <- "Mangrove"
fl.2000_sum <- aggregate(fl.2000[, 2:39], by = list(LULC = fl.2000$LABEL), FUN = sum)
fl.2000_sum <- as.data.frame(t(fl.2000_sum), stringsAsFactors = FALSE)
fl.2000_sum$year <- "2000"
rownames(fl.2000_sum) <- paste0(rownames(fl.2000_sum), "_", fl.2000_sum$year)
colnames(fl.2000_sum) <- c("UNK","grassland/wetland", "mangrove", "open water", "other veg", "urban", "year")
fl.2000_sum <- fl.2000_sum[-1, ]
fl.2000_sum <- fl.2000_sum[, -1]

fl.2005 <- fl.2005[, -1]
fl.2005$LABEL <- as.character(fl.2005$LABEL)
class(fl.2005$LABEL)
fl.2005[[2, 1]] <- "Mangrove"
fl.2005[[3, 1]] <- "Mangrove"
fl.2005[[8, 1]] <- "Mangrove"
fl.2005_sum <- aggregate(fl.2005[, 2:39], by = list(LULC = fl.2005$LABEL), FUN = sum)
fl.2005_sum <- as.data.frame(t(fl.2005_sum), stringsAsFactors = FALSE)
fl.2005_sum$year <- "2005"
rownames(fl.2005_sum) <- paste0(rownames(fl.2005_sum), "_", fl.2005_sum$year)
colnames(fl.2005_sum) <- c("UNK","grassland/wetland", "mangrove", "open water", "other veg", "urban", "year")
fl.2005_sum <- fl.2005_sum[-1, ]
fl.2005_sum <- fl.2005_sum[, -1]

fl.2010 <- fl.2010[, -1]
fl.2010$LABEL <- as.character(fl.2010$LABEL)
class(fl.2010$LABEL)
fl.2010[[2, 1]] <- "Mangrove"
fl.2010[[3, 1]] <- "Mangrove"
fl.2010[[8, 1]] <- "Mangrove"
fl.2010_sum <- aggregate(fl.2010[, 2:39], by = list(LULC = fl.2010$LABEL), FUN = sum)
fl.2010_sum <- as.data.frame(t(fl.2010_sum), stringsAsFactors = FALSE)
fl.2010_sum$year <- "2010"
rownames(fl.2010_sum) <- paste0(rownames(fl.2010_sum), "_", fl.2010_sum$year)
colnames(fl.2010_sum) <- c("UNK","grassland/wetland", "mangrove", "open water", "other veg", "urban", "year")
fl.2010_sum <- fl.2010_sum[-1, ]
fl.2010_sum <- fl.2010_sum[, -1]

fl.2015 <- fl.2015[, -13]
fl.2015[1, 3:40] <- 0
fl.2015 <- fl.2015[, -1]
fl.2015$LABEL <- as.character(fl.2015$LABEL)
class(fl.2015$LABEL)
fl.2015[[2, 1]] <- "Mangrove"
fl.2015[[3, 1]] <- "Mangrove"
fl.2015[[8, 1]] <- "Mangrove"
fl.2015_sum <- aggregate(fl.2015[, 2:39], by = list(LULC = fl.2015$LABEL), FUN = sum)
fl.2015_sum <- as.data.frame(t(fl.2015_sum), stringsAsFactors = FALSE)
fl.2015_sum$year <- "2015"
rownames(fl.2015_sum) <- paste0(rownames(fl.2015_sum), "_", fl.2015_sum$year)
colnames(fl.2015_sum) <- c("UNK","grassland/wetland", "mangrove", "open water", "other veg", "urban", "year")
fl.2015_sum <- fl.2015_sum[-1, ]
fl.2015_sum <- fl.2015_sum[, -1]

##Re-lable and rearrange to sum up all mangrove classes 
tx.1980 <- tx.1980[, -1]
tx.1980 <- tx.1980[1:7, ]
tx.1980$LABEL <- as.character(tx.1980$LABEL)
class(tx.1980$LABEL)
tx.1980[[2, 1]] <- "Mangrove"
tx.1980[[3, 1]] <- "Mangrove"
tx.1980$Year <- as.character(tx.1980$Year)
tx.1980_sum <- aggregate(tx.1980[, 2:12], by = list(LULC = tx.1980$LABEL), FUN = sum)
tx.1980_sum <- as.data.frame(t(tx.1980_sum), stringsAsFactors = FALSE)
tx.1980_sum$year <- "1980"
rownames(tx.1980_sum) <- paste0(rownames(tx.1980_sum), "_", tx.1980_sum$year)
colnames(tx.1980_sum) <- c("grassland/wetland", "mangrove", "open water", "other veg", "unclassified", "urban", "year")
tx.1980_sum <- tx.1980_sum[-1, ]
tx.1980_sum <- tx.1980_sum[, -5]

tx.1985 <- tx.1985[, -1]
tx.1985 <- tx.1985[1:7, ]
tx.1985$LABEL <- as.character(tx.1985$LABEL)
class(tx.1985$LABEL)
tx.1985[[2, 1]] <- "Mangrove"
tx.1985[[3, 1]] <- "Mangrove"
tx.1985$Year <- as.character(tx.1985$Year)
tx.1985_sum <- aggregate(tx.1985[, 2:12], by = list(LULC = tx.1985$LABEL), FUN = sum)
tx.1985_sum <- as.data.frame(t(tx.1985_sum), stringsAsFactors = FALSE)
tx.1985_sum$year <- "1985"
rownames(tx.1985_sum) <- paste0(rownames(tx.1985_sum), "_", tx.1985_sum$year)
colnames(tx.1985_sum) <-c("grassland/wetland", "mangrove", "open water", "other veg", "unclassified", "urban", "year")
tx.1985_sum <- tx.1985_sum[-1, ]
tx.1985_sum <- tx.1985_sum[, -5]


tx.1990 <- tx.1990[, -1]
tx.1990$LABEL <- as.character(tx.1990$LABEL)
class(tx.1990$LABEL)
tx.1990[[2, 1]] <- "Mangrove"
tx.1990[[3, 1]] <- "Mangrove"
tx.1990$Year <- as.character(tx.1990$Year)
tx.1990_sum <- aggregate(tx.1990[, 2:12], by = list(LULC = tx.1990$LABEL), FUN = sum)
tx.1990_sum <- as.data.frame(t(tx.1990_sum), stringsAsFactors = FALSE)
tx.1990_sum$year <- "1990"
rownames(tx.1990_sum) <- paste0(rownames(tx.1990_sum), "_", tx.1990_sum$year)
colnames(tx.1990_sum) <- c("grassland/wetland", "mangrove", "open water", "other veg", "unclassified", "urban", "year")
tx.1990_sum <- tx.1990_sum[-1, ]
tx.1990_sum <- tx.1990_sum[, -5]

tx.1995 <- tx.1995[, -1]
tx.1995$LABEL <- as.character(tx.1995$LABEL)
class(tx.1995$LABEL)
tx.1995[[2, 1]] <- "Mangrove"
tx.1995[[3, 1]] <- "Mangrove"
tx.1995$Year <- as.character(tx.1995$Year)
tx.1995_sum <- aggregate(tx.1995[, 2:12], by = list(LULC = tx.1995$LABEL), FUN = sum)
tx.1995_sum <- as.data.frame(t(tx.1995_sum), stringsAsFactors = FALSE)
tx.1995_sum$year <- "1995"
rownames(tx.1995_sum) <- paste0(rownames(tx.1995_sum), "_", tx.1995_sum$year)
colnames(tx.1995_sum) <- c("grassland/wetland", "mangrove", "open water", "other veg", "unclassified", "urban", "year")
tx.1995_sum <- tx.1995_sum[-1, ]
tx.1995_sum <- tx.1995_sum[, -5]

tx.2000 <- tx.2000[, -1]
tx.2000$LABEL <- as.character(tx.2000$LABEL)
class(tx.2000$LABEL)
tx.2000[[2, 1]] <- "Mangrove"
tx.2000[[3, 1]] <- "Mangrove"
tx.2000$Year <- as.character(tx.2000$Year)
tx.2000_sum <- aggregate(tx.2000[, 2:12], by = list(LULC = tx.2000$LABEL), FUN = sum)
tx.2000_sum <- as.data.frame(t(tx.2000_sum), stringsAsFactors = FALSE)
tx.2000_sum$year <- "2000"
rownames(tx.2000_sum) <- paste0(rownames(tx.2000_sum), "_", tx.2000_sum$year)
colnames(tx.2000_sum) <- c("grassland/wetland", "mangrove", "open water", "other veg", "unclassified", "urban", "year")
tx.2000_sum <- tx.2000_sum[-1, ]
tx.2000_sum <- tx.2000_sum[, -5]

tx.2005 <- tx.2005[, -1]
tx.2005$LABEL <- as.character(tx.2005$LABEL)
class(tx.2005$LABEL)
tx.2005[[2, 1]] <- "Mangrove"
tx.2005[[3, 1]] <- "Mangrove"
tx.2005$Year <- as.character(tx.2005$Year)
tx.2005_sum <- aggregate(tx.2005[, 2:12], by = list(LULC = tx.2005$LABEL), FUN = sum)
tx.2005_sum <- as.data.frame(t(tx.2005_sum), stringsAsFactors = FALSE)
tx.2005_sum$year <- "2005"
rownames(tx.2005_sum) <- paste0(rownames(tx.2005_sum), "_", tx.2005_sum$year)
colnames(tx.2005_sum) <- c("grassland/wetland", "mangrove", "open water", "other veg", "unclassified", "urban", "year")
tx.2005_sum <- tx.2005_sum[-1, ]
tx.2005_sum <- tx.2005_sum[, -5]

tx.2010 <- tx.2010[, -1]
tx.2010$LABEL <- as.character(tx.2010$LABEL)
class(tx.2010$LABEL)
tx.2010[[2, 1]] <- "Mangrove"
tx.2010[[3, 1]] <- "Mangrove"
tx.2010$Year <- as.character(tx.2010$Year)
tx.2010_sum <- aggregate(tx.2010[, 2:12], by = list(LULC = tx.2010$LABEL), FUN = sum)
tx.2010_sum <- as.data.frame(t(tx.2010_sum), stringsAsFactors = FALSE)
tx.2010_sum$year <- "2010"
rownames(tx.2010_sum) <- paste0(rownames(tx.2010_sum), "_", tx.2010_sum$year)
colnames(tx.2010_sum) <- c("grassland/wetland", "mangrove", "open water", "other veg", "unclassified", "urban", "year")
tx.2010_sum <- tx.2010_sum[-1, ]
tx.2010_sum <- tx.2010_sum[, -5]

tx.2015 <- tx.2015[, -1]
tx.2015$LABEL <- as.character(tx.2015$LABEL)
class(tx.2015$LABEL)
tx.2015[[2, 1]] <- "Mangrove"
tx.2015[[3, 1]] <- "Mangrove"
tx.2015$Year <- as.character(tx.2015$Year)
tx.2015_sum <- aggregate(tx.2015[, 2:12], by = list(LULC = tx.2015$LABEL), FUN = sum)
tx.2015_sum <- as.data.frame(t(tx.2015_sum), stringsAsFactors = FALSE)
tx.2015_sum$year <- "2015"
rownames(tx.2015_sum) <- paste0(rownames(tx.2015_sum), "_", tx.2015_sum$year)
colnames(tx.2015_sum) <- c("grassland/wetland", "mangrove", "open water", "other veg", "unclassified", "urban", "year")
tx.2015_sum <- tx.2015_sum[-1, ]
tx.2015_sum <- tx.2015_sum[, -5]


###LA
la.1980 <- la.1980[, -1]
la.1980 <- la.1980[1:4, ]
la.1980$LABEL <- as.character(la.1980$LABEL)
class(la.1980$LABEL)
la.1980 <- as.data.frame(t(la.1980), stringsAsFactors = FALSE)
la.1980$year <- "1980"
rownames(la.1980) <- paste0(rownames(la.1980), "_", la.1980$year)
colnames(la.1980) <- c("unclassified", "mangrove", "terrestrial", "water", "year")
la.1980 <- la.1980[-1, ]
la.1980 <- la.1980[-5, ]

la.1985 <- la.1985[, -1]
la.1985 <- la.1985[1:4, ]
la.1985$LABEL <- as.character(la.1985$LABEL)
class(la.1985$LABEL)
la.1985 <- as.data.frame(t(la.1985), stringsAsFactors = FALSE)
la.1985$year <- "1985"
rownames(la.1985) <- paste0(rownames(la.1985), "_", la.1985$year)
colnames(la.1985) <- c("unclassified", "mangrove", "terrestrial", "water", "year")
la.1985 <- la.1985[-1, ]
la.1985 <- la.1985[-5, ]

la.1990 <- la.1990[, -1]
la.1990 <- la.1990[1:4, ]
la.1990$LABEL <- as.character(la.1990$LABEL)
class(la.1990$LABEL)
la.1990 <- as.data.frame(t(la.1990), stringsAsFactors = FALSE)
la.1990$year <- "1990"
rownames(la.1990) <- paste0(rownames(la.1990), "_", la.1990$year)
colnames(la.1990) <- c("unclassified", "mangrove", "terrestrial", "water", "year")
la.1990 <- la.1990[-1, ]
la.1990 <- la.1990[-5, ]

la.1995 <- la.1995[, -1]
la.1995 <- la.1995[1:4, ]
la.1995$LABEL <- as.character(la.1995$LABEL)
class(la.1995$LABEL)
la.1995 <- as.data.frame(t(la.1995), stringsAsFactors = FALSE)
la.1995$year <- "1995"
rownames(la.1995) <- paste0(rownames(la.1995), "_", la.1995$year)
colnames(la.1995) <- c("unclassified", "mangrove", "terrestrial", "water", "year")
la.1995 <- la.1995[-1, ]
la.1995 <- la.1995[-5, ]

la.2000 <- la.2000[, -1]
la.2000 <- la.2000[1:4, ]
la.2000$LABEL <- as.character(la.2000$LABEL)
class(la.2000$LABEL)
la.2000 <- as.data.frame(t(la.2000), stringsAsFactors = FALSE)
la.2000$year <- "2000"
rownames(la.2000) <- paste0(rownames(la.2000), "_", la.2000$year)
colnames(la.2000) <- c("unclassified", "mangrove", "terrestrial", "water", "year")
la.2000 <- la.2000[-1, ]
la.2000 <- la.2000[-5, ]

la.2005 <- la.2005[, -1]
la.2005 <- la.2005[1:4, ]
la.2005$LABEL <- as.character(la.2005$LABEL)
class(la.2005$LABEL)
la.2005 <- as.data.frame(t(la.2005), stringsAsFactors = FALSE)
la.2005$year <- "2005"
rownames(la.2005) <- paste0(rownames(la.2005), "_", la.2005$year)
colnames(la.2005) <- c("unclassified", "mangrove", "terrestrial", "water", "year")
la.2005 <- la.2005[-1, ]
la.2005 <- la.2005[-5, ]

la.2010 <- la.2010[, -1]
la.2010 <- la.2010[1:4, ]
la.2010$LABEL <- as.character(la.2010$LABEL)
class(la.2010$LABEL)
la.2010 <- as.data.frame(t(la.2010), stringsAsFactors = FALSE)
la.2010$year <- "2010"
rownames(la.2010) <- paste0(rownames(la.2010), "_", la.2010$year)
colnames(la.2010) <- c("unclassified", "mangrove", "terrestrial", "water", "year")
la.2010 <- la.2010[-1, ]
la.2010 <- la.2010[-5, ]

la.2015 <- la.2015[, -1]
la.2015$LABEL <- as.character(la.2015$LABEL)
class(la.2015$LABEL)
la.2015[[1, 1]] <- "Mangrove"
la.2015[[2, 1]] <- "Mangrove"
la.2015[[7, 1]] <- "Mangrove"
la.2015$Year <- as.character(la.2015$Year)
la.2015_sum <- aggregate(la.2015[, 2:3], by = list(LULC = la.2015$LABEL), FUN = sum)
la.2015_sum <- as.data.frame(t(la.2015_sum), stringsAsFactors = FALSE)
la.2015_sum$year <- "2015"
rownames(la.2015_sum) <- paste0(rownames(la.2015_sum), "_", la.2015_sum$year)
colnames(la.2015_sum) <- c("unclassified", "grassland_wetland", "mangrove", "water", "other_veg", "urban", "year")
la.2015_sum <- la.2015_sum[-1, ]
la.2015_sum$urban <- as.numeric(la.2015_sum$urban)
la.2015_sum$mangrove <- as.numeric(la.2015_sum$mangrove)
la.2015_sum$water <- as.numeric(la.2015_sum$water)
la.2015_sum$other_veg <- as.numeric(la.2015_sum$other_veg)
la.2015_sum$grassland_wetland <- as.numeric(la.2015_sum$grassland_wetland)
la.2015_sum$unclassified <- as.numeric(la.2015_sum$unclassified)
la.2015_mut <- mutate(la.2015_sum, terrestrial = grassland_wetland + urban + other_veg)
la.2015_mut <- la.2015_mut[, c("unclassified", "mangrove", "terrestrial", "water", "year")]
rownames(la.2015_mut) <- rownames(la.2015_sum)

#Can bind the datasets using the unique col name for site + years 
# mangrove_bind <- rbind(fl.1980_sum, fl.1985_sum, fl.1990_sum, fl.1995_sum, fl.2000_sum, fl.2005_sum, fl.2010_sum, fl.2015_sum, tx.1980_sum, tx.1985_sum, tx.1990_sum, tx.1995_sum, tx.2000_sum, tx.2005_sum, tx.2010_sum, tx.2015_sum)
# mangrove_bind <- rownames_to_column(mangrove_bind, var = "unique_id")
# mangrove_cbc <- merge(cbc, mangrove_bind, by = "unique_id", all = T)


fl.bind <- rbind(fl.1980_sum, fl.1985_sum, fl.1990_sum, fl.1995_sum, fl.2000_sum, fl.2005_sum, fl.2010_sum, fl.2015_sum)
tx.bind <- rbind(tx.1980_sum, tx.1985_sum, tx.1990_sum, tx.1995_sum, tx.2000_sum, tx.2005_sum, tx.2010_sum, tx.2015_sum) 
la.bind <- rbind(la.1980, la.1985, la.1990, la.1995, la.2000, la.2005, la.2010, la.2015_mut)
fl.m <- fl.bind[, c("mangrove", "year")]
tx.m <- tx.bind[, c("mangrove", "year")]
la.m <- la.bind[, c("mangrove", "year")]

mangrove.bind <- rbind(fl.m, la.m, tx.m)
mangrove.bind <- rownames_to_column(mangrove.bind, var = "unique_id")
rownames(mangrove.bind) <- mangrove.bind$unique_id
# fl.bind$UNK <- as.numeric(fl.bind$UNK)

fl.bind$`grassland/wetland` <- as.numeric(fl.bind$`grassland/wetland`)

fl.bind$mangrove <- as.numeric(fl.bind$mangrove)

fl.bind$`open water` <- as.numeric(fl.bind$`open water`)

fl.bind$`other veg` <- as.numeric(fl.bind$`other veg`)

fl.bind$urban <- as.numeric(fl.bind$urban)

colnames(fl.bind)[1] <- "grassland_wetland"
colnames(fl.bind)[3] <- "open_water"
colnames(fl.bind)[4] <- "other_veg"

# fl.bind <- rbind(fl.1980_sum, fl.1985_sum, fl.1990_sum, fl.1995_sum, fl.2000_sum, fl.2005_sum, fl.2010_sum, fl.2015_sum)

# fl.bind$UNK <- as.numeric(fl.bind$UNK)
# colnames(mangrove_bind)[2] <- "grassland_wetland"
# colnames(mangrove_bind)[4] <- "open_water"
# colnames(mangrove_bind)[5] <- "other_veg"
# 
# mangrove_bind$grassland_wetland <- as.numeric(mangrove_bind$grassland_wetland)
# 
# mangrove_bind$mangrove <- as.numeric(mangrove_bind$mangrove)
# 
# mangrove_bind$open_water <- as.numeric(mangrove_bind$open_water)
# 
# mangrove_bind$other_veg<- as.numeric(mangrove_bind$other_veg)
# 
# mangrove_bind$urban <- as.numeric(mangrove_bind$urban)



fl.bind$percent_mangrove <- 0
fl.bind <- tibble::rownames_to_column(fl.bind, var = "Sites")
fl_reduced <- tidyr::separate(fl.bind, Sites, c("Abbrev", "yr"))
fl_reduced$unqiue_id <- paste0(fl_reduced$Abbrev, "_", fl_reduced$yr)
# fl.bind$total <- as.data.frame(rowSums(fl.bind))
# fl.bind$percent_mangrove <- fl.bind$mangrove / fl.bind$total
# fl.bind$percent_mangrove <- fl.bind$percent_mangrove * 100
# options(scipen = 999999)
# fl.bind$percent_mangrove <- round(fl.bind$percent_mangrove, digits = 2)


fl_reduced <- mutate(fl_reduced, 
              lulc_tot = grassland_wetland + mangrove + open_water + other_veg + urban) %>%
              mutate(percent_mangrove = (mangrove / lulc_tot) * 100)
options(scipen = 999)

##Have the total mangrove percentage based on Chandra's data
##Generate a %change column
##Plot the site beta diversity against %change?

# fl.bind.sort <- fl.bind[order(fl.bind$Sites),]
# fl.bind.sort$yr <- seq(1980, 2015, 5)

##Entire mangrove dataset 
# mangrove_bind <- mutate(mangrove_bind, tot_lulc = grassland_wetland + mangrove + open_water + other_veg + urban) %>%
#   mutate(pct_mangrove = (mangrove / tot_lulc) * 100)
# options(scipen = 999)
# mangrove_reduced <- mangrove_bind[, c("unique_id", "tot_lulc", "pct_mangrove")]



##Binning mangrove years to merge for site bins 
# years <- data.frame(yr = seq(1980, 2015, 5), yr_bin = NA)
# years$yr_bin[years$yr == 1980] <- 1
# years$yr_bin[years$yr == 1985] <- 2
# years$yr_bin[years$yr == 1990] <- 3
# years$yr_bin[years$yr == 1995] <- 4
# years$yr_bin[years$yr == 2000] <- 5
# years$yr_bin[years$yr == 2005] <- 6
# years$yr_bin[years$yr == 2010] <- 7
# years$yr_bin[years$yr == 2015] <- 8
# 
# 
# rownames(fl.bind) <- fl.bind$Sites
# fl.bind$year <- as.numeric(fl.bind$year)
# fl.binned <- merge(fl.bind, years, by.x = "year", by.y = "yr")
# rownames(fl.binned) <- fl.binned$Sites


# years <- data.frame(yr = seq(1980, 2015, 5), yr_bin = NA)
mangrove.bind$yr_bin[mangrove.bind$year == 1980] <- 1
mangrove.bind$yr_bin[mangrove.bind$year == 1985] <- 2
mangrove.bind$yr_bin[mangrove.bind$year == 1990] <- 3
mangrove.bind$yr_bin[mangrove.bind$year == 1995] <- 4
mangrove.bind$yr_bin[mangrove.bind$year == 2000] <- 5
mangrove.bind$yr_bin[mangrove.bind$year == 2005] <- 6
mangrove.bind$yr_bin[mangrove.bind$year == 2010] <- 7
mangrove.bind$yr_bin[mangrove.bind$year == 2015] <- 8

mangrove_binned  <- mangrove.bind

mangrove_binned <- tidyr::separate(mangrove_binned, unique_id, c("Abbrev", "yr"))
mangrove_binned$unique_id <- paste0(mangrove_binned$Abbrev, "_", mangrove_binned$yr)

# fl.binned <- tidyr::separate(fl.binned, Sites, c("Abbrev", "yr"))
# fl.binned$unique_id <- paste0(fl.binned$Abbrev, "_", fl.binned$yr)

#dplyr work-around 
# mangrove_binned <- mangrove_binned[, -9:-10]
mangrove_binned <- tibble::as_tibble(mangrove_binned)
mangrove_binned$mangrove <- as.numeric(mangrove_binned$mangrove)
# mangrove_only <- mangrove_binned[,c("Abbrev", "unique_id", "yr_bin", "yr", "mangrove")]
mangrove_binned <- mangrove_binned %>% mutate(pct.mangrove = (mangrove / 516047) * 100) %>% group_by(Abbrev) %>% mutate(lag = lag(mangrove, order_by = unique_id)) %>% 
  mutate(pct.change = ((mangrove - lag) / lag) *100) %>%  
  group_by(Abbrev) %>% 
  mutate(diff.from.first = mangrove - first(mangrove), 
         base.pct.diff = abs(diff.from.first)/((mangrove + first(mangrove) / 2)) * 100, 
         base.pct.nod = (diff.from.first/first(mangrove)) *100,
         change = pct.mangrove - first(pct.mangrove)) 
mangrove_binned[is.na(mangrove_binned)] <- 0

#Prepare for merging with CBC data
mangrove_binned$unqiue_bin <- paste0(mangrove_binned$Abbrev, "_", mangrove_binned$yr_bin)
mangrove_beta <- merge(site_bin_agg, mangrove_binned, by.x = "unique_id", by.y = "unqiue_bin")
mangrove_beta <- mutate(mangrove_beta, change_sq = change^2)
#Plots
plot(mangrove_beta$beta.bc ~ mangrove_beta$base.pct.diff)
plot(mangrove_beta$beta.jac ~ mangrove_beta$base.pct.diff)
sub_national <- cbc_df[, c("abbrev", "subnational_code")]
sub_national$abbrev <- as.character(sub_national$abbrev)
sub_national$subnational_code <- as.character(sub_national$subnational_code)
sub_national <- unique(sub_national)
latitude <- cbc_df[, c("abbrev", "latitude")]
latitude <- unique(latitude)
mangrove_beta <- merge(mangrove_beta, sub_national, by.x = "Abbrev.x", by.y = "abbrev")
mangrove_beta <- separate(mangrove_beta, subnational_code, c("country", "state"), sep = "-")
mangrove_beta <- merge(mangrove_beta, latitude, by.x = "Abbrev.x", by.y = "abbrev")
mangrove_beta$latitude_bin[mangrove_beta$latitude >= 27.50] <- "Northern Gulf" 
mangrove_beta$latitude_bin[mangrove_beta$latitude <= 27.49] <- "Southern Gulf"
colnames(mangrove_beta)[colnames(mangrove_beta) == "latitude_bin"] <- "Latitude"
plot1 <- ggplot(data = mangrove_beta, aes(x = base.pct.diff, y = beta.bc)) + geom_point(aes(color = Latitude), pch = 19, size = 8) + 
  labs(x = "Percent Difference in Mangrove Cover", y = expression(paste("Taxonomic ", beta, "-diversity"))) + scale_color_manual(values =  c("#88CCEE", "#CC6677")) +
  theme(legend.text = element_text(size = 18)) + theme(legend.title = element_blank()) +  theme(text = element_text(size=18)) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18)) +
  theme(legend.position = c(.65, .85))
plot1
plot2 <- ggplot(data = mangrove_beta, aes(x = base.pct.nod, y = beta.bc)) + geom_point(aes(color = Latitude), pch = 19, size = 8) + 
  labs(x = "Percent Change in Mangrove Cover", y = expression(paste("Taxnomic ", beta, "-diversity"))) + scale_color_manual(values =  c("#88CCEE", "#CC6677")) +
  theme(legend.text = element_text(size = 18)) + theme(legend.title = element_blank()) +  theme(text = element_text(size=18)) + theme(axis.text.x = element_text(size = 18), axis.text.y = element_text(size = 18))
plot2
#Florida only dataset 
fl.binned <- mutate(fl.binned, lag = lag(mangrove, order_by = Abbrev)) %>% 
  mutate(pct.change = ((mangrove - lag) / lag) * 100) %>% 
  group_by(Abbrev) %>% 
  mutate(diff.from.first = mangrove - first(mangrove), 
         base.pct.diff = diff.from.first/first(mangrove))
mod <- lmer(data = mangrove_beta, beta.bc ~ base.pct.diff + yr + (1|Abbrev.x))
summary(mod)
coeffs <- coef(summary(mod))
p <- pnorm(abs(coeffs[, "t value"]), lower.tail = FALSE) * 2
cbind(coeffs, "p value" = round(p,3))
coeff_abbrev <- coef(mod)$Abbrev
# mangrove_only <- merge(mangrove_only, fl_reduced, by = "unique_id")
#These calculations are likely not necessary. Raw % change is likley most appropriate 
master_mangrove <- mutate(mangrove_binned, lag2 = lag(pct_mangrove, order_by = Abbrev)) %>%
  mutate(pct.change.mangrove = (pct_mangrove - lag2) / lag2) %>% 
  group_by(Abbrev) %>% 
  mutate(pct.diff.from.first = pct_mangrove - first(pct_mangrove),
         final.base.diff = pct.diff.from.first/first(pct_mangrove))

# is.nan.data.frame <- function(x)
#   do.call(cbind, lapply(x, is.nan))
# 
# mangrove_binned[is.nan(mangrove_binned)] <- 0

# mangrove_filter <- filter(mangrove_only, yr %in% c(1980, 2015))
#master_mangrove$yr <- as.numeric(master_mangrove$yr)

#Calculate the slope of change in percent mangrove 
n.sites.mangrove <- length(unique(mangrove_binned$Abbrev))
site.list.mangrove <- as.character(unique(mangrove_binned$Abbrev))

slopes_sites_mangrove <- data.frame(site.list.mangrove, slope = NA)

for (i in 1:n.sites.mangrove){
  site.temp.mangrove <- site.list.mangrove[i]
  mangrove.binned.temp <- mangrove_binned[mangrove_binned$Abbrev == site.temp.mangrove,]
  if (nrow(mangrove.binned.temp) > 1){
    lm.temp.mangrove <- lm(mangrove.binned.temp$base.pct.diff ~ mangrove.binned.temp$yr)
    slope.temp.mangrove <- summary(lm.temp.mangrove)$coefficients[2,1]
    slopes_sites_mangrove[i,2] <- slope.temp.mangrove
  }
}

slopes <- merge(slopes_sites_mangrove, slopes_sites_beta, by.x = "site.list.mangrove", by.y = "site.list")
plot(slopes$beta.slope ~ slopes$slope)
# colnames(slopes_sites)[1] <- "Abbrev"
# 
# mangrove_test <- merge(fl_reduced, slopes_sites, by = "Abbrev")
# mangrove_test <- filter(mangrove_test, yr == 1980)
# 
# colnames(mangrove_test)[1] <- "site.list"
# 
# mangrove_test <- merge(mangrove_test, slopes_sites_beta, by = "site.list")
# mangrove_test <- mangrove_test[complete.cases(mangrove_test), ]

fl.binned$site_bin <- paste0(fl.binned$Abbrev, "_", fl.binned$yr_bin)
fl.sub <- fl.binned[, c("Abbrev", "base.pct.diff", "site_bin", "mangrove", "pct.change", "yr_bin")]
fl.merge <- merge(fl.sub, fl_beta_sub, by = "site_bin")

mod <- lmer(data = fl.merge, beta_bin ~ base.pct.diff + yr_bin + (1|Abbrev))
summary(mod)
coeffs <- coef(summary(mod))
p <- pnorm(abs(coeffs[, "t value"]), lower.tail = FALSE) * 2
cbind(coeffs, "p value" = round(p,3))
coeff_abbrev <- coef(mod)$Abbrev

ggplot(fl.binned, aes(x = unique_id, y = pct.change)) + geom_smooth(method = lm, se = FALSE)

#plotting
plot_mangrove <- ggplot(fl.binned, aes(x = yr, y = pct.change, group = Abbrev)) + 
                                           #colour = factor(subnational_code, 
                                                          # labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                 #geom_point(size = 3) +
                 #geom_line()+
                 geom_smooth(method = lm, se = FALSE, aes(x = yr, y = pct.change, group = Abbrev)) + 
                                                         # colour = factor(subnational_code, 
                                                          #                labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")))) +
                 xlab("Year") +
                 ylab("Mangrove Cover") +
                 #labs(colour = "States") +
                 #scale_colour_manual(name = "States", values = colour, labels = c("Alabama", "Florida", "Louisiana", "Mississippi", "Texas")) +
                 #scale_x_continuous(limits = c(1980,2017), expand = c(0, 0), 
                 #breaks = c(2006:2016), labels = c("2006", "", "2008", "", "2010", "","2012", "", "2014", "", "2016" )) + 
                 #scale_y_continuous(limits = c(0,1), expand = c(0, 0), breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
                 theme_bw() +
                 theme(axis.line = element_line(colour = "black", size =1.2),
                       axis.text.x = element_text(size = 14),
                       axis.text.y = element_text(size = 14),
                       axis.title.x = element_text(vjust = -1, size = 14),
                       axis.title.y = element_text(vjust = 1.5, size = 14),
                       axis.ticks = element_line(size = 1.2),
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(),
                       plot.margin = unit(c(1,1,2,2), "lines"),
                       text = element_text(size=14))

ggdraw() + 
  draw_plot(plot_mangrove + theme(legend.justification = "top"), 
            x = 0, y = 0, width = .9, height = 1) +
  #draw_plot(gghist, x = 0.77, y = .025, width = .2, height = .75, scale = 1) 


  
  
  
#Attempted loop, ask Lauren for curiosity sake
# mangrove.loop <- as.data.frame(mangrove_binned[, c("unique_id", "mangrove")])
# rownames(mangrove.loop) <- mangrove.loop[,1]
# mangrove.loop[,1] <- NULL

# mangrove_cast <- dcast(mangrove_binned, Abbrev~ yr, value.var = "mangrove")
# rownames(mangrove_cast) <- mangrove_cast[,1]
# mangrove_cast <- mangrove_cast[, -1]

# site.list <- unique(mangrove_binned[,c("Abbrev", "unique_id")])
# 
# sites <- as.character(unique(site.list$Abbrev))
# 
# #create data.frame to hold beta output
# mangrove.out <- data.frame(second_bin = site.list$unique_id, change.from.baseline = NA)

# # Loop through to compare the baseline mangrove cover to the new mangrove cover 
# pct <- function(u){
#   u/lag(u)
# }
# 
# for (u in 1:length(sites)){
#   site.temp <- sites[u]
#   surveys.temp <- site.list$unique_id[site.list$Abbrev == site.temp]
#   
#   mangrove_temp <- as.matrix(mangrove.loop[row.names(mangrove.loop) %in% surveys.temp, ])
#   
#   
#   change.temp <- pct(mangrove_temp)
#   
#   for (y in 1:length(change.temp)){
#     temp.unique.id <- row.names(surveys.temp)[change.temp]
#     mangrove.out[mangrove.out$second_bin == temp.unique.id,"change.from.baseline"] <- change.temp[y,1]
#   }
# }





