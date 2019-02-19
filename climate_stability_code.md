Climate stability code
================

<br/><br/> **Fordham, D., Brown, S., Wigley, T., Rahbek, C., 2019. Cradles of diversity: unlikely relics of regional climate stability. *Current Biology* (in press)** <br/><br/>

### R Code

Code used to identify stable cells and perform the sensitivity analysis

``` r
library(tidyverse)
library(gtools)
library(raster)
library(sf)
library(doParallel)

# load in the trace-21ka windows considered extremely unstable
load("./RObjects/unstableWindows.Rdata", verbose = TRUE)

# Calculate the amount of time that a cell is considered "stable" inside unstable windows

# Directory to trends
palDir <- "./path/to/trace/trends/"

# Trend
trendFiles <- mixedsort(
  list.files(palDir, "_TrendStand.grd$", full.names = TRUE),
  decreasing = TRUE
)
head(trendFiles)
tail(trendFiles)

# Stack. Files are single band and identical so can use stack(x, quick = TRUE)
trendStack <- raster::stack(trendFiles, quick = TRUE)
trendStack <- setZ(trendStack,
                   z = seq(-21000, by = 1, length.out = nlayers(trendStack)),
                   name = "Years BP")
names(trendStack) <- paste0("BP_", sub("-", "", getZ(trendStack)))
trendStack

# subset the trend stack to the extreme centuries of interest
trendStack <- trendStack[[unstableCents]]

# Variability
varFiles <- mixedsort(
  list.files(palDir, "Variability_.*grd$", full.names = TRUE),
  decreasing = TRUE
)
head(varFiles)
tail(varFiles)

varStack <- raster::stack(varFiles, quick = TRUE)
varStack <- setZ(varStack,
                   z = seq(-21000, by = 1, length.out = nlayers(varStack)),
                   name = "Years BP")
names(varStack) <- paste0("BP_", sub("-", "", getZ(varStack)))
varStack

# subset the trend stack to the extreme centuries of interest
varStack <- varStack[[unstableCents]]
compareRaster(trendStack, varStack)

# Paleo masks
## Land/Sea masks
masks <- mixedsort(list.files("path/to/paleo/masks/",
                              "\\.grd$", full.names = TRUE), decreasing = TRUE)
masks <- masks[-length(masks)]
paleoMasks <- stack(masks, quick = TRUE)

# Define the output directory
outDir <- "./OutputRasters/sensAnalysis/"

# Percentiles of interest
# 5% either side of P10 and P25 for sensitivity testing
timePercP10 <- seq(10 - (0.1*10), 10 + (0.1*10), by = 0.1)/100
timePercP25 <- seq(25 - (0.1*25), 25 + (0.1*25), length.out = length(timePercP10))/100
timePerc <- c(timePercP10, timePercP25)
timePerc

# Loop through the percentile and layers
for (pVal in timePerc) {
  pValtxt <- paste0("P", as.character(pVal*100))
  message("Processing ", pValtxt)
  cl <- detectCores()
  cl <- cl - 1
  cl <- makeCluster(cl)
  registerDoParallel(cl)
  foreach(
    i = 1:nlayers(trendStack), .inorder = TRUE,
    .packages = c("tidyverse", "raster"),
    .verbose = FALSE
  ) %dopar% {
    # Determine nearest century and grab correct mask
    timeCent <- plyr::round_any(x = getZ(trendStack)[i], 100, f = round)
    centMask <- paleoMasks[[paste0(sub("-", "BP_", timeCent), "_mask")]]
    # Remove Antarctica
    palCover <- raster(res = 2.5, xmn = -180, xmx = 180, ymn = -90, ymx = -60)
    palAnt <- crop(centMask, palCover)
    palAnt[palAnt == 1] <- 3
    palAnt[palAnt != 3] <- NA
    palAnt <- extend(palAnt, centMask)
    coverAnt <- cover(palAnt, centMask)
    centMask <- coverAnt
    rm(palCover, palAnt, coverAnt)
    trend <- trendStack[[i]]
    var <- varStack[[i]]
    trendYear <- as.numeric(sub("BP_", "-", names(trend)))
    varYear <- as.numeric(sub("BP_", "-", names(var)))
    # Check years are identical
    if (!identical(trendYear, varYear)) {
      stop(paste0("Trend year = ", trendYear, "\n", "Var year = ", varYear))
    }
    # Stack and rename
    rasStack <- stack(trend, var)
    names(rasStack) <- c("Trend", "Var")
    # Mask
    landMask <- centMask
    landMask[landMask != 1] <- NA
    oceanMask <- centMask
    oceanMask[oceanMask != 0] <- NA
    rasTerr <- mask(rasStack, landMask)
    rasOcean <- mask(rasStack, oceanMask)
    # Quantiles
    lowTrendTerr <- quantile(rasTerr[["Trend"]],
      prob = pVal,
      na.rm = TRUE, type = 8
    )
    lowVarTerr <- quantile(rasTerr[["Var"]],
      prob = pVal,
      na.rm = TRUE, type = 8
    )
    lowTrendOcean <- quantile(rasOcean[["Trend"]],
      prob = pVal,
      na.rm = TRUE, type = 8
    )
    lowVarOcean <- quantile(rasOcean[["Var"]],
      prob = pVal,
      na.rm = TRUE, type = 8
    )
    # create binary layers for trend, variability
    # 1 <= lowTrend/Var; 0 > lowTrend/Var
    rasTerr[["Trend"]] <- ifelse(
      values(rasTerr[["Trend"]]) <= lowTrendTerr, 1, 0
    )
    rasTerr[["Var"]] <- ifelse(
      values(rasTerr[["Var"]]) <= lowVarTerr, 1, 0
    )
    rasOcean[["Trend"]] <- ifelse(
      values(rasOcean[["Trend"]]) <= lowTrendOcean, 1, 0
    )
    rasOcean[["Var"]] <- ifelse(
      values(rasOcean[["Var"]]) <= lowVarOcean, 1, 0
    )
    # Combine trend and variability to create a stable climate raster for the given window
    stableTV <- stackApply(stack(
      rasTerr[["Trend"]], rasTerr[["Var"]],
      rasOcean[["Trend"]], rasOcean[["Var"]]
    ),
    indices = c(1, 1, 2, 2), fun = sum, na.rm = TRUE
    )
    stableTV[] <- ifelse(values(stableTV) != 2, NA, 1)
    stableTV <- cover(stableTV[[1]], stableTV[[2]])
    stableTV[is.na(stableTV)] <- 0
    names(stableTV) <- paste0("StableClim", sub("-", "BP_", trendYear))
    rasOut <- paste0(outDir, pValtxt, "/")
    if (!dir.exists(rasOut)) {
        dir.create(rasOut, recursive = TRUE)
    }
    writeRaster(stableTV,
      format = "raster", datatype = "LOG1S",
      file = paste0(rasOut, "stableTV_",
        sub("-", "BP_", trendYear)
      )
    )
    rm(
      timeCent, centMask, trend, var, trendYear, varYear, rasStack, lowTrendTerr,
      lowTrendOcean, lowVarTerr, lowVarOcean, rasTerr, rasOcean, stableTV,
      landMask, oceanMask, rasOut
    )
  }
  stopCluster(cl)
}

# Subset the directories based on these values
# return 21 positions - 10 either side of the median.
# length(P10) = 21, so a subset of P25 needs to equal 21
P25Sub <- paste0("P", sub("\\.0", "", as.character(
  timePercP25[c(which(timePercP25 == median(timePercP25))-10):c(which(timePercP25 == median(timePercP25))+10)]*100)))
P25Sub
P10Sub <- paste0("P", sub("\\.0", "", as.character(timePercP10*100)))
P10Sub

# load in the directories
rDirs <- mixedsort(list.dirs("./OutputRasters/sensAnalysis/", full.names = TRUE),
                   decreasing = FALSE)[-1]
rDirsP25 <- rDirs[str_detect(rDirs,
                             pattern = paste(P25Sub, collapse = "|"))]
rDirsP10 <- rDirs[str_detect(rDirs,
                         pattern = paste(P10Sub, collapse = "|"))]
rDirsSub <- c(rDirsP10, rDirsP25)
rDirsSub

for (d in rDirsSub) {
  message("Directory: ", d)
  p <- tail(unlist(str_split(d, "P")), 1)
  message("Percentile: ", p)
  # Load in all the binary stability files
  stableFiles <- mixedsort(list.files(d,
    pattern = "\\.grd$",
    recursive = TRUE,
    full.names = TRUE
  ), decreasing = TRUE)
  stableStack <- readAll(stack(stableFiles, quick = TRUE))
  # calculate % time in stable conditions
  calc(stableStack,
    fun = function(x, ...) {
      sum(x == 1, na.rm = TRUE) / length(x)
    },
    na.rm = TRUE,
    progress = "text",
    filename = paste0(d, "/paleo_prop_ts_stable_p", p, ".tif"),
    options = c("COMPRESS=NONE"), overwrite = TRUE,
    format = "GTiff", datatype = "FLT4S"
  )
}

stableFiles <- mixedsort(
  list.files("./OutputRasters/sensAnalysis/",
    pattern = "\\.tif$",
    recursive = TRUE,
    full.names = TRUE
  ),
  decreasing = TRUE
)

# Create lists for P10 and P25 files
stable <- stableFiles[str_detect(stableFiles, "prop_ts_stable_p")]
stable
stableP10 <- mixedsort(stable[str_detect(stable,
                               paste(sub("P", "p", P10Sub), collapse = "|"))])
stableP10
stableP25 <- mixedsort(stable[str_detect(stable,
                               paste(sub("P", "p", P25Sub), collapse = "|"))])
stableP25

# Create raster stacks
## P10
stableStackP10 <- readAll(stack(stableP10, quick = TRUE))
plot(stableStackP10)
stableStackP10[stableStackP10 == 0] <- NA # convert 0 to NA
plot(stableStackP10)
## P25
stableStackP25 <- readAll(stack(stableP25, quick = TRUE))
stableStackP25
plot(stableStackP25)
stableStackP25[stableStackP25 == 0] <- NA
plot(stableStackP25)

# Determine the mean and SD using the expanded percentiles
P10 <- stableStackP10[["paleo_prop_ts_stable_p10"]]
meanP10 <- calc(stableStackP10, fun = mean, na.rm = TRUE)
sdP10 <- calc(stableStackP10, fun = sd, na.rm = TRUE)

P25 <- stableStackP25[["paleo_prop_ts_stable_p25"]]
meanP25 <- calc(stableStackP25, fun = mean, na.rm = TRUE)
sdP25 <- calc(stableStackP25, fun = sd, na.rm = TRUE)

# Convert to dataframe and plot
P10DF <- stack(P10, meanP10, sdP10) %>%
  as.data.frame(., xy = TRUE) %>%
  tbl_df() %>%
  rename("Time" = 3, "Mean" = layer.1, "SD" = layer.2) %>%
  gather(., key = Layer, value = Value, -c(x, y)) %>%
  mutate(Thresh = "P10")
P10DF

P25DF <- stack(P25, meanP25, sdP25) %>%
  as.data.frame(., xy = TRUE) %>%
  tbl_df() %>%
  rename("Time" = 3, "Mean" = layer.1, "SD" = layer.2) %>%
  gather(., key = Layer, value = Value, -c(x, y)) %>%
  mutate(Thresh = "P25")
P25DF

threshDF <- bind_rows(P10DF, P25DF)
threshDF
# Reordering P10DF$Layer
threshDF$Layer <- factor(P10DF$Layer,
                      levels = c("Time",
                                 "Mean", "SD"))
## Reordering threshDF$Thresh
threshDF$Thresh <- factor(threshDF$Thresh,
                          levels = c("P25", "P10"))

# Import a continuous colour ramp for the plot, and convert to RGB()
speedCols <- read_tsv("./plottingObjects/speedCols.txt") %>%
  arrange(ID) %>%
  dplyr::select(-ID, Alpha) %>%
  as.matrix(.) %>%
  rgb(., maxColorValue = 255)

conts <- read_sf("./plottingObjects/ne_110m_land.shp")
border <- read_sf("./plottingObjects/ne_50m_wgs84_bounding_box.shp")

sensPlot <- ggplot(threshDF %>% filter(Layer == "Mean")) +
  facet_wrap(~Thresh, nrow = 1,
             labeller = as_labeller(c("P25" = "a", "P10" = "b"))) +
  geom_raster(data = threshDF %>% filter(Layer == "Mean"),
              aes(x = x, y = y, fill = Value)) +
  scale_fill_gradientn(
    limits = c(0, 0.92),
    breaks = seq(0, 0.90, length.out = 10),
    labels = c(">0", "10", "20", "30", "40", "50", "60", "70", "80", "90+"),
    na.value = "transparent",
    colours = speedCols,
    guide =
      guide_legend(
        title = "Time spent in stable conditions",
        title.position = "top",
        title.theme = element_text(size = 8, colour = "black"),
        title.hjust = 0.5,
        label.position = "bottom",
        label.theme = element_text(size = 8, colour = "black"),
        label.hjust = 0.5,
        keywidth = unit(0.5, "cm"),
        keyheight = unit(0.5, "cm"),
        direction = "horizontal",
        nrow = 1)) +
  geom_sf(data = conts, fill = NA, colour = "black") +
  geom_sf(data = border, fill = NA, colour = "black") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.spacing.x = unit(0, "cm"),
        strip.text = element_text(size = 12,
                                  face = "bold",
                                  colour = "black", hjust = 0, vjust = 1)) +
  coord_sf(expand = FALSE, datum = NA) +
  scale_x_continuous(limits = c(-180,180)) +
  scale_y_continuous(limits = c(-90,90))
sensPlot
sdPlot <- ggplot(threshDF %>% filter(Layer == "SD")) +
  facet_wrap(~Thresh, nrow = 1,
             labeller = as_labeller(c("P25" = "c", "P10" = "d"))) +
  geom_raster(data = threshDF %>% filter(Layer == "SD"),
              aes(x = x, y = y, fill = Value)) +
  scale_fill_gradientn(
    limits = c(0, 0.161),
    breaks = seq(0, 0.16, by = 0.04),
    labels = seq(0, 0.16, by = 0.04),
    na.value = "transparent",
    colours = colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(101),
    guide =
      guide_colourbar(
        title = "S.D.",
        title.position = "top",
        title.theme = element_text(size = 8, colour = "black"),
        title.hjust = 0.5,
        label.position = "bottom",
        label.theme = element_text(size = 8, colour = "black"),
        label.hjust = 0.5,
        barwidth = unit(10, "cm"),
        barheight = unit(0.5, "cm"),
        nbin = 101,
        frame.colour = "black",
        frame.linewidth = 0.5,
        ticks.colour = "black",
        ticks.linewidth = 0.5,
        draw.ulim = TRUE,
        draw.llim = TRUE,
        direction = "horizontal",
        nrow = 1)
    ) +
  geom_sf(data = conts, fill = NA, colour = "black") +
  geom_sf(data = border, fill = NA, colour = "black") +
  theme_void() +
  theme(legend.position = "bottom",
        legend.spacing.x = unit(0, "cm"),
        strip.text = element_text(size = 12,
                                  face = "bold",
                                  colour = "black", hjust = 0, vjust = 1)) +
  coord_sf(expand = FALSE, datum = NA) +
  scale_x_continuous(limits = c(-180,180)) +
  scale_y_continuous(limits = c(-90,90))
sdPlot
sensComb <- cowplot::plot_grid(sensPlot, sdPlot,
                               align = "hv",
                               axis = "tblr",
                               nrow = 2,
                               labels = NULL)
sensComb
```

<br/><br/>

Code used to determine the % of cells that are considered stable/unstable

``` r
library(raster)
library(gtools)

# Load in the paleo masks
masks <- mixedsort(list.files("./path/to/paleo/masks/",
                              "\\.grd$",  full.names = TRUE),
                   decreasing = TRUE)
masks <- masks[-length(masks)]

# Calculate the proportion of time a cell is land
# Land == 1; Sea == 0
paleoMask <- stack(masks, quick = TRUE)
paleoMask <- calc(paleoMask, fun = function(x, ...) {
    sum(x == 1, na.rm = TRUE) / length(x)
  }, na.rm = TRUE)

# If a cell spends <= 50% of the time as land
# set as NA
paleoMask[paleoMask < 0.5] <- NA

# Change any non-NA cell to 1, NA cells to 0
paleoMask[!is.na(paleoMask)] <- 1
paleoMask[is.na(paleoMask)] <- 0

# Reclassify antarctica so it's not counted in % area calc
palCover <- raster(res = 2.5, xmn = -180, xmx = 180, ymn = -90, ymx = -60)
palAnt <- crop(paleoMask, palCover)
palAnt[palAnt == 1] <- 3
palAnt[palAnt != 3] <- NA
palAnt <- extend(palAnt, paleoMask)
coverAnt <- cover(palAnt, paleoMask)
paleoMask <- coverAnt
plot(paleoMask)

# Define the tropics as 20°N to 20°S
tropics <- extent(-180, 180, -20, 20)

# Tropical mask
tropMask <- crop(paleoMask, tropics)

# Load in the classified raster
stableRas <- readAll(raster("./OutputASCII/trace21_classified.asc"))
stableRas[] <- ifelse(stableRas[] == 5 | stableRas[] == 7, 1, NA)
crs(stableRas) <- crs(raster())

# Crop the stable raster to the tropics
stableTrop <- crop(stableRas, tropics)

# Determine % of cells considered stable for land and sea:
# Global
landPerStable <- round(cellStats(
  mask(stableRas, paleoMask, maskvalue = 0),
  function(i, ...) sum(i == 1, na.rm = TRUE)
) /
  cellStats(paleoMask, function(j, ...) sum(j == 1, na.rm = TRUE)), 3)
landPerStable * 100 # 22 % of land cells are considered stable globally

oceanPerStable <- round(cellStats(
  mask(stableRas, paleoMask, maskvalue = 1),
  function(i, ...) sum(i == 1, na.rm = TRUE)
) /
  cellStats(paleoMask, function(j, ...) sum(j == 0, na.rm = TRUE)), 3)
oceanPerStable * 100 # 19 % of ocean cells are considered stable in the tropics

# limit to the tropics
# 1) % of tropical terrestrial stable
# 2) % of tropical marine stable
landPerStable <- round(cellStats(
  mask(stableTrop, tropMask, maskvalue = 0),
  function(i, ...) sum(i == 1, na.rm = TRUE)
) /
  cellStats(tropMask, function(j, ...) sum(j == 1, na.rm = TRUE)), 3)
landPerStable * 100 # 68% of land cells are considered stable in the tropics

oceanPerStable <- round(cellStats(
  mask(stableTrop, tropMask, maskvalue = 1),
  function(i, ...) sum(i == 1, na.rm = TRUE)
) /
  cellStats(tropMask, function(j, ...) sum(j == 0, na.rm = TRUE)), 3)
oceanPerStable * 100 # 56% of ocean cells are considered stable in the tropics

# Determine min and max % of time spent in stable condtions
# 1) for land
# 2) for ocean

# Load in the p25 raster
p25 <- readAll(raster("./OutputASCII/paleo_prop_ts_stable_p25.asc"))
crs(p25) <- crs(raster())
# If p25 == 0, set NA
p25[p25 == 0] <- NA

stableLand <- mask(stableRas, paleoMask, maskvalue = 0)
stableOcean <- mask(stableRas, paleoMask, maskvalue = 1)

# plot(p25, col = rev(heat.colors(101)))

# Using all cells
landTimeStable <- c(minValue(mask(p25, paleoMask, maskvalue = 0)),
                    maxValue(mask(p25, paleoMask, maskvalue = 0)))
landTimeStable # min 0.0005; max 0.9218

oceanTimeStable <-  c(minValue(mask(p25, paleoMask, maskvalue = 1)),
                      maxValue(mask(p25, paleoMask, maskvalue = 1)))
oceanTimeStable # min 0.0005; max 0.6794

# Using only cells defined as <= P25 from Figure 1A
landTimeStable <- c(minValue(mask(p25, stableLand)),
                    maxValue(mask(p25, stableLand)))*100
landTimeStable # ~20% to 92%

oceanTimeStable <-  c(minValue(mask(p25, stableOcean)),
                      maxValue(mask(p25, stableOcean)))*100
oceanTimeStable # ~11% to 68%
```

<br/><br/>

Some plotting code

``` r
library(tidyverse)
library(raster)
library(rgdal)
library(classInt)
library(RColorBrewer)

# Define the Mollweide Projection
prjCode <- CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs")

# Import natural earth vector and raster data for plotting.
# Can be downloaded using the 'rnaturalearth' library, or directly
# from the natural earth website
conts <- readOGR("./plottingObjects/ne_50m_land.shp") %>%
  spTransform(., prjCode)
boundBox <- readOGR("./plottingObjects/ne_50m_wgs84_bounding_box.shp") %>%
  spTransform(., prjCode)

# Define classes for the classified TraCE-21ka trends and variability
ClassDict <- tibble(ID = seq(1:8),
                    Class = c(
                      "Low-Low", "Low-High", "High-Low", "High-High",
                      "P25",  "P75", "P10", "P90"))

# Import the geotiff and convert to categorical raster
class_pal <- raster("./plottingObjects/paleo_classified.tif") %>%
  ratify(.)

# Define the classes
rat <- levels(class_pal)[[1]] %>%
  left_join(., ClassDict, by = "ID")
levels(class_pal) <- rat

# Convert to dataframe
class_pal <- class_pal %>%
  as.data.frame(., xy = TRUE) %>%
  na.omit(.) %>%
  tbl_df()

# re-order the factors in the categorical raster
class_pal$paleo_classified_Class <- as.factor(class_pal$paleo_classified_Class)
levels(class_pal$paleo_classified_Class)
class_pal$paleo_classified_Class <- factor(
  class_pal$paleo_classified_Class,
  levels(class_pal$paleo_classified_Class)[c(5, 6, 4, 3, 2, 1, 7, 8)]
)
levels(class_pal$paleo_classified_Class)

# Resample, and reproject the terrain raster (basemap) to ~0.50° and turn to dataframe
terrain05deg <- raster("./plottingObjects/GRAY_50M_SR_OB.tif") %>%
  resample(., y = raster(res = 0.50)) %>%
  projectRaster(., crs = prjCode, method = "bilinear",
                res = 50100,
                over = TRUE) %>%
  as.data.frame(., xy = TRUE) %>%
  na.omit(.) %>%
  tbl_df(.)

# Import and resample the TraCE-21ka mask
traceMask <- raster("./plottingObjects/paleo_mask.tif") %>%
  projectRaster(., crs = prjCode, method = "ngb",
                res = 150000,
                over = TRUE) %>%
  extend(., y = raster(res = 150000, crs = prjCode,
                              ext = extent(terrain05deg)))

# Import a continuous colour ramp for the plot, and convert to RGB()
speedCols <- read_tsv("./plottingObjects/speedCols.txt") %>%
  arrange(ID) %>%
  dplyr::select(-ID, Alpha) %>%
  as.matrix(.) %>%
  rgb(., maxColorValue = 255)

# Import the raster of time spent in stable conditions
timeStableP10 <- raster("./OutputRasters/paleo_prop_ts_stable_p10.tif") %>%
  projectRaster(., crs = prjCode, method = "bilinear",
           res = 150000,
           over = TRUE) %>%
  as.data.frame(., xy = TRUE) %>%
  rename("value" = paleo_prop_ts_stable_p10) %>%
  # When resampling during projection some values will be < 0
  # replace with 0
  mutate(., value = ifelse(value < 0, 0, round(value, 2))) %>%
  # convert 0 to NA
  mutate(., value = ifelse(value == 0, NA, round(value, 2))) %>%
  # remove NA
  na.omit(.) %>%
  # Add an ID variable
  mutate(., layer = "P10") %>%
  tbl_df(.)

# Import the raster of time spent in stable conditions
timeStableP25 <- raster("./OutputRasters/paleo_prop_ts_stable_p25.tif") %>%
  projectRaster(., crs = prjCode, method = "bilinear",
                res = 150000,
                over = TRUE) %>%
  as.data.frame(., xy = TRUE) %>%
  rename("value" = paleo_prop_ts_stable_p25) %>%
  # When resampling during projection some values will be < 0
  # replace with 0
  mutate(., value = ifelse(value < 0, 0, round(value, 2))) %>%
  # convert 0 to NA
  mutate(., value = ifelse(value == 0, NA, round(value, 2))) %>%
  # remove NA
  na.omit(.) %>%
  # Add an ID variable
  mutate(., layer = "P25") %>%
  tbl_df(.)

timeStable <- bind_rows(timeStableP10, timeStableP25)

# Create a raster of cells in either P10 or P25 to use as a mask
stableMask <- class_pal %>%
  mutate("paleo_classified_Class" = as.character(.$paleo_classified_Class)) %>%
  mutate("P10" = ifelse(.$paleo_classified_Class == "P10", 1, NA),
         "P25" = ifelse(.$paleo_classified_Class == "P25", 1, NA)) %>%
  filter(paleo_classified_Class %in% c("P10", "P25")) %>%
  dplyr::select(-"paleo_classified_Class") %>%
  rasterFromXYZ(., crs = prjCode) %>%
  resample(., y = raster(res = 150000, crs = prjCode,
                         ext = extent(traceMask)),
           method = "ngb") %>%
  extend(., y = raster(res = 150000, crs = prjCode,
                       ext = extent(traceMask)))

# Define number of classes for discrete colours
no_classes <- 11

# Find the breaks
breaks <- classIntervals(
  var = timeStable$value,
  n = no_classes,
  intervalClosure = "right",
  style = "fixed",
  fixedBreaks = seq(0, 1, length.out = no_classes)
)

# Add the discrete classes as a variable
timeStable$classInt <- base::cut(timeStable$value,
  breaks = breaks$brks,
  right = TRUE,
  labels = as.character(breaks$brks)[-11]
)

# Define the breaks and labels
brks_scale <- levels(timeStable$classInt)
labels_scale <- rev(as.numeric(brks_scale) * 100)
labels_scale[1] <- ">90%"
labels_scale[10] <- ">0%"

# Break the continuous colour ramp into discrete classes
discreteSpeedCols <- colorRampPalette(speedCols)(no_classes)

# Custom ggplot2 theme
theme_map <- function(...) {
  theme_minimal(base_size = 8) +
    theme(
      text = element_text(color = "black", size = 8),
      axis.line = element_blank(),
      axis.text.x = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      panel.grid.major = element_line(color = "#ebebe5", size = 0.5),
      panel.grid.minor = element_blank(),
      plot.background = element_blank(),
      panel.background = element_rect(fill = NA, color = "black", size = 0.5),
      legend.background = element_blank(),
      legend.text = element_text(size = 8, colour = "black"),
      legend.title = element_text(size = 8, colour = "black"),
      panel.border = element_blank(),
      ...
    )
}

percentLabs <- c("P10" = "P10", "P25" = "P25")

#### Produce maps ####
# Panel A #

# Define the colours for the classes
classCols <- rev(brewer.pal(8, "Spectral"))

classPlot <- ggplot() +
  geom_raster(data = terrain05deg,
              aes(x = x, y = y, alpha = GRAY_50M_SR_OB),
              fill = "grey20") +
  scale_alpha(range = c(0.5, 0), guide = FALSE) +
  geom_tile(data = class_pal,
            aes(x = x, y = y, fill = paleo_classified_Class)) +
  scale_fill_manual(
    values = classCols,
    name = "Stability class",
    drop = TRUE,
    labels = c("P10", "P25", "Low-Low", "Low-High", "High-Low", "High-High",
               "P75", "P90"),
    na.value = "transparent",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5,
      label.vjust = 0.5,
      label.theme = element_text(colour = "black", angle = 45, size = 8),
      nrow = 1,
      byrow = TRUE,
      label.position = "bottom"
    )
  ) +
  geom_polygon(
    data = conts, aes(x = long, y = lat, group = group),
    fill = NA, colour = "black", size = 0.50, show.legend = FALSE
  ) +
  geom_polygon(
    data = boundBox, aes(x = long, y = lat, group = group),
    fill = NA, colour = "black", size = 0.50, show.legend = FALSE
  ) +
  coord_equal(expand = FALSE) +
  theme_map() +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  )
classPlot

# Panels B & C #

# Maps of time spent in globally stable conditions
p10Stable <- ggplot() +
  geom_raster(data = terrain05deg,
            aes(x = x, y = y, alpha = GRAY_50M_SR_OB), fill = "grey20") +
  scale_alpha(range = c(0.5, 0), guide = FALSE) +
  geom_tile(data = timeStable %>%
              dplyr::filter(layer == "P10"),
                 aes(x = x, y = y, fill = classInt)) +
  scale_fill_manual(
    values = discreteSpeedCols,
    breaks = rev(brks_scale),
    name = "Time spent in stable conditions (%)",
    drop = FALSE,
    labels = labels_scale,
    na.value = "transparent",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      reverse = TRUE,
      keywidth = 1,
      keyheight = 1,
      label.position = "bottom",
      label.theme = element_text(size = 8, colour = "black")
    )
  ) +
  geom_polygon(
    data = conts, aes(x = long, y = lat, group = group),
    fill = NA, colour = "black", size = 0.50, show.legend = FALSE
  ) +
  geom_polygon(
    data = boundBox, aes(x = long, y = lat, group = group),
    fill = NA, colour = "black", size = 0.50, show.legend = FALSE
  ) +
  coord_equal(expand = FALSE) +
  theme_map() +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  )

p25Stable <- ggplot() +
  geom_raster(data = terrain05deg,
              aes(x = x, y = y, alpha = GRAY_50M_SR_OB), fill = "grey20") +
  scale_alpha(range = c(0.5, 0), guide = FALSE) +
  geom_tile(data = timeStable %>%
              dplyr::filter(layer == "P25"),
            aes(x = x, y = y, fill = classInt)) +
  scale_fill_manual(
    values = discreteSpeedCols,
    breaks = rev(brks_scale),
    name = "Time spent in stable conditions (%)",
    drop = FALSE,
    labels = labels_scale,
    na.value = "transparent",
    guide = guide_legend(
      direction = "horizontal",
      title.position = "top",
      title.hjust = 0.5,
      label.hjust = 0.5,
      nrow = 1,
      byrow = TRUE,
      reverse = TRUE,
      keywidth = 1,
      keyheight = 1,
      label.position = "bottom",
      label.theme = element_text(size = 8, colour = "black")
    )
  ) +
  geom_polygon(
    data = conts, aes(x = long, y = lat, group = group),
    fill = NA, colour = "black", size = 0.50, show.legend = FALSE
  ) +
  geom_polygon(
    data = boundBox, aes(x = long, y = lat, group = group),
    fill = NA, colour = "black", size = 0.50, show.legend = FALSE
  ) +
  coord_equal(expand = FALSE) +
  theme_map() +
  theme(
    legend.direction = "horizontal",
    legend.position = "bottom",
    legend.text = element_text(colour = "black", size = 8),
    panel.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank()
  )
p10Stable
p25Stable

# Panel D #

# Import the % time spent in stable conditions
timeStableP10 <- raster("./OutputASCII/paleo_prop_ts_stable_p10.asc",
                        crs = crs(raster()))
timeStableP25 <- readAll(raster("./OutputASCII/paleo_prop_ts_stable_p25.asc",
                        crs = crs(raster())))

# Create a raster of cells in either P10 or P25 to use as a mask
stableMask <- readAll(raster("./OutputASCII/trace21_classified.asc",
                             crs = crs(raster())))
stableMask25 <- stableMask
# If P25 or P10, return 1, else NA
stableMask25[] <- ifelse(stableMask25[] == 5 | stableMask25[] == 7, 1, NA)
stableMask10 <- stableMask
# If P10 return 1, else NA
stableMask10[] <- ifelse(stableMask10[] == 7, 1, NA)

# Mask the % time stable to the breaks from the classification
p10r <- mask(timeStableP10, stableMask10)
p25r <- mask(timeStableP25, stableMask25)

# Import the TraCE-21ka mask
traceMask2 <- raster("./plottingObjects/paleo_mask.tif")

# Convert to a dataframe
tR <- stack(p10r, p25r) %>%
  mask(., traceMask2) %>%
  addLayer(., mask(stack(p10r, p25r), traceMask2, inverse = TRUE)) %>%
  as.data.frame(., xy = TRUE) %>%
  tbl_df(.) %>%
  rename("LandP10" = 3, "LandP25" = 4,
         "OceanP10" = 5, "OceanP25" = 6) %>%
  gather(., key = Layer, value = PerTime, na.rm = TRUE,
         LandP10, LandP25, OceanP10, OceanP25) %>%
  mutate(Realm = ifelse(grepl(pattern = "Land", .$Layer), "Land", "Ocean"),
         Perc = ifelse(grepl(pattern = "P10", .$Layer), "P10", "P25"))

# Convert the wrapping variable to a factor and re-order
tR$Perc <- as.factor(tR$Perc)
levels(tR$Perc)
tR$Perc <- factor(
  tR$Perc,
  levels(tR$Perc)[c(2, 1)]
)
levels(tR$Perc)

# Histogram of % time in 25th/10th for land & sea
histStable <- ggplot(tR, aes(x = PerTime*100)) +
  facet_wrap(~Perc, ncol = 1, scales = "free_y",
             labeller = as_labeller(percentLabs),
             strip.position = "right") +
  geom_density(aes(x = PerTime*100, fill = Realm,
                   y = ..density..), trim = TRUE,
               size = 0.5, alpha = 0.5) +
  scale_x_continuous(limits = c(0, 100.1),
                     breaks = seq(0, 100, 10)) +
  scale_fill_manual(values = c("#33A02C", "#1F78B4"),
                    guide = guide_legend(
                      keywidth = 1,
                      keyheight = 1,
                      ncol = 1,
                      direction = "vertical",
                      title.position = "top",
                      title.hjust = 0.5,
                      label.hjust = 0.5,
                      byrow = TRUE,
                      label.position = "right",
                      label.theme = element_text(colour = "black", size = 8)
                    )
                    ) +
  scale_colour_manual(values = c("#33A02C", "#1F78B4")) +
  theme_minimal(base_size = 8) +
  coord_cartesian(expand = FALSE) +
  theme(legend.position = c(0.85, 0.20),
        legend.direction = "vertical",
        legend.title = element_text(colour = "black", size = 8),
        axis.text.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 8),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(size = 8, colour = "black", angle = 0)) +
  labs(x = "Time spent in stable conditions (%)",
       y = "Density")
histStable
```
