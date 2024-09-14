# Ingrid Slette
# started 13 September 2024

# getting MSWEP (http://www.gloh2o.org/mswep/) precipitation data for NutNet sites
# request access to MSWEP data, which is shared via Google Drive
# install and set up rclone to acess MSWEP Google Drive

# I copied all years of monthly MSWEP data into a file on my Desktop using the R terminal using the following code
# MacBook-Pro-169:~ ingridslette$ rclone  copy --drive-shared-with-me gdrive:MSWEP_V280/Past/Monthly /Users/ingridslette/Desktop/mswep
# MacBook-Pro-169:~ ingridslette$ rclone  copy --drive-shared-with-me gdrive:MSWEP_V280/NRT/Monthly /Users/ingridslette/Desktop/mswep

library(terra)
library(lubridate)
library(tidyverse)

# read in list of all NutNet site codes and coordinates 
sites <- read.csv("/Users/ingridslette/Desktop/NutNet/NutNet-site-coords.csv")

# make that a SpatVector
site <- sites %>% vect(geom = c("longitude", "latitude"), crs = "EPSG:4326")

# list all of the monthly mswep precip data files
r_paths <- list.files("/Users/ingridslette/Desktop/mswep",
                      full.names = TRUE) %>% 
  sort()

# make that a SpatRaster
r <- rast(r_paths)

# extract monthly precip data for each site
ppt_monthly <- terra::extract(r, site, bind = TRUE)

df <- as.data.frame(ppt_monthly)

names(df) <- c("site_code", paste0("precip_", time(r)))

out <- pivot_longer(df, -"site_code", names_to = "date",
                    values_to = "precip") %>% 
  mutate(date = str_replace(date, "^precip_", ""))

out$date <- substr(out$date, 1, 7)





