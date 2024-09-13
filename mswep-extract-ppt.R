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

sites <- read.csv("/Users/ingridslette/Desktop/NutNet/NutNet-site-coords.csv")

sites <- sites %>% vect(geom = c("longitude", "latitude"), crs = "EPSG:4326")

r_paths <- list.files("/Users/ingridslette/Desktop/mswep",
                      full.names = TRUE) %>% 
  sort()

years <- basename(r_paths) %>% 
  str_extract("^\\d{4}") %>% 
  unique() %>% 
  sort()                       

# extract clim data from mswep --------------------------------------------
#' extract mswep ppt data
#' @param r spatraster (containing daily ppt 'rasters')
#' @param site spatvector containing coordinates of ide sites
#' @return dataframe with site_code, date, and precip columns

extract_mswep <- function(r, site) {
  r_extracted <- terra::extract(r, site, bind = TRUE)
  
  # one row for each site, columns 
  # are site_code 
  df <- as.data.frame(r_extracted)
  # first cold should be site_code
  # additional cols should be precipitation for a given day
  stopifnot(names(df[, 1]) == "site_code",
            str_detect(names(df[, 2:ncol(df)]), "^precipitation"))
  df2 <- df
  names(df2) <- c("site_code", paste0("precip_", time(r)))
  out <- pivot_longer(df2, -"site_code", names_to = "date",
                      values_to = "precip") %>% 
    mutate(date = str_replace(date, "^precip_", ""))
  out
}

extract_mswep <- function(r, site) {
  r_extracted <- terra::extract(r, site, bind = TRUE)
  
  # one row for each site, columns 
  # are site_code 
  df <- as.data.frame(r_extracted)
  # first cold should be site_code
  # additional cols should be precipitation for a given day
  df2 <- df
  names(df2) <- c("site_code", paste0("precip_", time(r)))
  out <- pivot_longer(df2, -"site_code", names_to = "date",
                      values_to = "precip") %>% 
    mutate(date = str_replace(date, "^precip_", ""))
  out
}

dummy <- extract_mswep(rast(r_paths[1]), sites[1, ])

dummy <- dummy[c(), ]

for (yr in years) {
  
  path_names <- names(r_paths) %>% 
    str_subset(paste0("^", yr))
  
  paths_yr <- r_paths[path_names]
  
  r <- rast(paths_yr) # read in data just for the given year
  
  ppt_monthly <- extract_mswep(r, sites)
  write_csv(ppt_monthly, file = "/Users/ingridslette/Desktop/NutNet/mswep-monthly-ppt.csv", append = TRUE)
  message(yr, " complete")
}

ppt_monthly <- extract_mswep(r_paths, sites)




