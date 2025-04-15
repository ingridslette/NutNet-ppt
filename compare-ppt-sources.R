# Ingrid Slette
# started 14 September 2024

library(tidyverse)

mswep <- read.csv("/Users/ingridslette/Desktop/NutNet/ppt_annual_gs_only.csv")

cru <- read.csv("/Users/ingridslette/Dropbox/NutNet data/climate/CRU/CRU-annual-gs-pre-1901-2023.csv")

stations <- read.csv("/Users/ingridslette/Desktop/Bharath et al papers/Data_S1/Precip_data_20210310.csv")

mswep <- mswep %>%
  rename(year = gs_year)

mswep <- mswep %>%
  rename(mswep_ppt = precip)

cru <- cru %>%
  rename(cru_ppt = gs_ppt)

stations <- stations %>%
  rename(station_ppt = ppt_gs)

unique(cru$site_code)
unique(mswep$site_code)

mswep_cru <- left_join(mswep, cru, by = c("site_code", "year"))

unique(mswep_cru$site_code)

ggplot(mswep_cru, aes(x=mswep_ppt, y=cru_ppt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("mswep") + ylab("cru") +
  theme_bw()

ggplot(mswep_cru, aes(x=mswep_ppt, y=cru_ppt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("mswep") + ylab("cru") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

mswep_cru_station <- left_join(mswep_cru, stations, by = c("site_code", "year"))

ggplot(mswep_cru_station, aes(x=mswep_ppt, y=station_ppt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("mswep") + ylab("stations") +
  theme_bw()

ggplot(mswep_cru_station, aes(x=cru_ppt, y=station_ppt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("cru") + ylab("stations") +
  theme_bw()

ggplot(mswep_cru_station, aes(x=mswep_ppt, y=station_ppt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("mswep") + ylab("stations") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()

ggplot(mswep_cru_station, aes(x=cru_ppt, y=station_ppt)) +
  geom_point() + geom_smooth(method = lm) +
  xlab("cru") + ylab("stations") +
  facet_wrap(vars(site_code), scales = "free") +
  theme_bw()
