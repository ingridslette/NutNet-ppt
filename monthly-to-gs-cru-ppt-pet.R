library(tidyverse)

local <- read.csv("/Users/ingridslette/Desktop/Weather_monthly_20220615.csv")

unique(local$site_code)

local <- local %>% 
  filter(site_code %in% c("bayr.de", "comp.pt", "cowi.ca", "ethamc.au", "ethass.au", "frue.ch",
                       "hero.uk", "lancaster.uk", "look.us", "marc.ar", "mcla.us", "sevi.us",
                       "ukul.za", "valm.ch"))

local$local <- local$pre 

local <- local %>% 
  mutate(local = 1)


ghcn <- read.csv("/Users/ingridslette/Desktop/NutNet/nutnet-ghcn-ppt-data.csv")

unique(ghcn$STATION)

ghcn <- ghcn %>% 
  mutate(site_code = case_when(
    STATION == "USS0020K04S" ~ "sage.us",
    STATION == "ASN00010044" ~ "mtca.au",
    STATION == "USC00350652" ~ "bnch.us",
    STATION == "USC00130203" ~ "cbgb.us",
    STATION == "USC00252020" ~ "doane.us",
    STATION == "USC00212881" ~ "cdcr.us",
    STATION == "USC00353692" ~ "hart.us",
    STATION == "GME00125074" ~ "badlau.de",
    STATION == "ASN00083084" ~ "bogong.au",
    STATION == "USC00157049" ~ "hall.us",
    STATION == "USW00003131" ~ "elliot.us",
    STATION == "USC00080236" ~ "arch.us",
    STATION == "USW00023275" ~ "hopl.us",
    STATION == "GM000004204" ~ "jena.de",
    STATION == "USC00153194" ~ "spin.us",
    STATION == "USC00145628" ~ "saline.us",
    STATION == "ASN00082165" ~ "nilla.au",
    STATION == "USW00094074" ~ "sgs.us",
    STATION == "USC00203504" ~ "kbs.us",
    STATION == "ASN00080002" ~ "kiny.au",
    STATION == "ASN00014847" ~ "kidman.au",
    STATION == "USC00418646" ~ "temple.us",
    STATION == "USC00247894" ~ "msla.us",
    STATION == "ASN00010626" ~ "ping.au",
    STATION == "USC00111743" ~ "trel.us",
    STATION == "ASN00067105" ~ "yarra.au",
    STATION == "NLE00152497" ~ "veluwe.nl",
    STATION == "CAW00064757" ~ "koffler.ca"
  ))

unique(ghcn$site_code)

ghcn_monthly <- ghcn %>% 
  group_by(site_code, month, year) %>% 
  summarise(PRCP_monthly = sum(PRCP))

ghcn_monthly$ghcn <- ghcn_monthly$PRCP_monthly

ghcn_monthly <- ghcn_monthly %>% 
  mutate(ghcn = 1)


cru <- read.csv('/Users/ingridslette/Desktop/CRU/CRU-monthly-pre-pet-1901-2024.csv')

unique(cru$site_code)

cru$cru <- cru$pre_mm.month

cru <- cru %>% 
  mutate(cru = 1)

cru <- cru %>%
  rename(month_abb = month)

cru$month <- match(cru$month_abb, month.abb)

cru <- cru %>% 
  dplyr::select(site_code, month, year, plotdate, pre_mm.month, pet_mm.month, pre_pet, cru)

cru_ghcn <- left_join(cru, ghcn_monthly, by = c("site_code", "month", "year"))

ppt_data <- left_join(cru_ghcn, local, by = c("site_code", "month", "year"))

ppt_data <- ppt_data %>% 
  mutate(ppt = case_when(
    !is.na(pre) ~ pre,
    is.na(pre) & !is.na(PRCP_monthly) ~ PRCP_monthly,
    is.na(pre) & is.na(PRCP_monthly) ~ pre_mm.month
  ))

# Define a "not in" operator if not already defined
`%notin%` <- Negate(`%in%`)

# Define site groups
local_sites <- c("bayr.de", "comp.pt", "cowi.ca", "ethamc.au", "ethass.au", "frue.ch",
                 "hero.uk", "lancaster.uk", "look.us", "marc.ar", "mcla.us", "sevi.us",
                 "ukul.za", "valm.ch")

ghcn_sites <- c("cbgb.us", "doane.us", "cdcr.us", "hart.us", "badlau.de", "bogong.au",
                "hall.us", "elliot.us", "arch.us", "hopl.us", "jena.de", "spin.us",
                "saline.us", "nilla.au", "sgs.us", "kbs.us", "kiny.au", "kidman.au",
                "temple.us", "msla.us", "ping.au", "trel.us", "yarra.au", "veluwe.nl",
                "koffler.ca", "sage.us", "mtca.au", "bnch.us")

# All listed sites
non_cru_sites <- c(local_sites, ghcn_sites)

# Add the 'keep' column
ppt_data <- ppt_data %>%
  mutate(keep = case_when(
    site_code %in% local_sites & !is.na(local) ~ 1,
    site_code %in% ghcn_sites & !is.na(ghcn) ~ 1,
    site_code %notin% non_cru_sites & !is.na(cru) ~ 1,
    TRUE ~ 0  # default for everything else
  ))


# for sites where the growing season spans multiple calendar years:
ppt_data$gs_year <- ppt_data$year 

str(ppt_data)
ppt_data$gs_year <- as.numeric(ppt_data$gs_year)
ppt_data$year <- as.numeric(ppt_data$year)
ppt_data$month <- as.numeric(ppt_data$month)

unique(ppt_data$site_code)

ppt_data <- ppt_data %>%
  mutate(gs_year = case_when(
    site_code == "bogong.au" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "burrawan.au" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "chilcas.ar" & month %in% c(08, 09, 10, 11, 12) ~ year + 1,
    site_code == "comp.pt" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "elliot.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "ethamc.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "ethass.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "gilb.za" & month %in% c(9, 10, 11, 12) ~ year + 1,
    site_code == "hart.us" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "hopl.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "kidman.au" & month %in% c(11, 12) ~ year + 1,
    site_code == "lagoas.br" & month %in% c(7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "mcla.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "nilla.au" & month %in% c(2, 3 ,4, 5, 6, 7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "sedg.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "sier.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "smith.us" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "spin.us" & month %in% c(6, 7, 8, 9, 10, 11, 12) ~ year + 1,
    site_code == "ukul.za" & month %in% c(9, 10, 11, 12) ~ year + 1,
    site_code == "yarra.au" & month %in% c(9, 10, 11, 12) ~ year + 1,
    TRUE ~ year
  ))

# keep only growing season months at each site (there is definitely a better way to do this...)
ppt_gs_only <- filter(ppt_data, site_code =="arch.us" & month %in% c(5, 6, 7, 8, 9, 10) |
                        site_code =="badlau.de" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                        site_code =="bayr.de" & month %in% c(3, 4, 5, 6, 7, 8, 9) |
                        site_code =="bldr.us" & month %in% c(3, 4, 5, 6, 7, 8) |
                        site_code =="bnbt.us" & month %in% c(5, 6, 7, 8, 9) |
                        site_code =="bnch.us" & month %in% c(4, 5, 6, 7, 8) |
                        site_code =="bogong.au" & month %in% c(10, 11, 12, 1) |
                        site_code =="burrawan.au" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5) |
                        site_code =="cbgb.us" & month %in% c(5, 6, 7, 8, 9, 10) |
                        site_code =="cdcr.us" & month %in% c(4, 5, 6, 7, 8) |
                        site_code =="cedr.us" & month %in% c(4, 5, 6, 7, 8) |
                        site_code =="cdpt.us" & month %in% c(4, 5, 6, 7) |
                        site_code =="chilcas.ar" & month %in% c(8, 9, 10, 11, 12, 1, 2, 3) |
                        site_code =="comp.pt" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5) |
                        site_code =="cowi.ca" & month %in% c(4, 5, 6, 7) |
                        site_code =="doane.us" & month %in% c(5, 6, 7, 8, 9, 10, 11) |
                        site_code =="elliot.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                        site_code =="ethamc.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4) |
                        site_code =="ethass.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4) |
                        site_code =="frue.ch" & month %in% c(4, 5, 6, 7, 8, 9) |
                        site_code =="hall.us" & month %in% c(4, 5, 6, 7, 8, 9) |
                        site_code =="hero.uk" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                        site_code =="hopl.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                        site_code =="jena.de" & month %in% c(3, 4, 5, 6, 7, 8, 9, 10) |
                        site_code =="kbs.us" & month %in% c(4, 5, 6, 7, 8, 9) |
                        site_code =="kidman.au" & month %in% c(11, 12, 1, 2, 3, 4) |
                        site_code =="kilp.fi" & month %in% c(6, 7, 8) |
                        site_code =="kiny.au" & month %in% c(5, 6, 7, 8, 9, 10) |
                        site_code =="koffler.ca" & month %in% c(4, 5, 6, 7, 8) |
                        site_code =="konz.us" & month %in% c(5, 6, 7, 8, 9) |
                        site_code =="lagoas.br" & month %in% c(7, 8, 9, 10, 11, 12, 1, 2) |
                        site_code =="lancaster.uk" & month %in% c(3, 4, 5, 6, 7, 8) |
                        site_code =="look.us" & month %in% c(3, 4, 5, 6, 7, 8) | 
                        site_code =="lubb.us" & month %in% c(3, 4, 5, 6, 7, 8, 9, 10) |
                        site_code =="marc.ar" & month %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) |
                        site_code =="mcla.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                        site_code =="msla_3.us" & month %in% c(4, 5, 6, 7) |
                        site_code =="msla.us" & month %in% c(4, 5, 6, 7) |
                        site_code =="msla_2.us" & month %in% c(4, 5, 6, 7) |
                        site_code =="mtca.au" & month %in% c(8, 9, 10) |
                        site_code =="nilla.au" & month %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) |
                        site_code =="pape.de" & month %in% c(6, 7, 8, 9) |
                        site_code =="ping.au" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                        site_code =="potrok.ar" & month %in% c(1, 2, 3, 4) |
                        site_code =="saana.fi" & month %in% c(6, 7, 8) |
                        site_code =="sage.us" & month %in% c(4, 5, 6, 7) |
                        site_code =="saline.us" & month %in% c(5, 6, 7, 8, 9) |
                        site_code =="sedg.us" & month %in% c(11, 12, 1, 2, 3, 4, 5, 6, 7) |
                        site_code =="sevi.us" & month %in% c(4, 5, 6, 7, 8, 9, 10, 11) |
                        site_code =="sgs.us" & month %in% c(4, 5, 6, 7, 8) |
                        site_code =="sier.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                        site_code =="smith.us" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5, 6) |
                        site_code =="spin.us" & month %in% c(1, 2, 3, 4, 5) |
                        site_code =="temple.us" & month %in% c(3, 4, 5, 6, 7, 8, 9, 10) |
                        site_code =="trel.us" & month %in% c(1, 2, 3, 4, 5, 6, 7, 8, 9) |
                        site_code =="ukul.za" & month %in% c(9, 10, 11, 12, 1, 2, 3, 4) |
                        site_code =="valm.ch" & month %in% c(6, 7, 8) |
                        site_code =="veluwe.nl" & month %in% c(3, 4, 5, 6, 7, 8) |
                        site_code =="yarra.au" & month %in% c(9, 10, 11, 12, 1, 2, 3)
                      
)

unique(ppt_gs_only$site_code)

# sum to get total growing season precip per year per site

ppt_annual_gs_only <- ppt_gs_only %>% 
  group_by(site_code, gs_year) %>% 
  summarise(ags_ppt = sum(ppt), ags_pet = sum(pet_mm.month))

ppt_annual_gs_only <- ppt_annual_gs_only %>%
  rename(year = gs_year, ppt = ags_ppt, pet = ags_pet)

ppt_annual_gs_only <- ppt_annual_gs_only %>% 
  mutate(ppt_pet = ppt-pet)

unique(ppt_annual_gs_only$site_code)

write.csv(ppt_annual_gs_only, file = "/Users/ingridslette/Desktop/NutNet/ppt_annual_gs_only_2025-10-07.csv", row.names=FALSE)




