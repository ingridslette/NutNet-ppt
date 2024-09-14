library(tidyverse)

# for sites where the growing season spans multiple calendar years:
out$gs_year <- out$year 

str(out)
out$gs_year <- as.numeric(out$gs_year)
out$year <- as.numeric(out$year)
out$month <- as.numeric(out$month)

out <- out %>%
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
    site_code == "sedg.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "sier.us" & month %in% c(11, 12) ~ year + 1,
    site_code == "smith.us" & month %in% c(10, 11, 12) ~ year + 1,
    site_code == "ukul.za" & month %in% c(9, 10, 11, 12) ~ year + 1,
    site_code == "yarra.au" & month %in% c(9, 10, 11, 12) ~ year + 1,
    TRUE ~ year
  ))

# keep only growing season months at each site (there is definitely a better way to do this...)
ppt_monthly_gs_only <- filter(out, 
                         site_code == "arch.us" & month %in% c(5, 6, 7, 8, 9, 10) |
                           site_code =="badlau.de" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                           site_code =="bayr.de" & month %in% c(3, 4, 5, 6, 7, 8, 9) |
                           site_code =="bldr.us" & month %in% c(3, 4, 5, 6, 7, 8) |
                           site_code =="bnbt.us" & month %in% c(5, 6, 7, 8, 9) |
                           site_code =="bnch.us" & month %in% c(4, 5, 6, 7, 8) |
                           site_code =="bogong.au" & month %in% c(10, 11, 12, 1) |
                           site_code =="burrawan.au" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5) |
                           site_code =="burren.ie" & month %in% c(2, 3, 4, 5, 6, 7, 8) |
                           site_code =="cbgb.us" & month %in% c(5, 6, 7, 8, 9, 10) |
                           site_code =="cdcr.us" & month %in% c(4, 5, 6, 7, 8) |
                           site_code =="cdpt.us" & month %in% c(4, 5, 6, 7) |
                           site_code =="cereep.fr" & month %in% c(6, 7, 8, 9, 10, 11) |
                           site_code =="chilcas.ar" & month %in% c(8, 9, 10, 11, 12, 1, 2, 3) |
                           site_code =="comp.pt" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5) |
                           site_code =="cowi.ca" & month %in% c(4, 5, 6, 7) |
                           site_code =="doane.us" & month %in% c(5, 6, 7, 8, 9, 10, 11) |
                           site_code =="elliot.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                           site_code =="ethamc.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5) |
                           site_code =="ethass.au" & month %in% c(5, 6, 7, 8, 9, 10, 11, 12, 1, 2, 3, 4, 5) |
                           site_code =="frue.ch" & month %in% c(4, 5, 6, 7, 8, 9) |
                           site_code =="gilb.za" & month %in% c(9, 10, 11, 12, 1, 2, 3, 4) |
                           site_code =="hall.us" & month %in% c(4, 5, 6, 7, 8, 9) |
                           site_code =="hart.us" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5, 6, 7) |
                           site_code =="hero.uk" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                           site_code =="hopl.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                           site_code =="jena.de" & month %in% c(3, 4, 5, 6, 7, 8, 9, 10) |
                           site_code =="kbs.us" & month %in% c(4, 5, 6, 7, 8, 9) |
                           site_code =="kibber.in" & month %in% c(5, 6, 7, 8) |
                           site_code =="kidman.au" & month %in% c(11, 12, 1, 2, 3, 4) |
                           site_code =="kilp.fi" & month %in% c(6, 7, 8) |
                           site_code =="kiny.au" & month %in% c(5, 6, 7, 8, 9, 10) |
                           site_code =="koffler.ca" & month %in% c(4, 5, 6, 7, 8) |
                           site_code =="konz.us" & month %in% c(5, 6, 7, 8, 9) |
                           site_code =="lagoas.br" & month %in% c(7, 8, 9, 10, 11, 12, 1, 2) |
                           site_code =="lancaster.uk" & month %in% c(3, 4, 5, 6, 7, 8) |
                           site_code =="look.us" & month %in% c(3, 4, 5, 6, 7, 8) |
                           site_code =="marc.ar" & month %in% c(4, 5, 6, 7, 8, 9, 10, 11, 12) |
                           site_code =="mcla.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                           site_code =="msum.us" & month %in% c(4, 5, 6, 7, 8) |
                           site_code =="mtca.au" & month %in% c(8, 9, 10) |
                           site_code =="nilla.au" & month %in% c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12) |
                           site_code =="pape.de" & month %in% c(6, 7, 8, 9) |
                           site_code =="ping.au" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                           site_code =="pinj.au" & month %in% c(1, 2, 3, 4) |
                           site_code =="potrok.ar" & month %in% c(1, 2, 3, 4) |
                           site_code =="rook.uk" & month %in% c(4, 5, 6, 7, 8, 9, 10) |
                           site_code =="saana.fi" & month %in% c(6, 7, 8) |
                           site_code =="sage.us" & month %in% c(4, 5, 6, 7) |
                           site_code =="saline.us" & month %in% c(5, 6, 7, 8, 9) |
                           site_code =="sava.us" & month %in% c(3, 4, 5, 6, 7, 8, 9, 10) |
                           site_code =="sedg.us" & month %in% c(11, 12, 1, 2, 3, 4, 5, 6, 7) |
                           site_code =="sevi.us" & month %in% c(4, 5, 6, 7, 8, 9, 10, 11) |
                           site_code =="sgs.us" & month %in% c(4, 5, 6, 7, 8) |
                           site_code =="shps.us" & month %in% c(4, 5, 6, 7, 8, 9) |
                           site_code =="sier.us" & month %in% c(11, 12, 1, 2, 3, 4) |
                           site_code =="smith.us" & month %in% c(10, 11, 12, 1, 2, 3, 4, 5, 6) |
                           site_code =="spin.us" & month %in% c(3, 4, 5) |
                           site_code =="temple.us" & month %in% c(3, 4, 5, 6, 7, 8, 9, 10) |
                           site_code =="trel.us" & month %in% c(4, 5, 6, 7, 8, 9) |
                           site_code =="ukul.za" & month %in% c(9, 10, 11, 12, 1, 2, 3, 4) |
                           site_code =="valm.ch" & month %in% c(6, 7, 8) |
                           site_code =="veluwe.nl" & month %in% c(3, 4, 5, 6, 7, 8) |
                           site_code =="yarra.au" & month %in% c(9, 10, 11, 12, 1, 2, 3)
)

# sum over all growing season months to get total growing season precip per year per site
ppt_annual_gs_only <- aggregate(precip ~ site_code + gs_year, data = ppt_monthly_gs_only, sum)

write.csv(ppt_annual_gs_only, file = "/Users/ingridslette/Desktop/NutNet/ppt_annual_gs_only.csv")




