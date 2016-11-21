library(rfishbase)

spp_for_fb_lookup <- 
c('Caranx melampygus', 'Carangoides orthogrammus', 'Aphareus furca', 
  'Aprion virescens', 'Lutjanus bohar', 'Lutjanus kasmira', 'Cephalopholis argus', 
  'Cephalopholis urodeta', 'Epinephelus hexagonatus', 'Epinephelus maculatus', 
  'Epinephelus spilotoceps', 'Epinephelus tauvina', 'Variola louti', 
  'Paracirrhites arcatus', 'Monotaxis grandoculis', 'Parupeneus insularis', 
  'Caesio teres', 'Pterocaesio tile', 'Chromis vanderbilti', 
  'Pseudanthias bartlettorum', 'Pseudanthias dispar', 'Pseudanthias olivaceus', 
  'Acanthurus nigricans', 'Acanthurus olivaceus', 'Centropyge flavissima', 
  'Chlorurus sordidus', 'Scarus frenatus', 'Scarus rubroviolaceus', 
  'Chaetodon ornatissimus')

fb_fish_data <- 
  species(spp_for_fb_lookup) %>%
    select(Genus, Species, BodyShapeI, Length, LTypeMaxM, LengthFemale, 
           LTypeMaxF, MaxLengthRef, CommonLength, LTypeComM) %>% 
    as.data.frame

readr::write_csv(fb_fish_data, 'fishbase_max_lengths.csv')