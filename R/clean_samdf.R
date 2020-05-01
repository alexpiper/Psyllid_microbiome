
# Cleaning sample data ----------------------------------------------------

samdf <- read_csv("sample_data/Sample_info.csv")  %>%
  dplyr::rename_all(funs( stringr::str_replace_all(., '\\ ', '_')) ) %>%
  mutate(psyllid_spp = psyllid_spp %>% str_remove("\\.") %>%
           str_replace_all(" ", "_") %>%
           str_replace_all("Casuarinicola_sp", "Triozid_sp")) %>%
  mutate(hostplant_spp = hostplant_spp %>% str_remove("\\.") %>%
           str_replace_all(" ", "_") %>%
           str_to_sentence() %>%
           str_replace_all("Kamahi", "Weinmannia_racemosa") %>%
           str_replace_all("Oleaeria_odorata", "Olearia_odorata"))

# Flag replicated samples
samdf <- samdf %>%
  dplyr::mutate(replicated = Sample_Name %in% (samdf %>% group_by(Sample_Name) %>% 
                                                 add_count() %>%
                                                 filter(n >1) %>% 
                                                 pull(Sample_Name))) %>%
  separate(Collection, into= c("Lat", "Long"), sep=" ", remove = TRUE) %>%
  separate(Lat, paste("lat",c("d","m","s"), sep="_")) %>%
  separate(Long, paste("long",c("d","m","s"), sep="_" )) %>%
  mutate_at(vars(starts_with("lat")), .funs=as.numeric) %>%
  mutate_at(vars(starts_with("long")), .funs=as.numeric) %>%
  mutate(lat=-(lat_d + lat_m/60 + lat_s/60^2),
         long=long_d + long_m/60 + long_s/60^2) %>%
  dplyr::select(c("SampleID", "Sample_Name", "seqrun", "psyllid_spp", "psyllid_genus", "psyllid_family", "hostplant_spp","Collection_Date","replicated", "lat", "long")) %>%
  as.data.frame(stringsAsFactors=FALSE) %>%
  magrittr::set_rownames(.$SampleID) #Collection

write_csv(samdf, "sample_data/Sample_info2.csv")
