library(dplyr)
library(tidyr)

# 1. Carregar els fitxers segons la teva llista de 2009
pcr <- read.csv("qPCR.csv")
hchar <- read.csv("hchar_h.csv")
symp <- read.csv("symp_d.csv")

# 2. Reformatar hchar_h per tenir els dies reals de cada visita (v1, v2, v3)
visites_dies <- hchar %>%
  select(hhID, v1_day, v2_day, v3_day) %>%
  pivot_longer(cols = starts_with("v"), 
               names_to = "visit", 
               values_to = "day_since_index") %>%
  mutate(visit = as.integer(gsub("[^0-9]", "", visit)))

# 3. Unir els resultats de laboratori amb els dies reals
pcr_amb_dies <- pcr %>%
  inner_join(visites_dies, by = c("hhID", "visit"))

# 4. Identificar l'ARI onset per als contactes secundaris (member != 0)
# Criteri Lau et al. (2010): >= 2 símptomes de la llista de 7
ari_onset_2009 <- symp %>%
  filter(member != 0) %>%
  mutate(symptom_count = headache + sthroat + cough + pmuscle + rnose + phlegm + (bodytemp >= 37.8)) %>%
  filter(symptom_count >= 2) %>%
  group_by(hhID, member) %>%
  summarize(day_onset = min(day, na.rm = TRUE), .groups = 'drop')

# 5. Calcular la durada de la infecció (des de l'ARI onset fins a l'últim qPCR positiu)
# Nota: qPCR > 900 copies/mL segons diccionari
durada_2009 <- pcr_amb_dies %>%
  filter(member != 0 & qPCR > 900) %>%
  inner_join(ari_onset_2009, by = c("hhID", "member")) %>%
  mutate(days_since_onset = day_since_index - day_onset) %>%
  group_by(hhID, member) %>%
  summarize(max_duration = max(days_since_onset), .groups = 'drop') %>%
  filter(max_duration >= 0)

# 6. ESTIMACIÓ DE r, d i q_i (mantenim r=1 per comparativa)
r_fix <- 1
max_d_2009 <- max(durada_2009$max_duration, na.rm = TRUE)

resultats_qi_2009 <- data.frame(i = 0:(max_d_2009 - r_fix)) %>%
  rowwise() %>%
  mutate(q_i_2009 = sum(durada_2009$max_duration >= (r_fix + i)) / nrow(durada_2009))

print(resultats_qi_2009)
