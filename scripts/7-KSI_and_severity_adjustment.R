library(tidyverse)
library(lubridate)
library(here)

# ------------------------------------------------------------
# Load matched RTI data
# -------------------------------
injuries <- read_rds(here("data", "processed", "injuries_matched_OA.rds"))
# -------------------------------
#  severity adjustment 
# --------------------------------

# STATS19 coding:
# 1 = Fatal
# 2 = Serious
# 3 = Slight
# KSI = Fatal + Serious

injuries <- injuries %>%
  mutate(
    # ----------------------------
    # Unadjusted indicators
    # ----------------------------
    KSI_unadj = as.numeric(casualty_severity %in% c(1, 2)),
    Slight_unadj = as.numeric(casualty_severity == 3),
    
    # ----------------------------
    # Adjusted (proper KSI)
    # ----------------------------
    KSI_adj = case_when(
      casualty_severity == 1 ~ 1,  # Fatal
      casualty_severity == 2 ~ casualty_adjusted_severity_serious,
      casualty_severity == 3 ~ 0
    ),
    
    Slight_adj = case_when(
      casualty_severity == 3 ~ casualty_adjusted_severity_slight,
      TRUE ~ 0
    )
  )
# ------------------------------------------------------------
#checks
# ------------------------------------------------------------

stopifnot(
  all(injuries$KSI_adj >= 0 & injuries$KSI_adj <= 1),
  all(injuries$Slight_adj >= 0 & injuries$Slight_adj <= 1)
)

# :  comparison summary
national_totals <- injuries %>%
  summarise(
    KSI_unadj_total = sum(KSI_unadj, na.rm = TRUE),
    KSI_adj_total   = sum(KSI_adj, na.rm = TRUE),
    Slight_unadj_total = sum(Slight_unadj, na.rm = TRUE),
    Slight_adj_total   = sum(Slight_adj, na.rm = TRUE)
  )

print(national_totals)


injuries %>%
  filter(casualty_severity == 1) %>%
  summarise(mean_KSI_adj = mean(KSI_adj))

# # # #Adjusted KSI is very close to unadjusted KSI (difference = 20). 
##Slight decreases a lot after adjustment - -  ie reweighted toward Serious.



write_rds(
  injuries,
  here("data", "processed", "injuries_matched_final.rds")
)


