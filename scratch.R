# https://bioconnector.github.io/workshops/r-survival.html

library(RTCGA.clinical)
?clinical

dim(LUSC.clinical)
head(names(LUSC.clinical))

# Create the clinical data
clin <- survivalTCGA(LUSC.clinical, 
                     extract.cols="admin.disease_code")
# Show the first few lines
head(clin)

# subset the data to get only ids in jacki factor analysis
clin <- clin[clin$bcr_patient_barcode %in% LUSC_jacki_factor$Tumor_Sample_ID,]

clin <- clin %>% mutate(cdr3_score = ifelse(clin$bcr_patient_barcode %in% top50$Tumor_Sample_ID, "top", "bottom"))

# How many samples of each type?
# table(clin$admin.disease_code)
table(clin$cdr3_score)

# Tabulate by outcome
# xtabs(~admin.disease_code+patient.vital_status, data=clin) %>% addmargins()
xtabs(~cdr3_score+patient.vital_status, data=clin) %>% addmargins()

# sfit <- survfit(Surv(times, patient.vital_status)~admin.disease_code, data=clin)
sfit <- survfit(Surv(times, patient.vital_status)~cdr3_score, data=clin)

summary(sfit, times=seq(0,365*5,365))

ggsurvplot(sfit, conf.int=TRUE, pval=TRUE)
