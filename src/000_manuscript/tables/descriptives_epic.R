rm(list=ls())
set.seed(821)

# environment ====
library(dplyr)
library(openxlsx)
source("scripts/000_functions_epic-somalogic-data.R")

# phenofile ====
data_phenofile <- fread("analysis/001_phenofile/cancer_NormalizedSoma_PlateCorrAccountforDisease_Log10_SMP.txt", header = T, sep = "\t")
phenofile_cancer <- table(data_phenofile$cncr_mal_anyc)
phenofile_cancer_overall <- table(data_phenofile$cncr_mal_clrt)
phenofile_cancer_colon <- table(data_phenofile$cncr_mal_clrt_colon)
phenofile_cancer_rectum <- table(data_phenofile$cncr_mal_clrt_rectum)

# processed phenofile ====
data_phenofile_processed <- epic_somalogic_format_data(data_phenofile)
data_phenofile_processed <- epic_somalogic_create_subtypes(data_phenofile_processed)
data_phenofile_processed <- data_phenofile_processed[1:3]
data_phenofile_processed <- lapply(data_phenofile_processed, epic_somalogic_create_sex)
data_phenofile_processed <- epic_somalogic_create_analysis(data_phenofile_processed)

# descriptive table ====
result_list <- list()
for (i in names(data_phenofile_processed)){
  
  for (j in names(data_phenofile_processed[[i]])){
    result_list_j <- list()
    
    for (k in names(data_phenofile_processed[[i]][[j]])){
      
      pheno_cancer <- i
      pheno_sex <- j
      
      # age
      pheno_age_median <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(age) %>%
        median(na.rm = TRUE)
      pheno_age_sd <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(age) %>%
        sd(na.rm = TRUE)
      
      pheno_age_median_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(age) %>%
        median(na.rm = TRUE)
      pheno_age_sd_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(age) %>%
        sd(na.rm = TRUE)
      
      pheno_age_median_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(age) %>%
        median(na.rm = TRUE)
      pheno_age_sd_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(age) %>%
        sd(na.rm = TRUE)
      
      # followup
      pheno_followup_median <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(followup) %>%
        median(na.rm = TRUE)
      pheno_followup_sd <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(followup) %>%
        sd(na.rm = TRUE)
      
      pheno_followup_median_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(followup) %>%
        median(na.rm = TRUE)
      pheno_followup_sd_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(followup) %>%
        sd(na.rm = TRUE)
      
      pheno_followup_median_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(followup) %>%
        median(na.rm = TRUE)
      pheno_followup_sd_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(followup) %>%
        sd(na.rm = TRUE)
      
      # BMI
      pheno_bmi_median <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(bmi_c) %>%
        median(na.rm = TRUE)
      pheno_bmi_sd <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(bmi_c) %>%
        sd(na.rm = TRUE)
      
      pheno_bmi_median_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(bmi_c) %>%
        median(na.rm = TRUE)
      pheno_bmi_sd_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(bmi_c) %>%
        sd(na.rm = TRUE)
      
      pheno_bmi_median_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(bmi_c) %>%
        median(na.rm = TRUE)
      pheno_bmi_sd_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(bmi_c) %>%
        sd(na.rm = TRUE)
      
      # alcohol
      pheno_alcohol_median <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(alc_re) %>%
        median(na.rm = TRUE)
      pheno_alcohol_sd <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(alc_re) %>%
        sd(na.rm = TRUE)
      
      pheno_alcohol_median_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(alc_re) %>%
        median(na.rm = TRUE)
      pheno_alcohol_sd_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(alc_re) %>%
        sd(na.rm = TRUE)
      
      pheno_alcohol_median_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(alc_re) %>%
        median(na.rm = TRUE)
      pheno_alcohol_sd_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(alc_re) %>%
        sd(na.rm = TRUE)
      
      # smoke_stat
      pheno_smoke <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(smoke_stat) %>%
        table() %>%
        as.data.frame() %>%
        rename(smoke_stat = '.', pheno_smoke = Freq)
      pheno_smoke_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(smoke_stat) %>%
        table() %>%
        as.data.frame() %>%
        rename(smoke_stat = '.', pheno_smoke_control = Freq)
      pheno_smoke_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(smoke_stat) %>%
        table() %>%
        as.data.frame() %>%
        rename(smoke_stat = '.', pheno_smoke_control = Freq)
      
      # pa_index
      pheno_pa <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(pa_index) %>%
        table() %>%
        as.data.frame() %>%
        rename(pa_index = '.', pheno_pa = Freq)
      pheno_pa_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(pa_index) %>%
        table() %>%
        as.data.frame() %>%
        rename(pa_index = '.', pheno_pa_control = Freq)
      pheno_pa_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(pa_index) %>%
        table() %>%
        as.data.frame() %>%
        rename(pa_index = '.', pheno_pa_case = Freq)
      
      # l_school
      pheno_education <- data_phenofile_processed[[i]][[j]][[k]] %>%
        pull(l_school) %>%
        table() %>%
        as.data.frame() %>%
        rename(l_school = '.', pheno_education = Freq)
      pheno_education_control <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 0) %>%
        pull(l_school) %>%
        table() %>%
        as.data.frame() %>%
        rename(l_school = '.', pheno_education_control = Freq)
      pheno_education_case <- data_phenofile_processed[[i]][[j]][[k]] %>%
        filter(indevent == 1) %>%
        pull(l_school) %>%
        table() %>%
        as.data.frame() %>%
        rename(l_school = '.', pheno_education_case = Freq)
      
      # df_combined
      df_combined <- data.frame(
        cancer = pheno_cancer,
        sex = pheno_sex,
        case_status = "all",
        n = sum(table(data_phenofile_processed[[i]][[j]][[k]]$indevent)),
        age_median = round(pheno_age_median,2),
        age_sd = round(pheno_age_sd,2),
        followup_median = round(pheno_followup_median,2),
        followup_sd = round(pheno_followup_sd,2),
        bmi_median = round(pheno_bmi_median,2),
        bmi_sd = round(pheno_bmi_sd,2),
        alcohol_median = round(pheno_alcohol_median,2),
        alcohol_sd = round(pheno_alcohol_sd,2),
        
        pheno_smoke[1, 2],
        pheno_smoke[2, 2],
        pheno_smoke[3, 2],
        
        pheno_pa[1, 2],
        pheno_pa[2, 2],
        pheno_pa[3, 2],
        pheno_pa[4, 2],
        pheno_pa[5, 2],
        
        pheno_education[1, 2],
        pheno_education[2, 2],
        pheno_education[3, 2],
        pheno_education[4, 2],
        pheno_education[5, 2],
        pheno_education[6, 2]
        
      )
      
      names(df_combined)[13:15] <- paste0("smoke_", pheno_smoke[1:3, 1])
      names(df_combined)[16:20] <- paste0("pa_", pheno_pa[1:5, 1])
      names(df_combined)[21:26] <- paste0("education_", pheno_education[1:6, 1])
      
      # df_control
      df_control <- data.frame(
        cancer = pheno_cancer,
        sex = pheno_sex,
        case_status = "non-case",
        n = table(data_phenofile_processed[[i]][[j]][[k]]$indevent)[1],
        age_median = round(pheno_age_median_control,2),
        age_sd = round(pheno_age_sd_control,2),
        followup_median = round(pheno_followup_median_control,2),
        followup_sd = round(pheno_followup_sd_control,2),
        bmi_median = round(pheno_bmi_median_control,2),
        bmi_sd = round(pheno_bmi_sd_control,2),
        alcohol_median = round(pheno_alcohol_median_control,2),
        alcohol_sd = round(pheno_alcohol_sd_control,2),
        
        pheno_smoke_control[1, 2],
        pheno_smoke_control[2, 2],
        pheno_smoke_control[3, 2],
        
        pheno_pa_control[1, 2],
        pheno_pa_control[2, 2],
        pheno_pa_control[3, 2],
        pheno_pa_control[4, 2],
        pheno_pa_control[5, 2],
        
        pheno_education_control[1, 2],
        pheno_education_control[2, 2],
        pheno_education_control[3, 2],
        pheno_education_control[4, 2],
        pheno_education_control[5, 2],
        pheno_education_control[6, 2]
        
      )
      
      names(df_control)[13:15] <- paste0("smoke_", pheno_smoke_control[1:3, 1])
      names(df_control)[16:20] <- paste0("pa_", pheno_pa_control[1:5, 1])
      names(df_control)[21:26] <- paste0("education_", pheno_education_control[1:6, 1])
      
      # df_case
      df_case <- data.frame(
        cancer = pheno_cancer,
        sex = pheno_sex,
        case_status = "case",
        n = table(data_phenofile_processed[[i]][[j]][[k]]$indevent)[2],
        age_median = round(pheno_age_median_case,2),
        age_sd = round(pheno_age_sd_case,2),
        followup_median = round(pheno_followup_median_case,2),
        followup_sd = round(pheno_followup_sd_case,2),
        bmi_median = round(pheno_bmi_median_case,2),
        bmi_sd = round(pheno_bmi_sd_case,2),
        alcohol_median = round(pheno_alcohol_median_case,2),
        alcohol_sd = round(pheno_alcohol_sd_case,2),
        
        pheno_smoke_case[1, 2],
        pheno_smoke_case[2, 2],
        pheno_smoke_case[3, 2],
        
        pheno_pa_case[1, 2],
        pheno_pa_case[2, 2],
        pheno_pa_case[3, 2],
        pheno_pa_case[4, 2],
        pheno_pa_case[5, 2],
        
        pheno_education_case[1, 2],
        pheno_education_case[2, 2],
        pheno_education_case[3, 2],
        pheno_education_case[4, 2],
        pheno_education_case[5, 2],
        pheno_education_case[6, 2]
        
      )
      
      names(df_case)[13:15] <- paste0("smoke_", pheno_smoke_case[1:3, 1])
      names(df_case)[16:20] <- paste0("pa_", pheno_pa_case[1:5, 1])
      names(df_case)[21:26] <- paste0("education_", pheno_education_case[1:6, 1])
      
      # combine table
      table_temp <- as.data.frame(t(rbind(df_combined, df_control, df_case)))
      result_list_j[[k]] <- table_temp
    }
    result_list[[paste(i, j, sep = "_")]] <- result_list_j
    
  }
}

table_0year <- c(result_list[[1]]["analysis"],
                 result_list[[2]]["analysis"],
                 result_list[[3]]["analysis"],
                 result_list[[4]]["analysis"],
                 result_list[[5]]["analysis"],
                 result_list[[6]]["analysis"],
                 result_list[[7]]["analysis"],
                 result_list[[8]]["analysis"],
                 result_list[[9]]["analysis"])
table_0year <- bind_cols(table_0year)

table_2year <- c(result_list[[1]]["followup_excl2year"],
                 result_list[[2]]["followup_excl2year"],
                 result_list[[3]]["followup_excl2year"],
                 result_list[[4]]["followup_excl2year"],
                 result_list[[5]]["followup_excl2year"],
                 result_list[[6]]["followup_excl2year"],
                 result_list[[7]]["followup_excl2year"],
                 result_list[[8]]["followup_excl2year"],
                 result_list[[9]]["followup_excl2year"])
table_2year <- bind_cols(table_2year)

table_5year <- c(result_list[[1]]["followup_excl5year"],
                 result_list[[2]]["followup_excl5year"],
                 result_list[[3]]["followup_excl5year"],
                 result_list[[4]]["followup_excl5year"],
                 result_list[[5]]["followup_excl5year"],
                 result_list[[6]]["followup_excl5year"],
                 result_list[[7]]["followup_excl5year"],
                 result_list[[8]]["followup_excl5year"],
                 result_list[[9]]["followup_excl5year"])
table_5year <- bind_cols(table_5year)


# save ====
write.table(table_0year, "analysis/tables/manuscript/descriptives-0year.txt", 
            row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(table_2year, "analysis/tables/manuscript/descriptives-2year.txt", 
            row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")
write.table(table_5year, "analysis/tables/manuscript/descriptives-5year.txt", 
            row.names = TRUE, col.names = FALSE, quote = FALSE, sep = "\t")

wb <- createWorkbook()
addWorksheet(wb = wb, sheetName = "all")
writeData(wb = wb, sheet = "all", table_0year, colNames = F, rowNames = T)

addWorksheet(wb = wb, sheetName = "2year-fup-exclusion")
writeData(wb = wb, sheet = "2year-fup-exclusion", table_2year, colNames = F, rowNames = T)

addWorksheet(wb = wb, sheetName = "5year-fup-exclusion")
writeData(wb = wb, sheet = "5year-fup-exclusion", table_5year, colNames = F, rowNames = T)

saveWorkbook(wb, "analysis/tables/manuscript/descriptives.xlsx", overwrite = TRUE)
