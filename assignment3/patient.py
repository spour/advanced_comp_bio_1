# q1: create a table that captures the VAF for every variant observed in that control samples, 
# e.g. the columns are VAFs and the rwos are the control samples. if the variant was not observed
# in a particular control sample, set the VAF to 0

controls = "/path/to/controls/files"
patients = "/path/to/patient/input/files"
mutations = "/path/to/mutations.csv"
one_table <- function(path) {
  list_files = list.files(path, full.names = T, pattern = "\\.csv$")
  cols = c("chr", "pos", "nucleotide")
  cols_rm = c(cols,'count', 'coverage_at_position' )
  table = read.table(list_files[1], header=T, sep =",")
  table$x <- apply( table[ , cols ] , 1 , paste , collapse = "_" )
  table <- table[ , !( names( table ) %in% cols_rm ) ]
  names(table)[names(table) == 'VAF'] <- tools::file_path_sans_ext(basename(list_files[1]))
  for (file in list_files[-1]){
    table2 = read.table(file, header=T, sep =",")
    table2$x <- apply( table2[ , cols ] , 1 , paste , collapse = "_" )
    table2 <- table2[ , !( names( table2 ) %in% cols_rm )]
    names(table2)[names(table2) == 'VAF'] <- tools::file_path_sans_ext(basename(file))
    table <- merge(table, table2, all = TRUE, by = 'x')
  }
  table$x <- gsub(" ", "", table$x)
  table[is.na(table)] <- 0
  return(table)
}

table = one_table(controls) # call the fucntion
# write.table(table, row.names = F, quote = F, sep = ",", file="~/Desktop/acb1/Sara Pour A3/Q1.csv")

# q2: for each variant that was observed at least once in the control sample, fit an exponential 
# distribution to the set of teh VAF values for that variant. the fitted rate parameter(lambda) for
# each variant captured by the model should be put in a new rate column. 
library(tidyr)
library(purrr)
library(broom)
library(reshape2)
require(data.table)
rate_exp <- function(table){
  table_t = data.frame(t(table))
  names(table_t) <- as.matrix(table_t[1, ])
  table_t <- table_t[-1, ]
  table_t[] <- lapply(table_t, function(x) type.convert(as.numeric(as.character(x))))
  df <- gather(table_t)
  df$value <- as.numeric(as.character(df$value))
  df2 <- df %>% group_by(key) %>% do(broom::tidy(MASS::fitdistr(.$value, "exponential")))
  return(df2)
  # write.table(format(table_t, digits=20), quote = F, sep = "\t", file = "~/Desktop/acb1/table_t.tsv")
  # table_t$`chr1-43814991-T` <- round(table_t$`chr1-43814991-T`,5)
  # fit1 <- MASS::fitdistr(table_t$`chr1-43814991-T`, "exponential")
  # ks.test(table_t$`chr1-43814991-T`, "pexp", fit1$estimate)
  # hist(table_t$`chr1-43814991-T`, freq = FALSE, breaks = 100, xlim = c(0, quantile(table_t$`chr1-43814991-T`, 0.99)))
  # curve(dexp(x, rate = fit1$estimate), from = 0, col = "red", add = TRUE)
}
rates <- rate_exp(table)
# rates$estimate = c(1170.2842824371633,
#                     1357.4678731186439,
#                     12242.431545605452,
#                     2123.0490902060346,
#                     5377.054569462006,
#                     7293.570516490116,
#                     16399.841133648202,
#                     1106.5102017019626,
#                     1050.4904750119156,
#                     42015.99999999998,
#                     5200.599214359727,
#                     1856.5833945973473,
#                     4476.745702447119,
#                     3841.895706164958,
#                     13267.031275728114,
#                     3891.834274142676,
#                     12667.932671386963,
#                     15923.750313738436,
#                     12009.570812352507,
#                     23031.01909636433,
#                     1585.6906264369022,
#                     2529.5575525027084,
#                     20633.108339484505,
#                     18791.746311482264,
#                     54419.29478511458,
#                     3073.7097194337857,
#                     1785.0865010711293,
#                     1747.490955400892)
table <- merge(x = table, y = rates[ , c("key", "estimate")], by.x = "x", by.y = "key", all=TRUE) #rate column is 'estiamte'
# table <-table %>% select(-estimate.y) %>% rename(estimate.x, estimate)
table = as.data.frame(apply(table,2,function(x)gsub('\\s+', '',x)))
rates = as.data.frame(apply(rates,2,function(x)gsub('\\s+', '',x)))

# q3: so for each of the investigated variants (all variants detected in any patient sample) and for each
# patient, assess the signficance of departure of the patient variant allele frequency from the null
# dist that was modeled from controls. 
# so you are calculating the probability (under null) of observing a VAF larger than the patient variant
# frequency e.g. (P(VAF >VAFp)) 
# patient variant combination:
# 3-5, your patient/variant combination is actually the union of all variants in all patient samples.
# For each variant, the null model is the exponential model fitted to the given variant's VAFs in control samples.
# you need to use the VAF of patient variants after treatment. 
# get all patient VAFs. 
patients_q3 = one_table(patients)
patients_q3[patients_q3 == 0] <- NA
# patients_q3 = merge(x = patients_q3, y = rates[ , c("key", "estimate")], by.x = "x", by.y = "key")
patients_q3 = data.frame(t(patients_q3))
colnames(patients_q3) <- as.matrix(patients_q3[1,])
patients_q3_vaf <- patients_q3[-1, ] 
patients_q3_vaf = data.frame(apply(patients_q3_vaf,2,function(x)gsub('\\s+', '',x)))
rates_ord <- rates[match(colnames(patients_q3), rates$key),]
library(plyr)
patients_q3_vaf_p <- ldply(1:(ncol(patients_q3_vaf)), function(i)
  pexp(as.numeric(as.character(patients_q3_vaf[, i])), as.numeric(as.character(rates_ord$estimate[i])), lower.tail = F))
rownames(patients_q3_vaf_p) = colnames(patients_q3_vaf)
colnames(patients_q3_vaf_p) = rownames(patients_q3_vaf)
# write.table(table, row.names = F, quote = F, sep = ",", file="~/Desktop/acb1/Sara Pour A3/Q1.csv")

# q4: correct p value for multiple testing using bonferroni, list variants that received a corrected
# p value lower or equal to 0.05. set number hypotheses tested equal to total number of patient
# variant combinations. 
num_patient_var_combo_q4 = 24*28-sum(is.na(patients_q3_vaf)) #24 patient, 28 variants, minus the variants that weren't recorded
# divide the p values above:
patients_q4_vaf_p_bf  <-patients_q3_vaf_p[]*num_patient_var_combo_q4
patients_q4_vaf_p_bf<-melt(as.matrix(patients_q4_vaf_p_bf), value.name = "p.corrected", varnames=c('variant', 'sample'))
patients_q4_vaf_p_bf = filter(patients_q4_vaf_p_bf, sample != "estimate")
# remove non-0.05 values
patients_q4_vaf_p_bf = filter(patients_q4_vaf_p_bf, p.corrected < 0.05)
# add uncorrected and rate values
unc_q3 <- melt(as.matrix(patients_q3_vaf_p))
patients_q4_vaf_p_bf= merge(patients_q4_vaf_p_bf, unc_q3, by.x = c('variant', 'sample'), by.y = c('Var1', 'Var2'), all.x = T)
patients_q4_vaf_p_bf= merge(patients_q4_vaf_p_bf, rates[,c('key', 'estimate')], by.x = c('variant'), by.y = c('key'), all.x = T)
names(patients_q4_vaf_p_bf)[names(patients_q4_vaf_p_bf) == 'value'] <- 'p'
names(patients_q4_vaf_p_bf)[names(patients_q4_vaf_p_bf) == 'estimate'] <- 'rate'
patients_q4_vaf_p_bf <-patients_q4_vaf_p_bf[,c(1,2,5,4,3)]
write.table(patients_q4_vaf_p_bf, file = "~/Desktop/acb1/Sara Pour A3/Q4_R_rates.csv", sep = ",", row.names = F, quote = F )
# output: csv listing varaints that received a corrected p value of less than or equal to 0.05
# col1: variant (genome psotion, and variant), col2: sample(patient sample), col3: rate parameters
# from step 2, col4= p values, col5 = corrected p values

# q5. Determine whether patient-variant pairs that were identified at diagnosis are more likely 
# to be significantly above background after treatment. 2-by-2 contingency table based on two 
# binary variables: 1) whether the VAF for a patient-variant pair is significantly above background 
# after treatment and 2) whether the patient-variant pair was in the list of patient-variant pairs 
# observed at diagnosis.
# create df of all patient variant pairs df
p_v_pairs_q5 = melt(as.matrix(patients_q3_vaf)) #505 combos
p_v_pairs_q5 = filter(p_v_pairs_q5, !is.na(value))
p_v_pairs_q5$p_x <-  apply( p_v_pairs_q5[ , c("Var1", "Var2") ] , 1 , paste , collapse = "_" )
# get patient variant combos for the sig values in q4
p_v_pairs_q4_q5 = melt(patients_q4_vaf_p_bf)
p_v_pairs_q4_q5 = p_v_pairs_q4_q5 %>% filter(variable == "p.corrected") #27 patient_variant pairs are enriched above bg
p_v_pairs_q4_q5$p_x <-  apply( p_v_pairs_q4_q5[ , c("sample", "variant") ] , 1 , paste , collapse = "_" )
# get patient variant combos that were in the list of pv pairs at diagnosis
list_mut = read.table(mutations, sep = ",", header = T)
list_mut$x <- apply( list_mut[ , c("chr", "pos", 'alt') ] , 1 , paste , collapse = "_" )
list_mut$p_x <-  apply( list_mut[ , c("sample", "x") ] , 1 , paste , collapse = "_" )
list_mut = as.data.frame(apply(list_mut,2,function(x)gsub('\\s+', '',x)))
# create table
above_bg = ifelse(is.na(match( p_v_pairs_q5$p_x, 
                               p_v_pairs_q4_q5$p_x)), 
                  "not above background", "above background")
in_list =ifelse(is.na(match(p_v_pairs_q5$p_x, 
                            list_mut$p_x)), 
                "not in list at diagnosis", "in list at diagnosis")
q5 = cbind(above_bg, in_list)
tbl <- table(in_list, above_bg)
fisher.test(as.matrix(tbl), alternative="greater") # < 2.2e-16
#write.table(tbl, file = "~/Desktop/acb1/Sara Pour A3/Q5_R_rates.csv", sep = ",", row.names = F, quote = F )

# null hypothesis that identification of a patient-variant pair at diagnosis is independent of whether 
# the patient-variant pair is above background after treatment.
# it is unlikely that the sample labels were randomly swapped because when you make the contingency table
# comparing if the VAF of a patient-variant pair is above background after treatment and if it was in the
# list observed at diagnosis you get a p value of above 0.5, meaning you can't reject the null. 

# q6. you now revisit the nominal p-values calculated in Step 3. Now we only need to consider patient-variant
# pairs where the variant was detected at diagnosis. Please perform the correction for multiple testing with 
# the updated set of tests, again using a Bonferroni-type correction.
p_v_pairs_q3_q6_nominal_pv = melt(as.matrix(patients_q3_vaf_p))
p_v_pairs_q3_q6_nominal_pv = p_v_pairs_q3_q6_nominal_pv %>% filter(!(is.na(value))) #28 unique mutations show up 
# only consider patient variant pairs where the variant was in the list at diagnosis
p_v_pairs_q3_q6_nominal_pv$p_x <-  apply( p_v_pairs_q3_q6_nominal_pv[ , c("Var2", "Var1") ] , 1 , paste , collapse = "_" )
p_v_pairs_q3_q6_nominal_pv = p_v_pairs_q3_q6_nominal_pv %>% filter(p_x %in% list_mut$p_x)
p_v_pairs_q3_q6_nominal_pv$p.corrected <-  p_v_pairs_q3_q6_nominal_pv$value*28
p_v_pairs_q3_q6_nominal_pv_still_sig <- p_v_pairs_q3_q6_nominal_pv %>% filter(p.corrected < 0.05) 
p_v_pairs_q3_q6_nominal_pv_non_sig <- p_v_pairs_q3_q6_nominal_pv %>% filter(p.corrected >= 0.05) 
# find dropouts, q7
still_sig =  apply( p_v_pairs_q3_q6_nominal_pv_still_sig[ , c("Var2", "Var1") ] , 1 , paste , collapse = "_" )
diff =setdiff(list_mut$p_x, still_sig) 
length(unique(gsub("(.+?_.+?)_.*" ,"\\1",  diff)))  
#write.table(diff, sep =",", quote = F, row.names = F, file="~/Desktop/acb1/Sara Pour A3/Q7_R_values.csv")
# IDS: "Patient_0298" "Patient_1815" "Patient_2662" "Patient_0531" "Patient_0872"
# add rate values
p_v_pairs_q3_q6_nominal_pv = merge(p_v_pairs_q3_q6_nominal_pv_still_sig, rates[,c('key', 'estimate')], by.x = c('Var1'), by.y = c('key'), all.x = T)
names(p_v_pairs_q3_q6_nominal_pv)[names(p_v_pairs_q3_q6_nominal_pv) == 'value'] <- 'p'
names(p_v_pairs_q3_q6_nominal_pv)[names(p_v_pairs_q3_q6_nominal_pv) == 'estimate'] <- 'rate'
names(p_v_pairs_q3_q6_nominal_pv)[names(p_v_pairs_q3_q6_nominal_pv) == 'Var1'] <- 'variant'
names(p_v_pairs_q3_q6_nominal_pv)[names(p_v_pairs_q3_q6_nominal_pv) == 'Var2'] <- 'sample'
p_v_pairs_q3_q6_nominal_pv <-p_v_pairs_q3_q6_nominal_pv[,c(1,2,6,3,5)]
#write.table(p_v_pairs_q3_q6_nominal_pv, sep =",", quote = F, row.names = F, file="~/Desktop/acb1/Sara Pour A3/Q6_R_values.csv")
# frac different:
# frac that were different when p-adjusting for controls:
(dim(patients_q4_vaf_p_bf)[1]-dim(p_v_pairs_q3_q6_nominal_pv_still_sig)[1])/dim(patients_q4_vaf_p_bf)[1] #the diff between those that are sig relative to the ones different from background when considering all mutation
# For the python rates, 0.1 were different, for the R one, 0.04 were different
