library(targets)
library(tarchetypes)
options(tidyverse.quiet = TRUE)
library(tidyverse)
source("R/functions.R")
library(future)
library(future.callr)
tar_option_set(packages = c("brms", "cowplot", "ggdag", "ggrepel", "haven", "lavaan",
                            "papaja", "readxl", "rmarkdown", "rnaturalearth", 
                            "semTools", "simr", "sjlabelled", "stringr", "tidyverse"))
plan(callr)
# full workflow
list(
  
  #### Study 1 - GPS ####
  
  # raw data files
  tar_target(fileGPS, "data/gps/individual_new.dta", format = "file"),
  tar_target(fileISO, "data/iso.csv", format = "file"),
  tar_target(fileRM1, "data/relationalMobility/rm.csv", format = "file"),
  tar_target(fileRM2, "data/relationalMobility/country-level_Ver15.sav", format = "file"),
  tar_target(fileGEO, "data/networks/1F Population Distance.xlsx", format = "file"),
  tar_target(fileLIN, "data/networks/2F Country Distance 1pml adj.xlsx", format = "file"),
  tar_target(fileISO2, "data/countries_codes_and_coordinates.csv", format = "file"),
  # load data
  tar_target(d1, loadData1(fileGPS, fileISO, fileRM1, fileRM2)),
  # plot globe
  tar_target(plotGPSWorld, plotGPSGlobe(d1, fileISO2)),
  # load covariance matrices
  tar_target(lingCov1, loadCovMat(d1, fileLIN, log = FALSE)),
  tar_target(geoCov1,  loadCovMat(d1, fileGEO, log = TRUE )),
  # fit CFAs and measurement invariance models
  tar_target(cfaListGPS, fitCFAModelCountryList_GPS(d1)),
  tar_target(configGPS, getMGCFA_GPS(d1, group.equal = "configural")),
  tar_target(metricGPS, getMGCFA_GPS(d1, group.equal = "loadings")),
  tar_target(scalarGPS, getMGCFA_GPS(d1, group.equal = c("loadings", "intercepts"))),
  tar_target(tableInvarianceGPS, createTableInvariance(configGPS, metricGPS, scalarGPS)),
  # fit bayesian multilevel models
  # altruism, trust, and positive reciprocity
  tar_target(m1.0, fitModel1.0(d1, lingCov1, geoCov1)),
  tar_target(m1.1, fitModel1.1(d1, lingCov1, geoCov1)),
  tar_target(m1.2, fitModel1.2(d1, lingCov1, geoCov1)),
  tar_target(m1.3, fitModel1.3(d1, lingCov1, geoCov1)),
  # altruism
  tar_target(m2.0, fitModel_IndItem_0(d1, lingCov1, geoCov1, itemName = "altruism")),
  tar_target(m2.1, fitModel_IndItem_1(d1, lingCov1, geoCov1, itemName = "altruism")),
  tar_target(m2.2, fitModel_IndItem_2(d1, lingCov1, geoCov1, itemName = "altruism")),
  # positive reciprocity
  tar_target(m3.0, fitModel_IndItem_0(d1, lingCov1, geoCov1, itemName = "posrecip")),
  tar_target(m3.1, fitModel_IndItem_1(d1, lingCov1, geoCov1, itemName = "posrecip")),
  tar_target(m3.2, fitModel_IndItem_2(d1, lingCov1, geoCov1, itemName = "posrecip")),
  # trust
  tar_target(m4.0, fitModel_IndItem_0(d1, lingCov1, geoCov1, itemName = "trust")),
  tar_target(m4.1, fitModel_IndItem_1(d1, lingCov1, geoCov1, itemName = "trust")),
  tar_target(m4.2, fitModel_IndItem_2(d1, lingCov1, geoCov1, itemName = "trust")),
  # loo
  tar_target(loo1.0, loo(m1.0)),
  tar_target(loo1.1, loo(m1.1)),
  tar_target(loo1.2, loo(m1.2)),
  tar_target(looCompare1a, loo_compare(loo1.0, loo1.1)),
  tar_target(looCompare1b, loo_compare(loo1.1, loo1.2)),
  # power analysis
  tar_target(power1, runPowerAnalysis1(d1, nsim = 500)),
  # posterior
  tar_target(post1.1, as_draws_array(m1.1, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "r_item"))),
  tar_target(post1.2, as_draws_array(m1.2, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "b_THREAT4", "b_rSUBSIST2", "r_item"))),
  tar_target(post1.3, as_draws_array(m1.3, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "bsp_ImeRMLRML_SEgrEQisoE2", 
                                                        "b_THREAT4", "b_rSUBSIST2", "r_item"))),
  # plots
  tar_target(plotGPS1, plotGPSResults(d1, m1.1, "figures/study1/plotGPSResults.pdf")),
  tar_target(plotGPS2, plotGPSResults(d1, m1.2, "figures/study1/plotGPSResultsWithControls.pdf")),
  tar_target(plotDAG, drawDAG()),
  
  #### Study 2 - WVS/EVS ####
  
  # raw data files
  tar_target(fileWVS, "data/wvs/EVS_WVS_Joint_v2_0_CLEAN.xlsx", format = "file"),
  # load and wrangle data
  tar_target(d2, loadData2(fileWVS, fileISO, fileRM1, fileRM2)),
  tar_target(d3, wrangleData1(d2)),
  tar_target(d4, wrangleData2(d2)),
  # plot globe
  tar_target(plotWVSWorld, plotWVSGlobe(d2, fileISO2)),
  # load covariance matrices
  tar_target(lingCov2, loadCovMat(d2, fileLIN, log = FALSE)),
  tar_target(geoCov2,  loadCovMat(d2, fileGEO, log = TRUE )),
  # fit CFAs and measurement invariance models
  # trust in different groups
  tar_target(cfaListWVS1, fitCFAModelCountryList_WVS1(d2)),
  tar_target(configWVS1, getMGCFA_WVS1(d2, group.equal = "configural")),
  tar_target(metricWVS1, getMGCFA_WVS1(d2, group.equal = "loadings")),
  tar_target(scalarWVS1, getMGCFA_WVS1(d2, group.equal = c("loadings", "intercepts"))),
  tar_target(tableInvarianceWVS1, createTableInvariance(configWVS1, metricWVS1, scalarWVS1)),
  # moral justifiability
  tar_target(cfaListWVS2, fitCFAModelCountryList_WVS2(d2)),
  tar_target(configWVS2, getMGCFA_WVS2(d2, group.equal = "configural")),
  tar_target(metricWVS2, getMGCFA_WVS2(d2, group.equal = "loadings")),
  tar_target(scalarWVS2, getMGCFA_WVS2(d2, group.equal = c("loadings", "intercepts"))),
  tar_target(tableInvarianceWVS2, createTableInvariance(configWVS2, metricWVS2, scalarWVS2)),
  # fit bayesian multilevel models
  # mentioned belonging to humanitarian or charitable organization - logistic regression
  tar_target(m5.0, fitModel5.0(d2, lingCov2, geoCov2)),
  tar_target(m5.1, fitModel5.1(d2, lingCov2, geoCov2)),
  tar_target(m5.2, fitModel5.2(d2, lingCov2, geoCov2)),
  tar_target(m5.3, fitModel5.3(d2, lingCov2, geoCov2)),
  # most people can be trusted - logistic regression
  tar_target(m6.0, fitModel6.0(d2, lingCov2, geoCov2)),
  tar_target(m6.1, fitModel6.1(d2, lingCov2, geoCov2)),
  tar_target(m6.2, fitModel6.2(d2, lingCov2, geoCov2)),
  tar_target(m6.3, fitModel6.3(d2, lingCov2, geoCov2)),
  # trusting different groups - ordinal regression
  tar_target(m7.0, fitModel7.0(d3, lingCov2, geoCov2)),
  tar_target(m7.1, fitModel7.1(d3, lingCov2, geoCov2)),
  tar_target(m7.2, fitModel7.2(d3, lingCov2, geoCov2)),
  tar_target(m7.3, fitModel7.3(d3, lingCov2, geoCov2)),
  # moral justifiability - ordinal regression
  tar_target(m8.0, fitModel8.0(d4, lingCov2, geoCov2)),
  tar_target(m8.1, fitModel8.1(d4, lingCov2, geoCov2)),
  tar_target(m8.2, fitModel8.2(d4, lingCov2, geoCov2)),
  tar_target(m8.3, fitModel8.3(d4, lingCov2, geoCov2)),
  # loo
  tar_target(loo5.0, loo(m5.0)),
  tar_target(loo5.1, loo(m5.1)),
  tar_target(loo6.0, loo(m6.0)),
  tar_target(loo6.1, loo(m6.1)),
  tar_target(loo8.0, loo(m8.0)),
  tar_target(loo8.1, loo(m8.1)),
  tar_target(looCompare5, loo_compare(loo5.0, loo5.1)),
  tar_target(looCompare6, loo_compare(loo6.0, loo6.1)),
  tar_target(looCompare8, loo_compare(loo8.0, loo8.1)),
  # power analysis
  tar_target(power5, runPowerAnalysis5(d2, nsim = 1000)),
  tar_target(power6, runPowerAnalysis6(d2, nsim = 1000)),
  tar_target(power7, runPowerAnalysis7(d3, nsim = 100)),
  tar_target(power8, runPowerAnalysis8(d4, nsim = 100)),
  # posterior
  tar_target(post5.1, as_draws_array(m5.1, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso"))),
  tar_target(post5.2, as_draws_array(m5.2, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso"))),
  tar_target(post5.3, as_draws_array(m5.3, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "bsp_ImeRMLRML_SEgrEQisoE2"))),
  tar_target(post6.1, as_draws_array(m6.1, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso"))),
  tar_target(post6.2, as_draws_array(m6.2, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso"))),
  tar_target(post6.3, as_draws_array(m6.3, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "bsp_ImeRMLRML_SEgrEQisoE2"))),
  tar_target(post7.1, as_draws_array(m7.1, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "r_group"))),
  tar_target(post7.2, as_draws_array(m7.2, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "r_group"))),
  tar_target(post7.3, as_draws_array(m7.3, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "bsp_ImeRMLRML_SEgrEQisoE2", "r_group"))),
  tar_target(post8.1, as_draws_array(m8.1, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "r_item"))),
  tar_target(post8.2, as_draws_array(m8.2, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "r_item"))),
  tar_target(post8.3, as_draws_array(m8.3, variable = c("b_Intercept", "bsp_meRMLRML_SEgrEQiso", "bsp_ImeRMLRML_SEgrEQisoE2", "r_item"))),
  # plots
  tar_target(plotWVS_5.1, plotWVSResults5(d2, m5.1, "figures/study2/plotCharitableOrg.pdf")),
  tar_target(plotWVS_5.2, plotWVSResults5(d2, m5.2, "figures/study2/plotCharitableOrgWithControls.pdf")),
  tar_target(plotWVS_6.1, plotWVSResults6(d2, m6.1, "figures/study2/plotGenTrust.pdf")),
  tar_target(plotWVS_6.2, plotWVSResults6(d2, m6.2, "figures/study2/plotGenTrustWithControls.pdf")),
  tar_target(plotWVS_7.1, plotWVSResults7(d3, m7.1, "figures/study2/plotTrustGroups.pdf")),
  tar_target(plotWVS_7.2, plotWVSResults7(d3, m7.2, "figures/study2/plotTrustGroupsWithControls.pdf")),
  tar_target(plotWVS_7.3, plotWVSResults7Quad(d3, post7.3, "figures/study2/plotTrustGroupsWithControlsAndQuadratic.pdf")),
  tar_target(plotWVS_8.1, plotWVSResults8(d4, m8.1, "figures/study2/plotMoralJust.pdf")),
  tar_target(plotWVS_8.2, plotWVSResults8(d4, m8.2, "figures/study2/plotMoralJustWithControls.pdf")),
  tar_target(plotWVS_8.3, plotWVSResults8Quad(d4, post8.3, "figures/study2/plotMoralJustWithControlsAndQuadratic.pdf")),
  
  #### Manuscript ####
  
  # tables for manuscript
  tar_target(tableRaw1, createTableRawStudy1(d1)),
  tar_target(tableRaw2, createTableRawStudy2(d2)),
  tar_target(tableQuadratic, createTableQuadratic(post1.3, post5.3, post6.3, post7.3, post8.3)),
  tar_target(tablePower, createTablePower(power1, power5, power6, power7, power8)),
  # render manuscript
  tar_render(manuscript, "manuscript.Rmd")
  
)