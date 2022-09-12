# custom functions

# load data
loadData1 <- function(fileGPS, fileISO, fileRM1, fileRM2) {
  out <-
    # load gps data
    read_dta(fileGPS) %>%
    # join 2-digit iso code
    left_join(read_csv(fileISO), by = "isocode") %>%
    # add unique id row and duplicate iso columns
    mutate(id = 1:nrow(.), iso2 = iso, iso3 = iso) %>%
    # join country-level rm data and standardised controls
    inner_join(read_csv(fileRM1), by = "iso") %>%
    left_join(read_sav(fileRM2) %>%
                mutate(COUNTRY = get_labels(COUNTRY, drop.unused = TRUE),
                       THREAT4 = as.numeric(scale(THREAT4)),
                       rSUBSIST2 = as.numeric(scale(rSUBSIST2))) %>%
                select(COUNTRY, THREAT4, rSUBSIST2), by = c("countryRM" = "COUNTRY")) %>%
    # pivot longer for cooperation items
    pivot_longer(cols = c(posrecip, altruism, trust),
                 names_to = "item",
                 values_to = "coop") %>%
    # drop NAs in outcome variable
    drop_na(coop)
  return(out)
}

# make global plot for gps
plotGPSGlobe <- function(d1, fileISO2) {
  # get counts, long, and lat for plotting
  counts <-
    d1 %>%
    group_by(iso) %>%
    summarise(count = length(unique(id))) %>%
    left_join(read_csv(fileISO2), by = c("iso" = "Alpha-2 code"))
  # world map with counts
  out <-
    ne_countries(scale = "medium", returnclass = "sf") %>%
    filter(continent %in% unique(continent)[1:6]) %>%
    left_join(counts, by = c("iso_a2" = "iso")) %>%
    rename(Continent = continent) %>%
    ggplot() +
    geom_sf(colour = NA) + 
    geom_point(aes(x = `Longitude (average)`, 
                   y = `Latitude (average)`, 
                   size = count,
                   colour = Continent),
               alpha = 0.6) +
    coord_sf(ylim = c(-50, 90), xlim = c(), datum = NA) +
    scale_size_area(breaks = c(5, seq(25,125,25))) +
    theme(panel.background = element_rect(fill = 'white'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  # save plot
  ggsave("figures/study1/plotGPSWorldMap.pdf", height = 4, width = 8)
  return(out)
}

# load covariance matrices
loadCovMat <- function(d, file, log) {
  # load distance matrix
  out <- 
    read_excel(file, na = "") %>%
    select(-ISO) %>%
    as.matrix()
  rownames(out) <- colnames(out)
  # keep only countries from data
  out <- out[rownames(out) %in% d$iso,
             colnames(out) %in% d$iso]
  # log distances?
  if (log) out <- log(out)
  # distances between 0 and 1
  out <- out / max(out)
  diag(out) <- 0
  # 1 - distance = proximity (covariance)
  out <- 1 - out
  return(out)
}

# basic cfa fitting function - gps
fitCFAModel_GPS <- function(d) {
  model <- "prosocial =~ altruism + trust + posrecip"
  out <- cfa(model, data = d)
  return(out)
}

# fit cfa model in each country separately - gps
fitCFAModelCountryList_GPS <- function(d) {
  d %>%
    pivot_wider(names_from = item, values_from = coop) %>%
    group_by(iso) %>%
    nest() %>%
    mutate(cfaModel = map(data, fitCFAModel_GPS)) %>%
    mutate(
      # which countries fail?
      converged  = map(cfaModel, function(x) lavInspect(x, "converged")),
      post.check = map(cfaModel, function(x) lavInspect(x, "post.check"))
    ) %>%
    select(-data) %>%
    unnest(c(converged, post.check)) %>%
    ungroup()
}

# multi-group confirmatory factor analysis - gps
getMGCFA_GPS <- function(d, group.equal, group.partial = "") {
  d <- 
    d %>%
    # pivot data wider
    pivot_wider(names_from = item, values_from = coop) %>%
    # remove japan, chile, and brazil (failed to converge in japan, did not find
    # factor structure in chile or brazil)
    filter(!(iso %in% c("JP", "CL", "BR")))
  # fit multi-group confirmatory factor analysis
  # https://m-clark.github.io/posts/2019-08-05-comparing-latent-variables/#supplemental-measurement-invariance
  out <- 
    measEq.syntax(configural.model = "prosocial =~ posrecip + altruism + trust", 
                  data = d, 
                  group = "iso",
                  group.equal = group.equal,
                  group.partial = group.partial,
                  return.fit = TRUE)
  return(out)
}

# measurement invariance results in table
createTableInvariance <- function(config, metric, scalar) {
  Round <- function(x) format(round(x, 2), nsmall = 2)
  tibble(
    Model = c("Configural invariance", "Metric invariance", "Scalar invariance"),
    RMSEA = c(
      Round(fitMeasures(config)["rmsea"]),
      Round(fitMeasures(metric)["rmsea"]),
      Round(fitMeasures(scalar)["rmsea"])
    ),
    CFI = c(
      Round(fitMeasures(config)["cfi"]),
      Round(fitMeasures(metric)["cfi"]),
      Round(fitMeasures(scalar)["cfi"])
    ),
    SRMR = c(
      Round(fitMeasures(config)["srmr"]),
      Round(fitMeasures(metric)["srmr"]),
      Round(fitMeasures(scalar)["srmr"])
    )
  )
}

# fit null model
fitModel1.0 <- function(d, lingCov, geoCov) {
  brm(formula = bf("coop ~ 0 + Intercept + (1 | item) + (1 | id) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = gaussian,
      prior = c(prior(normal(0, 0.1), class = b),
                prior(exponential(5), class = sd),
                prior(exponential(5), class = sigma)), 
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 4000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      sample_prior = "yes", seed = 2113)
}

# fit RML model without controls
fitModel1.1 <- function(d, lingCov, geoCov) {
  brm(formula = bf("coop ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + (1 + RML | item) + (1 | id) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = gaussian,
      prior = c(prior(normal(0, 0.1), class = b),
                prior(exponential(5), class = sd),
                prior(exponential(5), class = sigma),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)), 
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 4000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      sample_prior = "yes", seed = 2113)
}

# fit RML model with controls
fitModel1.2 <- function(d, lingCov, geoCov) {
  brm(formula = bf("coop ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + THREAT4 + rSUBSIST2 + (1 + RML + THREAT4 + rSUBSIST2 | item) + (1 | id) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = gaussian,
      prior = c(prior(normal(0, 0.1), class = b),
                prior(exponential(5), class = sd),
                prior(exponential(5), class = sigma),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)), 
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 4000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      sample_prior = "yes", seed = 2113)
}

# fit RML model with controls and quadratic effect
fitModel1.3 <- function(d, lingCov, geoCov) {
  brm(formula = bf('coop ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + I(me(RML, RML_SE, gr = iso)^2) + THREAT4 + rSUBSIST2 + 
                   (1 + RML + I(RML^2) + THREAT4 + rSUBSIST2 | item) + (1 | id) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)'), 
      data = d, family = gaussian,
      prior = c(prior(normal(0, 0.1), class = b),
                prior(exponential(5), class = sd),
                prior(exponential(5), class = sigma),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)), 
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 4000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      sample_prior = "yes", seed = 2113)
}

# each item alone - null model
fitModel_IndItem_0 <- function(d, lingCov, geoCov, itemName) {
  # get item alone
  d <- 
    d %>% 
    filter(item == itemName) %>%
    rename(!!itemName := coop)
  # fit model
  brm(formula = bf(paste0(itemName, " ~ 0 + Intercept + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)")), 
      data = d, family = gaussian,
      prior = c(prior(normal(0, 0.1), class = b),
                prior(exponential(5), class = sd),
                prior(exponential(5), class = sigma)), 
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 4000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      sample_prior = "yes", seed = 2113)
}

# each item alone - RML model without controls
fitModel_IndItem_1 <- function(d, lingCov, geoCov, itemName) {
  # get item alone
  d <- 
    d %>% 
    filter(item == itemName) %>%
    rename(!!itemName := coop)
  # fit model
  brm(formula = bf(paste0(itemName, " ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)")), 
      data = d, family = gaussian, 
      prior = c(prior(normal(0, 0.1), class = b),
                prior(exponential(5), class = sd),
                prior(exponential(5), class = sigma),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 4000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      sample_prior = "yes", seed = 2113)
}

# each item alone - RML model with controls
fitModel_IndItem_2 <- function(d, lingCov, geoCov, itemName) {
  # get item alone
  d <- 
    d %>% 
    filter(item == itemName) %>%
    rename(!!itemName := coop)
  # fit model
  brm(formula = bf(paste0(itemName, " ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + THREAT4 + rSUBSIST2 + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)")), 
      data = d, family = gaussian,
      prior = c(prior(normal(0, 0.1), class = b),
                prior(exponential(5), class = sd),
                prior(exponential(5), class = sigma),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 4000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99, max_treedepth = 15),
      sample_prior = "yes", seed = 2113)
}

# run power analysis #1
runPowerAnalysis1 <- function(d1, std = 0.28, nsim) {
  # standardise outcome
  d1$coop <- as.numeric(scale(d1$coop))
  # standardise RML predictor
  d1$RML <- as.numeric(scale(d1$RML))
  # listwise deletion for predictors
  d1 <- drop_na(d1, RML, THREAT4, rSUBSIST2)
  # model of interest
  model <- lmer(coop ~ 1 + RML + THREAT4 + rSUBSIST2 + 
                   (1 + RML + THREAT4 + rSUBSIST2 | item) + 
                  (1 | id) + (1 | iso), data = d1)
  # fix effect size for simulation
  fixef(model)['RML'] <- std
  # power sim function
  out <- powerSim(model, test = simr::fixed("RML", "z"), 
                  nsim = nsim, seed = 2113, progress = FALSE)
  return(out)
}

# plot gps results
plotGPSResults <- function(d, model, file) {
  # get overall conditional effects
  condOverall <- conditional_effects(model, re_formula = . ~ 1 + (1 + RML | item))
  # get conditional effects for individual items
  items <- data.frame(item = unique(d$item))
  rownames(items) <- unique(d$item)
  condItems <- conditional_effects(model, conditions = items, re_formula = . ~ 1 + (1 + RML | item))
  # coop plotting function
  pA <-
    condOverall$RML %>%
    ggplot() +
    geom_ribbon(aes(x = RML, ymax = upper__, ymin = lower__), fill = "grey", alpha = 0.5) +
    geom_line(aes(x = RML, y = estimate__)) +
    scale_y_continuous(name = "Overall prosociality", limits = c(-1.2, 0.7)) +
    scale_x_continuous(name = "Relational mobility", limits = c(-0.47, 0.47)) +
    theme_classic()
  # item plotting function
  plotFun <- function(var, ylab, col) {
    d %>%
      filter(item == var) %>%
      group_by(iso) %>%
      summarise(y = mean(coop), y_se = sd(coop) / sqrt(n()),
                x = mean(RML), x_se = mean(RML_SE)) %>%
      ggplot() +
      geom_errorbarh(aes(y = y, xmin = x - x_se, xmax = x + x_se), alpha = 0.3, colour = col) +
      geom_errorbar(aes(x = x, ymin = y - y_se, ymax = y + y_se), alpha = 0.3, colour = col) +
      geom_point(aes(x = x, y = y), size = 0.5, colour = col, fill = col) +
      geom_text_repel(aes(x = x, y = y, label = iso), size = 2) +
      geom_ribbon(data = condItems$RML[condItems$RML$item == var,], 
                  aes(x = RML, ymax = upper__, ymin = lower__), fill = col, alpha = 0.5) +
      geom_line(data = condItems$RML[condItems$RML$item == var,],
                aes(x = RML, y = estimate__)) +
      scale_y_continuous(name = ylab, limits = c(-1.2, 0.7)) +
      scale_x_continuous(name = "Relational mobility", limits = c(-0.47, 0.47)) +
      theme_classic()
  }
  # item plots
  pB <- plotFun("altruism", "Altruism", "#2E604A")
  pC <- plotFun("posrecip", "Positive reciprocity", "#DBB165")
  pD <- plotFun("trust", "Trust", "#D1362F")
  # put together
  out <- plot_grid(pA, pB, pC, pD, labels = letters[1:4], nrow = 2)
  ggsave(out, filename = file, height = 6.5, width = 6.5)
  return(out)
}

# draw causal model
drawDAG <- function() {
  # according to pnas paper: https://www.pnas.org/content/115/29/7521
  # main antecedents of relational mobility
  # - interdependent subsistence style ("Sub")
  # - geoclimate harshness ("Harsh")
  #
  # papers showing harsh environments also predict prosociality (independently of relational mobility)
  # - positive relationship: https://link.springer.com/chapter/10.1007/978-3-030-15800-2_4
  # - negative relationship: https://www.biorxiv.org/content/10.1101/663518v1
  # - null relationship: https://www.sciencedirect.com/science/article/pii/S1090513816300630
  #
  # theory linking subsistence-type to cooperation and interdependence (independently of relational mobility)
  # - https://science.sciencemag.org/content/344/6184/603.full
  #
  # and unmeasured confounds captured by geographic and linguistic distance
  #
  dag_coords <-
    tibble(name = c("RM", "Pro", "Sub", "Harsh", "U", "Geo", "Lin"),
           x    = c(1, 2, 1.5, 1.5, 0, 0, 0),
           y    = c(0, 0, 1, -1, 0, 1, -1))
  out <-
    dagify(Pro ~ RM + Harsh + Sub + U,
           RM ~ Harsh + Sub + U,
           Sub ~ U,
           Harsh ~ U,
           U ~ Geo + Lin,
           exposure = "RM",
           outcome = "Pro",
           latent = "U",
           coords = dag_coords) %>%
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point(
      data = function(x) filter(x, name == "U"),
      alpha = 1/2, size = 20, show.legend = F, colour = "lightgrey"
      ) +
    geom_dag_text(color = "black") +
    geom_dag_edges_link(
      data = function(x) filter(x, !(name == "U" & to == "Pro"))
    ) +
    geom_dag_edges_arc(
      data = function(x) filter(x, name == "U" & to == "Pro"),
      curvature = 0.3
    ) +
    theme(axis.line = element_blank(),
          axis.ticks = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank()) +
    theme_void()
  # save
  ggsave(out, filename = "figures/plotDAG.pdf", height = 4, width = 6)
  return(out)
}

# load data wvs evs
loadData2 <- function(fileWVS, fileISO, fileRM1, fileRM2) {
  out <-
    # load wvs / evs data
    readxl::read_xlsx(fileWVS, sheet = 1) %>%
    # add unique id row and duplicate iso columns
    mutate(iso = cntry_AN, iso2 = iso, iso3 = iso) %>%
    # join country labels
    left_join(read_csv(fileISO), by = "iso") %>%
    # join country-level rm data and standardised controls
    inner_join(read_csv(fileRM1), by = "iso") %>%
    left_join(read_sav(fileRM2) %>% 
                mutate(COUNTRY = get_labels(COUNTRY, drop.unused = TRUE),
                       THREAT4 = as.numeric(scale(THREAT4)),
                       rSUBSIST2 = as.numeric(scale(rSUBSIST2))) %>%
                select(COUNTRY, THREAT4, rSUBSIST2), by = c("countryRM" = "COUNTRY"))
  return(out)
}

# modify data for models 7.0 - 7.2
wrangleData1 <- function(d) {
  d %>%
    # pivot longer
    pivot_longer(cols = c(D001_B, G007_18_B, G007_33_B, 
                          G007_34_B, G007_35_B, G007_36_B),
                 names_to = "group",
                 values_to = "trust") %>%
    # group names
    mutate(group = ifelse(group == "D001_B",    "Your family",
                          ifelse(group == "G007_18_B", "Your neighbourhood",
                                 ifelse(group == "G007_33_B", "People you know personally",
                                        ifelse(group == "G007_34_B", "People you meet for the first time",
                                               ifelse(group == "G007_35_B", "People of another religion",
                                                      ifelse(group == "G007_36_B", "People of another nationality", NA))))))) %>%
    # listwise deletion of non-response
    filter(trust %in% 1:4) %>%
    # reorder trust dv (4 = trust completely, 3 = trust somewhat, 2 = do not trust very much, 1 = do not trust at all)
    mutate(trust = 5 - trust)
}

# modify data for models 8.0 - 8.2
wrangleData2 <- function(d) {
  d2 <- d %>%
    # pivot longer
    pivot_longer(cols = c(F114A, F115, F116, F117),
                 names_to = "item",
                 values_to = "justify") %>%
    # group names
    mutate(item = ifelse(item == "F114A", "Claiming government benefits to which you are not entitled",
                         ifelse(item == "F115",  "Avoiding a fare on public transport",
                                ifelse(item == "F116",  "Cheating on taxes",
                                       ifelse(item == "F117",  "Someone accepting a bribe", NA))))) %>%
    # listwise deletion of non-response
    filter(justify %in% 1:10) %>%
    # reorder justify dv (1 = always justifiable, 10 = never justifiable)
    mutate(justify = 11 - justify)
}

# make global plot for wvs
plotWVSGlobe <- function(d2, fileISO2) {
  # get counts, long, and lat for plotting
  counts <-
    d2 %>%
    group_by(iso) %>%
    summarise(count = n()) %>%
    left_join(read_csv(fileISO2), by = c("iso" = "Alpha-2 code"))
  # world map with counts
  out <-
    ne_countries(scale = "medium", returnclass = "sf") %>%
    filter(continent %in% unique(continent)[1:6]) %>%
    left_join(counts, by = c("iso_a2" = "iso")) %>%
    rename(Continent = continent) %>%
    ggplot() +
    geom_sf(colour = NA) + 
    geom_point(aes(x = `Longitude (average)`, 
                   y = `Latitude (average)`, 
                   size = count,
                   colour = Continent),
               alpha = 0.6) +
    coord_sf(ylim = c(-50, 90), xlim = c(), datum = NA) +
    scale_size_area(breaks = c(5, seq(25,125,25))) +
    theme(panel.background = element_rect(fill = 'white'),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
  # save plot
  ggsave("figures/study2/plotWVSWorldMap.pdf", height = 4, width = 8)
  return(out)
}

# basic cfa fitting function - wvs1
fitCFAModel_WVS1 <- function(d) {
  model <- "trustClose =~ D001_B + G007_18_B + G007_33_B\ntrustStranger =~ G007_34_B + G007_35_B + G007_36_B"
  out <- cfa(model, data = d)
  return(out)
}

# fit cfa model in each country separately - wvs1
fitCFAModelCountryList_WVS1 <- function(d) {
  d %>%
    # listwise deletion
    filter(
      D001_B    %in% 1:4 &
      G007_18_B %in% 1:4 &
      G007_33_B %in% 1:4 &
      G007_34_B %in% 1:4 &
      G007_35_B %in% 1:4 &
      G007_36_B %in% 1:4
    ) %>%
    # fit models
    group_by(iso) %>%
    nest() %>%
    mutate(cfaModel = map(data, fitCFAModel_WVS1)) %>%
    mutate(
      # which countries fail?
      converged  = map(cfaModel, function(x) lavInspect(x, "converged")),
      post.check = map(cfaModel, function(x) lavInspect(x, "post.check"))
    ) %>%
    select(-data) %>%
    unnest(c(converged, post.check)) %>%
    ungroup()
}

# multi-group confirmatory factor analysis - wvs 1
getMGCFA_WVS1 <- function(d, group.equal, group.partial = "") {
  # listwise deletion
  d <-
    filter(d, 
      D001_B    %in% 1:4 &
      G007_18_B %in% 1:4 &
      G007_33_B %in% 1:4 &
      G007_34_B %in% 1:4 &
      G007_35_B %in% 1:4 &
      G007_36_B %in% 1:4
    )
  # fit multi-group confirmatory factor analysis
  # https://m-clark.github.io/posts/2019-08-05-comparing-latent-variables/#supplemental-measurement-invariance
  out <- 
    measEq.syntax(configural.model = "trustClose =~ D001_B + G007_18_B + G007_33_B\ntrustStranger =~ G007_34_B + G007_35_B + G007_36_B", 
                  data = d, 
                  group = "iso",
                  group.equal = group.equal,
                  group.partial = group.partial,
                  return.fit = TRUE)
  return(out)
}

# basic cfa fitting function - wvs2
fitCFAModel_WVS2 <- function(d) {
  model <- "moral =~ F114A + F115 + F116 + F117"
  out <- cfa(model, data = d)
  return(out)
}

# fit cfa model in each country separately - wvs2
fitCFAModelCountryList_WVS2 <- function(d) {
  d %>%
    # listwise deletion
    filter(
      F114A %in% 1:10 &
      F115  %in% 1:10 &
      F116  %in% 1:10 &
      F117  %in% 1:10
    ) %>%
    # fit models
    group_by(iso) %>%
    nest() %>%
    mutate(cfaModel = map(data, fitCFAModel_WVS2)) %>%
    mutate(
      # which countries fail?
      converged  = map(cfaModel, function(x) lavInspect(x, "converged")),
      post.check = map(cfaModel, function(x) lavInspect(x, "post.check"))
    ) %>%
    select(-data) %>%
    unnest(c(converged, post.check)) %>%
    ungroup()
}

# multi-group confirmatory factor analysis - wvs 2
getMGCFA_WVS2 <- function(d, group.equal, group.partial = "") {
  # listwise deletion
  d <-
    filter(d, 
      F114A %in% 1:10 &
      F115  %in% 1:10 &
      F116  %in% 1:10 &
      F117  %in% 1:10
    )
  # fit multi-group confirmatory factor analysis
  # https://m-clark.github.io/posts/2019-08-05-comparing-latent-variables/#supplemental-measurement-invariance
  out <- 
    measEq.syntax(configural.model = "moral =~ F114A + F115 + F116 + F117", 
                  data = d, 
                  group = "iso",
                  group.equal = group.equal,
                  group.partial = group.partial,
                  return.fit = TRUE)
  return(out)
}

# fit model 5.0
fitModel5.0 <- function(d, lingCov, geoCov) {
  # listwise deletion of non-response
  d <- filter(d, A080_01 %in% 0:1)
  # fit model
  brm(formula = bf("A080_01 ~ 0 + Intercept + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(2), class = sd)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99),
      sample_prior = "yes", seed = 2113)
}

# fit model 5.1
fitModel5.1 <- function(d, lingCov, geoCov) {
  # listwise deletion of non-response
  d <- filter(d, A080_01 %in% 0:1)
  # fit model
  brm(formula = bf("A080_01 ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(2), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99),
      sample_prior = "yes", seed = 2113)
}

# fit model 5.2
fitModel5.2 <- function(d, lingCov, geoCov) {
  # listwise deletion of non-response
  d <- filter(d, A080_01 %in% 0:1)
  # fit model
  brm(formula = bf("A080_01 ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + THREAT4 + rSUBSIST2 + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(2), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99),
      sample_prior = "yes", seed = 2113)
}

# fit model 5.3
fitModel5.3 <- function(d, lingCov, geoCov) {
  # listwise deletion of non-response
  d <- filter(d, A080_01 %in% 0:1)
  # fit model
  brm(formula = bf('A080_01 ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + I(me(RML, RML_SE, gr = iso)^2) + 
                   THREAT4 + rSUBSIST2 + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)'), 
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(2), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99),
      sample_prior = "yes", seed = 2113)
}

# fit model 6.0
fitModel6.0 <- function(d, lingCov, geoCov) {
  d <- 
    d %>%
    # listwise deletion of non-response
    filter(A165 %in% 1:2) %>%
    # 0 = can't be too careful, 1 = most people can be trusted
    mutate(A165 = ifelse(A165 == 2, 0, A165))
  # fit model
  brm(formula = bf("A165 ~ 0 + Intercept + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(2), class = sd)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99),
      sample_prior = "yes", seed = 2113)
}

# fit model 6.1
fitModel6.1 <- function(d, lingCov, geoCov) {
  d <- 
    d %>%
    # listwise deletion of non-response
    filter(A165 %in% 1:2) %>%
    # 0 = can't be too careful, 1 = most people can be trusted
    mutate(A165 = ifelse(A165 == 2, 0, A165))
  # fit model
  brm(formula = bf("A165 ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(2), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99),
      sample_prior = "yes", seed = 2113)
}

# fit model 6.2
fitModel6.2 <- function(d, lingCov, geoCov) {
  d <- 
    d %>%
    # listwise deletion of non-response
    filter(A165 %in% 1:2) %>%
    # 0 = can't be too careful, 1 = most people can be trusted
    mutate(A165 = ifelse(A165 == 2, 0, A165))
  # fit model
  brm(formula = bf("A165 ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + THREAT4 + rSUBSIST2 + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(2), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99),
      sample_prior = "yes", seed = 2113)
}

# fit model 6.3
fitModel6.3 <- function(d, lingCov, geoCov) {
  d <- 
    d %>%
    # listwise deletion of non-response
    filter(A165 %in% 1:2) %>%
    # 0 = can't be too careful, 1 = most people can be trusted
    mutate(A165 = ifelse(A165 == 2, 0, A165))
  # fit model
  brm(formula = bf('A165 ~ 0 + Intercept + me(RML, RML_SE, gr = iso) + I(me(RML, RML_SE, gr = iso)^2) + 
                   THREAT4 + rSUBSIST2 + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)'), 
      data = d, family = bernoulli,
      prior = c(prior(normal(0, 1), class = b),
                prior(exponential(2), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      control = list(adapt_delta = 0.99),
      sample_prior = "yes", seed = 2113)
}

# fit model 7.0
fitModel7.0 <- function(d, lingCov, geoCov) {
  brm(formula = bf("trust ~ 1 + (1 | group) + (1 | uniqid) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = cumulative,
      prior = c(prior(normal(0, 2), class = Intercept),
                prior(exponential(4), class = sd)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 2113)
}

# fit model 7.1
fitModel7.1 <- function(d, lingCov, geoCov) {
  brm(formula = bf("trust ~ 1 + me(RML, RML_SE, gr = iso) + (1 + RML | group) + (1 | uniqid) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = cumulative,
      prior = c(prior(normal(0, 2), class = Intercept),
                prior(normal(0, 0.5), class = b),
                prior(exponential(4), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 2113)
}

# fit model 7.2
fitModel7.2 <- function(d, lingCov, geoCov) {
  brm(formula = bf("trust ~ 1 + me(RML, RML_SE, gr = iso) + THREAT4 + rSUBSIST2 + (1 + RML + THREAT4 + rSUBSIST2 | group) + (1 | uniqid) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = cumulative,
      prior = c(prior(normal(0, 2), class = Intercept),
                prior(normal(0, 0.5), class = b),
                prior(exponential(4), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 2113,
      control = list(adapt_delta = 0.9, max_treedepth = 15))
}

# fit model 7.3
fitModel7.3 <- function(d, lingCov, geoCov) {
  brm(formula = bf('trust ~ 1 + me(RML, RML_SE, gr = iso) + I(me(RML, RML_SE, gr = iso)^2) + 
                   THREAT4 + rSUBSIST2 + (1 + RML + I(RML^2) + THREAT4 + rSUBSIST2 | group) + (1 | uniqid) + 
                   (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)'), 
      data = d, family = cumulative,
      prior = c(prior(normal(0, 2), class = Intercept),
                prior(normal(0, 0.5), class = b),
                prior(exponential(4), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 2113,
      control = list(adapt_delta = 0.9, max_treedepth = 15))
}

# fit model 8.0
fitModel8.0 <- function(d, lingCov, geoCov) {
  brm(formula = bf("justify ~ 1 + (1 | item) + (1 | uniqid) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = cumulative,
      prior = c(prior(normal(0, 2), class = Intercept),
                prior(exponential(4), class = sd)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 2113,
      control = list(adapt_delta = 0.85))
}

# fit model 8.1
fitModel8.1 <- function(d, lingCov, geoCov) {
  brm(formula = bf("justify ~ 1 + me(RML, RML_SE, gr = iso) + (1 + RML | item) + (1 | uniqid) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = cumulative,
      prior = c(prior(normal(0, 2), class = Intercept),
                prior(normal(0, 0.5), class = b),
                prior(exponential(4), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 2113,
      control = list(adapt_delta = 0.85))
}

# fit model 8.2
fitModel8.2 <- function(d, lingCov, geoCov) {
  brm(formula = bf("justify ~ 1 + me(RML, RML_SE, gr = iso) + THREAT4 + rSUBSIST2 + (1 + RML + THREAT4 + rSUBSIST2 | item) + (1 | uniqid) + (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)"), 
      data = d, family = cumulative,
      prior = c(prior(normal(0, 2), class = Intercept),
                prior(normal(0, 0.5), class = b),
                prior(exponential(4), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 2113,
      control = list(adapt_delta = 0.85))
}

# fit model 8.3
fitModel8.3 <- function(d, lingCov, geoCov) {
  brm(formula = bf('justify ~ 1 + me(RML, RML_SE, gr = iso) + I(me(RML, RML_SE, gr = iso)^2) + 
                   THREAT4 + rSUBSIST2 + (1 + RML + I(RML^2) + THREAT4 + rSUBSIST2 | item) + (1 | uniqid) + 
                   (1 | gr(iso, cov = lingCov)) + (1 | gr(iso2, cov = geoCov)) + (1 | iso3)'), 
      data = d, family = cumulative,
      prior = c(prior(normal(0, 2), class = Intercept),
                prior(normal(0, 0.5), class = b),
                prior(exponential(4), class = sd),
                prior(normal(0, 0.1), class = meanme),
                prior(exponential(5), class = sdme)),
      data2 = list(lingCov = lingCov, geoCov = geoCov),
      iter = 2000, chains = 4, cores = 4,
      sample_prior = "yes", seed = 2113,
      control = list(adapt_delta = 0.85))
}

# run power analysis #5
# for log odds slope effect size thresholds, see: https://www.tandfonline.com/doi/full/10.1080/03610911003650383
# and see: https://easystats.github.io/effectsize/articles/interpret.html
runPowerAnalysis5 <- function(d2, stdLO = log(1.80), nsim) {
  # listwise deletion of non-response
  d2 <- filter(d2, A080_01 %in% 0:1)
  # standardise RML predictor
  d2$RML <- as.numeric(scale(d2$RML))
  # listwise deletion for predictors
  d2 <- drop_na(d2, RML, THREAT4, rSUBSIST2)
  # model of interest
  model <- glmer(A080_01 ~ 1 + RML + THREAT4 + rSUBSIST2 + (1 | iso), 
                 data = d2, family = "binomial")
  # fix effect size (log odds slope)
  fixef(model)['RML'] <- stdLO
  # power sim function
  out <- powerSim(model, test = simr::fixed("RML", "z"), 
                  nsim = nsim, seed = 2113, progress = FALSE)
  return(out)
}

# run power analysis #6
# for log odds slope effect size thresholds, see: https://www.tandfonline.com/doi/full/10.1080/03610911003650383
# and see: https://easystats.github.io/effectsize/articles/interpret.html
runPowerAnalysis6 <- function(d2, stdLO = log(1.78), nsim) {
  # amend data for modelling
  d2 <- 
    d2 %>%
    # listwise deletion of non-response
    filter(A165 %in% 1:2) %>%
    # 0 = can't be too careful, 1 = most people can be trusted
    mutate(A165 = ifelse(A165 == 2, 0, A165))
  # standardise RML predictor
  d2$RML <- as.numeric(scale(d2$RML))
  # listwise deletion for predictors
  d2 <- drop_na(d2, RML, THREAT4, rSUBSIST2)
  # model of interest
  model <- glmer(A165 ~ 1 + RML + THREAT4 + rSUBSIST2 + (1 | iso), 
                 data = d2, family = "binomial")
  # fix effect size (log odds slope)
  fixef(model)['RML'] <- stdLO
  # power sim function
  out <- powerSim(model, test = simr::fixed("RML", "z"), 
                  nsim = nsim, seed = 2113, progress = FALSE)
  return(out)
}

# run power analysis #7
runPowerAnalysis7 <- function(d3, std = 0.20, nsim) {
  # standardise outcome variable
  d3$trust <- as.numeric(scale(d3$trust))
  # standardise RML predictor
  d3$RML <- as.numeric(scale(d3$RML))
  # listwise deletion for predictors
  d3 <- drop_na(d3, RML, THREAT4, rSUBSIST2)
  # model of interest
  model <- lmer(trust ~ 1 + RML + THREAT4 + rSUBSIST2 + 
                  (1 + RML + THREAT4 + rSUBSIST2 | group) + 
                  (1 | uniqid) + (1 | iso), data = d3)
  # fix effect size (log odds slope)
  fixef(model)['RML'] <- std
  # power sim function
  out <- powerSim(model, test = simr::fixed("RML", "z"), 
                  nsim = nsim, seed = 2113, progress = FALSE)
  return(out)
}

# run power analysis #8
runPowerAnalysis8 <- function(d4, std = 0.25, nsim) {
  # standardise outcome
  d4$justify <- as.numeric(scale(d4$justify))
  # standardise RML predictor
  d4$RML <- as.numeric(scale(d4$RML))
  # listwise deletion for predictors
  d4 <- drop_na(d4, RML, THREAT4, rSUBSIST2)
  # model of interest
  model <- lmer(justify ~ 1 + RML + THREAT4 + rSUBSIST2 + 
                  (1 + RML + THREAT4 + rSUBSIST2 | item) + 
                  (1 | uniqid) + (1 | iso),
                data = d4)
  # fix effect size (log odds slope)
  fixef(model)['RML'] <- std
  # power sim function
  out <- powerSim(model, test = simr::fixed("RML", "z"), 
                  nsim = nsim, seed = 2113, progress = FALSE)
  return(out)
}

# plot wvs results - charitable organisations
plotWVSResults5 <- function(d, model, file) {
  # remove missing data and summarise
  d <-
    d %>% 
    filter(A080_01 %in% c(0, 1)) %>% 
    group_by(iso2) %>% 
    summarise(prop = mean(A080_01), RML = mean(RML), RML_SE = mean(RML_SE))
  # plot
  out <-
    conditional_effects(model)$RML %>%
    ggplot() +
    geom_errorbarh(data = d, aes(y = prop, xmin = RML - RML_SE, xmax = RML + RML_SE), alpha = 0.3) +
    geom_point(data = d, aes(x = RML, y = prop), colour = "darkgrey") +
    geom_text_repel(data = d, aes(x = RML, y = prop, label = iso2), size = 2) +
    geom_ribbon(aes(x = RML, ymax = upper__, ymin = lower__), fill = "grey", alpha = 0.5) +
    geom_line(aes(x = RML, y = estimate__)) +
    scale_y_continuous(name = "Probability of belonging to humanitarian\nor charitable organisation", 
                       limits = c(0, 1)) +
    scale_x_continuous(name = "Relational mobility", limits = c(-0.47, 0.47)) +
    theme_classic()
  # save
  ggsave(out, filename = file, height = 3.5, width = 5)
  return(out)
}

# plot wvs results - most people can be trusted
plotWVSResults6 <- function(d, model, file) {
  d <- 
    d %>%
    # listwise deletion of non-response
    filter(A165 %in% 1:2) %>%
    # 0 = can't be too careful, 1 = most people can be trusted
    mutate(A165 = ifelse(A165 == 2, 0, A165)) %>%
    group_by(iso2) %>% 
    summarise(prop = mean(A165), RML = mean(RML), RML_SE = mean(RML_SE))
  # plot
  out <-
    conditional_effects(model)$RML %>%
    ggplot() +
    geom_errorbarh(data = d, aes(y = prop, xmin = RML - RML_SE, xmax = RML + RML_SE), alpha = 0.3) +
    geom_point(data = d, aes(x = RML, y = prop), colour = "darkgrey") +
    geom_text_repel(data = d, aes(x = RML, y = prop, label = iso2), size = 2) +
    geom_ribbon(aes(x = RML, ymax = upper__, ymin = lower__), fill = "grey", alpha = 0.5) +
    geom_line(aes(x = RML, y = estimate__)) +
    scale_y_continuous(name = "Probability of agreeing that\n'most people can be trusted'", 
                       limits = c(0, 1)) +
    scale_x_continuous(name = "Relational mobility", limits = c(-0.47, 0.47)) +
    theme_classic()
  # save
  ggsave(out, filename = file, height = 3.5, width = 5)
  return(out)
}

# plot wvs results - trust specific groups
plotWVSResults7 <- function(d, model, file) {
  # get conditional effects for individual groups
  groups <- data.frame(group = unique(d$group))
  rownames(groups) <- unique(d$group)
  condGroups <- conditional_effects(model, conditions = groups, 
                                    re_formula = . ~ 1 + (1 + RML | group), categorical = FALSE)
  # group plotting function
  plotFun <- function(var) {
    d %>%
      filter(group == var) %>%
      group_by(iso) %>%
      summarise(y = mean(trust), y_se = sd(trust) / sqrt(n()),
                x = mean(RML), x_se = mean(RML_SE)) %>%
      ggplot() +
      geom_errorbarh(aes(y = y, xmin = x - x_se, xmax = x + x_se), alpha = 0.3) +
      geom_errorbar(aes(x = x, ymin = y - y_se, ymax = y + y_se), alpha = 0.3) +
      geom_point(aes(x = x, y = y), size = 0.5) +
      geom_text_repel(aes(x = x, y = y, label = iso), size = 2) +
      geom_ribbon(data = condGroups$RML[condGroups$RML$group == var,], 
                  aes(x = RML, ymax = upper__, ymin = lower__), fill = "grey", alpha = 0.5) +
      geom_line(data = condGroups$RML[condGroups$RML$group == var,],
                aes(x = RML, y = estimate__)) +
      scale_y_continuous(name = var, limits = c(0.8, 4.2)) +
      scale_x_continuous(name = "Relational mobility", limits = c(-0.47, 0.47)) +
      theme_classic()
  }
  # group plots
  pA <- plotFun(groups$group[1])
  pB <- plotFun(groups$group[2])
  pC <- plotFun(groups$group[3])
  pD <- plotFun(groups$group[4])
  pE <- plotFun(groups$group[5])
  pF <- plotFun(groups$group[6])
  # put together
  out <- plot_grid(pA, pB, pC, pD, pE, pF, labels = letters[1:6], nrow = 2)
  ggsave(out, filename = file, height = 6.5, width = 8)
  return(out)
}

# plot wvs results - trust specific groups - quadratic effects
plotWVSResults7Quad <- function(d, post, file) {
  # get cumulative probabilities given some predictor value
  getCumProb <- function(i, x, re) {
    brms::inv_logit_scaled(
      post[,,paste0("b_Intercept[", i, "]")] - 
        (post[,,paste0("r_group[", re, ",Intercept]")] +
           (post[,,"bsp_meRMLRML_SEgrEQiso"] + post[,,paste0("r_group[", re, ",RML]")])*x +
           (post[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post[,,paste0("r_group[", re, ",IRMLE2]")])*(x^2))
    )
  }
  # get posterior average given some predictor value
  getAverage <- function(x, re) {
    tibble(
      prob1 = as.vector(getCumProb(1, x, re)),
      prob2 = as.vector(getCumProb(2, x, re)),
      prob3 = as.vector(getCumProb(3, x, re)),
      prob4 = rep(1, times = 4000)
    ) %>%
      mutate(
        prob4 = prob4 - prob3,
        prob3 = prob3 - prob2,
        prob2 = prob2 - prob1,
        avg = prob1*1 + prob2*2 + prob3*3 + prob4*4
      ) %>%
      pull(avg)
  }
  # plot for given random effect
  plotAvg <- function(re, ylab) {
    pred <-
      tibble(x = seq(min(d$RML), max(d$RML), length.out = 1000), re = re) %>%
      mutate(
        y = map2(x, re, getAverage),
        med = map(y, median),
        low = map(y, function(x) quantile(x, 0.025)),
        upp = map(y, function(x) quantile(x, 0.975))
      ) %>%
      select(!y) %>%
      unnest(c(med, low, upp))
    d %>%
      filter(group == str_replace_all(re, "\\.", " ")) %>%
      group_by(iso) %>%
      summarise(y = mean(trust), y_se = sd(trust) / sqrt(n()),
                x = mean(RML), x_se = mean(RML_SE)) %>%
      ggplot() +
      geom_errorbarh(aes(y = y, xmin = x - x_se, xmax = x + x_se), alpha = 0.3) +
      geom_errorbar(aes(x = x, ymin = y - y_se, ymax = y + y_se), alpha = 0.3) +
      geom_point(aes(x = x, y = y), size = 0.5) +
      geom_text_repel(aes(x = x, y = y, label = iso), size = 2) +
      geom_ribbon(data = pred, aes(x = x, ymin = low, ymax = upp), fill = "grey", alpha = 0.5) +
      geom_line(data = pred, aes(x = x, y = med)) +
      scale_y_continuous(name = ylab, limits = c(1, 4)) +
      scale_x_continuous(name = "Relational mobility", limits = c(-0.47, 0.47)) +
      theme_classic()
  }
  # group plots
  pA <- plotAvg(re = "Your.family",                        ylab = "Your family")
  pB <- plotAvg(re = "Your.neighbourhood",                 ylab = "Your neighbourhood")
  pC <- plotAvg(re = "People.you.know.personally",         ylab = "People you know personally")
  pD <- plotAvg(re = "People.you.meet.for.the.first.time", ylab = "People you meet for the first time")
  pE <- plotAvg(re = "People.of.another.religion",         ylab = "People of another religion")
  pF <- plotAvg(re = "People.of.another.nationality",      ylab = "People of another nationality")
  # put together
  out <- plot_grid(pA, pB, pC, pD, pE, pF, labels = letters[1:6], nrow = 2)
  ggsave(out, filename = file, height = 6.5, width = 8)
  return(out)
}

# plot wvs results - moral justifiability
plotWVSResults8 <- function(d, model, file) {
  # get conditional effects for individual items
  items <- data.frame(item = unique(d$item))
  rownames(items) <- unique(d$item)
  condItems <- conditional_effects(model, conditions = items, 
                                    re_formula = . ~ 1 + (1 + RML | item), categorical = FALSE)
  # item plotting function
  plotFun <- function(var, ylab) {
    d %>%
      filter(item == var) %>%
      group_by(iso) %>%
      summarise(y = mean(justify), y_se = sd(justify) / sqrt(n()),
                x = mean(RML), x_se = mean(RML_SE)) %>%
      ggplot() +
      geom_errorbarh(aes(y = y, xmin = x - x_se, xmax = x + x_se), alpha = 0.3) +
      geom_errorbar(aes(x = x, ymin = y - y_se, ymax = y + y_se), alpha = 0.3) +
      geom_point(aes(x = x, y = y), size = 0.5) +
      geom_text_repel(aes(x = x, y = y, label = iso), size = 2) +
      geom_ribbon(data = condItems$RML[condItems$RML$item == var,], 
                  aes(x = RML, ymax = upper__, ymin = lower__), fill = "grey", alpha = 0.5) +
      geom_line(data = condItems$RML[condItems$RML$item == var,],
                aes(x = RML, y = estimate__)) +
      scale_y_continuous(name = ylab, limits = c(0.8, 10.2), breaks = c(2,4,6,8,10)) +
      scale_x_continuous(name = "Relational mobility", limits = c(-0.47, 0.47)) +
      theme_classic()
  }
  # group plots
  pA <- plotFun(items$item[1], "Claiming government benefits\nto which you are not entitled")
  pB <- plotFun(items$item[2], "Avoiding a fare\nonpublic transport")
  pC <- plotFun(items$item[3], " \nCheating on taxes")
  pD <- plotFun(items$item[4], " \nSomeone accepting a bribe")
  # put together
  out <- plot_grid(pA, pB, pC, pD, labels = letters[1:4], nrow = 2)
  ggsave(out, filename = file, height = 6, width = 6)
  return(out)
}

# plot wvs results - moral justifiability - quadratic effects
plotWVSResults8Quad <- function(d, post, file) {
  # get cumulative probabilities given some predictor value
  getCumProb <- function(i, x, re) {
    brms::inv_logit_scaled(
      post[,,paste0("b_Intercept[", i, "]")] - 
        (post[,,paste0("r_item[", re, ",Intercept]")] +
           (post[,,"bsp_meRMLRML_SEgrEQiso"] + post[,,paste0("r_item[", re, ",RML]")])*x +
           (post[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post[,,paste0("r_item[", re, ",IRMLE2]")])*(x^2))
    )
  }
  # get posterior average given some predictor value
  getAverage <- function(x, re) {
    tibble(
      prob01 = as.vector(getCumProb(1, x, re)),
      prob02 = as.vector(getCumProb(2, x, re)),
      prob03 = as.vector(getCumProb(3, x, re)),
      prob04 = as.vector(getCumProb(4, x, re)),
      prob05 = as.vector(getCumProb(5, x, re)),
      prob06 = as.vector(getCumProb(6, x, re)),
      prob07 = as.vector(getCumProb(7, x, re)),
      prob08 = as.vector(getCumProb(8, x, re)),
      prob09 = as.vector(getCumProb(9, x, re)),
      prob10 = rep(1, times = 4000)
    ) %>%
      mutate(
        prob10 = prob10 - prob09,
        prob09 = prob09 - prob08,
        prob08 = prob08 - prob07,
        prob07 = prob07 - prob06,
        prob06 = prob06 - prob05,
        prob05 = prob05 - prob04,
        prob04 = prob04 - prob03,
        prob03 = prob03 - prob02,
        prob02 = prob02 - prob01,
        avg = prob01*1 + prob02*2 + prob03*3 + prob04*4 + prob05*5 +
          prob06*6 + prob07*7 + prob08*8 + prob09*9 + prob10*10
      ) %>%
      pull(avg)
  }
  # plot for given random effect
  plotAvg <- function(re, ylab) {
    pred <-
      tibble(x = seq(min(d$RML), max(d$RML), length.out = 1000), re = re) %>%
      mutate(
        y = map2(x, re, getAverage),
        med = map(y, median),
        low = map(y, function(x) quantile(x, 0.025)),
        upp = map(y, function(x) quantile(x, 0.975))
      ) %>%
      select(!y) %>%
      unnest(c(med, low, upp))
    d %>%
      filter(item == str_replace_all(re, "\\.", " ")) %>%
      group_by(iso) %>%
      summarise(y = mean(justify), y_se = sd(justify) / sqrt(n()),
                x = mean(RML), x_se = mean(RML_SE)) %>%
      ggplot() +
      geom_errorbarh(aes(y = y, xmin = x - x_se, xmax = x + x_se), alpha = 0.3) +
      geom_errorbar(aes(x = x, ymin = y - y_se, ymax = y + y_se), alpha = 0.3) +
      geom_point(aes(x = x, y = y), size = 0.5) +
      geom_text_repel(aes(x = x, y = y, label = iso), size = 2) +
      geom_ribbon(data = pred, aes(x = x, ymin = low, ymax = upp), fill = "grey", alpha = 0.5) +
      geom_line(data = pred, aes(x = x, y = med)) +
      scale_y_continuous(name = ylab, limits = c(1, 10), breaks = c(2,4,6,8,10)) +
      scale_x_continuous(name = "Relational mobility", limits = c(-0.47, 0.47)) +
      theme_classic()
  }
  # group plots
  pA <- plotAvg(re = "Claiming.government.benefits.to.which.you.are.not.entitled",
                ylab = "Claiming government benefits\nto which you are not entitled")
  pB <- plotAvg(re = "Avoiding.a.fare.on.public.transport",
                ylab = "Avoiding a fare\nonpublic transport")
  pC <- plotAvg(re = "Cheating.on.taxes", ylab = " \nCheating on taxes")
  pD <- plotAvg(re = "Someone.accepting.a.bribe", ylab = " \nSomeone accepting a bribe")
  # put together
  out <- plot_grid(pA, pB, pC, pD, labels = letters[1:4], nrow = 2)
  ggsave(out, filename = file, height = 6, width = 6)
  return(out)
}

# create raw data table - study 1
createTableRawStudy1 <- function(d1) {
  d1 %>%
    pivot_wider(names_from = "item", values_from = "coop") %>%
    rename(Country = countryRM) %>%
    mutate(Country = ifelse(Country == "United Kingdom", "UK", Country)) %>%
    group_by(Country) %>%
    summarise(
      `Positive reciprocity` = mean(posrecip, na.rm = TRUE),
      `Trust` = mean(trust, na.rm = TRUE),
      `Altruism` = mean(altruism, na.rm = TRUE),
      `Relational mobility` = mean(RML, na.rm = TRUE),
      SE = mean(RML_SE, na.rm = TRUE)
    ) %>%
    mutate_if(is.numeric, function(x) format(round(x, 2), nsmall = 2))
}

# create raw data table - study 2
createTableRawStudy2 <- function(d2) {
  d2 %>%
    # non-response
    mutate(
      # NAs for non-response 
      A080_01 = ifelse(A080_01 %in% 0:1, A080_01, NA),
      A165 = ifelse(A165 %in% 1:2, A165, NA),
      D001_B = ifelse(D001_B %in% 1:4, D001_B, NA),
      G007_18_B = ifelse(G007_18_B %in% 1:4, G007_18_B, NA),
      G007_33_B = ifelse(G007_33_B %in% 1:4, G007_33_B, NA),
      G007_34_B = ifelse(G007_34_B %in% 1:4, G007_34_B, NA),
      G007_35_B = ifelse(G007_35_B %in% 1:4, G007_35_B, NA),
      G007_36_B = ifelse(G007_36_B %in% 1:4, G007_36_B, NA),
      F114A = ifelse(F114A %in% 1:10, F114A, NA),
      F115 = ifelse(F115 %in% 1:10, F115, NA),
      F116 = ifelse(F116 %in% 1:10, F116, NA),
      F117 = ifelse(F117 %in% 1:10, F117, NA),
      # A165: 0 = can't be too careful, 1 = most people can be trusted
      A165 = ifelse(A165 == 2, 0, A165),
      # trust: 1-4
      D001_B = 5 - D001_B,
      G007_18_B = 5 - G007_18_B, 
      G007_33_B = 5 - G007_33_B, 
      G007_34_B = 5 - G007_34_B, 
      G007_35_B = 5 - G007_35_B, 
      G007_36_B = 5 - G007_36_B,
      # justify: 1-10
      F114A = 11 - F114A,
      F115 = 11 - F115,
      F116 = 11 - F116,
      F117 = 11 - F117
    ) %>%
    rename(Country = countryRM) %>%
    mutate(Country = ifelse(Country == "United Kingdom", "UK", Country)) %>%
    group_by(Country) %>%
    summarise(
      CharOrg   = mean(A080_01, na.rm = TRUE),
      Trust     = mean(A165, na.rm = TRUE),
      TruFam    = mean(D001_B, na.rm = TRUE),
      TruNeigh  = mean(G007_18_B, na.rm = TRUE),
      TruKnow   = mean(G007_33_B, na.rm = TRUE),
      TruMeet   = mean(G007_34_B, na.rm = TRUE),
      TruRel    = mean(G007_35_B, na.rm = TRUE),
      TruNat    = mean(G007_36_B, na.rm = TRUE),
      JusGovBen = mean(F114A, na.rm = TRUE),
      JusFare   = mean(F115, na.rm = TRUE),
      JusTax    = mean(F116, na.rm = TRUE),
      JusBribe  = mean(F117, na.rm = TRUE),
      RelMob    = mean(RML, na.rm = TRUE),
      SE        = mean(RML_SE, na.rm = TRUE)
    ) %>%
    mutate_if(is.numeric, function(x) ifelse(is.na(x), "", format(round(x, 2), nsmall = 2)))
}

# create table of quadratic results
createTableQuadratic <- function(post1.3, post5.3, post6.3, post7.3, post8.3) {
  # Round function
  Round <- function(x) {
    x <- format(round(x, 2), nsmall = 2)
    out <- ifelse(x < 0, as.character(x), paste0(" ", as.character(x)))
    return(out)
  }
  # extract linear slopes
  slopeLin1.3a <- post1.3[,,"bsp_meRMLRML_SEgrEQiso"]
  slopeLin1.3b <- post1.3[,,"bsp_meRMLRML_SEgrEQiso"] + post1.3[,,"r_item[altruism,RML]"]
  slopeLin1.3c <- post1.3[,,"bsp_meRMLRML_SEgrEQiso"] + post1.3[,,"r_item[posrecip,RML]"]
  slopeLin1.3d <- post1.3[,,"bsp_meRMLRML_SEgrEQiso"] + post1.3[,,"r_item[trust,RML]"]
  slopeLin5.3a <- post5.3[,,"bsp_meRMLRML_SEgrEQiso"]
  slopeLin6.3a <- post6.3[,,"bsp_meRMLRML_SEgrEQiso"]
  slopeLin7.3a <- post7.3[,,"bsp_meRMLRML_SEgrEQiso"]
  slopeLin7.3b <- post7.3[,,"bsp_meRMLRML_SEgrEQiso"] + post7.3[,,"r_group[People.of.another.nationality,RML]"]
  slopeLin7.3c <- post7.3[,,"bsp_meRMLRML_SEgrEQiso"] + post7.3[,,"r_group[People.of.another.religion,RML]"]
  slopeLin7.3d <- post7.3[,,"bsp_meRMLRML_SEgrEQiso"] + post7.3[,,"r_group[People.you.know.personally,RML]"]
  slopeLin7.3e <- post7.3[,,"bsp_meRMLRML_SEgrEQiso"] + post7.3[,,"r_group[People.you.meet.for.the.first.time,RML]"]
  slopeLin7.3f <- post7.3[,,"bsp_meRMLRML_SEgrEQiso"] + post7.3[,,"r_group[Your.family,RML]"]
  slopeLin7.3g <- post7.3[,,"bsp_meRMLRML_SEgrEQiso"] + post7.3[,,"r_group[Your.neighbourhood,RML]"]
  slopeLin8.3a <- post8.3[,,"bsp_meRMLRML_SEgrEQiso"]
  slopeLin8.3b <- post8.3[,,"bsp_meRMLRML_SEgrEQiso"] + post8.3[,,"r_item[Avoiding.a.fare.on.public.transport,RML]"]
  slopeLin8.3c <- post8.3[,,"bsp_meRMLRML_SEgrEQiso"] + post8.3[,,"r_item[Cheating.on.taxes,RML]"]
  slopeLin8.3d <- post8.3[,,"bsp_meRMLRML_SEgrEQiso"] + post8.3[,,"r_item[Claiming.government.benefits.to.which.you.are.not.entitled,RML]"]
  slopeLin8.3e <- post8.3[,,"bsp_meRMLRML_SEgrEQiso"] + post8.3[,,"r_item[Someone.accepting.a.bribe,RML]"]
  # extract quadratic slopes
  slopeQuad1.3a <- post1.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"]
  slopeQuad1.3b <- post1.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post1.3[,,"r_item[altruism,IRMLE2]"]
  slopeQuad1.3c <- post1.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post1.3[,,"r_item[posrecip,IRMLE2]"]
  slopeQuad1.3d <- post1.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post1.3[,,"r_item[trust,IRMLE2]"]
  slopeQuad5.3a <- post5.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"]
  slopeQuad6.3a <- post6.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"]
  slopeQuad7.3a <- post7.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"]
  slopeQuad7.3b <- post7.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post7.3[,,"r_group[People.of.another.nationality,IRMLE2]"]
  slopeQuad7.3c <- post7.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post7.3[,,"r_group[People.of.another.religion,IRMLE2]"]
  slopeQuad7.3d <- post7.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post7.3[,,"r_group[People.you.know.personally,IRMLE2]"]
  slopeQuad7.3e <- post7.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post7.3[,,"r_group[People.you.meet.for.the.first.time,IRMLE2]"]
  slopeQuad7.3f <- post7.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post7.3[,,"r_group[Your.family,IRMLE2]"]
  slopeQuad7.3g <- post7.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post7.3[,,"r_group[Your.neighbourhood,IRMLE2]"]
  slopeQuad8.3a <- post8.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"]
  slopeQuad8.3b <- post8.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post8.3[,,"r_item[Avoiding.a.fare.on.public.transport,IRMLE2]"]
  slopeQuad8.3c <- post8.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post8.3[,,"r_item[Cheating.on.taxes,IRMLE2]"]
  slopeQuad8.3d <- post8.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post8.3[,,"r_item[Claiming.government.benefits.to.which.you.are.not.entitled,IRMLE2]"]
  slopeQuad8.3e <- post8.3[,,"bsp_ImeRMLRML_SEgrEQisoE2"] + post8.3[,,"r_item[Someone.accepting.a.bribe,IRMLE2]"]
  # create table
  tibble(
    Outcome = c("GPS Prosociality", "", "", "", "WVS Charitable", "WVS Trust", 
                "WVS Trust Groups", "", "", "", "", "", "", "WVS Justify", "", "", "", ""),
    Parameter = c("Population-level", "RE: Altruism", "RE: Positive reciprocity",
                  "RE: Trust", "Population-level", "Population-level", 
                  "Population-level", "RE: Another nationality", "RE: Another religion",
                  "RE: Know personally", "RE: Meet first time", "RE: Family",
                  "RE: Neighbourhood", "Population-level", "RE: Public transport",
                  "RE: Cheat taxes", "RE: Gov benefits", "RE: Accept bribe"),
    `Linear slope` = c(
      paste0("b = ", Round(median(slopeLin1.3a)), ", 95% CI [", Round(quantile(slopeLin1.3a, 0.025)), ", ", Round(quantile(slopeLin1.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin1.3b)), ", 95% CI [", Round(quantile(slopeLin1.3b, 0.025)), ", ", Round(quantile(slopeLin1.3b, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin1.3c)), ", 95% CI [", Round(quantile(slopeLin1.3c, 0.025)), ", ", Round(quantile(slopeLin1.3c, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin1.3d)), ", 95% CI [", Round(quantile(slopeLin1.3d, 0.025)), ", ", Round(quantile(slopeLin1.3d, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin5.3a)), ", 95% CI [", Round(quantile(slopeLin5.3a, 0.025)), ", ", Round(quantile(slopeLin5.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin6.3a)), ", 95% CI [", Round(quantile(slopeLin6.3a, 0.025)), ", ", Round(quantile(slopeLin6.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin7.3a)), ", 95% CI [", Round(quantile(slopeLin7.3a, 0.025)), ", ", Round(quantile(slopeLin7.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin7.3b)), ", 95% CI [", Round(quantile(slopeLin7.3b, 0.025)), ", ", Round(quantile(slopeLin7.3b, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin7.3c)), ", 95% CI [", Round(quantile(slopeLin7.3c, 0.025)), ", ", Round(quantile(slopeLin7.3c, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin7.3d)), ", 95% CI [", Round(quantile(slopeLin7.3d, 0.025)), ", ", Round(quantile(slopeLin7.3d, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin7.3e)), ", 95% CI [", Round(quantile(slopeLin7.3e, 0.025)), ", ", Round(quantile(slopeLin7.3e, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin7.3f)), ", 95% CI [", Round(quantile(slopeLin7.3f, 0.025)), ", ", Round(quantile(slopeLin7.3f, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin7.3g)), ", 95% CI [", Round(quantile(slopeLin7.3g, 0.025)), ", ", Round(quantile(slopeLin7.3g, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin8.3a)), ", 95% CI [", Round(quantile(slopeLin8.3a, 0.025)), ", ", Round(quantile(slopeLin8.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin8.3b)), ", 95% CI [", Round(quantile(slopeLin8.3b, 0.025)), ", ", Round(quantile(slopeLin8.3b, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin8.3c)), ", 95% CI [", Round(quantile(slopeLin8.3c, 0.025)), ", ", Round(quantile(slopeLin8.3c, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin8.3d)), ", 95% CI [", Round(quantile(slopeLin8.3d, 0.025)), ", ", Round(quantile(slopeLin8.3d, 0.975)), "]"),
      paste0("b = ", Round(median(slopeLin8.3e)), ", 95% CI [", Round(quantile(slopeLin8.3e, 0.025)), ", ", Round(quantile(slopeLin8.3e, 0.975)), "]")
    ),
    `Quadratic slope` = c(
      paste0("b = ", Round(median(slopeQuad1.3a)), ", 95% CI [", Round(quantile(slopeQuad1.3a, 0.025)), ", ", Round(quantile(slopeQuad1.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad1.3b)), ", 95% CI [", Round(quantile(slopeQuad1.3b, 0.025)), ", ", Round(quantile(slopeQuad1.3b, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad1.3c)), ", 95% CI [", Round(quantile(slopeQuad1.3c, 0.025)), ", ", Round(quantile(slopeQuad1.3c, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad1.3d)), ", 95% CI [", Round(quantile(slopeQuad1.3d, 0.025)), ", ", Round(quantile(slopeQuad1.3d, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad5.3a)), ", 95% CI [", Round(quantile(slopeQuad5.3a, 0.025)), ", ", Round(quantile(slopeQuad5.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad6.3a)), ", 95% CI [", Round(quantile(slopeQuad6.3a, 0.025)), ", ", Round(quantile(slopeQuad6.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad7.3a)), ", 95% CI [", Round(quantile(slopeQuad7.3a, 0.025)), ", ", Round(quantile(slopeQuad7.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad7.3b)), ", 95% CI [", Round(quantile(slopeQuad7.3b, 0.025)), ", ", Round(quantile(slopeQuad7.3b, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad7.3c)), ", 95% CI [", Round(quantile(slopeQuad7.3c, 0.025)), ", ", Round(quantile(slopeQuad7.3c, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad7.3d)), ", 95% CI [", Round(quantile(slopeQuad7.3d, 0.025)), ", ", Round(quantile(slopeQuad7.3d, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad7.3e)), ", 95% CI [", Round(quantile(slopeQuad7.3e, 0.025)), ", ", Round(quantile(slopeQuad7.3e, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad7.3f)), ", 95% CI [", Round(quantile(slopeQuad7.3f, 0.025)), ", ", Round(quantile(slopeQuad7.3f, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad7.3g)), ", 95% CI [", Round(quantile(slopeQuad7.3g, 0.025)), ", ", Round(quantile(slopeQuad7.3g, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad8.3a)), ", 95% CI [", Round(quantile(slopeQuad8.3a, 0.025)), ", ", Round(quantile(slopeQuad8.3a, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad8.3b)), ", 95% CI [", Round(quantile(slopeQuad8.3b, 0.025)), ", ", Round(quantile(slopeQuad8.3b, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad8.3c)), ", 95% CI [", Round(quantile(slopeQuad8.3c, 0.025)), ", ", Round(quantile(slopeQuad8.3c, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad8.3d)), ", 95% CI [", Round(quantile(slopeQuad8.3d, 0.025)), ", ", Round(quantile(slopeQuad8.3d, 0.975)), "]"),
      paste0("b = ", Round(median(slopeQuad8.3e)), ", 95% CI [", Round(quantile(slopeQuad8.3e, 0.025)), ", ", Round(quantile(slopeQuad8.3e, 0.975)), "]")
    )
  )
}

# create table of power analysis results
createTablePower <- function(power1, power5, power6, power7, power8) {
  Round <- function(x) format(round(x, 2), nsmall = 2)
  tibble(
    Outcome = c("GPS Prosociality", "WVS Charitable", "WVS Trust", 
                "WVS Trust Groups", "WVS Justify"),
    Model = c("Multilevel regression", "Multilevel logistic regression",
              "Multilevel logistic regression", "Multilevel regression",
              "Multilevel regression"),
    Slope = c(
      str_sub(power1$description[2], -4),
      str_sub(power5$description[2], -4),
      str_sub(power6$description[2], -4),
      str_sub(power7$description[2], -4),
      str_sub(power8$description[2], -4)
    ),
    `Effect size` = c("Medium", "Small", "Small", "Small", "Medium"),
    Power = c(
      Round(summary(power1)$mean),
      Round(summary(power5)$mean),
      Round(summary(power6)$mean),
      Round(summary(power7)$mean),
      Round(summary(power8)$mean)
    ),
    `Lower 95%` = c(
      Round(summary(power1)$lower),
      Round(summary(power5)$lower),
      Round(summary(power6)$lower),
      Round(summary(power7)$lower),
      Round(summary(power8)$lower)
    ),
    `Upper 95%` = c(
      Round(summary(power1)$upper),
      Round(summary(power5)$upper),
      Round(summary(power6)$upper),
      Round(summary(power7)$upper),
      Round(summary(power8)$upper)
    )
  )
}