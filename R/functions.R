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
      group_by(country) %>%
      summarise(y = mean(coop), y_se = sd(coop) / sqrt(n()),
                x = mean(RML), x_se = mean(RML_SE)) %>%
      ggplot() +
      geom_errorbarh(aes(y = y, xmin = x - x_se, xmax = x + x_se), alpha = 0.3, colour = col) +
      geom_errorbar(aes(x = x, ymin = y - y_se, ymax = y + y_se), alpha = 0.3, colour = col) +
      geom_point(aes(x = x, y = y), size = 0.5, colour = col, fill = col) +
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
      group_by(countryRM) %>%
      summarise(y = mean(trust), y_se = sd(trust) / sqrt(n()),
                x = mean(RML), x_se = mean(RML_SE)) %>%
      ggplot() +
      geom_errorbarh(aes(y = y, xmin = x - x_se, xmax = x + x_se), alpha = 0.3) +
      geom_errorbar(aes(x = x, ymin = y - y_se, ymax = y + y_se), alpha = 0.3) +
      geom_point(aes(x = x, y = y), size = 0.5) +
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
      group_by(countryRM) %>%
      summarise(y = mean(justify), y_se = sd(justify) / sqrt(n()),
                x = mean(RML), x_se = mean(RML_SE)) %>%
      ggplot() +
      geom_errorbarh(aes(y = y, xmin = x - x_se, xmax = x + x_se), alpha = 0.3) +
      geom_errorbar(aes(x = x, ymin = y - y_se, ymax = y + y_se), alpha = 0.3) +
      geom_point(aes(x = x, y = y), size = 0.5) +
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
