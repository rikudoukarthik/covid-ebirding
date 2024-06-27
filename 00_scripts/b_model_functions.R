# expand by species -----------------------------------------------------------------

expandbyspecies = function(data, species) {
  
  require(tidyverse)
  
  # considers only complete lists
  
  occinfo <- data %>% 
    group_by(GROUP.ID, COMMON.NAME) %>% 
    reframe(REPORT = 1) %>% 
    ungroup()
  
  checklistinfo = data %>% 
    filter(ALL.SPECIES.REPORTED == 1) %>%
    distinct(GROUP.ID, M.YEAR, MONTH, DAY.M, URBAN, CELL.ID, SUBCELL.ID, NO.SP) %>%
    group_by(GROUP.ID) %>% 
    slice(1) %>% 
    ungroup()
  
  # expand data frame to include the bird species in every list
  expanded = checklistinfo %>% 
    # getting species of interest
    group_by(GROUP.ID) %>% 
    reframe(COMMON.NAME = species) %>% 
    left_join(occinfo, by = c("GROUP.ID", "COMMON.NAME")) %>%
    left_join(checklistinfo, by = "GROUP.ID") %>% 
    # for species not reported in lists, filling in NAs in COMMON.NAME and REPORT
    mutate(REPORT = replace_na(REPORT, 0))

  return(expanded)
  
}


# single species model --------------------------------------------------------------

singlespeciesmodel = function(data, species, specieslist, iter = NULL) {
  
  require(tidyverse)
  require(lme4)
  require(merTools)
  require(glue)


  # getting median list length for prediction later
  median_length <- data %>% 
    distinct(M.YEAR, MONTH, GROUP.ID, NO.SP) %>% 
    group_by(MONTH) %>% 
    reframe(NO.SP.MED = floor(median(NO.SP)))
  
  message("Completed preparations for modelling. Now starting modelling.")

  
  # "constraining data to expected range" (SoIB 2023)
  data_filt = data %>%
    filter(COMMON.NAME == species) %>%
    # using only CELL.ID-MONTH (space-time) combos in which species occurs
    # using 100 km to maximise denominator sample size in low-birded grids
    distinct(GRID.G3, MONTH) %>% 
    # joining median list length
    left_join(median_length, by = "MONTH") %>% 
    # join rest of dataset back
    left_join(data)

  # expand dataframe to include absences as well
  data_exp = expandbyspecies(data_filt, species) %>% 
    # join species categories
    left_join(specieslist %>% distinct(COMMON.NAME, SP.CATEGORY),
              by = "COMMON.NAME")
  
  
  # the model ---------------------------------------------------------------

  # for some species in KL and MH, using cloglog link is resulting in "PIRLS step-halvings
  # failed to reduce deviance in pwrssUpdate"
  
  if ((state_name == "Kerala" & species %in% fail_spec_KL) |
      (state_name == "Maharashtra" & species %in% fail_spec_MH)){
    
    model_spec <- glmer(REPORT ~ M.YEAR + MONTH + MONTH:log(NO.SP) + MONTH:M.YEAR +
                          (1|CELL.ID),
                        data = data_exp, family = binomial,
                        nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
    
  } else {
    
    model_spec <- glmer(REPORT ~ M.YEAR + MONTH + MONTH:log(NO.SP) + MONTH:M.YEAR +
                          (1|CELL.ID),
                        data = data_exp, family = binomial(link = "cloglog"),
                        nAGQ = 0, control = glmerControl(optimizer = "bobyqa"))
    
  }

  # predicting from model ---------------------------------------------------
  
  # prepare a new data file to predict
  birds_pred <- data_exp %>% 
    # selecting random CELL.ID because predictInterval needs input even if which == "fixed"
    mutate(CELL.ID = sample(unique(CELL.ID), 1)) %>% 
    distinct(MONTH, M.YEAR, CELL.ID, COMMON.NAME, SP.CATEGORY) %>% 
    # joining median list length
    left_join(median_length, by = "MONTH") %>% 
    rename(NO.SP = NO.SP.MED)
  
  
  pred = predictInterval(model_spec, newdata = birds_pred, which = "fixed",
                         level = 0.48, type = "linear.prediction")
  birds_pred$PRED.LINK = pred$fit
  birds_pred$SE.LINK = pred$fit - pred$lwr
  
  
  birds_pred = birds_pred %>%
    filter(!is.na(PRED.LINK) & !is.na(SE.LINK)) %>%
    group_by(COMMON.NAME, SP.CATEGORY, M.YEAR) %>% 
    reframe(PRED.LINK = mean(PRED.LINK), 
            # just averaging across months, so no need for propagation
            SE.LINK = mean(SE.LINK)) 
  
  # add iteration number as column so we can save all in one object
  if (!is.null(iter)) {
    birds_pred <- birds_pred %>% mutate(ITERATION = iter)
  }
  
  return(birds_pred)
  
}


# inverse link functions ------------------------------------------------------------

# choose inverse link function based on which model
inverse_link <- function(data, state, species) {
  
  if ((state == "Kerala" & species %in% fail_spec_KL) |
      (state == "Maharashtra" & species %in% fail_spec_MH)) {
    exp(data)/(1 + exp(data)) # inverse logit
  } else {
    clogloglink(data, inverse = TRUE) # inverse cloglog
  }
  
}


# bootstrapped error prop. for division ---------------------------------------------

# bootstrapped error calculation for division operation (% change)
simerrordiv = function(x1, x2, se1, se2, state, species)
{
  # takes untransformed (link) mean and SE values, and generates normal dist. from 1000 sims,
  # then transformed 
  # after the function, lower and upper quantiles are selected as limits of 95% CI
  tp = data.frame(num = inverse_link(rnorm(1000, x1, se1), state, species), 
                  den = inverse_link(rnorm(1000, x2, se2), state, species)) %>%
    reframe(rat = num/den, 
            val = num)
  
  return(tp)
}
# plotting function for bird model results -----------------

gg_b_model <- function(data, type, data_points) {

  require(tidyverse)
  require(patchwork)
  require(ggtext)
  
  if (type == "overall_ribbon") {

    lab_y <- "Change in abundance index<br>(from t~0~ = 2018\u201319)"
    lab_ribbons <- urbrur_palette %>% 
      mutate(HEX.LABEL = case_when(MONTHS.TYPE == "LD" ~ "Apr\u2013May\n(Peak impact)",
                                   MONTHS.TYPE == "NL" ~ "Jun\u2013Mar\n(Rest of year)")) %>% 
      mutate(HEX.LABEL = factor(HEX.LABEL, 
                                levels = c("Jun\u2013Mar\n(Rest of year)", 
                                           "Apr\u2013May\n(Peak impact)")))
    
    model_data <- data %>% 
      mutate(MONTHS.TYPE = factor(MONTHS.TYPE, levels = c("NL", "LD"))) %>%
      rename(PRED = PRED.PERC, SE = SE.PERC, CI.L = CI.L.PERC, CI.U = CI.U.PERC) %>% 
      # convert to + and - values
      mutate(PRED.LABEL = PRED - 100) %>% 
      mutate(PRED.LABEL = case_when(PRED.LABEL > 0 ~ glue("+{PRED.LABEL}%"),
                                    PRED == 100 ~ glue("0%"), 
                                    TRUE ~ glue("{PRED.LABEL}%"))) %>% 
      # joining colours
      left_join(lab_ribbons, by = c("SP.CATEGORY", "MONTHS.TYPE")) 
    
    model_plot <- model_data %>% 
      ggplot(mapping = aes(x = as.numeric(levels(M.YEAR))[M.YEAR])) +
      geom_line(mapping = aes(y = PRED, colour = HEX),
                linewidth = 1, lineend = "round") +
      geom_ribbon(mapping = aes(ymin = CI.L, ymax = CI.U, fill = HEX, colour = HEX),
                  alpha = 0.5, linewidth = 1) +
      geom_hline(yintercept = 100, linewidth = 1) +
      geom_vline(xintercept = c(2019, 2020, 2021), 
                 linetype = "dashed", linewidth = 0.5) +
      # annotate("text",
      #          y = 102, angle = 90, size = 5, 
      #          x = c(2018.92, 2019.92, 2020.92),
      #          label = rep(c("2019\u201320", "2020\u201321", "2021\u201322"), 2),
      #          colour = rep(c("red", "red", "black"), 2)) +
      scale_color_identity(guide = "legend", name = "Months of year", labels = lab_ribbons$HEX.LABEL) +
      scale_fill_identity(guide = "legend", name = "Months of year", labels = lab_ribbons$HEX.LABEL) +
      scale_x_continuous(labels = c("2018\u201319", "2019\u201320", "2020\u201321", "2021\u201322"),
                         limits = c(2018, 2021.003)) +
      labs(x = "Migratory Year", 
           y = lab_y,
           title = unique(data$STATE.NAME)) +
      coord_cartesian(expand = FALSE, ylim = c(74.2, 135)) +
      # guides(fill = guide_legend(byrow = TRUE, reverse = FALSE),
      #        colour = guide_legend(byrow = TRUE, reverse = FALSE)) +
      theme(legend.position = "right",
            axis.line.y = element_line(linewidth = 1, arrow = grid::arrow()),
            axis.title.y = element_markdown(lineheight = 1.2, margin = margin(0, 10, 0, 0),
                                            colour = "black"),
            axis.text.y = element_text(colour = "black"),
            axis.ticks.length.y.left = unit(0.3, "cm"),
            axis.line.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            # legend.key.size = unit(0.5, "cm"),
            # legend.spacing.y = unit(0.3, "cm"),
            plot.margin = margin(10, 20, 10, 20),
            plot.title = element_text(colour = "black", hjust = 0.5, vjust = 1, 
                                      margin = margin(0, 0, 10, 0))) +
      facet_wrap(~ MONTHS.TYPE, nrow = 2)
    
    return(model_plot)
    
  } else if (type == "overall") {
    
    plot_title <- glue("{state_name} state")
    plot_subtitle <- paste0(
      "Predicted reporting frequencies of ",
      n_distinct(data_occ$COMMON.NAME), " species (",
      n_distinct(filter(data_occ, SP.CATEGORY == "U")$COMMON.NAME), " urban, ",
      n_distinct(filter(data_occ, SP.CATEGORY == "R")$COMMON.NAME), " rural)",
      " in three separate monthwise models")
    
    # different y lims for different states
    if (state_name == "Karnataka") {
      plot_ylims <- c(0.13, 0.35)
    } else if (state_name == "Kerala") {
      plot_ylims <- c(0.05, 0.22)
    } else if (state_name == "Maharashtra") {
      plot_ylims <- c(0.175, 0.35)
    } else if (state_name == "Assam") {
      plot_ylims <- c(0.1, 0.5)
    }
    
    
    model_plot <- (ggplot(filter(data, MONTHS.TYPE == "LD"), 
                          aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = plot_ylims) +
                     labs(title = "For the months of April and May",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     geom_point(size = 1.75, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1, width = 0.15, position = position_dodge(0.5)) |
                     ggplot(filter(data, MONTHS.TYPE == "NL"), 
                            aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = plot_ylims) +
                     labs(title = "For other ten months",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     geom_point(size = 1.75, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1, width = 0.15, position = position_dodge(0.5)) |
                     ggplot(filter(data, MONTHS.TYPE == "ALL"), 
                            aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = plot_ylims) +
                     labs(title = "For all twelve months",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     geom_point(size = 1.75, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1, width = 0.15, position = position_dodge(0.5))) +
      plot_layout(guides = "collect") +
      plot_annotation(title = plot_title,
                      subtitle = plot_subtitle) 
    
    return(model_plot)
    
  } else if (type == "prolific") {
    
    plot_title <- "Change in bird species reporting from locations with consistent effort"
    plot_subtitle <- paste0(
      "Predicted reporting frequencies of ",
      n_distinct(data_prolific$COMMON.NAME), " species (",
      n_distinct(filter(data_prolific, SP.CATEGORY == "U")$COMMON.NAME), " urban, ",
      n_distinct(filter(data_prolific, SP.CATEGORY == "R")$COMMON.NAME), " rural)",
      " in three separate monthwise models")
    
    
    model_plot <- (ggplot(filter(data, MONTHS.TYPE == "LD"), 
                          aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = c(0.1, 1.0)) +
                     labs(title = "For the months of April and May",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     # data points
                     geom_point(data = filter(data_points, MONTHS.TYPE == "LD"), 
                                aes(group = PATH.GROUP),
                                size = 3, alpha = 0.25,
                                position = position_dodge(0.25)) + 
                     geom_path(data = filter(data_points, MONTHS.TYPE == "LD"), 
                               aes(x = as.numeric(M.YEAR), y = PRED, group = PATH.GROUP),
                               size = 1, alpha = 0.15, position = position_dodge(0.25)) + 
                     # main points
                     geom_point(size = 3, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1.5, width = 0.2, position = position_dodge(0.5)) |
                     ggplot(filter(data, MONTHS.TYPE == "NL"), 
                            aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = c(0.1, 1.0)) +
                     labs(title = "For other ten months",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     # data points
                     geom_point(data = filter(data_points, MONTHS.TYPE == "NL"), 
                                aes(group = PATH.GROUP),
                                size = 3, alpha = 0.25,
                                position = position_dodge(0.25)) + 
                     geom_path(data = filter(data_points, MONTHS.TYPE == "NL"), 
                               aes(x = as.numeric(M.YEAR), y = PRED, group = PATH.GROUP),
                               size = 1, alpha = 0.15, position = position_dodge(0.25)) + 
                     # main points
                     geom_point(size = 3, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1.5, width = 0.2, position = position_dodge(0.5)) |
                     ggplot(filter(data, MONTHS.TYPE == "ALL"), 
                            aes(M.YEAR, PRED, col = SP.CATEGORY)) +
                     scale_color_manual(values = c("#8F85C1", "#A3383C"),
                                        name = "Species\ncategory",
                                        labels = c("Rural", "Urban")) +
                     scale_y_continuous(limits = c(0.1, 1.0)) +
                     labs(title = "For all twelve months",
                          x = "Migratory year", y = "Predicted reporting frequency") +
                     # data points
                     geom_point(data = filter(data_points, MONTHS.TYPE == "ALL"), 
                                aes(group = PATH.GROUP),
                                size = 3, alpha = 0.25,
                                position = position_dodge(0.25)) + 
                     geom_path(data = filter(data_points, MONTHS.TYPE == "ALL"), 
                               aes(x = as.numeric(M.YEAR), y = PRED, group = PATH.GROUP),
                               size = 1, alpha = 0.15, position = position_dodge(0.25)) + 
                     # main points
                     geom_point(size = 3, position = position_dodge(0.5)) +
                     geom_errorbar(aes(ymin = CI.L, ymax = CI.U), 
                                   size = 1.5, width = 0.2, position = position_dodge(0.5))) +
      plot_layout(guides = "collect") +
      plot_annotation(title = plot_title,
                      subtitle = plot_subtitle) 
    
    return(model_plot)
    
  } else {print("Please select valid analysis type!")}
  
}
