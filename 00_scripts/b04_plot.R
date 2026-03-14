require(ggplot2)

theme_set(theme_classic(base_size = 16))

covid_palette <- c("#1B9E77", "#EF4050", "#E89005", "#9678B6")
# urbrur_palette <- tibble(SP.CATEGORY = c("R", "R", "U", "U"),
#                          MONTHS.TYPE = c("LD", "ALL", "LD", "ALL"),
#                          # HEX = c("#a6611a", "#dfc27d", "#018571", "#80cdc1"),
#                          HEX = c(rep("#d8b365", 2), rep("#5ab4ac", 2))) %>% 
#   mutate(MONTHS.TYPE = factor(MONTHS.TYPE, levels = c("LD", "ALL")))
urbrur_palette <- c("#d8b365", "#5ab4ac") # rur, urb

source("00_scripts/b03_model_functions.R")



# overall -----------------------------------------------------------------

# aggregating all states and months data 

data_base <- tibble(STATE = factor(c("Karnataka", "Kerala", "Maharashtra", "Assam"),
                                   levels = c("Karnataka", "Kerala", "Maharashtra", "Assam"))) %>% 
  group_by(STATE) %>% 
  reframe(MT = c("LD", "ALL")) 


data_all <- map2(data_base$STATE, data_base$MT, ~ {
  
  data_path <- glue("00_outputs/bird_models/{.x}/b04_results_{.y}.RData")
  
  if(file.exists(data_path)) {
    load(data_path)
    data_res_avg %>% mutate(STATE.NAME = .x)
  } 
  
}) %>% 
  list_rbind()


# ribbon fig 

create_bird_graph <- function() {
  
  data_all <- data %>% filter(MONTHS.TYPE = "ALL")
  data_ld <- data %>% filter(MONTHS.TYPE = "LD")
  
  data_all %>% 
    gg_b_model("overall_ribbon") +
    facet_wrap(~ STATE.NAME)
  
  
}


# prolific ----------------------------------------------------------------

load("00_outputs/bird_models/prolific.RData")

prolific_plot <- gg_b_model(prolific_pred, type = "prolific", prolific_points)

ggsave(filename = "02_bird_figs/prolific_model.png", plot = prolific_plot,
       dpi = 300, width = 11, height = 9, units = "in")

