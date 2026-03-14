require(ggplot2)
require(glue)

theme_set(theme_classic(base_size = 16))

covid_palette <- c("#1B9E77", "#EF4050", "#E89005", "#9678B6")
urbrur_palette_df <- tibble(SP.CATEGORY = c("R", "R", "U", "U"),
                         MONTHS.TYPE = c("LD", "ALL", "LD", "ALL"),
                         HEX = c("#a6611a", "#dfc27d", "#018571", "#80cdc1")) %>%
  mutate(MONTHS.TYPE = factor(MONTHS.TYPE, levels = c("LD", "ALL")))
urbrur_palette <- c("#d8b365", "#5ab4ac") # rur, urb

source("00_scripts/b03_model_functions.R")


# overall: ribbon ---------------------------------------------------------

# data_all %>% 
#   gg_b_model("overall_ribbon") +
#   facet_grid(SP.CATEGORY ~ STATE.NAME)


# overall: point-error --------------------------------------------------------

# aggregating all states and months data 

to_iter <- tibble(STATE = factor(spat_iter, levels = spat_iter)) %>% 
  group_by(STATE) %>% 
  reframe(MT = factor(c("LD", "ALL"), levels = c("LD", "ALL"))) 

# load and combine all data
data_all <- map2(to_iter$STATE, to_iter$MT, ~ {
  
  data_path <- glue("00_outputs/bird_models/{.x}/b04_results_{.y}.RData")
  
  if(file.exists(data_path)) {
    load(data_path)
    data_res_avg %>% mutate(STATE.NAME = .x)
  } else stop(glue("Missing data file: {.x}, {.y}"))
  
}) %>% 
  list_rbind()


# iterate plotting+saving across all spatial units
walk(unique(data_all$STATE.NAME), ~{
  
  overall_plot <- data_all |> 
    filter(STATE.NAME == .x) |> 
    gg_b_model("overall")
  
  ggsave(filename = glue("02_bird_figs/overall_{.x}.png"), 
         plot = overall_plot,
         dpi = 300, width = 11, height = 9, units = "in")
  
})


# prolific ----------------------------------------------------------------

load("00_outputs/bird_models/prolific.RData")

prolific_plot <- gg_b_model(prolific_pred, type = "prolific", prolific_points)

ggsave(filename = "02_bird_figs/prolific.png", plot = prolific_plot,
       dpi = 300, width = 11, height = 9, units = "in")

