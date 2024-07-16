# aggregating all states and months data --------------------------------------------

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


# ribbon fig ------------------------------------------------------------------------

create_bird_graph <- function() {
  
  data_all <- data %>% filter(MONTHS.TYPE = "ALL")
  data_ld <- data %>% filter(MONTHS.TYPE = "LD")
  
  data_all %>% 
    gg_b_model("overall_ribbon") +
    facet_wrap(~ STATE.NAME)
  
  
}