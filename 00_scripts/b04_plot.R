require(ggplot2)

theme_set(theme_classic())
covid_palette <- c("#1B9E77", "#EF4050", "#E89005", "#9678B6")

# prolific ----------------------------------------------------------------

prolific_plot <- gg_b_model(prolific_pred, type = "prolific", prolific_points)

ggsave(filename = glue("03_wrap_figs/b_prolific_model.png"), plot = prolific_plot,
       dpi = 300, width = 21, height = 6, units = "in")

