

boot_se <- function(samplecol, fn = median, B = 100){
  sd(unlist(replicate(B, 
                      do.call(fn, list(sample(samplecol, n(), replace = T))), 
                      simplify = F)))
}

boot_se2 <- function(samplecol, fn = median, B = 100){
  replicate(B,
            do.call("fn", list(sample(samplecol, n(), replace = T))),
            simplify = F) %>% 
    unlist() %>% 
    sd()
}


x <- left_join(temp1, temp2) %>% 
  mutate(NO.SLISTS = ifelse(is.na(NO.SLISTS), 0, NO.SLISTS))

set.seed(420)
y <- x %>% ungroup() %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  summarise(NO = median(NO.SLISTS),
            NO.SE = boot_se2(NO.SLISTS, B = 100))

boot_fn = function(x, fn = median, B = 100) {
  1:B %>%
    # For each iteration, generate a sample of x with replacement
    map(~ x[sample(1:length(x), replace = TRUE)]) %>%
    # Obtain the fn estimate for each bootstrap sample
    map_dbl(fn) %>%
    # Obtain the standard error
    sd()
}


boot_fn2 = function(x, fn = median, B = 100) {
  1:B %>%
    # For each iteration, generate a sample of x with replacement
    map_dbl(~ x[sample(1:length(x), replace = TRUE)]) %>%
    # Obtain the fn estimate for each bootstrap sample
    map_dbl(fn) %>%
    # Obtain the standard error
    sd()
}

set.seed(420)
z <- x %>% ungroup() %>% 
  group_by(COVID, YEAR, MONTH, STATE) %>% 
  summarise(NO = median(NO.SLISTS),
            NO.SE = boot_fn2(NO.SLISTS, B = 100))


set.seed(57575)
rbind(system.time(x %>% ungroup() %>% 
                    group_by(COVID, YEAR, MONTH, STATE) %>% 
                    summarise(NO = median(NO.SLISTS),
                              NO.SE = boot_se2(NO.SLISTS, B = 100))),
      system.time(x %>% ungroup() %>% 
                    group_by(COVID, YEAR, MONTH, STATE) %>% 
                    summarise(NO = median(NO.SLISTS),
                              NO.SE = boot_fn(NO.SLISTS, B = 100))))


data1 <- mtcars %>% 
  select(gear, hp) %>% 
  group_by(gear)

data2 <- data %>% 
  summarise(hpmed = median(hp),
            hpse = boot_se2(hp))

data3 <- data %>% 
  summarise(hpmed = median(hp),
            hpse = boot_fn(hp))

library(microbenchmark)
microbenchmark(data2, data3, times = 1000)
