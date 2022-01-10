

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



