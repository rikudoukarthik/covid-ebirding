# global setup steps

library(tidyverse)
source("00_scripts/functions.R")

# spatial units to iterate over
spat_iter <- c("India", "Karnataka", "Kerala", "Maharashtra", "Assam")

# species for analysis --------------------------------------------------------------

# Guwahati lacks PBFl, WCBa (but BTBa), InRo, GrFr, PiBu, BWKi, GHSw (in good numbers)
# Ernakulam lacks GrFr, SmMi, BWKi, LTSh, PiBu (in good numbers)
# different from list for BLR and PUN because 5 non-urban species from that list do not occur in sufficient numbers in EKM

bird_anal_spec_list <- function() {
  
  species_bang <- data.frame(COMMON.NAME = c(
    "House Crow", "Common Myna", "Jungle Myna", "Rock Pigeon", "Black Kite", 
    "Common Tailorbird", "Rose-ringed Parakeet", "Asian Koel", "Purple Sunbird", 
    "Pale-billed Flowerpecker", "White-cheeked Barbet", 
    
    "Little Cormorant", "Gray-headed Swamphen", "Indian Pond-Heron",
    
    "Green Bee-eater", "Black Drongo", "Pied Bushchat", "Plain Prinia", "Indian Robin", 
    "House Sparrow", "Spotted Dove", "Coppersmith Barbet", "Black-winged Kite", 
    "Small Minivet", "Long-tailed Shrike", "Gray Francolin"
  ),
  SP.CATEGORY = c(rep("U", 11), rep("W", 3), rep("R", 12)))
  
  # removing Black Drongo for models to run
  species_pune <- data.frame(COMMON.NAME = c(
    "House Crow", "Common Myna", "Jungle Myna", "Rock Pigeon", "Black Kite", 
    "Common Tailorbird", "Rose-ringed Parakeet", "Asian Koel", "Purple Sunbird", 
    "Pale-billed Flowerpecker", "White-cheeked Barbet", 
    
    "Little Cormorant", "Gray-headed Swamphen", "Indian Pond-Heron",
    
    "Green Bee-eater", "Pied Bushchat", "Plain Prinia", "Indian Robin", 
    "House Sparrow", "Spotted Dove", "Coppersmith Barbet", "Black-winged Kite", 
    "Small Minivet", "Long-tailed Shrike", "Gray Francolin"
  ),
  SP.CATEGORY = c(rep("U", 11), rep("W", 3), rep("R", 11)))
  
  species_koch <- data.frame(COMMON.NAME = c(
    "House Crow", "Common Myna", "Jungle Myna", "Rock Pigeon", "Black Kite", 
    "Common Tailorbird", "Rose-ringed Parakeet", "Asian Koel", "Purple Sunbird",
    "Pale-billed Flowerpecker", "White-cheeked Barbet", 
    
    "Little Cormorant", "Gray-headed Swamphen", "Indian Pond-Heron",
    
    "Green Bee-eater", "Black Drongo", "Plain Prinia", "Indian Robin", 
    "House Sparrow", "Spotted Dove", "Coppersmith Barbet", "White-throated Kingfisher", 
    "Black-rumped Flameback", "Rufous Treepie"
  ),
  SP.CATEGORY = c(rep("U", 11), rep("W", 3), rep("R", 10)))
  
  species_guwa <- data.frame(COMMON.NAME = c(
    "House Crow", "Common Myna", "Jungle Myna", "Rock Pigeon", "Black Kite", "Common Tailorbird", 
    "Rose-ringed Parakeet", "Asian Koel", "Purple Sunbird", "Blue-throated Barbet",
    
    "Little Cormorant", "Indian Pond-Heron",
    
    "Green Bee-eater", "Black Drongo", "Plain Prinia", "House Sparrow", "Spotted Dove",
    "Coppersmith Barbet", "White-throated Kingfisher", "Rufous Treepie", "Long-tailed Shrike"
  ),
  SP.CATEGORY = c(rep("U", 10), rep("W", 2), rep("R", 9)))
  
  
  species_list <- bind_rows("Karnataka" = species_bang,
                            "Maharashtra" = species_pune,
                            "Kerala" = species_koch,
                            "Assam" = species_guwa,
                            .id = "STATE") %>% 
    filter(!(SP.CATEGORY %in% "W")) # removing wetland spp. cos we want to compare UNU
  
  # adding species list for full country analysis
  species_list <- bind_rows(species_list,
                            species_list |> 
                              distinct(COMMON.NAME, SP.CATEGORY) |> 
                              mutate(STATE = "India") |> 
                              relocate(STATE, COMMON.NAME, SP.CATEGORY))
  
  return(species_list)
  
}

