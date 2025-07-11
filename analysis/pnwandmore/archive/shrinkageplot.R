
library(ggplot2)
library(ggforce)

test <- data.frame(
  species = c(1,1,2,2,3,3),
  q10 = c(0.4,0.35,0.15,0.18,0.07,0.16),
  q50 = c(0.5,0.4,0.2,0.2,0.1,0.19),
  pooling = rep(c('no', 'partial'), 3)
)
test$q90 <- 2*test$q50-test$q10


ggplot(data = test) + 
  geom_ribbon(aes(x = pooling, ymin = q10, ymax = q90, group = species),
              position = position_dodge(width = 0.1), alpha = 0.1) +
  geom_pointrange(aes(x = pooling, y = q50, ymin = q10, ymax = q90, group = species, color = as.character(species)),
                  position = position_dodge(width = 0.1)) +
  theme_classic()


