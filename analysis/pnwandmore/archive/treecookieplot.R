
library(ggplot2)

test <- data.frame(x = 1, y = 0, ring = 0)
sigma <- 0.095/50
npoints <- 359
rw <- c(1e-10, 0.2,0.1, 0.2)
rwacc <- 0
for(i in 1:length(rw)){
  gen <- rnorm(npoints, rw[i], sigma)
  test <- rbind(
    test,
    data.frame(x = 1:360, y = c(rw[i]+rwacc, gen+rwacc), ring = i)
  )
  rwacc <- rwacc + rw[i]
}







ggplot(data = test) +
  geom_line(aes(x=x, y=y, group = ring)) +
  coord_polar()
