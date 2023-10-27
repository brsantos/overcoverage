## creating hex sticker

library(hexSticker)

library(ggplot2)
library(dplyr)

set.seed(123)
p <- ggplot(expand.grid(x = 1:60,
                  y = 1:60) %>%
         mutate(detected = sample(c(0, 1),
                                  size = 60 * 60,
                                  replace = TRUE,
                                  prob = c(0.95, 0.05)))) +
  geom_tile(aes(x, y, fill = factor(detected)),
            color = "grey90") +
  theme_void() +
  scale_fill_manual(values = c("royalblue", "white")) +
  theme(legend.position = "none") +
  coord_fixed() +
  theme_transparent()

library(showtext)
## Loading Google fonts (http://www.google.com/fonts)
font_add_google("Oswald")
## Automatically use showtext to render text for future devices
showtext_auto()

sticker(
  p,
  package = "overcoverage",
  p_size = 24,
  s_x = 1,
  s_y = 1,
  p_x = 1,
  p_y = 1.15,
  s_width = 2.5,
  s_height = 2.5,
  filename = "man/figures/overcoverage.png",
  p_family = "Oswald",
  p_color = "grey25",
  h_color = "grey75",
  white_around_sticker = TRUE,
  dpi = 300
)





