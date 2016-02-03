#===============================================================================
# Allometry graph for predictions
#===============================================================================
int   <- rep(0, 3)
slope <- c(1, 2, 3)
text  <- rep(NA, 3)
func  <- rep(NA, 3)
lines <- data.frame(int=int, slope=slope, text=text, fun = func)

gg_colour_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# calculates ggplot default colour values for 'n' colours
gg_colour_hue(3)

for (i in (1:3)) {
  df  <- lines
  l   <- list(b = df$slope[i])
  eq  <- substitute(italic(y) == italic(x)^b, l)
  eqn <- as.character(as.expression(eq))
  lines[i,3] <- eqn
}
lines
lines$colour <- c("#56B4E9", "#009E73", "#D55E00")

# Colour blind-friendly palette
#               Blue       Green     Dark Blue   Red
cbPalette <- c("#56B4E9", "#009E73", "#0072B2", "#D55E00")

dev.new(height = 5.5, width = 6)

# Single line for isometry
ggplot() +
  scale_x_continuous(limits = c(-1, 30)) +
  scale_y_continuous(limits = c(-4, 100)) +
  geom_abline(data=lines[2, ], mapping=aes(slope=slope, intercept=int, colour = text),
              size = 2) +
  scale_color_manual(values=c("#009E73")) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  theme(plot.title = element_text(hjust= -0.16, vjust= -1)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text = element_blank()) +
  geom_text(label = c("10", "20", "30"), aes(x = c(10, 20, 30), y = -4), size = 5) +
  geom_text(label = c("25", "50", "75", "100"), aes(x = c(-0.5), y = c(25, 50, 75, 100)), size = 5, hjust = 1) +
  geom_segment(aes(x = c(10, 20, 30), xend = c(10, 20, 30), y = -1, yend = 1)) +
  geom_segment(aes(x = -0.3, xend = 0.3, y = c(25, 50, 75, 100), yend = c(25, 50, 75, 100))) +
  theme(axis.title = element_blank()) 

# Two lines, isometric and negatively allometric
ggplot() +
  scale_x_continuous(limits = c(-1, 30)) +
  scale_y_continuous(limits = c(-4, 100)) +
  geom_abline(data=lines[1:2, ], mapping=aes(slope=slope, intercept=int, colour = text),
              size = 2) +
  scale_color_manual(values=c("#56B4E9", "#009E73", "#D55E00")) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  theme(plot.title = element_text(hjust= -0.16, vjust= -1)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text = element_blank()) +
  geom_text(label = c("10", "20", "30"), aes(x = c(10, 20, 30), y = -4), size = 5) +
  geom_text(label = c("25", "50", "75", "100"), aes(x = c(-0.5), y = c(25, 50, 75, 100)), size = 5, hjust = 1) +
  geom_segment(aes(x = c(10, 20, 30), xend = c(10, 20, 30), y = -1, yend = 1)) +
  geom_segment(aes(x = -0.3, xend = 0.3, y = c(25, 50, 75, 100), yend = c(25, 50, 75, 100))) +
  theme(axis.title = element_blank()) 

# All three lines, isometry, negative allometry, and positive allometry
ggplot() +
  scale_x_continuous(limits = c(-1, 30)) +
  scale_y_continuous(limits = c(-4, 100)) +
  geom_abline(data=lines, mapping=aes(slope=slope, intercept=int, colour = text),
              size = 2) +
  scale_color_manual(values=c("#56B4E9", "#009E73", "#D55E00")) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  theme(plot.title = element_text(hjust= -0.16, vjust= -1)) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.ticks = element_blank()) +
  theme(axis.text = element_blank()) +
  geom_text(label = c("10", "20", "30"), aes(x = c(10, 20, 30), y = -4), size = 5) +
  geom_text(label = c("25", "50", "75", "100"), aes(x = c(-0.5), y = c(25, 50, 75, 100)), size = 5, hjust = 1) +
  geom_segment(aes(x = c(10, 20, 30), xend = c(10, 20, 30), y = -1, yend = 1)) +
  geom_segment(aes(x = -0.3, xend = 0.3, y = c(25, 50, 75, 100), yend = c(25, 50, 75, 100))) +
  theme(axis.title = element_blank()) 