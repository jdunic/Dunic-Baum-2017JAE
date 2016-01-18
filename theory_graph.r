# Creates log-log power function graph  ####

library('directlabels')

# Possibly more biologically relevant values:
int   <- rep(1, 8)
slope <- c(0, 0.67, 0.75, 0.9, 1, 1.1, 1.25, 1.33)
text <- rep(NA, 8)
lines <- data.frame(int=int, slope=slope, text=text)

for (i in (1:8)) {
  df  <- lines
  l   <- list(b = df$slope[i])
  eq  <- substitute(italic(y) == italic(x)^b, l)
  eqn <- as.character(as.expression(eq))
  lines[i,3] <- eqn
}
lines

# But these look better...
int   <- rep(1, 7)
slope <- c(0, 0.5, 0.75, 1, 1.25, 1.5, 2)
text  <- rep(NA, 7)
func  <- rep(NA, 7)
lines <- data.frame(int=int, slope=slope, text=text, fun = func)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# calculates ggplot default colour values for 'n' colours
gg_color_hue(7)

for (i in (1:7)) {
  df  <- lines
  l   <- list(b = df$slope[i])
  eq  <- substitute(italic(y) == italic(x)^b, l)
  eqn <- as.character(as.expression(eq))
  lines[i,3] <- eqn
}
lines

transformed <- ggplot() +
  scale_x_continuous(name = "log(x)", limits = c(-1, 5)) +
  scale_y_continuous(name = "log(y)", limits = c(-1, 10)) +
  geom_abline(data=lines, mapping=aes(slope=slope, intercept=int, colour=text)) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  ggtitle("b)") +
  theme(plot.title = element_text(hjust= -0.16, vjust= -1))

#transformed

raw <- ggplot(data.frame(x = c(-1, 5)), aes(x)) +
  scale_x_continuous(name = "x", limits = c(-1, 5)) +
  scale_y_continuous(name = "y", limits = c(-1, 10)) +
  #geom_abline(data=lines, mapping=aes(slope=slope, intercept=int, colour=text)) +
  geom_vline(xintercept=0) +
  geom_hline(yintercept=0) +
  stat_function(fun = function(x) { x^0    }, colour = "#F8766D") +
  stat_function(fun = function(x) { x^0.5  }, colour = "#C49A00") +
  stat_function(fun = function(x) { x^0.75 }, colour = "#53B400") +
  stat_function(fun = function(x) { x^1    }, colour = "#00C094") +
  stat_function(fun = function(x) { x^1.25 }, colour = "#00B6EB") +
  stat_function(fun = function(x) { x^1.5  }, colour = "#A58AFF") +
  stat_function(fun = function(x) { x^2    }, colour = "#FB61D7") +
  ggtitle("a)") +
  theme(plot.title = element_text(hjust= -0.16, vjust= -1))

#raw



# Hypothesis graph ####
# But these look better...
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

ggplot() +
  scale_x_continuous(name = "log(standard length)", limits = c(-1, 30)) +
  scale_y_continuous(name = "log(gape area)", limits = c(-4, 100)) +
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
  geom_text(label = c("10", "20", "30"), aes(x = c(10, 20, 30), y = -4), size = 4) +
  geom_text(label = c("25", "50", "75", "100"), aes(x = c(-1), y = c(25, 50, 75, 100)), size = 4) +
  geom_segment(aes(x = c(10, 20, 30), xend = c(10, 20, 30), y = -1, yend = 1)) +
  geom_segment(aes(x = -0.3, xend = 0.3, y = c(25, 50, 75, 100), yend = c(25, 50, 75, 100))) 





ggplot() +
  #scale_x_continuous(name = "log(x)", limits = c(-1, 30)) +
  #scale_y_continuous(name = "log(y)", limits = c(-1, 100)) +
  scale_x_log10(limits = c(0.1, 10000)) +
  scale_y_log10(limits = c(0.1, 1000)) +
  geom_abline(data=lines, mapping=aes(slope=slope, intercept=int, colour=text),
              size = 2) +
  geom_vline(xintercept=1) +
  geom_hline(yintercept=1) +
  theme_bw() +
  theme(panel.border = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  #theme(axis.line = element_line(color = 'black')) +
  theme(plot.title = element_text(hjust= -0.16, vjust= -1))