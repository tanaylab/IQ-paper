theme_arial_bw <- function(size = 8) {
    theme_arial <- ggplot2::theme_bw() %+replace% theme(
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(
            colour = "black", size = size,
            family = "ArialMT"
        )
    )
    theme_arial <- theme_arial + theme(text = element_text(
        size = size,
        family = "ArialMT"
    ), axis.text = element_text(
        size = size,
        family = "ArialMT"
    ), legend.text = element_text(
        size = size,
        family = "ArialMT"
    ), plot.title = element_text(
        size = size,
        family = "ArialMT"
    ))
    return(theme_arial)
}

theme_arial_void <- function(size = 8) {
    theme_arial <- ggplot2::theme_void() %+replace% theme(
        strip.text = element_text(
            colour = "black", size = size,
            family = "ArialMT"
        )
    )
    theme_arial <- theme_arial + theme(
        text = element_text(
            size = size,
            family = "ArialMT"
        ), legend.text = element_text(
            size = size,
            family = "ArialMT"
        ), plot.title = element_text(
            size = size,
            family = "ArialMT"
        ),
        plot.caption = element_text(
            size = size / 2,
            family = "ArialMT"
        )
    )
    return(theme_arial)
}

theme_arial_classic <- function(size = 8) {
    theme_arial <- ggplot2::theme_classic() %+replace% theme(
        panel.background = element_blank(),
        panel.grid.minor = element_blank(), strip.background = element_blank(),
        strip.text = element_text(
            colour = "black", size = size,
            family = "ArialMT"
        )
    )
    theme_arial <- theme_arial + theme(text = element_text(
        size = size,
        family = "ArialMT"
    ), axis.text = element_text(
        size = size,
        family = "ArialMT"
    ), legend.text = element_text(
        size = size,
        family = "ArialMT"
    ), plot.title = element_text(
        size = size,
        family = "ArialMT"
    ))
    return(theme_arial)
}