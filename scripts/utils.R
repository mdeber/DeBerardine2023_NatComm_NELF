# Misc global parameters --------------------------------------------------

# KD colors
color_LacZ <- rgb(0, 0, 0.5)
color_NELF <- rgb(0.6, 0, 0)
color_DSIF <- rgb(0.5, 0.75, 0)
color_Dual <- rgb(0.3, 0.1, 0.3)
KD_colors <- c(color_LacZ, color_NELF, color_DSIF, color_Dual)

# "change" colors
color_down <- cet_pal(10, "d1a")[1]
color_unchanged <- "#BEBEBE" # just "gray"
color_up <- cet_pal(10, "d1a")[9]
change_colors <- c(color_down, color_unchanged, color_up)

# FP TC colors
fptc_colors <- c("#000000", "#FCA729", "#F0471B", "#96042A")


# ggplot themes
theme_md_classic <- function() {
    theme_classic(
        base_size = 10, 
        base_family = "Helvetica",
        base_line_size = 0.25,
        base_rect_size = 0.25
    ) %+replace% theme(
        axis.text = element_text(color = "black", size = 8),
        axis.ticks = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        complete = TRUE
    )
}

theme_md_bw <- function() {
    theme_bw(
        base_size = 10, 
        base_family = "Helvetica",
        base_line_size = 0.25,
        base_rect_size = 0.25
    ) %+replace% theme(
        axis.text = element_text(color = "black", size = 8),
        axis.ticks = element_line(color = "black"),
        strip.background = element_blank(),
        strip.text = element_text(size = 10),
        strip.text.y = element_text(angle = 0),
        plot.title = element_text(hjust = 0.5, face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(fill = NA, color = "black"),
        complete = TRUE
    )
}

# Misc Functions ----------------------------------------------------------

pseudo_log10 <- function(x) {
    # pseudo-log transformation, see:
    # https://win-vector.com/2012/03/01/modeling-trick-the-signed-pseudo-logarithm/
    asinh(x/2) / log(10)
}

pseudotrans_rgb <- function(x, alpha) {
    # calculate a (non-transparent) output color based on what the input color
    # would look like on a white background with the given alpha level
    require(dplyr)
    if (length(x) > 1) {
        sapply(x, pseudotrans_rgb, alpha, USE.NAMES = FALSE)
    } else {
        # From transparency function: y = x*alpha + 255*(1-alpha)
        col2rgb(x) %>% 
            apply(1, . %>% "*"(alpha) %>% "+"(255*(1-alpha))) %>% 
            lapply("/", 255) %>% 
            do.call(rgb, .)
    }
}

dfList2df <- function(df.list, col_name = "sample") {
    # convert a list of dataframes into a single dataframe, with a column
    # for the original list element name
    df.list <- lapply(seq_along(df.list), function(i) {
        df <- df.list[[i]]
        df[[col_name]] <- names(df.list)[i]
        df
    })
    do.call(rbind, c(df.list, make.row.names = FALSE))
}

binVector <- function(x, binsize = NULL, nbins = NULL, FUN = sum) {
    # superimpose evenly spaced bins over a vector, and perform FUN on them
    #   (output length = nbins); this is faster than any R implementations
    #   I've seen (e.g. methods using findInterval)
    # (all conditions for safety)
    if (is.null(nbins)) {
        binsize <- as.integer(binsize)
        nbins <- as.integer(length(x) / binsize)
    } else if (is.null(binsize)) {
        nbins <- as.integer(nbins)
        binsize <- as.integer(length(x) / nbins)
    } else {
        nbins <- as.integer(nbins)
        binsize <- as.integer(binsize)
    }
    mat <- matrix(x[seq_len(nbins*binsize)], nrow = binsize)
    apply(mat, 2, FUN)
}

format_prettysci = function(x, parse = TRUE) {
    # for plot labels, formats text as pretty scientific format, e.g.
    # 3000 is written 3 x 10^3 (with superscripted 3)
    out <- format(x, scientific = TRUE) %>%
        gsub("0e\\+00","0", .) %>%
        gsub("^(.*)e", "'\\1'e", .) %>%
        gsub("\\+", "", .) %>%
        gsub("e", "%*%10^", .)
    
    if (parse) parse(text = out) else out
}
