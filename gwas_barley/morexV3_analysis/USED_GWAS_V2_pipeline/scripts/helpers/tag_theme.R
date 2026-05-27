# tag_theme.R
# TAG (Theoretical and Applied Genetics, Springer) figure-spec constants and helpers.
# Sourced by every plot-producing R script in this pipeline.

# Always set TMPDIR before anything that might allocate (data.table, ggplot, etc.).
Sys.setenv(TMPDIR = "/mnt/data/shahar/.tmp")

# --- Springer / TAG figure spec ---
TAG_WIDTH_MM       <- 174    # double-column width
TAG_SINGLECOL_MM   <- 84     # single-column width
TAG_FONT_PT        <- 10     # body text size
TAG_LINE_LWD       <- 0.6    # line width (axes, reference lines)

mm_to_in <- function(mm) mm / 25.4

# Vector-PDF device with TAG defaults; size in inches (default 4x4 per panel)
tag_cairo_pdf <- function(file, width_in = 4, height_in = 4) {
  cairo_pdf(file, width = width_in, height = height_in,
            family = "sans", pointsize = TAG_FONT_PT)
}

# Raster-PNG device with TAG defaults at 300 dpi
tag_png <- function(file, width_in = 4, height_in = 4, res = 300) {
  png(file, width = width_in, height = height_in, units = "in", res = res,
      family = "sans", pointsize = TAG_FONT_PT)
}

# Open BOTH PDF and PNG for the same plot, run the given drawing function once
# per device. Pattern:
#   tag_plot_dual(file_pdf, file_png, function() { plot(...); abline(...); })
tag_plot_dual <- function(file_pdf, file_png, draw_fn,
                          width_in = 4, height_in = 4, res = 300) {
  tag_cairo_pdf(file_pdf, width_in, height_in); draw_fn(); dev.off()
  tag_png(file_png, width_in, height_in, res);  draw_fn(); dev.off()
}
