% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.r
\name{plot.sreg}
\alias{plot.sreg}
\title{Plot Method for `sreg` Objects}
\usage{
\method{plot}{sreg}(
  x,
  treatment_labels = NULL,
  title = "Estimated ATEs with Confidence Intervals",
  bar_fill = NULL,
  point_shape = 23,
  point_size = 3,
  point_fill = "white",
  point_stroke = 1.2,
  point_color = "black",
  label_color = "black",
  label_size = 4,
  bg_color = NULL,
  grid = TRUE,
  zero_line = TRUE,
  y_axis_title = NULL,
  x_axis_title = NULL,
  ...
)
}
\arguments{
\item{x}{An object of class \code{sreg}.}

\item{treatment_labels}{Optional vector of treatment labels to display on the y-axis. If \code{NULL}, default labels like "Treatment 1", "Treatment 2", etc., are used.}

\item{title}{Optional plot title. Defaults to "Estimated ATEs with Confidence Intervals".}

\item{bar_fill}{Optional fill color(s) for the confidence interval bars. Can be \code{NULL} (default viridis scale), a single color, or a vector of two colors for a gradient.}

\item{point_shape}{Optional shape of the point used to mark the estimated ATE. Default is 23 (a diamond).}

\item{point_size}{Optional size of the point marking the ATE.}

\item{point_fill}{Optional fill color of the ATE point shape.}

\item{point_stroke}{Optional stroke (border) thickness of the ATE point shape.}

\item{point_color}{Optional outline color of the ATE point.}

\item{label_color}{Optional color of the text label displaying the estimate and standard error.}

\item{label_size}{Optional size of the text label displaying the estimate and standard error.}

\item{bg_color}{Optional background color of the plot panel. If \code{NULL}, the default theme background is used.}

\item{grid}{Optional logical flag. If \code{TRUE} (default), grid lines are shown; if \code{FALSE}, they are removed.}

\item{zero_line}{Optional logical flag. If \code{TRUE} (default), a vertical dashed line at 0 is added for reference.}

\item{y_axis_title}{Optional title of the y-axis. If \code{NULL}, no y-axis label is added.}

\item{x_axis_title}{Optional title of the x-axis. If \code{NULL}, no x-axis label is added.}

\item{...}{Additional arguments passed to other methods.}
}
\value{
Invisibly returns the ggplot object. Called for its side effects (i.e., generating a plot).
}
\description{
Visualize estimated ATEs and confidence intervals for objects of class \code{sreg}.
}
