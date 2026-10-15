# ============================================================
# Lab Stats Explorer
# A Shiny app for flexible statistics on biological lab-experiment
# data (e.g. sorted cell populations, cell assays, enzyme kinetics,
# survival studies, dose-response/titration experiments).
#
# Analysis modes:
#  1. Compare groups
#       - One grouping factor:
#           parametric    -> t-test (2 groups) / one-way ANOVA (3+)
#                             multiple comparisons: Tukey HSD or
#                             Bonferroni pairwise t-test
#           non-parametric-> Wilcoxon (2 groups) / Kruskal-Wallis (3+)
#                             multiple comparisons: Dunn's test (BH),
#                             pairwise Wilcoxon (BH or Bonferroni)
#           -> plot shows p-value brackets for 2-4 groups
#           -> optional: plot ALL numeric columns at once, faceted
#       - Two grouping factors (Two-way ANOVA):
#           parametric only -> main effects + interaction + Tukey HSD
#  2. Correlation & regression -> Pearson/Spearman + simple linear
#     regression
#  3. Survival analysis -> Kaplan-Meier curves + log-rank test
#  4. Titration / dose-response -> four-parameter logistic curve fit
#     (EC50, Hill-type slope) per group
#
# Plot customization: color palette, x/y axis label override,
# x-axis label angle, text size.
#
# HOW TO RUN:
#   install.packages(c("shiny","bslib","readxl","dplyr","tidyr",
#                       "ggplot2","DT","ggsignif","dunn.test","survival"))
#   shiny::runApp("app.R")
# ============================================================

library(shiny)
library(bslib)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)
library(ggsignif)
library(dunn.test)
library(survival)
library(gridExtra)
library(grid)
library(colourpicker)
library(scales)
library(svglite)
library(httr2)

options(shiny.maxRequestSize = 25 * 1024^2)  # 25 MB upload limit

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || a == "") b else a

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, digits = 3, format = "f"))
}

# Same as format_p() but used for embedding inline in a sentence (just an
# alias kept separate in case the two ever need to diverge in wording).
fmt_p_clause <- function(p) format_p(p)

posthoc_label <- function(method) {
  switch(method,
    "tukey"       = "Tukey HSD",
    "bonf"        = "Bonferroni-corrected pairwise t-tests",
    "dunn_bh"     = "Dunn's test, BH-adjusted",
    "wilcox_bh"   = "pairwise Wilcoxon tests, BH-adjusted",
    "wilcox_bonf" = "pairwise Wilcoxon tests, Bonferroni-adjusted",
    "the selected post-hoc method"
  )
}

get_fill_scale <- function(palette, custom_colors = NULL) {
  if (palette == "Custom" && !is.null(custom_colors)) return(scale_fill_manual(values = custom_colors))
  switch(palette,
    "Viridis"  = scale_fill_viridis_d(),
    "Set1"     = scale_fill_brewer(palette = "Set1"),
    "Set2"     = scale_fill_brewer(palette = "Set2"),
    "Dark2"    = scale_fill_brewer(palette = "Dark2"),
    "Paired"   = scale_fill_brewer(palette = "Paired"),
    "Pastel1"  = scale_fill_brewer(palette = "Pastel1"),
    NULL
  )
}

get_color_scale <- function(palette, custom_colors = NULL) {
  if (palette == "Custom" && !is.null(custom_colors)) return(scale_color_manual(values = custom_colors))
  switch(palette,
    "Viridis"  = scale_color_viridis_d(),
    "Set1"     = scale_color_brewer(palette = "Set1"),
    "Set2"     = scale_color_brewer(palette = "Set2"),
    "Dark2"    = scale_color_brewer(palette = "Dark2"),
    "Paired"   = scale_color_brewer(palette = "Paired"),
    "Pastel1"  = scale_color_brewer(palette = "Pastel1"),
    NULL
  )
}

axis_angle_theme <- function(angle) {
  if (angle == 0) theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  else theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1))
}

# Axis line/text color + optional gridline removal, applied on top of the
# base theme for every plot.
extra_style_theme <- function(axis_color, hide_gridlines) {
  thm <- theme(
    axis.text  = element_text(color = axis_color),
    axis.title = element_text(color = axis_color),
    axis.line  = element_line(color = axis_color),
    axis.ticks = element_line(color = axis_color)
  )
  if (isTRUE(hide_gridlines)) thm <- thm + theme(panel.grid = element_blank())
  thm
}

# Log10-transforms the y-axis if requested. Returns NULL (a no-op layer)
# otherwise, since ggplot objects can have NULL added to them safely.
log_y_scale <- function(apply_log) {
  if (isTRUE(apply_log)) scale_y_log10() else NULL
}

# Builds a label with a small (?) icon that shows `tip` text on hover.
# Use this in place of a plain string wherever an input's `label = ` argument
# is set, to give the person a short explanation of what that option does.
tooltip_label <- function(text, tip, placement = "right") {
  tagList(
    text,
    tooltip(
      trigger = icon("circle-question",
                      style = "font-size:0.8em; color:#888; margin-left:5px; cursor:help;"),
      tip,
      placement = placement
    )
  )
}

# Run a group comparison test (omnibus + chosen multiple-comparison method)
# on an arbitrary (group, value) data frame. Returns NULL on failure.
run_group_test <- function(ad, family, posthoc_method,
                            paired = FALSE, equal_var = FALSE, paired_np = FALSE) {
  ad$group <- droplevels(ad$group)
  n <- nlevels(ad$group)
  if (n < 2) return(NULL)

  if (n == 2) {
    if (family == "param") {
      res <- tryCatch(t.test(value ~ group, data = ad, paired = paired, var.equal = equal_var),
                       error = function(e) NULL)
      if (is.null(res)) return(NULL)
      return(list(type = "t-test", result = res, posthoc = NULL, posthoc_kind = NULL))
    } else {
      if (isTRUE(paired_np)) {
        g <- levels(ad$group)
        v1 <- ad$value[ad$group == g[1]]; v2 <- ad$value[ad$group == g[2]]
        if (length(v1) != length(v2)) return(NULL)
        res <- tryCatch(wilcox.test(v1, v2, paired = TRUE), error = function(e) NULL)
      } else {
        res <- tryCatch(wilcox.test(value ~ group, data = ad), error = function(e) NULL)
      }
      if (is.null(res)) return(NULL)
      return(list(type = "wilcoxon", result = res, posthoc = NULL, posthoc_kind = NULL))
    }
  }

  # 3+ groups
  if (family == "param") {
    fit <- tryCatch(aov(value ~ group, data = ad), error = function(e) NULL)
    if (is.null(fit)) return(NULL)
    if (posthoc_method == "bonf") {
      ph <- tryCatch(pairwise.t.test(ad$value, ad$group, p.adjust.method = "bonferroni"),
                      error = function(e) NULL)
      list(type = "anova", result = fit, posthoc = ph, posthoc_kind = "pairwise_matrix")
    } else {
      ph <- tryCatch(TukeyHSD(fit), error = function(e) NULL)
      list(type = "anova", result = fit, posthoc = ph, posthoc_kind = "tukey")
    }
  } else {
    res <- tryCatch(kruskal.test(value ~ group, data = ad), error = function(e) NULL)
    if (is.null(res)) return(NULL)
    if (posthoc_method == "wilcox_bh") {
      ph <- tryCatch(pairwise.wilcox.test(ad$value, ad$group, p.adjust.method = "BH"),
                      error = function(e) NULL)
      list(type = "kruskal", result = res, posthoc = ph, posthoc_kind = "pairwise_matrix")
    } else if (posthoc_method == "wilcox_bonf") {
      ph <- tryCatch(pairwise.wilcox.test(ad$value, ad$group, p.adjust.method = "bonferroni"),
                      error = function(e) NULL)
      list(type = "kruskal", result = res, posthoc = ph, posthoc_kind = "pairwise_matrix")
    } else {
      dt <- tryCatch({
        invisible(capture.output(
          out <- dunn.test::dunn.test(ad$value, ad$group, method = "bh",
                                       table = FALSE, list = FALSE, kw = FALSE)
        ))
        out
      }, error = function(e) NULL)
      list(type = "kruskal", result = res, posthoc = dt, posthoc_kind = "dunn")
    }
  }
}

# Extract pairwise group1/group2/p data frame from a run_group_test() result.
extract_pairwise <- function(tr, group_levels) {
  empty <- data.frame(group1 = character(), group2 = character(), p = numeric())
  if (is.null(tr)) return(empty)

  if (tr$type %in% c("t-test", "wilcoxon")) {
    data.frame(group1 = group_levels[1], group2 = group_levels[2], p = tr$result$p.value,
               stringsAsFactors = FALSE)
  } else if (!is.null(tr$posthoc_kind) && tr$posthoc_kind == "tukey") {
    tk <- as.data.frame(tr$posthoc$group)
    comp <- strsplit(rownames(tk), "-")
    data.frame(group1 = sapply(comp, `[`, 2), group2 = sapply(comp, `[`, 1),
               p = tk[["p adj"]], stringsAsFactors = FALSE)
  } else if (!is.null(tr$posthoc_kind) && tr$posthoc_kind == "pairwise_matrix") {
    if (is.null(tr$posthoc)) return(empty)
    pm <- tr$posthoc$p.value
    out <- as.data.frame(as.table(pm), stringsAsFactors = FALSE)
    names(out) <- c("group2", "group1", "p")
    out[!is.na(out$p), c("group1", "group2", "p")]
  } else if (!is.null(tr$posthoc_kind) && tr$posthoc_kind == "dunn") {
    if (is.null(tr$posthoc)) return(empty)
    comp <- strsplit(tr$posthoc$comparisons, " - ")
    data.frame(group1 = trimws(sapply(comp, `[`, 1)), group2 = trimws(sapply(comp, `[`, 2)),
               p = tr$posthoc$P.adjusted, stringsAsFactors = FALSE)
  } else {
    empty
  }
}

# For the two-way ANOVA plot: finds pairwise comparisons from the
# factor1:factor2 Tukey HSD table that compare different factor2 levels
# *within the same factor1 level* (e.g. within "Kidney", Control vs Drug) --
# these are the comparisons that make sense to draw as a bracket between
# two dodged bars/boxes sitting at the same x position. Comparisons that
# cross factor1 levels (e.g. Kidney:Control vs Spleen:Drug) are skipped,
# since they don't correspond to two bars next to each other.
twoway_bracket_df <- function(ad, fit, dodge_width = 0.75) {
  f1_levels <- levels(ad$factor1); f2_levels <- levels(ad$factor2); n2 <- length(f2_levels)
  tk <- tryCatch(TukeyHSD(fit)[["factor1:factor2"]], error = function(e) NULL)
  if (is.null(tk)) return(data.frame())
  tk_df <- as.data.frame(tk); tk_df$comparison <- rownames(tk_df)

  offset <- function(j) dodge_width * ((j - 0.5) / n2 - 0.5)

  rows <- list()
  for (i in seq_len(nrow(tk_df))) {
    parts <- strsplit(tk_df$comparison[i], "-")[[1]]
    left  <- strsplit(parts[1], ":")[[1]]
    right <- strsplit(parts[2], ":")[[1]]
    if (left[1] != right[1]) next  # different factor1 level -- not a same-cluster comparison

    x_base <- match(left[1], f1_levels)
    j_left <- match(left[2], f2_levels); j_right <- match(right[2], f2_levels)
    xmin <- x_base + offset(j_left); xmax <- x_base + offset(j_right)
    if (xmin > xmax) { tmp <- xmin; xmin <- xmax; xmax <- tmp }

    base_y <- suppressWarnings(max(ad$value[ad$factor1 == left[1]], na.rm = TRUE))
    rows[[length(rows) + 1]] <- data.frame(
      factor1 = left[1], xmin = xmin, xmax = xmax, p = tk_df[["p adj"]][i], base_y = base_y
    )
  }
  if (length(rows) == 0) return(data.frame())

  out <- do.call(rbind, rows)
  rng <- suppressWarnings(max(ad$value, na.rm = TRUE) - min(ad$value, na.rm = TRUE))
  step <- if (rng > 0) rng * 0.1 else abs(max(ad$value, na.rm = TRUE)) * 0.1 + 0.01

  out <- out %>%
    group_by(factor1) %>%
    mutate(y_position = base_y + step * row_number()) %>%
    ungroup()
  out$label <- sapply(out$p, format_p)
  out
}

# ---------------------------------------------------------
# UI
# ---------------------------------------------------------
ui <- page_sidebar(
  title = "Lab Stats Explorer",
  theme = bs_theme(version = 5, primary = "#2C6E49", base_font = font_google("Inter"),
                    font_size_base = "0.85rem"),
  tags$head(tags$style(HTML(
    "pre { white-space: pre-wrap; word-wrap: break-word; overflow-wrap: break-word; max-width: 100%; }"
  ))),

  sidebar = sidebar(
    width = 400,
    fileInput("datafile",
              tooltip_label("Upload data (CSV or Excel)",
                             "Each row should be one sample/replicate; each column a variable (e.g. group, treatment, measurement)."),
              accept = c(".csv", ".xlsx", ".xls")),
    checkboxInput("header",
                  tooltip_label("File has column headers",
                                 "Check this if the first row of your file contains column names rather than data."),
                  value = TRUE),
    hr(),

    radioButtons(
      "analysis_mode",
      tooltip_label("Analysis type", "Choose what kind of statistical analysis to run on your data."),
      choices = c("Compare groups" = "groups",
                  "Correlation & regression" = "corr",
                  "Survival analysis" = "survival",
                  "Titration / dose-response" = "titration"),
      selected = "groups"
    ),
    hr(),

    conditionalPanel(
      "output.file_uploaded == 'yes'",
      uiOutput("mode_ui"),
      hr()
    ),
    numericInput("alpha",
                 tooltip_label("Significance level (alpha)",
                                "The p-value threshold below which a result is called statistically significant. 0.05 is the conventional default."),
                 value = 0.05, min = 0.001, max = 0.5, step = 0.01),
    helpText("Upload a file, choose an analysis type, then pick the ",
             "relevant columns below. All plot styling and sizing options ",
             "are on the right side of the Plot tab.")
  ),

  navset_card_tab(
    nav_panel("Data Preview", DTOutput("data_preview")),

    nav_panel(
      "Descriptive Stats",
      conditionalPanel("input.analysis_mode == 'groups' && input.group_design != 'two'",
                        DTOutput("desc_table")),
      conditionalPanel("input.analysis_mode == 'groups' && input.group_design == 'two'",
                        DTOutput("desc_table_2f")),
      conditionalPanel("input.analysis_mode == 'corr'", DTOutput("corr_desc_table")),
      conditionalPanel("input.analysis_mode == 'survival'", DTOutput("surv_summary_table")),
      conditionalPanel("input.analysis_mode == 'titration'", DTOutput("titration_fit_table"))
    ),

    nav_panel(
      "Normality Check",
      conditionalPanel(
        "input.analysis_mode == 'groups' && input.group_design != 'two'",
        helpText("Shapiro-Wilk test per group. p < 0.05 suggests the data ",
                 "deviate from a normal distribution -- consider a ",
                 "non-parametric test in that case."),
        DTOutput("normality_table")
      ),
      conditionalPanel(
        "input.analysis_mode == 'groups' && input.group_design == 'two'",
        helpText("Shapiro-Wilk test on the two-way ANOVA model residuals. ",
                 "p < 0.05 suggests the residuals deviate from normality, ",
                 "so ANOVA results should be interpreted with some caution."),
        verbatimTextOutput("normality_2f")
      ),
      conditionalPanel(
        "input.analysis_mode == 'corr' || input.analysis_mode == 'survival' || input.analysis_mode == 'titration'",
        helpText("Normality checking is not applicable to this analysis type.")
      )
    ),

    nav_panel(
      "Test Results",
      conditionalPanel(
        "input.analysis_mode == 'groups' && input.group_design != 'two'",
        verbatimTextOutput("test_summary"),
        conditionalPanel(
          "output.show_posthoc == 'yes'",
          h5("Multiple comparisons"),
          DTOutput("posthoc_table")
        )
      ),
      conditionalPanel(
        "input.analysis_mode == 'groups' && input.group_design == 'two'",
        verbatimTextOutput("test_summary_2f"),
        h5("Tukey HSD post-hoc (main effects + interaction)"),
        DTOutput("posthoc_table_2f")
      ),
      conditionalPanel("input.analysis_mode == 'corr'", verbatimTextOutput("corr_summary")),
      conditionalPanel("input.analysis_mode == 'survival'", verbatimTextOutput("surv_test_summary")),
      conditionalPanel("input.analysis_mode == 'titration'", verbatimTextOutput("titration_summary"))
    ),

    nav_panel(
      "Summary",
      h5("Auto-generated summary"),
      helpText("A plain-language summary written directly from the statistics above -- ",
               "every number in it comes straight from your results, nothing is estimated or guessed."),
      verbatimTextOutput("auto_summary"),
      conditionalPanel(
        "output.claude_key_missing == 'no'",
        hr(),
        actionButton("polish_btn", "Polish with Claude", icon = icon("wand-magic-sparkles")),
        hr(),
        h5("Polished version"),
        uiOutput("polished_summary_ui")
      )
    ),

    nav_panel(
      "Plot",
      tags$script(HTML("
        function copyPlotToClipboard() {
          var container = document.getElementById('main_plot');
          var img = container ? container.querySelector('img') : null;
          if (!img) { alert('Plot is not ready yet -- wait for it to render, then try again.'); return; }
          fetch(img.src)
            .then(function(res) { return res.blob(); })
            .then(function(blob) {
              return navigator.clipboard.write([new ClipboardItem({ [blob.type]: blob })]);
            })
            .then(function() { alert('Plot copied to clipboard!'); })
            .catch(function(err) { alert('Could not copy the plot: ' + err); });
        }
      ")),
      layout_sidebar(
        sidebar = sidebar(
          position = "right",
          width = 300,
          open = "always",

          h6("Export"),
          tags$div(
            class = "dropdown",
            tags$button(
              class = "btn btn-default dropdown-toggle", type = "button",
              `data-bs-toggle` = "dropdown", `aria-expanded` = "false",
              icon("download"), " Export plot"
            ),
            tags$ul(
              class = "dropdown-menu",
              tags$li(tags$a(class = "dropdown-item", href = "#",
                              onclick = "copyPlotToClipboard(); return false;",
                              icon("copy"), " Copy to clipboard")),
              tags$li(downloadLink("download_plot", label = tagList(icon("image"), " PNG"), class = "dropdown-item")),
              tags$li(downloadLink("download_plot_svg", label = tagList(icon("bezier-curve"), " SVG"), class = "dropdown-item")),
              tags$li(downloadLink("download_plot_eps", label = tagList(icon("file"), " EPS"), class = "dropdown-item")),
              tags$li(downloadLink("download_report", label = tagList(icon("file-pdf"), " Full report (PDF)"), class = "dropdown-item"))
            )
          ),
          hr(),

          h6("Size & spacing"),
          sliderInput("plot_width_px",
                      tooltip_label("Plot width (px)", "Overall width of the plot. Narrowing this (along with a smaller box/bar width) keeps things looking compact instead of scattered."),
                      min = 400, max = 1400, value = 800, step = 50),
          sliderInput("plot_height_px",
                      tooltip_label("Plot height (px)", "Overall height of the plot."),
                      min = 300, max = 900, value = 500, step = 50),
          sliderInput("point_size",
                      tooltip_label("Point size", "Size of the individual data points/dots shown on the plot."),
                      min = 0.5, max = 6, value = 2, step = 0.5),
          conditionalPanel(
            "input.analysis_mode == 'groups'",
            sliderInput("box_width",
                        tooltip_label("Box/Bar width", "Width of each individual box or bar."),
                        min = 0.2, max = 0.9, value = 0.6, step = 0.05),
            sliderInput("category_spacing",
                        tooltip_label("Space between groups", "Controls the gap left on either side of the groups and between them, independent of box/bar width."),
                        min = 0, max = 0.4, value = 0.1, step = 0.05)
          ),
          hr(),

          h6("Style"),
          selectInput(
            "palette",
            tooltip_label("Color palette", "The color scheme used for groups/lines. Choose \"Custom\" to pick each group's color yourself."),
            choices = c("Default", "Viridis", "Set1", "Set2", "Dark2", "Paired", "Pastel1", "Custom"),
            selected = "Default"
          ),
          conditionalPanel("input.palette == 'Custom'", uiOutput("custom_color_ui")),
          conditionalPanel(
            "input.analysis_mode == 'groups'",
            selectInput("chart_type",
                        tooltip_label("Chart type", "Boxplot shows the full distribution. Bar plot shows mean + SE. Violin plot shows density shape. Scatter/dot plot shows only raw points. Density plot overlays distribution curves per group."),
                        choices = c("Boxplot" = "box", "Bar plot (mean + SE)" = "bar",
                                    "Violin plot" = "violin", "Scatter / dot plot" = "scatter",
                                    "Density plot" = "density"),
                        selected = "box")
          ),
          colourpicker::colourInput("axis_color",
                                     tooltip_label("Axis color", "Color of the axis lines, ticks, and text."),
                                     value = "#000000"),
          checkboxInput("hide_gridlines",
                        tooltip_label("Remove gridlines", "Removes the background gridlines for a cleaner look."),
                        value = FALSE),
          conditionalPanel(
            "input.analysis_mode != 'survival'",
            checkboxInput("log_y_axis",
                          tooltip_label("Log-transform Y axis", "Switches the y-axis to log10 scale. Values must be greater than 0."),
                          value = FALSE)
          ),
          sliderInput("axis_angle",
                      tooltip_label("X-axis label angle", "Rotate the x-axis text -- helpful for long category names."),
                      min = 0, max = 90, value = 0, step = 15),
          sliderInput("text_size",
                      tooltip_label("Plot text size", "Font size for titles, axis labels, and legends."),
                      min = 8, max = 24, value = 14, step = 1),
          textInput("x_axis_label",
                    tooltip_label("X-axis label (override)", "Leave blank to auto-use the column name."),
                    value = "", placeholder = "Column name by default"),
          textInput("y_axis_label",
                    tooltip_label("Y-axis label (override)", "Leave blank to auto-use the column name."),
                    value = "", placeholder = "Column name by default"),
          conditionalPanel(
            "input.analysis_mode == 'groups' && input.group_design != 'two'",
            checkboxInput("show_pvalues",
                          tooltip_label("Show p-value brackets (2-4 groups)", "Draws brackets with the p-value for every pairwise comparison."),
                          value = TRUE),
            checkboxInput("plot_all_cols",
                          tooltip_label("Plot ALL numeric columns (faceted)", "Shows every other numeric column as its own mini comparison plot."),
                          value = FALSE)
          )
        ),

        plotOutput("main_plot")
      )
    )
  )
)

# ---------------------------------------------------------
# SERVER
# ---------------------------------------------------------
server <- function(input, output, session) {

  # ---- Load data ----
  raw_data <- reactive({
    req(input$datafile)
    path <- input$datafile$datapath
    ext  <- tolower(tools::file_ext(input$datafile$name))

    df <- switch(
      ext,
      csv  = read.csv(path, header = input$header, stringsAsFactors = FALSE),
      xlsx = as.data.frame(read_excel(path, col_names = input$header)),
      xls  = as.data.frame(read_excel(path, col_names = input$header)),
      {
        showNotification("Unsupported file type. Please upload a CSV or Excel file.",
                         type = "error")
        NULL
      }
    )
    df
  })

  numeric_cols <- reactive({
    df <- raw_data()
    req(df)
    names(df)[sapply(df, function(x) suppressWarnings(!all(is.na(as.numeric(x)))))]
  })

  output$file_uploaded <- reactive({ if (!is.null(raw_data())) "yes" else "no" })
  outputOptions(output, "file_uploaded", suspendWhenHidden = FALSE)

  # The set of categories currently being colored/filled in the plot --
  # used to generate one color picker per category when "Custom" palette
  # is selected.
  current_group_levels <- reactive({
    tryCatch({
      if (input$analysis_mode == "groups") {
        if (input$group_design == "two") levels(analysis_data_2f()$factor2)
        else levels(analysis_data()$group)
      } else if (input$analysis_mode == "survival") {
        levels(surv_data()$group)
      } else if (input$analysis_mode == "titration") {
        levels(titration_data()$group)
      } else {
        character(0)
      }
    }, error = function(e) character(0))
  })

  output$custom_color_ui <- renderUI({
    levels_vec <- current_group_levels()
    if (length(levels_vec) == 0) return(helpText("Select your data/columns above to choose colors per group."))
    default_cols <- scales::hue_pal()(length(levels_vec))
    tagList(
      lapply(seq_along(levels_vec), function(i) {
        colourpicker::colourInput(
          paste0("custom_color_", make.names(levels_vec[i])),
          label = levels_vec[i],
          value = default_cols[i]
        )
      })
    )
  })

  # Builds the named color vector (level -> hex) from the current color
  # picker inputs, for use with scale_fill_manual()/scale_color_manual().
  get_custom_colors <- function(levels_vec) {
    if (length(levels_vec) == 0) return(NULL)
    cols <- sapply(levels_vec, function(lv) {
      id <- paste0("custom_color_", make.names(lv))
      val <- input[[id]]
      if (is.null(val) || !nzchar(val)) "#888888" else val
    })
    names(cols) <- levels_vec
    cols
  }

  output$data_preview <- renderDT({
    df <- raw_data()
    req(df)
    datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- Sidebar controls that depend on analysis mode ----
  output$mode_ui <- renderUI({
    df <- raw_data()
    if (is.null(df)) return(NULL)

    if (input$analysis_mode == "groups") {
      tagList(
        radioButtons(
          "group_design",
          tooltip_label("Design",
                         "One factor compares groups along a single variable. Two factors tests two variables plus whether they interact (two-way ANOVA)."),
          choices = c("One grouping factor" = "one",
                      "Two grouping factors (Two-way ANOVA)" = "two"),
          selected = "one"
        ),

        conditionalPanel(
          "input.group_design == 'one'",
          selectInput("group_col",
                      tooltip_label("Grouping column", "The column that defines which group each row belongs to (e.g. Treatment, Tissue, Genotype)."),
                      choices = names(df)),
          selectInput("value_col",
                      tooltip_label("Numeric value column", "The measurement you want to compare across groups."),
                      choices = numeric_cols()),
          radioButtons(
            "test_family",
            tooltip_label("Test type",
                           "Parametric tests (t-test/ANOVA) assume normally distributed data. Non-parametric tests make no such assumption -- check the Normality Check tab if unsure which to use."),
            choices = c("Parametric (t-test / ANOVA)" = "param",
                        "Non-parametric (Wilcoxon / Kruskal-Wallis)" = "nonparam"),
            selected = "param"
          ),
          conditionalPanel(
            "output.is_two_group == 'yes' && input.test_family == 'param'",
            checkboxInput("paired",
                          tooltip_label("Paired samples", "Check this if the same subjects were measured in both groups (e.g. before/after treatment)."),
                          value = FALSE),
            checkboxInput("equal_var",
                          tooltip_label("Assume equal variances", "Check this if you believe both groups have similar spread/variability. Leave unchecked for the safer Welch's t-test."),
                          value = FALSE)
          ),
          conditionalPanel(
            "output.is_two_group == 'yes' && input.test_family == 'nonparam'",
            checkboxInput("paired_np",
                          tooltip_label("Paired samples (signed-rank test)", "Check this if the same subjects were measured in both groups."),
                          value = FALSE)
          ),
          conditionalPanel(
            "output.show_posthoc == 'yes' && input.test_family == 'param'",
            selectInput("posthoc_method_param",
                        tooltip_label("Multiple comparison method",
                                       "How to adjust for testing multiple group pairs at once, to control false positives. Tukey HSD is the standard default."),
                        choices = c("Tukey HSD" = "tukey",
                                    "Pairwise t-test (Bonferroni)" = "bonf"),
                        selected = "tukey")
          ),
          conditionalPanel(
            "output.show_posthoc == 'yes' && input.test_family == 'nonparam'",
            selectInput("posthoc_method_nonparam",
                        tooltip_label("Multiple comparison method",
                                       "How to adjust for testing multiple group pairs at once. Dunn's test is the standard non-parametric choice."),
                        choices = c("Dunn's test (BH-adjusted)" = "dunn_bh",
                                    "Pairwise Wilcoxon (BH-adjusted)" = "wilcox_bh",
                                    "Pairwise Wilcoxon (Bonferroni)" = "wilcox_bonf"),
                        selected = "dunn_bh")
          )
        ),

        conditionalPanel(
          "input.group_design == 'two'",
          selectInput("factor1_col",
                      tooltip_label("Grouping factor 1", "The first categorical variable to test (e.g. Tissue)."),
                      choices = names(df)),
          selectInput("factor2_col",
                      tooltip_label("Grouping factor 2", "The second categorical variable to test (e.g. Treatment). The app also tests whether factor 1 and factor 2 interact."),
                      choices = names(df)),
          selectInput("value_col_2f",
                      tooltip_label("Numeric value column", "The measurement you want to compare across both factors."),
                      choices = numeric_cols()),
          helpText("Two-way ANOVA tests the main effect of each factor plus ",
                   "whether they interact. Parametric only (Tukey HSD post-hoc).")
        )
      )
    } else if (input$analysis_mode == "corr") {
      tagList(
        selectInput("x_col",
                    tooltip_label("X variable (numeric)", "The variable plotted on the horizontal axis."),
                    choices = numeric_cols()),
        selectInput("y_col",
                    tooltip_label("Y variable (numeric)", "The variable plotted on the vertical axis."),
                    choices = numeric_cols()),
        radioButtons(
          "corr_method",
          tooltip_label("Correlation method",
                         "Pearson measures linear correlation and assumes roughly normal data. Spearman measures rank-based (monotonic) correlation and is more robust to outliers/non-linear trends."),
          choices = c("Pearson (linear, parametric)" = "pearson",
                      "Spearman (rank-based, non-parametric)" = "spearman"),
          selected = "pearson"
        ),
        checkboxInput("show_regression",
                      tooltip_label("Fit simple linear regression line", "Fits a straight line through the data and reports its slope, intercept, and R-squared."),
                      value = TRUE)
      )
    } else if (input$analysis_mode == "survival") {
      tagList(
        selectInput("surv_time_col",
                    tooltip_label("Time column", "How long each subject was followed, in whatever time unit your data uses (e.g. days)."),
                    choices = numeric_cols()),
        selectInput("surv_status_col",
                    tooltip_label("Event/status column", "Indicates whether the event (e.g. death) occurred, or the subject was censored (survived to the end of follow-up / left the study)."),
                    choices = names(df)),
        helpText("Status should be coded 1 = event occurred, 0 = censored."),
        selectInput("surv_group_col",
                    tooltip_label("Grouping column (optional)", "If provided, plots a separate survival curve per group and runs a log-rank test comparing them."),
                    choices = c("None", names(df)), selected = "None")
      )
    } else {
      tagList(
        selectInput("dose_col",
                    tooltip_label("Dose / concentration column (X)", "The dose, concentration, or titration variable."),
                    choices = numeric_cols()),
        selectInput("response_col",
                    tooltip_label("Response column (Y)", "The measured response at each dose (e.g. % viability, signal intensity)."),
                    choices = numeric_cols()),
        checkboxInput("log_dose",
                      tooltip_label("Log10-transform dose (recommended)", "Standard practice for titrations where doses span orders of magnitude -- makes the curve shape and EC50 estimate more reliable."),
                      value = TRUE),
        selectInput("titration_group_col",
                    tooltip_label("Grouping column (optional, separate curve per group)", "If provided, fits and plots a separate dose-response curve for each group (e.g. comparing two drugs)."),
                    choices = c("None", names(df)), selected = "None")
      )
    }
  })

  # =========================================================
  # ONE-FACTOR group comparisons
  # =========================================================
  analysis_data <- reactive({
    req(input$analysis_mode == "groups", input$group_design == "one")
    df <- raw_data()
    req(df, input$group_col, input$value_col)

    out <- df %>%
      transmute(
        group = as.factor(.data[[input$group_col]]),
        value = suppressWarnings(as.numeric(.data[[input$value_col]]))
      ) %>%
      filter(!is.na(group), !is.na(value))

    validate(need(nrow(out) > 0,
                  "No valid numeric data found in the selected value column."))
    validate(need(nlevels(droplevels(out$group)) >= 2,
                  "Need at least 2 groups with data to run a comparison."))
    out$group <- droplevels(out$group)
    out
  })

  n_groups <- reactive({ nlevels(analysis_data()$group) })

  output$is_two_group <- reactive({ if (n_groups() == 2) "yes" else "no" })
  outputOptions(output, "is_two_group", suspendWhenHidden = FALSE)

  output$show_posthoc <- reactive({ if (isTRUE(n_groups() > 2)) "yes" else "no" })
  outputOptions(output, "show_posthoc", suspendWhenHidden = FALSE)

  current_posthoc_method <- reactive({
    if ((input$test_family %||% "param") == "param") input$posthoc_method_param %||% "tukey"
    else input$posthoc_method_nonparam %||% "dunn_bh"
  })

  desc_stats <- reactive({
    analysis_data() %>%
      group_by(group) %>%
      summarise(
        n = n(), mean = mean(value), sd = sd(value), se = sd / sqrt(n),
        median = median(value), min = min(value), max = max(value),
        .groups = "drop"
      ) %>%
      mutate(across(where(is.numeric) & !matches("^n$"), ~ round(.x, 4)))
  })

  output$desc_table <- renderDT({
    req(input$analysis_mode == "groups", input$group_design == "one")
    datatable(desc_stats(), options = list(dom = "t", paging = FALSE))
  })

  output$normality_table <- renderDT({
    req(input$analysis_mode == "groups", input$group_design == "one")
    ad <- analysis_data()
    res <- ad %>%
      group_by(group) %>%
      summarise(
        n = n(),
        shapiro_p = if (n() >= 3 && n() <= 5000) {
          tryCatch(shapiro.test(value)$p.value, error = function(e) NA_real_)
        } else NA_real_,
        .groups = "drop"
      ) %>%
      mutate(
        shapiro_p = round(shapiro_p, 4),
        interpretation = case_when(
          is.na(shapiro_p) ~ "Not enough data (need n >= 3)",
          shapiro_p < input$alpha ~ "Deviates from normal -> consider non-parametric",
          TRUE ~ "Consistent with normal distribution"
        )
      )
    datatable(res, options = list(dom = "t", paging = FALSE))
  })

  test_result <- reactive({
    ad <- analysis_data()
    family <- input$test_family %||% "param"
    tr <- run_group_test(ad, family, current_posthoc_method(),
                          isTRUE(input$paired), isTRUE(input$equal_var), isTRUE(input$paired_np))
    validate(need(!is.null(tr), "Could not run the test on this column (check group sizes)."))
    tr
  })

  output$test_summary <- renderPrint({
    tr <- test_result(); alpha <- input$alpha
    report_p <- function(p, label) {
      cat("\n---\n")
      if (p < alpha) cat(sprintf("Result: p = %.4g < alpha (%.3g) -> statistically significant %s.\n", p, alpha, label))
      else cat(sprintf("Result: p = %.4g >= alpha (%.3g) -> no statistically significant %s.\n", p, alpha, label))
    }
    if (tr$type == "t-test") { print(tr$result); report_p(tr$result$p.value, "difference between groups") }
    else if (tr$type == "wilcoxon") { print(tr$result); report_p(tr$result$p.value, "difference between groups") }
    else if (tr$type == "anova") {
      print(summary(tr$result))
      report_p(summary(tr$result)[[1]][["Pr(>F)"]][1], "difference among groups (see multiple comparisons tab)")
    } else if (tr$type == "kruskal") {
      print(tr$result)
      report_p(tr$result$p.value, "difference among groups (see multiple comparisons tab)")
    }
  })

  output$posthoc_table <- renderDT({
    tr <- test_result()
    pw <- extract_pairwise(tr, levels(analysis_data()$group))
    req(nrow(pw) > 0)
    pw$p <- round(pw$p, 4)
    names(pw) <- c("group1", "group2", "p_adj")
    datatable(pw, options = list(dom = "t", paging = FALSE))
  })

  # Pairwise p-values for the single-column plot's brackets
  pairwise_pvalues <- reactive({
    tr <- test_result()
    extract_pairwise(tr, levels(analysis_data()$group))
  })

  # =========================================================
  # TWO-FACTOR group comparisons (Two-way ANOVA)
  # =========================================================
  analysis_data_2f <- reactive({
    req(input$analysis_mode == "groups", input$group_design == "two")
    df <- raw_data()
    req(df, input$factor1_col, input$factor2_col, input$value_col_2f)

    out <- df %>%
      transmute(
        factor1 = as.factor(.data[[input$factor1_col]]),
        factor2 = as.factor(.data[[input$factor2_col]]),
        value   = suppressWarnings(as.numeric(.data[[input$value_col_2f]]))
      ) %>%
      filter(!is.na(factor1), !is.na(factor2), !is.na(value))

    validate(need(nrow(out) > 0, "No valid numeric data found in the selected value column."))
    out$factor1 <- droplevels(out$factor1)
    out$factor2 <- droplevels(out$factor2)
    validate(need(nlevels(out$factor1) >= 2 && nlevels(out$factor2) >= 2,
                  "Each grouping factor needs at least 2 levels for a two-way ANOVA."))
    out
  })

  two_way_fit <- reactive({
    ad <- analysis_data_2f()
    aov(value ~ factor1 * factor2, data = ad)
  })

  output$desc_table_2f <- renderDT({
    ad <- analysis_data_2f()
    out <- ad %>%
      group_by(factor1, factor2) %>%
      summarise(n = n(), mean = mean(value), sd = sd(value), se = sd / sqrt(n),
                median = median(value), min = min(value), max = max(value),
                .groups = "drop") %>%
      mutate(across(where(is.numeric) & !matches("^n$"), ~ round(.x, 4)))
    datatable(out, options = list(dom = "t", paging = FALSE))
  })

  output$normality_2f <- renderPrint({
    fit <- two_way_fit()
    resid_vals <- residuals(fit)
    if (length(resid_vals) >= 3 && length(resid_vals) <= 5000) {
      st <- shapiro.test(resid_vals)
      print(st)
      cat("\n---\n")
      if (st$p.value < input$alpha) {
        cat("Residuals deviate from normality -- interpret the ANOVA with some caution ",
            "(consider transforming the data, e.g. log-transform).\n")
      } else {
        cat("Residuals are consistent with a normal distribution.\n")
      }
    } else {
      cat("Not enough data points to run a normality check (need 3-5000 residuals).\n")
    }
  })

  output$test_summary_2f <- renderPrint({
    fit <- two_way_fit()
    s <- summary(fit)
    print(s)

    tbl <- s[[1]]
    pvals <- tbl[["Pr(>F)"]]
    terms <- trimws(rownames(tbl))
    alpha <- input$alpha

    cat("\n---\n")
    for (i in seq_along(terms)) {
      if (terms[i] == "Residuals" || is.na(pvals[i])) next
      label <- switch(terms[i],
        "factor1" = paste("Main effect of", input$factor1_col),
        "factor2" = paste("Main effect of", input$factor2_col),
        "factor1:factor2" = paste("Interaction between", input$factor1_col, "and", input$factor2_col),
        terms[i]
      )
      if (pvals[i] < alpha) {
        cat(sprintf("%s: p = %.4g < alpha (%.3g) -> statistically significant.\n", label, pvals[i], alpha))
      } else {
        cat(sprintf("%s: p = %.4g >= alpha (%.3g) -> not statistically significant.\n", label, pvals[i], alpha))
      }
    }
    cat("\nIf the interaction is significant, interpret the two main effects with caution --\n")
    cat("the effect of one factor depends on the level of the other. Check the Plot tab\n")
    cat("to visualize the interaction.\n")
  })

  output$posthoc_table_2f <- renderDT({
    fit <- two_way_fit()
    tk <- TukeyHSD(fit)

    all_terms <- lapply(names(tk), function(term_name) {
      d <- as.data.frame(tk[[term_name]])
      d$comparison <- rownames(d)
      d$term <- term_name
      d[, c("term", "comparison", "diff", "lwr", "upr", "p adj")]
    })
    out <- do.call(rbind, all_terms)
    out[, c("diff", "lwr", "upr", "p adj")] <- round(out[, c("diff", "lwr", "upr", "p adj")], 4)
    datatable(out, options = list(pageLength = 15, scrollX = TRUE))
  })

  # =========================================================
  # Correlation & regression
  # =========================================================
  corr_data <- reactive({
    req(input$analysis_mode == "corr", input$x_col, input$y_col)
    df <- raw_data()
    x <- suppressWarnings(as.numeric(df[[input$x_col]]))
    y <- suppressWarnings(as.numeric(df[[input$y_col]]))
    out <- data.frame(x = x, y = y)
    out <- out[complete.cases(out), ]
    validate(need(nrow(out) >= 3, "Need at least 3 complete numeric pairs to run a correlation."))
    out
  })

  output$corr_desc_table <- renderDT({
    req(input$analysis_mode == "corr", input$x_col, input$y_col)
    df <- raw_data()
    x <- suppressWarnings(as.numeric(df[[input$x_col]]))
    y <- suppressWarnings(as.numeric(df[[input$y_col]]))
    keep <- !is.na(x) & !is.na(y); x <- x[keep]; y <- y[keep]
    out <- data.frame(
      variable = c(input$x_col, input$y_col), n = c(length(x), length(y)),
      mean = round(c(mean(x), mean(y)), 4), sd = round(c(sd(x), sd(y)), 4),
      median = round(c(median(x), median(y)), 4),
      min = round(c(min(x), min(y)), 4), max = round(c(max(x), max(y)), 4)
    )
    datatable(out, options = list(dom = "t", paging = FALSE))
  })

  output$corr_summary <- renderPrint({
    cd <- corr_data(); method <- input$corr_method %||% "pearson"
    ct <- cor.test(cd$x, cd$y, method = method)
    print(ct)
    cat("\n---\n")
    if (ct$p.value < input$alpha) {
      cat(sprintf("Result: p = %.4g < alpha (%.3g) -> statistically significant %s correlation.\n",
                  ct$p.value, input$alpha, method))
    } else {
      cat(sprintf("Result: p = %.4g >= alpha (%.3g) -> no statistically significant %s correlation.\n",
                  ct$p.value, input$alpha, method))
    }
    if (isTRUE(input$show_regression)) {
      fit <- lm(y ~ x, data = cd); s <- summary(fit)
      cat("\n--- Simple linear regression (y ~ x) ---\n")
      cat(sprintf("Slope     : %.4f\n", coef(fit)[2]))
      cat(sprintf("Intercept : %.4f\n", coef(fit)[1]))
      cat(sprintf("R-squared : %.4f\n", s$r.squared))
      cat(sprintf("Model p-value : %.4g\n", pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)))
    }
  })

  # =========================================================
  # Survival analysis
  # =========================================================
  surv_data <- reactive({
    req(input$analysis_mode == "survival")
    df <- raw_data()
    req(df, input$surv_time_col, input$surv_status_col)

    time <- suppressWarnings(as.numeric(df[[input$surv_time_col]]))
    status <- suppressWarnings(as.numeric(df[[input$surv_status_col]]))
    group <- if (!is.null(input$surv_group_col) && input$surv_group_col != "None") {
      as.factor(df[[input$surv_group_col]])
    } else {
      factor(rep("All", nrow(df)))
    }

    out <- data.frame(time = time, status = status, group = group)
    out <- out[complete.cases(out[, c("time", "status")]), ]
    out$group <- droplevels(out$group)
    validate(need(nrow(out) > 0, "No valid survival data (check time/status columns)."))
    validate(need(all(out$status %in% c(0, 1)),
                  "Status column must be coded 0 (censored) or 1 (event)."))
    out
  })

  surv_fit <- reactive({
    sd <- surv_data()
    survival::survfit(survival::Surv(time, status) ~ group, data = sd)
  })

  output$surv_summary_table <- renderDT({
    sd <- surv_data()
    out <- sd %>%
      group_by(group) %>%
      group_modify(~ {
        fit_g <- survival::survfit(survival::Surv(time, status) ~ 1, data = .x)
        s <- summary(fit_g)$table
        med <- if ("median" %in% names(s)) unname(s["median"]) else NA
        data.frame(n = nrow(.x), events = sum(.x$status == 1), median_survival = round(med, 4))
      }) %>%
      ungroup()
    datatable(out, options = list(dom = "t", paging = FALSE))
  })

  output$surv_test_summary <- renderPrint({
    sd <- surv_data()
    fit <- surv_fit()
    cat("--- Kaplan-Meier fit summary table ---\n")
    print(summary(fit)$table)

    if (nlevels(sd$group) > 1) {
      test <- survival::survdiff(survival::Surv(time, status) ~ group, data = sd)
      cat("\n--- Log-rank test (comparing survival curves across groups) ---\n")
      print(test)
      p_val <- 1 - pchisq(test$chisq, length(test$n) - 1)
      cat("\n---\n")
      if (p_val < input$alpha) {
        cat(sprintf("Result: p = %.4g < alpha (%.3g) -> statistically significant difference in survival between groups.\n",
                    p_val, input$alpha))
      } else {
        cat(sprintf("Result: p = %.4g >= alpha (%.3g) -> no statistically significant difference in survival between groups.\n",
                    p_val, input$alpha))
      }
    } else {
      cat("\nNo grouping column selected -- showing a single survival curve (no comparison test run).\n")
    }
  })

  make_surv_plot <- function() {
    fit <- surv_fit()
    if (!is.null(fit$strata)) {
      strata_names <- rep(names(fit$strata), fit$strata)
      strata_names <- gsub("^group=", "", strata_names)
    } else {
      strata_names <- rep("All", length(fit$time))
    }
    sfit_df <- data.frame(time = fit$time, surv = fit$surv, strata = strata_names)
    starts <- sfit_df %>% group_by(strata) %>% slice(1) %>% mutate(time = 0, surv = 1)
    sfit_df <- bind_rows(starts, sfit_df) %>% arrange(strata, time)

    x_lab <- if (nzchar(input$x_axis_label)) input$x_axis_label else (input$surv_time_col %||% "Time")
    y_lab <- if (nzchar(input$y_axis_label)) input$y_axis_label else "Survival probability"

    surv_caption <- tryCatch({
      sd <- surv_data()
      if (nlevels(sd$group) > 1) {
        test <- survival::survdiff(survival::Surv(time, status) ~ group, data = sd)
        p <- 1 - pchisq(test$chisq, length(test$n) - 1)
        paste("Log-rank test:", format_p(p))
      } else ""
    }, error = function(e) "")

    ggplot(sfit_df, aes(x = time, y = surv, color = strata)) +
      geom_step(linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = x_lab, y = y_lab, color = input$surv_group_col %||% "Group",
           title = "Kaplan-Meier survival curve", caption = surv_caption) +
      theme_minimal(base_size = input$text_size) +
      get_color_scale(input$palette, get_custom_colors(unique(sfit_df$strata))) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines) +
      theme(plot.caption = element_text(size = rel(1.1), face = "bold", color = input$axis_color %||% "#000000"))
  }

  # =========================================================
  # Titration / dose-response
  # =========================================================
  titration_data <- reactive({
    req(input$analysis_mode == "titration")
    df <- raw_data()
    req(df, input$dose_col, input$response_col)

    dose <- suppressWarnings(as.numeric(df[[input$dose_col]]))
    response <- suppressWarnings(as.numeric(df[[input$response_col]]))
    group <- if (!is.null(input$titration_group_col) && input$titration_group_col != "None") {
      as.factor(df[[input$titration_group_col]])
    } else {
      factor(rep("All", nrow(df)))
    }

    out <- data.frame(dose = dose, response = response, group = group)
    out <- out[complete.cases(out[, c("dose", "response")]), ]

    if (isTRUE(input$log_dose)) {
      validate(need(all(out$dose > 0),
                    "Log10 transform requires all doses > 0. Uncheck the log option or remove non-positive doses."))
      out$x <- log10(out$dose)
    } else {
      out$x <- out$dose
    }
    out$group <- droplevels(out$group)
    validate(need(nrow(out) >= 5, "Need at least 5 valid data points to fit a dose-response curve."))
    out
  })

  fit_one_curve <- function(sub) {
    tryCatch(nls(response ~ SSfpl(x, A, B, xmid, scal), data = sub), error = function(e) NULL)
  }

  titration_fits <- reactive({
    td <- titration_data()
    td %>%
      group_by(group) %>%
      group_modify(~ {
        sub <- .x
        fit <- fit_one_curve(sub)
        if (is.null(fit)) {
          return(data.frame(A = NA, B = NA, xmid = NA, scal = NA, EC50 = NA, R2 = NA,
                             note = "Fit failed - check data"))
        }
        co <- coef(fit)
        pred <- predict(fit)
        ss_res <- sum((sub$response - pred)^2)
        ss_tot <- sum((sub$response - mean(sub$response))^2)
        r2 <- 1 - ss_res / ss_tot
        ec50 <- if (isTRUE(input$log_dose)) 10 ^ co["xmid"] else co["xmid"]
        data.frame(A = round(co["A"], 4), B = round(co["B"], 4), xmid = round(co["xmid"], 4),
                   scal = round(co["scal"], 4), EC50 = round(unname(ec50), 4), R2 = round(r2, 4),
                   note = "OK")
      }) %>%
      ungroup()
  })

  output$titration_fit_table <- renderDT({
    datatable(titration_fits(), options = list(dom = "t", paging = FALSE))
  })

  output$titration_summary <- renderPrint({
    cat("Four-parameter logistic fit: response = A + (B - A) / (1 + exp((xmid - x) / scal))\n")
    cat("  A = bottom asymptote, B = top asymptote\n")
    cat("  xmid = inflection point", if (isTRUE(input$log_dose)) "(in log10(dose) units)" else "(in dose units)", "\n")
    cat("  scal = slope factor (smaller = steeper curve)\n")
    cat("  EC50 = dose at the inflection point (back-transformed to original dose units)\n\n")
    print(titration_fits())
  })

  make_titration_plot <- function() {
    td <- titration_data()

    pred_list <- td %>%
      group_by(group) %>%
      group_modify(~ {
        sub <- .x
        fit <- fit_one_curve(sub)
        if (is.null(fit)) return(data.frame(x = numeric(), response = numeric()))
        xs <- seq(min(sub$x), max(sub$x), length.out = 100)
        data.frame(x = xs, response = as.numeric(predict(fit, newdata = data.frame(x = xs))))
      }) %>%
      ungroup()

    x_lab <- if (nzchar(input$x_axis_label)) input$x_axis_label else
      (if (isTRUE(input$log_dose)) paste0("log10(", input$dose_col, ")") else (input$dose_col %||% "Dose"))
    y_lab <- if (nzchar(input$y_axis_label)) input$y_axis_label else (input$response_col %||% "Response")

    ggplot(td, aes(x = x, y = response, color = group)) +
      geom_point(size = input$point_size %||% 2, alpha = 0.6) +
      geom_line(data = pred_list, aes(x = x, y = response, color = group), linewidth = 1) +
      labs(x = x_lab, y = y_lab, color = input$titration_group_col %||% "Group",
           title = "Titration / dose-response curve") +
      theme_minimal(base_size = input$text_size) +
      get_color_scale(input$palette, get_custom_colors(levels(td$group))) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines) +
      log_y_scale(input$log_y_axis)
  }

  # =========================================================
  # Multi-column faceted plot (all numeric columns at once)
  # =========================================================
  multi_col_data <- reactive({
    req(input$analysis_mode == "groups", input$group_design == "one",
        isTRUE(input$plot_all_cols))
    df <- raw_data()
    req(df, input$group_col)

    value_cols <- setdiff(numeric_cols(), input$group_col)
    validate(need(length(value_cols) > 0, "No other numeric columns found to plot."))

    out <- df %>%
      select(all_of(c(input$group_col, value_cols))) %>%
      rename(group = all_of(input$group_col)) %>%
      mutate(group = as.factor(group)) %>%
      pivot_longer(cols = -group, names_to = "colname", values_to = "value") %>%
      mutate(value = suppressWarnings(as.numeric(value))) %>%
      filter(!is.na(group), !is.na(value))

    validate(need(nrow(out) > 0, "No valid numeric data found across the other columns."))
    out
  })

  multi_col_annotations <- reactive({
    df <- multi_col_data()
    family <- input$test_family %||% "param"
    method <- current_posthoc_method()

    df %>%
      group_by(colname) %>%
      group_modify(~ {
        sub <- .x
        sub$group <- droplevels(sub$group)
        tr <- tryCatch(
          run_group_test(sub, family, method, isTRUE(input$paired), isTRUE(input$equal_var), isTRUE(input$paired_np)),
          error = function(e) NULL
        )
        pw <- extract_pairwise(tr, levels(sub$group))
        if (nrow(pw) == 0) return(data.frame())

        rng  <- max(sub$value, na.rm = TRUE) - min(sub$value, na.rm = TRUE)
        base <- max(sub$value, na.rm = TRUE)
        step <- if (rng > 0) rng * 0.12 else abs(base) * 0.12 + 0.01

        pw$xmin  <- match(pw$group1, levels(sub$group))
        pw$xmax  <- match(pw$group2, levels(sub$group))
        pw$y_position <- base + step * seq_len(nrow(pw))
        pw$label <- sapply(pw$p, format_p)
        pw
      }) %>%
      ungroup()
  })

  make_multi_col_plot <- function() {
    df <- multi_col_data()
    x_lab <- if (nzchar(input$x_axis_label)) input$x_axis_label else (input$group_col %||% "Group")
    y_lab <- if (nzchar(input$y_axis_label)) input$y_axis_label else "Value"
    pt_size <- input$point_size %||% 1.4
    bw <- input$box_width %||% 0.6

    chart_type <- input$chart_type %||% "box"
    ds_multi <- df %>% group_by(colname, group) %>%
      summarise(mean = mean(value), se = sd(value) / sqrt(n()), .groups = "drop")

    if (chart_type == "bar") {
      p <- ggplot(df, aes(x = group, y = value, fill = group)) +
        geom_col(data = ds_multi, aes(x = group, y = mean, fill = group), alpha = 0.7, width = bw) +
        geom_errorbar(data = ds_multi, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                      inherit.aes = FALSE, width = 0.15, color = "#C0392B", linewidth = 0.6) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.6, color = "#2C3E50") +
        facet_wrap(~ colname, scales = "free_y") +
        labs(x = x_lab, y = y_lab, title = "All numeric columns, compared by group")
    } else if (chart_type == "violin") {
      p <- ggplot(df, aes(x = group, y = value, fill = group)) +
        geom_violin(alpha = 0.5, width = bw, trim = FALSE) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.6, color = "#2C3E50") +
        facet_wrap(~ colname, scales = "free_y") +
        labs(x = x_lab, y = y_lab, title = "All numeric columns, compared by group")
    } else if (chart_type == "scatter") {
      p <- ggplot(df, aes(x = group, y = value, fill = group, color = group)) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.8) +
        geom_errorbar(data = ds_multi, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                      inherit.aes = FALSE, width = 0.1, color = "#C0392B", linewidth = 0.6) +
        facet_wrap(~ colname, scales = "free_y") +
        labs(x = x_lab, y = y_lab, title = "All numeric columns, compared by group")
    } else if (chart_type == "density") {
      p <- ggplot(df, aes(x = value, fill = group, color = group)) +
        geom_density(alpha = 0.35, linewidth = 0.7) +
        facet_wrap(~ colname, scales = "free") +
        labs(x = y_lab, y = "Density", fill = x_lab, color = x_lab,
             title = "All numeric columns, density by group")
    } else {
      p <- ggplot(df, aes(x = group, y = value, fill = group)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA, width = bw) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.6, color = "#2C3E50") +
        facet_wrap(~ colname, scales = "free_y") +
        labs(x = x_lab, y = y_lab, title = "All numeric columns, compared by group")
    }

    p <- p +
      theme_minimal(base_size = input$text_size) +
      get_fill_scale(input$palette, get_custom_colors(levels(droplevels(df$group)))) +
      get_color_scale(input$palette, get_custom_colors(levels(droplevels(df$group)))) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines)

    if (chart_type == "density") {
      p <- p + theme(legend.position = "right")
      if (isTRUE(input$log_y_axis)) p <- p + scale_x_log10()
    } else {
      p <- p +
        theme(legend.position = "none") +
        scale_x_discrete(expand = expansion(mult = input$category_spacing %||% 0.1)) +
        log_y_scale(input$log_y_axis)
    }

    n <- n_groups()
    if (chart_type != "density" && isTRUE(input$show_pvalues) && n >= 2 && n <= 4) {
      ann <- multi_col_annotations()
      if (nrow(ann) > 0) {
        p <- p + geom_signif(
          data = ann,
          aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y_position),
          manual = TRUE, tip_length = 0.01,
          textsize = max(2.2, input$text_size / 4.5), vjust = -0.2
        )
      }
    }
    p
  }

  # =========================================================
  # One-factor group comparison plot
  # =========================================================
  make_groups_plot <- function() {
    ad <- analysis_data(); ds <- desc_stats()
    y_lab <- if (nzchar(input$y_axis_label)) input$y_axis_label else (input$value_col %||% "Value")
    x_lab <- if (nzchar(input$x_axis_label)) input$x_axis_label else (input$group_col %||% "Group")
    pt_size <- input$point_size %||% 2
    bw <- input$box_width %||% 0.6

    chart_type <- input$chart_type %||% "box"

    overall_p_caption <- tryCatch({
      tr0 <- test_result()
      if (tr0$type == "t-test") paste("t-test:", format_p(tr0$result$p.value))
      else if (tr0$type == "wilcoxon") paste("Wilcoxon:", format_p(tr0$result$p.value))
      else if (tr0$type == "anova") paste("ANOVA (overall):", format_p(summary(tr0$result)[[1]][["Pr(>F)"]][1]))
      else if (tr0$type == "kruskal") paste("Kruskal-Wallis (overall):", format_p(tr0$result$p.value))
      else ""
    }, error = function(e) "")

    if (chart_type == "bar") {
      p <- ggplot(ad, aes(x = group, y = value, fill = group)) +
        geom_col(data = ds, aes(x = group, y = mean, fill = group), inherit.aes = FALSE, alpha = 0.7, width = bw) +
        geom_errorbar(data = ds, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                      inherit.aes = FALSE, width = 0.15, color = "#C0392B", linewidth = 0.8) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.7, color = "#2C3E50") +
        labs(x = x_lab, y = y_lab, title = "Group comparison",
             subtitle = "Bar height = mean, error bar = SE, points = individual data")
    } else if (chart_type == "violin") {
      p <- ggplot(ad, aes(x = group, y = value, fill = group)) +
        geom_violin(alpha = 0.5, width = bw, trim = FALSE) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.7, color = "#2C3E50") +
        geom_errorbar(data = ds, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                      inherit.aes = FALSE, width = 0.1, color = "#C0392B", linewidth = 0.8) +
        geom_point(data = ds, aes(x = group, y = mean), inherit.aes = FALSE,
                   color = "#C0392B", size = 3, shape = 18) +
        labs(x = x_lab, y = y_lab, title = "Group comparison",
             subtitle = "Violin plot (shows data density) with individual points; red diamond = mean +/- SE")
    } else if (chart_type == "scatter") {
      p <- ggplot(ad, aes(x = group, y = value, fill = group, color = group)) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.8) +
        geom_errorbar(data = ds, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                      inherit.aes = FALSE, width = 0.1, color = "#C0392B", linewidth = 0.8) +
        geom_point(data = ds, aes(x = group, y = mean), inherit.aes = FALSE,
                   color = "#C0392B", size = 3, shape = 18) +
        labs(x = x_lab, y = y_lab, title = "Group comparison",
             subtitle = "Scatter/dot plot of individual data points; red diamond = mean +/- SE")
    } else if (chart_type == "density") {
      # Density plots put the measured value on the x-axis and overlay a
      # distribution curve per group -- a fundamentally different layout
      # from the box/bar/violin/scatter options above (which put the
      # group on the x-axis).
      p <- ggplot(ad, aes(x = value, fill = group, color = group)) +
        geom_density(alpha = 0.35, linewidth = 0.8) +
        labs(x = y_lab, y = "Density", fill = x_lab, color = x_lab,
             title = "Group comparison", subtitle = "Overlaid distribution density per group")
    } else {
      p <- ggplot(ad, aes(x = group, y = value, fill = group)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA, width = bw) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.7, color = "#2C3E50") +
        geom_errorbar(data = ds, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                      inherit.aes = FALSE, width = 0.15, color = "#C0392B", linewidth = 0.8) +
        geom_point(data = ds, aes(x = group, y = mean), inherit.aes = FALSE,
                   color = "#C0392B", size = 3, shape = 18) +
        labs(x = x_lab, y = y_lab, title = "Group comparison",
             subtitle = "Boxplot with individual points; red diamond = mean +/- SE")
    }

    p <- p +
      labs(caption = overall_p_caption) +
      theme_minimal(base_size = input$text_size) +
      get_fill_scale(input$palette, get_custom_colors(levels(ad$group))) +
      get_color_scale(input$palette, get_custom_colors(levels(ad$group))) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines) +
      theme(plot.caption = element_text(size = rel(1.05), face = "bold", color = input$axis_color %||% "#000000"))

    if (chart_type == "density") {
      p <- p + theme(legend.position = "right")
      if (isTRUE(input$log_y_axis)) p <- p + scale_x_log10()  # log the measured value, now on the x-axis
    } else {
      p <- p +
        theme(legend.position = "none") +
        scale_x_discrete(expand = expansion(mult = input$category_spacing %||% 0.1)) +
        log_y_scale(input$log_y_axis)
    }

    n <- n_groups()
    if (chart_type != "density" && isTRUE(input$show_pvalues) && n >= 2 && n <= 4) {
      pw <- pairwise_pvalues()
      if (nrow(pw) > 0) {
        comparisons <- Map(c, pw$group1, pw$group2)
        annotations <- sapply(pw$p, format_p)

        data_range <- max(ad$value, na.rm = TRUE) - min(ad$value, na.rm = TRUE)
        base_y <- max(ad$value, na.rm = TRUE)
        step <- data_range * 0.12
        y_positions <- base_y + step * seq_len(nrow(pw))

        p <- p + geom_signif(
          comparisons = comparisons,
          annotations = annotations,
          y_position  = y_positions,
          tip_length  = 0.01,
          textsize    = max(2.2, input$text_size / 4),
          vjust       = -0.2
        )
      }
    }
    p
  }

  make_groups_plot_2f <- function() {
    ad <- analysis_data_2f()
    y_lab <- if (nzchar(input$y_axis_label)) input$y_axis_label else (input$value_col_2f %||% "Value")
    x_lab <- if (nzchar(input$x_axis_label)) input$x_axis_label else (input$factor1_col %||% "Factor 1")
    pt_size <- input$point_size %||% 2
    bw <- input$box_width %||% 0.6

    chart_type <- input$chart_type %||% "box"
    ds2 <- ad %>% group_by(factor1, factor2) %>%
      summarise(mean = mean(value), se = sd(value) / sqrt(n()), .groups = "drop")

    twoway_p_caption <- tryCatch({
      s2 <- summary(two_way_fit())[[1]]
      sprintf("%s: %s | %s: %s | Interaction: %s",
              input$factor1_col %||% "Factor 1", format_p(s2[["Pr(>F)"]][1]),
              input$factor2_col %||% "Factor 2", format_p(s2[["Pr(>F)"]][2]),
              format_p(s2[["Pr(>F)"]][3]))
    }, error = function(e) "")

    if (chart_type == "bar") {
      p <- ggplot(ad, aes(x = factor1, y = value, fill = factor2)) +
        geom_col(data = ds2, aes(x = factor1, y = mean, fill = factor2),
                 position = position_dodge(width = 0.75), width = bw, alpha = 0.7) +
        geom_errorbar(data = ds2, aes(x = factor1, y = mean, ymin = mean - se, ymax = mean + se, group = factor2),
                      position = position_dodge(width = 0.75), width = 0.15, color = "#C0392B", linewidth = 0.7) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                   size = pt_size, alpha = 0.6, color = "#2C3E50") +
        labs(x = x_lab, y = y_lab, fill = input$factor2_col %||% "Factor 2",
             title = "Two-way comparison", subtitle = "Bar height = mean, error bar = SE")
    } else if (chart_type == "violin") {
      p <- ggplot(ad, aes(x = factor1, y = value, fill = factor2)) +
        geom_violin(alpha = 0.5, position = position_dodge(width = 0.75), width = bw, trim = FALSE) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                   size = pt_size, alpha = 0.6, color = "#2C3E50") +
        labs(x = x_lab, y = y_lab, fill = input$factor2_col %||% "Factor 2",
             title = "Two-way comparison", subtitle = "Violin plot grouped by factor 1, colored by factor 2")
    } else if (chart_type == "scatter") {
      p <- ggplot(ad, aes(x = factor1, y = value, fill = factor2, color = factor2)) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                   size = pt_size, alpha = 0.8) +
        geom_errorbar(data = ds2, aes(x = factor1, y = mean, ymin = mean - se, ymax = mean + se, group = factor2),
                      position = position_dodge(width = 0.75), width = 0.1, color = "#C0392B", linewidth = 0.7) +
        labs(x = x_lab, y = y_lab, fill = input$factor2_col %||% "Factor 2", color = input$factor2_col %||% "Factor 2",
             title = "Two-way comparison", subtitle = "Scatter/dot plot grouped by factor 1, colored by factor 2")
    } else if (chart_type == "density") {
      p <- ggplot(ad, aes(x = value, fill = factor2, color = factor2)) +
        geom_density(alpha = 0.35, linewidth = 0.8) +
        facet_wrap(~ factor1, scales = "free") +
        labs(x = y_lab, y = "Density", fill = input$factor2_col %||% "Factor 2",
             color = input$factor2_col %||% "Factor 2",
             title = "Two-way comparison", subtitle = "Density per factor 2 level, faceted by factor 1")
    } else {
      p <- ggplot(ad, aes(x = factor1, y = value, fill = factor2)) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75), width = bw) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                   size = pt_size, alpha = 0.6, color = "#2C3E50") +
        labs(x = x_lab, y = y_lab, fill = input$factor2_col %||% "Factor 2",
             title = "Two-way comparison", subtitle = "Grouped by factor 1, colored by factor 2")
    }

    p <- p +
      labs(caption = twoway_p_caption) +
      theme_minimal(base_size = input$text_size) +
      get_fill_scale(input$palette, get_custom_colors(levels(ad$factor2))) +
      get_color_scale(input$palette, get_custom_colors(levels(ad$factor2))) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines) +
      theme(plot.caption = element_text(size = rel(1.0), face = "bold", color = input$axis_color %||% "#000000"))

    if (chart_type == "density") {
      if (isTRUE(input$log_y_axis)) p <- p + scale_x_log10()
      p
    } else {
      p <- p +
        scale_x_discrete(expand = expansion(mult = input$category_spacing %||% 0.1)) +
        log_y_scale(input$log_y_axis)

      n2 <- nlevels(ad$factor2)
      if (isTRUE(input$show_pvalues) && n2 >= 2 && n2 <= 4) {
        ann2 <- twoway_bracket_df(ad, two_way_fit())
        if (nrow(ann2) > 0) {
          p <- p + geom_signif(
            data = ann2,
            aes(xmin = xmin, xmax = xmax, annotations = label, y_position = y_position),
            manual = TRUE, inherit.aes = FALSE, tip_length = 0.01,
            textsize = max(2.2, input$text_size / 4), vjust = -0.2
          )
        }
      }
      p
    }
  }

  make_corr_plot <- function() {
    cd <- corr_data()
    y_lab <- if (nzchar(input$y_axis_label)) input$y_axis_label else (input$y_col %||% "Y")
    x_lab <- if (nzchar(input$x_axis_label)) input$x_axis_label else (input$x_col %||% "X")

    corr_caption <- tryCatch({
      ct <- cor.test(cd$x, cd$y, method = input$corr_method %||% "pearson")
      sprintf("r = %.2f, %s", unname(ct$estimate), format_p(ct$p.value))
    }, error = function(e) "")

    p <- ggplot(cd, aes(x = x, y = y)) +
      geom_point(size = input$point_size %||% 2.5, alpha = 0.7, color = "#2C6E49") +
      labs(x = x_lab, y = y_lab, title = "Correlation / regression", caption = corr_caption) +
      theme_minimal(base_size = input$text_size) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines) +
      theme(plot.caption = element_text(size = rel(1.1), face = "bold", color = input$axis_color %||% "#000000")) +
      log_y_scale(input$log_y_axis)
    if (isTRUE(input$show_regression)) p <- p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#C0392B")
    p
  }

  make_plot <- function() {
    if (input$analysis_mode == "groups") {
      if (input$group_design == "two") make_groups_plot_2f()
      else if (isTRUE(input$plot_all_cols)) make_multi_col_plot()
      else make_groups_plot()
    } else if (input$analysis_mode == "corr") {
      make_corr_plot()
    } else if (input$analysis_mode == "survival") {
      make_surv_plot()
    } else {
      make_titration_plot()
    }
  }

  # =========================================================
  # Auto-generated plain-language summary (rule-based, no AI)
  # =========================================================
  generate_summary_groups_one <- function() {
    ad <- analysis_data(); tr <- test_result(); ds <- desc_stats(); alpha <- input$alpha
    val_col <- input$value_col %||% "the measurement"
    group_desc <- paste(sprintf("%s (n=%d)", ds$group, ds$n), collapse = ", ")
    n_grp <- n_groups()

    if (tr$type == "t-test") {
      test_name <- if (isTRUE(input$paired)) "paired t-test"
                   else if (isTRUE(input$equal_var)) "Student's t-test"
                   else "Welch's t-test"
      p <- tr$result$p.value
      sig <- if (p < alpha) "a statistically significant difference" else "no statistically significant difference"
      sprintf("A %s compared %s between %s. There was %s (%s).",
              test_name, val_col, group_desc, sig, fmt_p_clause(p))
    } else if (tr$type == "wilcoxon") {
      p <- tr$result$p.value
      sig <- if (p < alpha) "a statistically significant difference" else "no statistically significant difference"
      sprintf("A Wilcoxon rank-sum test compared %s between %s. There was %s (%s).",
              val_col, group_desc, sig, fmt_p_clause(p))
    } else if (tr$type == "anova") {
      fs <- summary(tr$result)[[1]]
      f_val <- fs[["F value"]][1]; p <- fs[["Pr(>F)"]][1]
      sig <- if (p < alpha) "a statistically significant difference" else "no statistically significant difference"
      base <- sprintf("A one-way ANOVA compared %s across %d groups (%s). There was %s among groups (F = %.2f, %s).",
                       val_col, n_grp, group_desc, sig, f_val, fmt_p_clause(p))
      if (p < alpha) {
        pw <- extract_pairwise(tr, levels(ad$group))
        sig_pairs <- pw[pw$p < alpha, ]
        method_lbl <- posthoc_label(current_posthoc_method())
        if (nrow(sig_pairs) > 0) {
          pair_strs <- sprintf("%s vs %s (%s)", sig_pairs$group1, sig_pairs$group2, sapply(sig_pairs$p, fmt_p_clause))
          base <- paste0(base, " Post-hoc comparisons (", method_lbl, ") showed significant differences for ",
                         paste(pair_strs, collapse = "; "), ".")
        } else {
          base <- paste0(base, " However, no individual pairwise comparisons reached significance after correction (",
                         method_lbl, ").")
        }
      }
      base
    } else if (tr$type == "kruskal") {
      p <- tr$result$p.value; stat <- unname(tr$result$statistic)
      sig <- if (p < alpha) "a statistically significant difference" else "no statistically significant difference"
      base <- sprintf("A Kruskal-Wallis test compared %s across %d groups (%s). There was %s among groups (chi-squared = %.2f, %s).",
                       val_col, n_grp, group_desc, sig, stat, fmt_p_clause(p))
      if (p < alpha) {
        pw <- extract_pairwise(tr, levels(ad$group))
        sig_pairs <- pw[pw$p < alpha, ]
        method_lbl <- posthoc_label(current_posthoc_method())
        if (nrow(sig_pairs) > 0) {
          pair_strs <- sprintf("%s vs %s (%s)", sig_pairs$group1, sig_pairs$group2, sapply(sig_pairs$p, fmt_p_clause))
          base <- paste0(base, " Post-hoc comparisons (", method_lbl, ") showed significant differences for ",
                         paste(pair_strs, collapse = "; "), ".")
        } else {
          base <- paste0(base, " However, no individual pairwise comparisons reached significance after correction (",
                         method_lbl, ").")
        }
      }
      base
    } else {
      "Could not generate a summary for this result."
    }
  }

  generate_summary_groups_two <- function() {
    fit <- two_way_fit(); s <- summary(fit)[[1]]; alpha <- input$alpha
    f1 <- input$factor1_col %||% "Factor 1"; f2 <- input$factor2_col %||% "Factor 2"
    val_col <- input$value_col_2f %||% "the measurement"
    p_f1 <- s[["Pr(>F)"]][1]; p_f2 <- s[["Pr(>F)"]][2]; p_int <- s[["Pr(>F)"]][3]
    sig_txt <- function(p) if (p < alpha) "a statistically significant" else "no statistically significant"

    base <- sprintf(
      "A two-way ANOVA examined the effects of %s and %s on %s. There was %s main effect of %s (%s), %s main effect of %s (%s), and %s interaction between %s and %s (%s).",
      f1, f2, val_col,
      sig_txt(p_f1), f1, fmt_p_clause(p_f1),
      sig_txt(p_f2), f2, fmt_p_clause(p_f2),
      sig_txt(p_int), f1, f2, fmt_p_clause(p_int)
    )
    if (p_int < alpha) {
      base <- paste0(base, sprintf(
        " Because the interaction was significant, the main effects of %s and %s should be interpreted with caution, as the effect of one factor depends on the level of the other.",
        f1, f2))
    }
    base
  }

  generate_summary_corr <- function() {
    cd <- corr_data(); method <- input$corr_method %||% "pearson"; alpha <- input$alpha
    ct <- cor.test(cd$x, cd$y, method = method)
    r <- unname(ct$estimate); p <- ct$p.value
    direction <- if (r > 0) "positive" else "negative"
    strength <- if (abs(r) < 0.3) "weak" else if (abs(r) < 0.6) "moderate" else "strong"
    sig_txt <- if (p < alpha) "statistically significant" else "not statistically significant"
    method_lbl <- if (method == "pearson") "Pearson" else "Spearman"

    base <- sprintf(
      "A %s correlation test examined the relationship between %s and %s (n = %d). There was a %s %s correlation (r = %.2f, %s), which was %s.",
      method_lbl, input$x_col %||% "X", input$y_col %||% "Y", nrow(cd), strength, direction, r, fmt_p_clause(p), sig_txt
    )
    if (isTRUE(input$show_regression)) {
      fit <- lm(y ~ x, data = cd); s <- summary(fit)
      base <- paste0(base, sprintf(
        " A linear regression indicated that %s explained %.1f%% of the variance in %s (R-squared = %.3f).",
        input$x_col %||% "X", s$r.squared * 100, input$y_col %||% "Y", s$r.squared))
    }
    base
  }

  generate_summary_survival <- function() {
    sd <- surv_data(); alpha <- input$alpha
    n_total <- nrow(sd); n_events <- sum(sd$status == 1)

    if (nlevels(sd$group) > 1) {
      test <- survival::survdiff(survival::Surv(time, status) ~ group, data = sd)
      p <- 1 - pchisq(test$chisq, length(test$n) - 1)
      sig_txt <- if (p < alpha) "a statistically significant" else "no statistically significant"
      grp_summary <- sd %>% group_by(group) %>% summarise(n = n(), events = sum(status == 1), .groups = "drop")
      grp_txt <- paste(sprintf("%s (n=%d, %d events)", grp_summary$group, grp_summary$n, grp_summary$events), collapse = ", ")
      sprintf("A Kaplan-Meier survival analysis compared %s. A log-rank test found %s difference in survival between groups (chi-squared = %.2f, %s).",
              grp_txt, sig_txt, test$chisq, fmt_p_clause(p))
    } else {
      sprintf("A Kaplan-Meier survival analysis was performed on %d subjects (%d events, %d censored).",
              n_total, n_events, n_total - n_events)
    }
  }

  generate_summary_titration <- function() {
    tf <- titration_fits()
    lines <- apply(tf, 1, function(row) {
      sprintf("%s (EC50 = %s, R-squared = %s)", row[["group"]], row[["EC50"]], row[["R2"]])
    })
    paste0("A four-parameter logistic curve was fit to the dose-response data for: ", paste(lines, collapse = "; "), ".")
  }

  generate_auto_summary <- function() {
    if (input$analysis_mode == "groups") {
      if (input$group_design == "two") generate_summary_groups_two() else generate_summary_groups_one()
    } else if (input$analysis_mode == "corr") {
      generate_summary_corr()
    } else if (input$analysis_mode == "survival") {
      generate_summary_survival()
    } else {
      generate_summary_titration()
    }
  }

  output$auto_summary <- renderText({ generate_auto_summary() })

  # =========================================================
  # Optional AI polish via the Claude API (requires an
  # ANTHROPIC_API_KEY environment variable set on the server --
  # never hardcode a key in this file, since it is in a public repo)
  # =========================================================
  output$claude_key_missing <- reactive({ if (nzchar(Sys.getenv("ANTHROPIC_API_KEY"))) "no" else "yes" })
  outputOptions(output, "claude_key_missing", suspendWhenHidden = FALSE)

  polished_text  <- reactiveVal(NULL)
  polish_error   <- reactiveVal(NULL)

  observeEvent(input$polish_btn, {
    polished_text(NULL)
    polish_error(NULL)

    api_key <- Sys.getenv("ANTHROPIC_API_KEY")
    if (!nzchar(api_key)) {
      polish_error("No Claude API key is configured for this app.")
      return()
    }

    summary_text <- tryCatch(generate_auto_summary(), error = function(e) NULL)
    if (is.null(summary_text)) {
      polish_error("Could not generate a summary to polish -- check the Test Results tab for errors first.")
      return()
    }

    result <- tryCatch({
      req <- httr2::request("https://api.anthropic.com/v1/messages") %>%
        httr2::req_headers(
          `x-api-key` = api_key,
          `anthropic-version` = "2023-06-01",
          `content-type` = "application/json"
        ) %>%
        httr2::req_body_json(list(
          model = "claude-haiku-5-5",
          max_tokens = 300,
          messages = list(list(
            role = "user",
            content = paste0(
              "Rewrite the following statistical results summary in clear, natural, polished prose suitable for ",
              "a lab notebook or a manuscript results section. Do not change, add, or remove any numbers, ",
              "p-values, or conclusions -- only improve the phrasing and flow. Keep it to one short paragraph ",
              "and return only the rewritten paragraph, nothing else.\n\n", summary_text
            )
          ))
        ))
      resp <- httr2::req_perform(req)
      body <- httr2::resp_body_json(resp)
      body$content[[1]]$text
    }, error = function(e) {
      polish_error(paste("Could not reach Claude:", conditionMessage(e)))
      NULL
    })

    if (!is.null(result)) polished_text(result)
  })

  output$polished_summary_ui <- renderUI({
    if (!is.null(polish_error())) {
      return(div(style = "color:#a33;", polish_error()))
    }
    if (!is.null(polished_text())) {
      return(div(polished_text()))
    }
    helpText("Click \"Polish with Claude\" above to generate a more natural-language version of the summary.")
  })

  output$main_plot <- renderPlot(
    { make_plot() },
    width  = function() input$plot_width_px %||% 800,
    height = function() input$plot_height_px %||% 500
  )

  # Converts the on-screen pixel width/height (at ~96 px/inch) to inches for
  # ggsave, so downloaded files match what's shown on screen. The "wide"
  # multi-column view gets extra width/height on top of that, since it has
  # many small facets that need the room regardless of the slider settings.
  export_dimensions <- function() {
    wide <- isTRUE(input$plot_all_cols) && input$analysis_mode == "groups" && input$group_design == "one"
    w <- (input$plot_width_px %||% 800) / 96 * (if (wide) 1.5 else 1)
    h <- (input$plot_height_px %||% 500) / 96 * (if (wide) 1.5 else 1)
    list(w = w, h = h)
  }

  output$download_plot <- downloadHandler(
    filename = function() "lab_stats_plot.png",
    content  = function(file) {
      d <- export_dimensions()
      ggsave(file, plot = make_plot(), width = d$w, height = d$h, dpi = 300)
    }
  )

  output$download_plot_svg <- downloadHandler(
    filename = function() "lab_stats_plot.svg",
    content  = function(file) {
      d <- export_dimensions()
      ggsave(file, plot = make_plot(), width = d$w, height = d$h, device = "svg")
    }
  )

  output$download_plot_eps <- downloadHandler(
    filename = function() "lab_stats_plot.eps",
    content  = function(file) {
      d <- export_dimensions()
      ggsave(file, plot = make_plot(), width = d$w, height = d$h, device = "eps")
    }
  )

  # =========================================================
  # PDF report (plot + descriptive stats + test results)
  # Reports the currently-selected single analysis; does not cover the
  # "plot all columns" faceted view (that one is PNG-only).
  # =========================================================
  report_title <- function() {
    if (input$analysis_mode == "groups") {
      if (input$group_design == "two") "Two-Way ANOVA Report" else "Group Comparison Report"
    } else if (input$analysis_mode == "corr") "Correlation & Regression Report"
    else if (input$analysis_mode == "survival") "Survival Analysis Report"
    else "Dose-Response Report"
  }

  report_desc_table <- function() {
    if (input$analysis_mode == "groups") {
      if (input$group_design == "two") as.data.frame(analysis_data_2f() %>%
          group_by(factor1, factor2) %>%
          summarise(n = n(), mean = round(mean(value), 4), sd = round(sd(value), 4),
                    .groups = "drop"))
      else as.data.frame(desc_stats())
    } else if (input$analysis_mode == "corr") {
      cd <- corr_data()
      data.frame(
        variable = c(input$x_col, input$y_col), n = c(nrow(cd), nrow(cd)),
        mean = round(c(mean(cd$x), mean(cd$y)), 4), sd = round(c(sd(cd$x), sd(cd$y)), 4)
      )
    } else if (input$analysis_mode == "survival") {
      sd <- surv_data()
      as.data.frame(sd %>% group_by(group) %>% group_modify(~ {
        fit_g <- survival::survfit(survival::Surv(time, status) ~ 1, data = .x)
        s <- summary(fit_g)$table
        med <- if ("median" %in% names(s)) unname(s["median"]) else NA
        data.frame(n = nrow(.x), events = sum(.x$status == 1), median_survival = round(med, 4))
      }) %>% ungroup())
    } else {
      as.data.frame(titration_fits())
    }
  }

  report_text_lines <- function() {
    capture.output({
      if (input$analysis_mode == "groups" && input$group_design == "one") {
        tr <- test_result(); alpha <- input$alpha
        if (tr$type == "t-test" || tr$type == "wilcoxon") print(tr$result)
        else if (tr$type == "anova") print(summary(tr$result))
        else if (tr$type == "kruskal") print(tr$result)
        pw <- extract_pairwise(tr, levels(analysis_data()$group))
        if (nrow(pw) > 0) {
          cat("\n--- Multiple comparisons ---\n")
          pw$p <- round(pw$p, 4)
          print(pw)
        }
      } else if (input$analysis_mode == "groups" && input$group_design == "two") {
        fit <- two_way_fit()
        print(summary(fit))
      } else if (input$analysis_mode == "corr") {
        cd <- corr_data(); method <- input$corr_method %||% "pearson"
        print(cor.test(cd$x, cd$y, method = method))
        if (isTRUE(input$show_regression)) {
          fit <- lm(y ~ x, data = cd); print(summary(fit))
        }
      } else if (input$analysis_mode == "survival") {
        sd <- surv_data(); fit <- surv_fit()
        print(summary(fit)$table)
        if (nlevels(sd$group) > 1) print(survival::survdiff(survival::Surv(time, status) ~ group, data = sd))
      } else {
        print(titration_fits())
      }
    })
  }

  output$download_report <- downloadHandler(
    filename = function() "lab_stats_report.pdf",
    content  = function(file) {
      pdf(file, width = 8.5, height = 11)

      # Page 1: the plot (opening pdf() already starts on a fresh page,
      # so no grid.newpage() is needed here -- calling it would skip to a
      # second blank page before drawing anything)
      print(make_plot())

      # Page 2: descriptive stats table (grid.arrange() advances to a new
      # page on its own by default, so no manual grid.newpage() here either)
      desc_df <- tryCatch(report_desc_table(), error = function(e) NULL)
      if (!is.null(desc_df) && nrow(desc_df) > 0) {
        grid.arrange(tableGrob(desc_df, rows = NULL), top = "Descriptive Statistics")
      }

      # Page 3+: test results text, paginated
      lines <- tryCatch(report_text_lines(), error = function(e) character())
      if (length(lines) > 0) {
        lines_per_page <- 55
        chunks <- split(lines, ceiling(seq_along(lines) / lines_per_page))
        for (chunk in chunks) {
          grid.newpage()
          grid.text(paste(chunk, collapse = "\n"),
                     x = unit(0.04, "npc"), y = unit(0.97, "npc"),
                     just = c("left", "top"),
                     gp = gpar(fontfamily = "mono", fontsize = 8))
        }
      }

      dev.off()
    }
  )
}

shinyApp(ui, server)
