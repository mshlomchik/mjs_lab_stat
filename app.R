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

options(shiny.maxRequestSize = 25 * 1024^2)  # 25 MB upload limit

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || a == "") b else a

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, digits = 3, format = "f"))
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

# ---------------------------------------------------------
# UI
# ---------------------------------------------------------
ui <- page_sidebar(
  title = "Lab Stats Explorer",
  theme = bs_theme(version = 5, primary = "#2C6E49", base_font = font_google("Inter"),
                    font_size_base = "0.85rem"),

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
    hr(),

    h5("Plot options"),
    selectInput(
      "palette",
      tooltip_label("Color palette", "The color scheme used for groups/lines in the plot. Choose \"Custom\" to pick each group's color yourself."),
      choices = c("Default", "Viridis", "Set1", "Set2", "Dark2", "Paired", "Pastel1", "Custom"),
      selected = "Default"
    ),
    conditionalPanel(
      "input.palette == 'Custom'",
      uiOutput("custom_color_ui")
    ),
    conditionalPanel(
      "input.analysis_mode == 'groups'",
      radioButtons("chart_type",
                   tooltip_label("Chart type", "Boxplot shows the full distribution (median, quartiles, outliers). Bar plot shows just the mean with an SE error bar -- a more traditional presentation style."),
                   choices = c("Boxplot" = "box", "Bar plot (mean + SE)" = "bar"),
                   selected = "box")
    ),
    sliderInput("point_size",
                tooltip_label("Point size", "Size of the individual data points/dots shown on the plot."),
                min = 0.5, max = 6, value = 2, step = 0.5),
    conditionalPanel(
      "input.analysis_mode == 'groups'",
      sliderInput("box_width",
                  tooltip_label("Box/Bar width", "Width of each box or bar. Smaller values add more gap between groups."),
                  min = 0.2, max = 0.9, value = 0.6, step = 0.05)
    ),
    colourpicker::colourInput("axis_color",
                               tooltip_label("Axis color", "Color of the axis lines, ticks, and text."),
                               value = "#000000"),
    checkboxInput("hide_gridlines",
                  tooltip_label("Remove gridlines", "Removes the background gridlines for a cleaner, more minimal look."),
                  value = FALSE),
    conditionalPanel(
      "input.analysis_mode != 'survival'",
      checkboxInput("log_y_axis",
                    tooltip_label("Log-transform Y axis", "Switches the y-axis to a log10 scale -- useful when your data spans a wide range (e.g. several orders of magnitude). Values must be greater than 0."),
                    value = FALSE)
    ),
    sliderInput("axis_angle",
                tooltip_label("X-axis label angle",
                               "Rotate the x-axis text -- helpful when category names are long and overlap each other."),
                min = 0, max = 90, value = 0, step = 15),
    sliderInput("text_size",
                tooltip_label("Plot text size", "Font size for titles, axis labels, and legends in the plot."),
                min = 8, max = 24, value = 14, step = 1),
    textInput("x_axis_label",
              tooltip_label("X-axis label (optional override)",
                             "Type your own x-axis label. Leave blank to auto-use the column name."),
              value = "", placeholder = "Leave blank to use the column name"),
    textInput("y_axis_label",
              tooltip_label("Y-axis label (optional override)",
                             "Type your own y-axis label. Leave blank to auto-use the column name."),
              value = "", placeholder = "Leave blank to use the column name"),
    conditionalPanel(
      "input.analysis_mode == 'groups' && input.group_design != 'two'",
      checkboxInput("show_pvalues",
                    tooltip_label("Show p-value brackets on plot (2-4 groups)",
                                   "Draws brackets with the p-value for every pairwise group comparison directly on the plot."),
                    value = TRUE),
      checkboxInput("plot_all_cols",
                    tooltip_label("Plot ALL numeric columns at once (faceted)",
                                   "Instead of one chosen column, shows every other numeric column in your file as its own mini comparison plot."),
                    value = FALSE)
    ),

    helpText("Upload a file, choose an analysis type, then pick the ",
             "relevant columns below.")
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
      "Plot",
      plotOutput("main_plot", height = "600px"),
      downloadButton("download_plot", "Download plot (PNG)"),
      downloadButton("download_plot_svg", "Download plot (SVG)"),
      downloadButton("download_report", "Download Report (PDF)")
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

    ggplot(sfit_df, aes(x = time, y = surv, color = strata)) +
      geom_step(linewidth = 1) +
      scale_y_continuous(limits = c(0, 1)) +
      labs(x = x_lab, y = y_lab, color = input$surv_group_col %||% "Group",
           title = "Kaplan-Meier survival curve") +
      theme_minimal(base_size = input$text_size) +
      get_color_scale(input$palette, get_custom_colors(unique(sfit_df$strata))) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines)
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

    if ((input$chart_type %||% "box") == "bar") {
      ds_multi <- df %>% group_by(colname, group) %>%
        summarise(mean = mean(value), se = sd(value) / sqrt(n()), .groups = "drop")
      p <- ggplot(df, aes(x = group, y = value, fill = group)) +
        geom_col(data = ds_multi, aes(x = group, y = mean, fill = group), alpha = 0.7, width = bw) +
        geom_errorbar(data = ds_multi, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                      inherit.aes = FALSE, width = 0.15, color = "#C0392B", linewidth = 0.6) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.6, color = "#2C3E50") +
        facet_wrap(~ colname, scales = "free_y") +
        labs(x = x_lab, y = y_lab, title = "All numeric columns, compared by group")
    } else {
      p <- ggplot(df, aes(x = group, y = value, fill = group)) +
        geom_boxplot(alpha = 0.5, outlier.shape = NA, width = bw) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.6, color = "#2C3E50") +
        facet_wrap(~ colname, scales = "free_y") +
        labs(x = x_lab, y = y_lab, title = "All numeric columns, compared by group")
    }

    p <- p +
      theme_minimal(base_size = input$text_size) +
      theme(legend.position = "none") +
      get_fill_scale(input$palette, get_custom_colors(levels(droplevels(df$group)))) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines) +
      log_y_scale(input$log_y_axis)

    n <- n_groups()
    if (isTRUE(input$show_pvalues) && n >= 2 && n <= 4) {
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

    if ((input$chart_type %||% "box") == "bar") {
      p <- ggplot(ad, aes(x = group, y = value, fill = group)) +
        geom_col(data = ds, aes(x = group, y = mean, fill = group), inherit.aes = FALSE, alpha = 0.7, width = bw) +
        geom_errorbar(data = ds, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                      inherit.aes = FALSE, width = 0.15, color = "#C0392B", linewidth = 0.8) +
        geom_jitter(width = 0.12, size = pt_size, alpha = 0.7, color = "#2C3E50") +
        labs(x = x_lab, y = y_lab, title = "Group comparison",
             subtitle = "Bar height = mean, error bar = SE, points = individual data")
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
      theme_minimal(base_size = input$text_size) +
      theme(legend.position = "none") +
      get_fill_scale(input$palette, get_custom_colors(levels(ad$group))) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines) +
      log_y_scale(input$log_y_axis)

    n <- n_groups()
    if (isTRUE(input$show_pvalues) && n >= 2 && n <= 4) {
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

    if ((input$chart_type %||% "box") == "bar") {
      ds2 <- ad %>% group_by(factor1, factor2) %>%
        summarise(mean = mean(value), se = sd(value) / sqrt(n()), .groups = "drop")
      p <- ggplot(ad, aes(x = factor1, y = value, fill = factor2)) +
        geom_col(data = ds2, aes(x = factor1, y = mean, fill = factor2),
                 position = position_dodge(width = 0.75), width = bw, alpha = 0.7) +
        geom_errorbar(data = ds2, aes(x = factor1, y = mean, ymin = mean - se, ymax = mean + se, group = factor2),
                      position = position_dodge(width = 0.75), width = 0.15, color = "#C0392B", linewidth = 0.7) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                   size = pt_size, alpha = 0.6, color = "#2C3E50") +
        labs(x = x_lab, y = y_lab, fill = input$factor2_col %||% "Factor 2",
             title = "Two-way comparison", subtitle = "Bar height = mean, error bar = SE")
    } else {
      p <- ggplot(ad, aes(x = factor1, y = value, fill = factor2)) +
        geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75), width = bw) +
        geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                   size = pt_size, alpha = 0.6, color = "#2C3E50") +
        labs(x = x_lab, y = y_lab, fill = input$factor2_col %||% "Factor 2",
             title = "Two-way comparison", subtitle = "Grouped by factor 1, colored by factor 2")
    }

    p +
      theme_minimal(base_size = input$text_size) +
      get_fill_scale(input$palette, get_custom_colors(levels(ad$factor2))) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines) +
      log_y_scale(input$log_y_axis)
  }

  make_corr_plot <- function() {
    cd <- corr_data()
    y_lab <- if (nzchar(input$y_axis_label)) input$y_axis_label else (input$y_col %||% "Y")
    x_lab <- if (nzchar(input$x_axis_label)) input$x_axis_label else (input$x_col %||% "X")

    p <- ggplot(cd, aes(x = x, y = y)) +
      geom_point(size = input$point_size %||% 2.5, alpha = 0.7, color = "#2C6E49") +
      labs(x = x_lab, y = y_lab, title = "Correlation / regression") +
      theme_minimal(base_size = input$text_size) +
      axis_angle_theme(input$axis_angle) +
      extra_style_theme(input$axis_color, input$hide_gridlines) +
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

  output$main_plot <- renderPlot({ make_plot() })

  output$download_plot <- downloadHandler(
    filename = function() "lab_stats_plot.png",
    content  = function(file) {
      wide <- isTRUE(input$plot_all_cols) && input$analysis_mode == "groups" && input$group_design == "one"
      w <- if (wide) 12 else 8
      h <- if (wide) 9 else 6
      ggsave(file, plot = make_plot(), width = w, height = h, dpi = 300)
    }
  )

  output$download_plot_svg <- downloadHandler(
    filename = function() "lab_stats_plot.svg",
    content  = function(file) {
      wide <- isTRUE(input$plot_all_cols) && input$analysis_mode == "groups" && input$group_design == "one"
      w <- if (wide) 12 else 8
      h <- if (wide) 9 else 6
      ggsave(file, plot = make_plot(), width = w, height = h, device = "svg")
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
