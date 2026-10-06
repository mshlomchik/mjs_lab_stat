# ============================================================
# Lab Stats Explorer
# A Shiny app for flexible statistics on biological lab-experiment
# data (e.g. sorted cell populations, cell assays, enzyme kinetics).
#
# Analysis modes:
#  1. Compare groups
#       - One grouping factor:
#           parametric    -> t-test (2 groups) / one-way ANOVA + Tukey (3+)
#           non-parametric-> Wilcoxon (2 groups) / Kruskal-Wallis +
#                             pairwise Wilcoxon post-hoc (3+)
#           -> plot automatically shows p-value brackets for 2, 3, or 4
#              groups (all pairwise comparisons)
#           -> optional: plot ALL numeric columns at once, faceted,
#              each with its own p-value brackets
#       - Two grouping factors (Two-way ANOVA):
#           parametric only -> main effects + interaction + Tukey HSD
#  2. Correlation & regression -> Pearson/Spearman correlation and
#     simple linear regression between two numeric variables
#
# Plot customization: color palette, x-axis label angle, text size,
# and an optional manual y-axis label override.
#
# HOW TO RUN:
#   install.packages(c("shiny","bslib","readxl","dplyr","tidyr",
#                       "ggplot2","DT","ggsignif"))
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

options(shiny.maxRequestSize = 25 * 1024^2)  # 25 MB upload limit

`%||%` <- function(a, b) if (is.null(a) || length(a) == 0 || a == "") b else a

format_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.001) return("p < 0.001")
  paste0("p = ", formatC(p, digits = 3, format = "f"))
}

get_fill_scale <- function(palette) {
  switch(palette,
    "Viridis"  = scale_fill_viridis_d(),
    "Set1"     = scale_fill_brewer(palette = "Set1"),
    "Set2"     = scale_fill_brewer(palette = "Set2"),
    "Dark2"    = scale_fill_brewer(palette = "Dark2"),
    "Paired"   = scale_fill_brewer(palette = "Paired"),
    "Pastel1"  = scale_fill_brewer(palette = "Pastel1"),
    NULL  # "Default" -> ggplot default colors
  )
}

axis_angle_theme <- function(angle) {
  if (angle == 0) {
    theme(axis.text.x = element_text(angle = 0, hjust = 0.5))
  } else {
    theme(axis.text.x = element_text(angle = angle, hjust = 1, vjust = 1))
  }
}

# Run a group comparison test on an arbitrary (group, value) data frame.
# Returns NULL if the test can't be run (too few/many groups, errors, etc).
run_group_test <- function(ad, family, paired = FALSE, equal_var = FALSE, paired_np = FALSE) {
  ad$group <- droplevels(ad$group)
  n <- nlevels(ad$group)
  if (n < 2 || n > 4) return(NULL)

  if (n == 2) {
    if (family == "param") {
      res <- tryCatch(t.test(value ~ group, data = ad, paired = paired, var.equal = equal_var),
                       error = function(e) NULL)
      if (is.null(res)) return(NULL)
      list(type = "t-test", result = res)
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
      list(type = "wilcoxon", result = res)
    }
  } else {
    if (family == "param") {
      fit <- tryCatch(aov(value ~ group, data = ad), error = function(e) NULL)
      if (is.null(fit)) return(NULL)
      list(type = "anova", result = fit, posthoc = tryCatch(TukeyHSD(fit), error = function(e) NULL))
    } else {
      res <- tryCatch(kruskal.test(value ~ group, data = ad), error = function(e) NULL)
      if (is.null(res)) return(NULL)
      ph <- tryCatch(pairwise.wilcox.test(ad$value, ad$group, p.adjust.method = "BH"),
                      error = function(e) NULL)
      list(type = "kruskal", result = res, posthoc = ph)
    }
  }
}

# Extract pairwise group1/group2/p data frame from a run_group_test() result.
extract_pairwise <- function(tr, group_levels) {
  if (is.null(tr)) return(data.frame(group1 = character(), group2 = character(), p = numeric()))

  if (tr$type %in% c("t-test", "wilcoxon")) {
    data.frame(group1 = group_levels[1], group2 = group_levels[2], p = tr$result$p.value,
               stringsAsFactors = FALSE)
  } else if (tr$type == "anova" && !is.null(tr$posthoc)) {
    tk <- as.data.frame(tr$posthoc$group)
    comp <- strsplit(rownames(tk), "-")
    data.frame(group1 = sapply(comp, `[`, 2), group2 = sapply(comp, `[`, 1),
               p = tk[["p adj"]], stringsAsFactors = FALSE)
  } else if (tr$type == "kruskal" && !is.null(tr$posthoc)) {
    pm <- tr$posthoc$p.value
    out <- as.data.frame(as.table(pm), stringsAsFactors = FALSE)
    names(out) <- c("group2", "group1", "p")
    out[!is.na(out$p), c("group1", "group2", "p")]
  } else {
    data.frame(group1 = character(), group2 = character(), p = numeric())
  }
}

# ---------------------------------------------------------
# UI
# ---------------------------------------------------------
ui <- page_sidebar(
  title = "Lab Stats Explorer",
  theme = bs_theme(version = 5, primary = "#2C6E49", base_font = font_google("Inter")),

  sidebar = sidebar(
    width = 380,
    fileInput("datafile", "Upload data (CSV or Excel)",
              accept = c(".csv", ".xlsx", ".xls")),
    checkboxInput("header", "File has column headers", value = TRUE),
    hr(),

    radioButtons(
      "analysis_mode", "Analysis type",
      choices = c("Compare groups" = "groups",
                  "Correlation & regression" = "corr"),
      selected = "groups"
    ),
    hr(),

    uiOutput("mode_ui"),
    hr(),
    numericInput("alpha", "Significance level (alpha)", value = 0.05,
                 min = 0.001, max = 0.5, step = 0.01),
    hr(),

    h5("Plot options"),
    selectInput(
      "palette", "Color palette",
      choices = c("Default", "Viridis", "Set1", "Set2", "Dark2", "Paired", "Pastel1"),
      selected = "Default"
    ),
    sliderInput("axis_angle", "X-axis label angle", min = 0, max = 90, value = 0, step = 15),
    sliderInput("text_size", "Plot text size", min = 8, max = 24, value = 14, step = 1),
    textInput("x_axis_label", "X-axis label (optional override)", value = "",
              placeholder = "Leave blank to use the column name"),
    textInput("y_axis_label", "Y-axis label (optional override)", value = "",
              placeholder = "Leave blank to use the column name"),
    conditionalPanel(
      "input.analysis_mode == 'groups' && input.group_design != 'two'",
      checkboxInput("show_pvalues", "Show p-value brackets on plot (2-4 groups)", value = TRUE),
      checkboxInput("plot_all_cols", "Plot ALL numeric columns at once (faceted)", value = FALSE)
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
      conditionalPanel("input.analysis_mode == 'corr'", DTOutput("corr_desc_table"))
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
      )
    ),

    nav_panel(
      "Test Results",
      conditionalPanel(
        "input.analysis_mode == 'groups' && input.group_design != 'two'",
        verbatimTextOutput("test_summary"),
        conditionalPanel(
          "output.show_posthoc == 'yes'",
          h5("Post-hoc pairwise comparisons"),
          DTOutput("posthoc_table")
        )
      ),
      conditionalPanel(
        "input.analysis_mode == 'groups' && input.group_design == 'two'",
        verbatimTextOutput("test_summary_2f"),
        h5("Tukey HSD post-hoc (main effects + interaction)"),
        DTOutput("posthoc_table_2f")
      ),
      conditionalPanel(
        "input.analysis_mode == 'corr'",
        verbatimTextOutput("corr_summary")
      )
    ),

    nav_panel(
      "Plot",
      plotOutput("main_plot", height = "600px"),
      downloadButton("download_plot", "Download plot (PNG)")
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

  # ---- Sidebar controls that depend on analysis mode ----
  output$mode_ui <- renderUI({
    df <- raw_data()
    if (is.null(df)) return(helpText("Upload a file to get started."))

    if (input$analysis_mode == "groups") {
      tagList(
        radioButtons(
          "group_design", "Design",
          choices = c("One grouping factor" = "one",
                      "Two grouping factors (Two-way ANOVA)" = "two"),
          selected = "one"
        ),

        conditionalPanel(
          "input.group_design == 'one'",
          selectInput("group_col", "Grouping column", choices = names(df)),
          selectInput("value_col", "Numeric value column", choices = numeric_cols()),
          radioButtons(
            "test_family", "Test type",
            choices = c("Parametric (t-test / ANOVA)" = "param",
                        "Non-parametric (Wilcoxon / Kruskal-Wallis)" = "nonparam"),
            selected = "param"
          ),
          conditionalPanel(
            "output.is_two_group == 'yes' && input.test_family == 'param'",
            checkboxInput("paired", "Paired samples", value = FALSE),
            checkboxInput("equal_var", "Assume equal variances", value = FALSE)
          ),
          conditionalPanel(
            "output.is_two_group == 'yes' && input.test_family == 'nonparam'",
            checkboxInput("paired_np", "Paired samples (signed-rank test)", value = FALSE)
          )
        ),

        conditionalPanel(
          "input.group_design == 'two'",
          selectInput("factor1_col", "Grouping factor 1", choices = names(df)),
          selectInput("factor2_col", "Grouping factor 2", choices = names(df)),
          selectInput("value_col_2f", "Numeric value column", choices = numeric_cols()),
          helpText("Two-way ANOVA tests the main effect of each factor plus ",
                   "whether they interact. Parametric only.")
        )
      )
    } else {
      tagList(
        selectInput("x_col", "X variable (numeric)", choices = numeric_cols()),
        selectInput("y_col", "Y variable (numeric)", choices = numeric_cols()),
        radioButtons(
          "corr_method", "Correlation method",
          choices = c("Pearson (linear, parametric)" = "pearson",
                      "Spearman (rank-based, non-parametric)" = "spearman"),
          selected = "pearson"
        ),
        checkboxInput("show_regression", "Fit simple linear regression line", value = TRUE)
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
    tr <- run_group_test(ad, family, isTRUE(input$paired), isTRUE(input$equal_var), isTRUE(input$paired_np))
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
      report_p(summary(tr$result)[[1]][["Pr(>F)"]][1], "difference among groups (see post-hoc tab)")
    } else if (tr$type == "kruskal") {
      print(tr$result)
      report_p(tr$result$p.value, "difference among groups (see post-hoc tab)")
    }
  })

  output$posthoc_table <- renderDT({
    tr <- test_result(); req(!is.null(tr$posthoc))
    if (tr$type == "anova") {
      tk <- as.data.frame(tr$posthoc$group)
      tk$comparison <- rownames(tk)
      tk <- tk[, c("comparison", "diff", "lwr", "upr", "p adj")]
      tk[, c("diff", "lwr", "upr", "p adj")] <- round(tk[, c("diff", "lwr", "upr", "p adj")], 4)
      datatable(tk, options = list(dom = "t", paging = FALSE))
    } else if (tr$type == "kruskal") {
      pm <- tr$posthoc$p.value
      out <- as.data.frame(as.table(pm))
      names(out) <- c("group1", "group2", "p_adj")
      out <- out[!is.na(out$p_adj), ]
      out$p_adj <- round(out$p_adj, 4)
      datatable(out, options = list(dom = "t", paging = FALSE))
    }
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

    df %>%
      group_by(colname) %>%
      group_modify(~ {
        sub <- .x
        sub$group <- droplevels(sub$group)
        tr <- tryCatch(
          run_group_test(sub, family, isTRUE(input$paired), isTRUE(input$equal_var), isTRUE(input$paired_np)),
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

    p <- ggplot(df, aes(x = group, y = value, fill = group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.12, size = 1.4, alpha = 0.6, color = "#2C3E50") +
      facet_wrap(~ colname, scales = "free_y") +
      labs(x = (if (nzchar(input$x_axis_label)) input$x_axis_label else (input$group_col %||% "Group")),
           y = (if (nzchar(input$y_axis_label)) input$y_axis_label else "Value"),
           title = "All numeric columns, compared by group") +
      theme_minimal(base_size = input$text_size) +
      theme(legend.position = "none") +
      get_fill_scale(input$palette) +
      axis_angle_theme(input$axis_angle)

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
  # Plots
  # =========================================================
  make_groups_plot <- function() {
    ad <- analysis_data(); ds <- desc_stats()
    y_lab <- if (nzchar(input$y_axis_label)) input$y_axis_label else (input$value_col %||% "Value")
    x_lab <- if (nzchar(input$x_axis_label)) input$x_axis_label else (input$group_col %||% "Group")

    p <- ggplot(ad, aes(x = group, y = value, fill = group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.12, size = 2, alpha = 0.7, color = "#2C3E50") +
      geom_errorbar(data = ds, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                    inherit.aes = FALSE, width = 0.15, color = "#C0392B", linewidth = 0.8) +
      geom_point(data = ds, aes(x = group, y = mean), inherit.aes = FALSE,
                 color = "#C0392B", size = 3, shape = 18) +
      labs(x = x_lab, y = y_lab,
           title = "Group comparison",
           subtitle = "Boxplot with individual points; red diamond = mean +/- SE") +
      theme_minimal(base_size = input$text_size) +
      theme(legend.position = "none") +
      get_fill_scale(input$palette) +
      axis_angle_theme(input$axis_angle)

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

    ggplot(ad, aes(x = factor1, y = value, fill = factor2)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75)) +
      geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                 size = 1.8, alpha = 0.6, color = "#2C3E50") +
      labs(x = x_lab, y = y_lab,
           fill = input$factor2_col %||% "Factor 2",
           title = "Two-way comparison",
           subtitle = "Grouped by factor 1, colored by factor 2") +
      theme_minimal(base_size = input$text_size) +
      get_fill_scale(input$palette) +
      axis_angle_theme(input$axis_angle)
  }

  make_corr_plot <- function() {
    cd <- corr_data()
    y_lab <- if (nzchar(input$y_axis_label)) input$y_axis_label else (input$y_col %||% "Y")
    x_lab <- if (nzchar(input$x_axis_label)) input$x_axis_label else (input$x_col %||% "X")

    p <- ggplot(cd, aes(x = x, y = y)) +
      geom_point(size = 2.5, alpha = 0.7, color = "#2C6E49") +
      labs(x = x_lab, y = y_lab, title = "Correlation / regression") +
      theme_minimal(base_size = input$text_size) +
      axis_angle_theme(input$axis_angle)
    if (isTRUE(input$show_regression)) p <- p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#C0392B")
    p
  }

  make_plot <- function() {
    if (input$analysis_mode == "groups") {
      if (input$group_design == "two") {
        make_groups_plot_2f()
      } else if (isTRUE(input$plot_all_cols)) {
        make_multi_col_plot()
      } else {
        make_groups_plot()
      }
    } else {
      make_corr_plot()
    }
  }

  output$main_plot <- renderPlot({ make_plot() })

  output$download_plot <- downloadHandler(
    filename = function() "lab_stats_plot.png",
    content  = function(file) {
      w <- if (isTRUE(input$plot_all_cols) && input$analysis_mode == "groups" && input$group_design == "one") 12 else 8
      h <- if (isTRUE(input$plot_all_cols) && input$analysis_mode == "groups" && input$group_design == "one") 9 else 6
      ggsave(file, plot = make_plot(), width = w, height = h, dpi = 300)
    }
  )
}

shinyApp(ui, server)
