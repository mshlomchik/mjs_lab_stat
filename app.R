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
#       - Two grouping factors (Two-way ANOVA):
#           parametric only -> main effects + interaction + Tukey HSD
#  2. Correlation & regression -> Pearson/Spearman correlation and
#     simple linear regression between two numeric variables
#
# Plot customization: color palette and x-axis label angle, available
# for any plot.
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
    conditionalPanel(
      "input.analysis_mode == 'groups' && input.group_design != 'two'",
      checkboxInput("show_pvalues", "Show p-value brackets on plot (2-4 groups)", value = TRUE)
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
      plotOutput("main_plot", height = "500px"),
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

    if (n_groups() == 2) {
      if (family == "param") {
        res <- t.test(value ~ group, data = ad,
                       paired = isTRUE(input$paired), var.equal = isTRUE(input$equal_var))
        list(type = "t-test", result = res)
      } else {
        if (isTRUE(input$paired_np)) {
          g <- levels(ad$group)
          v1 <- ad$value[ad$group == g[1]]; v2 <- ad$value[ad$group == g[2]]
          validate(need(length(v1) == length(v2),
                        "Paired Wilcoxon requires equal sample sizes in both groups."))
          res <- wilcox.test(v1, v2, paired = TRUE)
        } else {
          res <- wilcox.test(value ~ group, data = ad)
        }
        list(type = "wilcoxon", result = res)
      }
    } else {
      if (family == "param") {
        fit <- aov(value ~ group, data = ad)
        list(type = "anova", result = fit, posthoc = TukeyHSD(fit))
      } else {
        res <- kruskal.test(value ~ group, data = ad)
        ph  <- pairwise.wilcox.test(ad$value, ad$group, p.adjust.method = "BH")
        list(type = "kruskal", result = res, posthoc = ph)
      }
    }
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

  # Pairwise p-values for plot brackets (works for 2, 3, or 4 groups)
  pairwise_pvalues <- reactive({
    tr <- test_result()
    if (tr$type %in% c("t-test", "wilcoxon")) {
      g <- levels(analysis_data()$group)
      data.frame(group1 = g[1], group2 = g[2], p = tr$result$p.value, stringsAsFactors = FALSE)
    } else if (tr$type == "anova") {
      tk <- as.data.frame(tr$posthoc$group)
      comp <- strsplit(rownames(tk), "-")
      data.frame(
        group1 = sapply(comp, `[`, 2),
        group2 = sapply(comp, `[`, 1),
        p = tk[["p adj"]],
        stringsAsFactors = FALSE
      )
    } else if (tr$type == "kruskal") {
      pm <- tr$posthoc$p.value
      out <- as.data.frame(as.table(pm), stringsAsFactors = FALSE)
      names(out) <- c("group2", "group1", "p")
      out <- out[!is.na(out$p), c("group1", "group2", "p")]
      out
    }
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
  # Plots
  # =========================================================
  make_groups_plot <- function() {
    ad <- analysis_data(); ds <- desc_stats()

    p <- ggplot(ad, aes(x = group, y = value, fill = group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.12, size = 2, alpha = 0.7, color = "#2C3E50") +
      geom_errorbar(data = ds, aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
                    inherit.aes = FALSE, width = 0.15, color = "#C0392B", linewidth = 0.8) +
      geom_point(data = ds, aes(x = group, y = mean), inherit.aes = FALSE,
                 color = "#C0392B", size = 3, shape = 18) +
      labs(x = input$group_col %||% "Group", y = input$value_col %||% "Value",
           title = "Group comparison",
           subtitle = "Boxplot with individual points; red diamond = mean +/- SE") +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none") +
      get_fill_scale(input$palette) +
      axis_angle_theme(input$axis_angle)

    n <- n_groups()
    if (isTRUE(input$show_pvalues) && n >= 2 && n <= 4) {
      pw <- pairwise_pvalues()
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
        textsize    = 3.6,
        vjust       = -0.2
      )
    }
    p
  }

  make_groups_plot_2f <- function() {
    ad <- analysis_data_2f()
    ggplot(ad, aes(x = factor1, y = value, fill = factor2)) +
      geom_boxplot(alpha = 0.6, outlier.shape = NA, position = position_dodge(width = 0.75)) +
      geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75),
                 size = 1.8, alpha = 0.6, color = "#2C3E50") +
      labs(x = input$factor1_col %||% "Factor 1", y = input$value_col_2f %||% "Value",
           fill = input$factor2_col %||% "Factor 2",
           title = "Two-way comparison",
           subtitle = "Grouped by factor 1, colored by factor 2") +
      theme_minimal(base_size = 14) +
      get_fill_scale(input$palette) +
      axis_angle_theme(input$axis_angle)
  }

  make_corr_plot <- function() {
    cd <- corr_data()
    p <- ggplot(cd, aes(x = x, y = y)) +
      geom_point(size = 2.5, alpha = 0.7, color = "#2C6E49") +
      labs(x = input$x_col %||% "X", y = input$y_col %||% "Y", title = "Correlation / regression") +
      theme_minimal(base_size = 14) +
      axis_angle_theme(input$axis_angle)
    if (isTRUE(input$show_regression)) p <- p + geom_smooth(method = "lm", formula = y ~ x, se = TRUE, color = "#C0392B")
    p
  }

  make_plot <- function() {
    if (input$analysis_mode == "groups") {
      if (input$group_design == "two") make_groups_plot_2f() else make_groups_plot()
    } else {
      make_corr_plot()
    }
  }

  output$main_plot <- renderPlot({ make_plot() })

  output$download_plot <- downloadHandler(
    filename = function() "lab_stats_plot.png",
    content  = function(file) ggsave(file, plot = make_plot(), width = 8, height = 6, dpi = 300)
  )
}

shinyApp(ui, server)
