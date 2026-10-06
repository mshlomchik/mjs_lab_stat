# ============================================================
# Lab Stats Explorer
# A Shiny app for descriptive statistics and group comparisons
# (t-test / ANOVA) on biological lab-experiment data
# (e.g. cell assays, enzyme kinetics, dose-response readouts).
#
# HOW TO RUN:
#   1. install.packages(c("shiny","bslib","readxl","dplyr",
#                          "ggplot2","DT","tidyr","broom"))
#   2. shiny::runApp("app.R")
#
# INPUT DATA FORMAT (long format):
#   One column = grouping variable (e.g. "Treatment": Control, DrugA, DrugB)
#   One column = numeric measurement (e.g. "Absorbance", "Activity")
#   One row per replicate/sample.
# ============================================================

library(shiny)
library(bslib)
library(readxl)
library(dplyr)
library(tidyr)
library(ggplot2)
library(DT)

options(shiny.maxRequestSize = 25 * 1024^2)  # 25 MB upload limit

# ---------------------------------------------------------
# UI
# ---------------------------------------------------------
ui <- page_sidebar(
  title = "Lab Stats Explorer",
  theme = bs_theme(version = 5, primary = "#2C6E49", base_font = font_google("Inter")),

  sidebar = sidebar(
    width = 340,
    fileInput("datafile", "Upload data (CSV or Excel)",
              accept = c(".csv", ".xlsx", ".xls")),
    checkboxInput("header", "File has column headers", value = TRUE),
    hr(),
    uiOutput("col_selectors"),
    hr(),
    numericInput("alpha", "Significance level (alpha)", value = 0.05,
                 min = 0.001, max = 0.5, step = 0.01),
    conditionalPanel(
      condition = "output.is_two_group == 'yes'",
      checkboxInput("paired", "Paired samples", value = FALSE),
      checkboxInput("equal_var", "Assume equal variances", value = FALSE)
    ),
    hr(),
    helpText("Upload a file, then pick which column identifies your groups ",
             "(e.g. Treatment) and which column holds the numeric measurement ",
             "(e.g. Absorbance, Enzyme Activity).")
  ),

  navset_card_tab(
    nav_panel("Data Preview", DTOutput("data_preview")),
    nav_panel("Descriptive Stats", DTOutput("desc_table")),
    nav_panel(
      "Test Results",
      verbatimTextOutput("test_summary"),
      conditionalPanel(
        condition = "output.is_anova == 'yes'",
        h5("Tukey HSD post-hoc comparisons"),
        DTOutput("tukey_table")
      )
    ),
    nav_panel(
      "Plot",
      plotOutput("group_plot", height = "500px"),
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

  # ---- Dynamic column selectors ----
  output$col_selectors <- renderUI({
    df <- raw_data()
    if (is.null(df)) return(NULL)

    tagList(
      selectInput("group_col", "Grouping column", choices = names(df)),
      selectInput("value_col", "Numeric value column", choices = names(df))
    )
  })

  # ---- Cleaned analysis data: group column as factor, value column numeric ----
  analysis_data <- reactive({
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
    validate(need(nlevels(out$group) >= 2,
                  "Need at least 2 groups to run a comparison."))
    out
  })

  n_groups <- reactive({ nlevels(analysis_data()$group) })

  output$is_two_group <- reactive({ if (n_groups() == 2) "yes" else "no" })
  outputOptions(output, "is_two_group", suspendWhenHidden = FALSE)

  output$is_anova <- reactive({ if (n_groups() > 2) "yes" else "no" })
  outputOptions(output, "is_anova", suspendWhenHidden = FALSE)

  # ---- Data preview ----
  output$data_preview <- renderDT({
    df <- raw_data()
    req(df)
    datatable(df, options = list(pageLength = 10, scrollX = TRUE))
  })

  # ---- Descriptive statistics ----
  desc_stats <- reactive({
    analysis_data() %>%
      group_by(group) %>%
      summarise(
        n      = n(),
        mean   = mean(value),
        sd     = sd(value),
        se     = sd / sqrt(n),
        median = median(value),
        min    = min(value),
        max    = max(value),
        .groups = "drop"
      ) %>%
      mutate(across(where(is.numeric) & !matches("^n$"), ~ round(.x, 4)))
  })

  output$desc_table <- renderDT({
    datatable(desc_stats(), options = list(dom = "t", paging = FALSE))
  })

  # ---- Statistical test (t-test or ANOVA) ----
  test_result <- reactive({
    ad <- analysis_data()

    if (n_groups() == 2) {
      res <- t.test(
        value ~ group,
        data   = ad,
        paired = isTRUE(input$paired),
        var.equal = isTRUE(input$equal_var)
      )
      list(type = "t-test", result = res)
    } else {
      fit <- aov(value ~ group, data = ad)
      list(type = "anova", result = fit, tukey = TukeyHSD(fit))
    }
  })

  output$test_summary <- renderPrint({
    tr <- test_result()
    alpha <- input$alpha

    if (tr$type == "t-test") {
      print(tr$result)
      cat("\n---\n")
      if (tr$result$p.value < alpha) {
        cat(sprintf("Result: p = %.4g < alpha (%.3g) -> statistically significant difference between groups.\n",
                    tr$result$p.value, alpha))
      } else {
        cat(sprintf("Result: p = %.4g >= alpha (%.3g) -> no statistically significant difference detected.\n",
                    tr$result$p.value, alpha))
      }
    } else {
      print(summary(tr$result))
      cat("\n---\n")
      p_val <- summary(tr$result)[[1]][["Pr(>F)"]][1]
      if (p_val < alpha) {
        cat(sprintf("Result: p = %.4g < alpha (%.3g) -> at least one group differs. See Tukey HSD tab for pairwise comparisons.\n",
                    p_val, alpha))
      } else {
        cat(sprintf("Result: p = %.4g >= alpha (%.3g) -> no statistically significant difference detected among groups.\n",
                    p_val, alpha))
      }
    }
  })

  output$tukey_table <- renderDT({
    tr <- test_result()
    req(tr$type == "anova")
    tk <- as.data.frame(tr$tukey$group)
    tk$comparison <- rownames(tk)
    tk <- tk[, c("comparison", "diff", "lwr", "upr", "p adj")]
    tk[ , c("diff","lwr","upr","p adj")] <- round(tk[ , c("diff","lwr","upr","p adj")], 4)
    datatable(tk, options = list(dom = "t", paging = FALSE))
  })

  # ---- Plot ----
  make_plot <- function() {
    ad <- analysis_data()
    ds <- desc_stats()

    ggplot(ad, aes(x = group, y = value, fill = group)) +
      geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.6) +
      geom_jitter(width = 0.12, size = 2, alpha = 0.7, color = "#2C3E50") +
      geom_errorbar(
        data = ds,
        aes(x = group, y = mean, ymin = mean - se, ymax = mean + se),
        inherit.aes = FALSE, width = 0.15, color = "#C0392B", linewidth = 0.8
      ) +
      geom_point(
        data = ds, aes(x = group, y = mean),
        inherit.aes = FALSE, color = "#C0392B", size = 3, shape = 18
      ) +
      labs(
        x = input$group_col %||% "Group",
        y = input$value_col %||% "Value",
        title = "Group comparison",
        subtitle = "Boxplot with individual points; red diamond = mean +/- SE"
      ) +
      theme_minimal(base_size = 14) +
      theme(legend.position = "none")
  }

  `%||%` <- function(a, b) if (is.null(a)) b else a

  output$group_plot <- renderPlot({ make_plot() })

  output$download_plot <- downloadHandler(
    filename = function() "group_comparison_plot.png",
    content  = function(file) {
      ggsave(file, plot = make_plot(), width = 8, height = 6, dpi = 300)
    }
  )
}

shinyApp(ui, server)
