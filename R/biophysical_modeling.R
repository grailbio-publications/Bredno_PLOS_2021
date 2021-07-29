#' @file biophysics.R
#' Functionality used to create and validate biophysical models that explain
#' estimated cfDNA tumor fraction with clinical covariates.
#' This code is used for different cancer types, for example breast, lung, and
#' colon cancers.

requireNamespace("dplyr", quietly = TRUE)
requireNamespace("ggplot2", quietly = TRUE)
requireNamespace("pROC", quietly = TRUE)
requireNamespace("relaimpo", quietly = TRUE)
requireNamespace("stringr", quietly = TRUE)
requireNamespace("svglite", quietly = TRUE)

#' Impute tumor fraction from CCGA1 WGBS classifier score and create a figure
#' to visualize the imputation.
#' @param df_model: Data frame with at least the columns wgbs_classifier_score and
#'   tumor_fraction. tumor_fraction is expected to be NA for at least some rows.
#' @param replace_with_imputed: Logical. If set to TRUE, then missing tumor
#'   fraction information in the column tumor_fraction is replaced with the
#'   results of imputation
#' @param draw_or_fixed: Either "draw" or "fixed". When imputing tumor fractions,
#'   the look-up from m-score either uses a fixed set of fitting parameters or
#'   draws the fitting parameters randomly from their estimated distributions
#'   for each case.
#' @param score_frac_min, score_frac_max: Valid WGBS classifier score range to
#'   allow imputation of tumor fraction
#' @return list with two items:
#'   df_model: incoming df_model with additional columns imputed_ctf,
#'   ctf_source set to "imputed" or "measured", and optionally updated column
#'   ctf
#'   hPlot: ggplot2 figure plot that visualizes the imputation process
#'
#' This function fits a sigmoid function to wgbs_classifier_score over log(ctf).
#' The fitted function's value range is then used to identify a range of classifier
#' scores for which a lookup using the inverse sigmoid function can be used to
#' impute ctf. The range of valid imputation is
#' [min + score_frac_min * (max-min), min + score_frac_max * (max-min)]. For samples without
#' cTF measurement and a classifier score in this range, tumor fraction is imputed via
#' this fitted function. When the parameter draw_or_fixed is set to "draw", then for the imputation
#' of each cRF, the parameters of the sigmoid are drawn from individuall Gaussian distributions
#' with center and width estimated from the confidence interval of the fitting parameters.
impute_ctf <- function(
  df_model,
  replace_with_imputed = TRUE,
  draw_or_fixed = c("fixed", "draw"),
  score_frac_min = 0.1,
  score_frac_max = 0.9) {
  draw_or_fixed <- match.arg(draw_or_fixed)
  ctf_model <- nls(data = df_model,
                    wgbs_classifier_score ~ a / (1 + exp(-b * (log(ctf) - c))) + d,
                    start=list(a = 0.5, b = 1, c = log(5e-3), d = 0.5),
                    control = list(maxiter = 500))
  fp <- ctf_model$m$getPars() %>% as.list
  # For imputation with a random draw, the 95% confidence interval is requested
  # from the nls model fit for each parameter and assumed to span 2 * 1.96
  # standard deviations of that parameter for a random draw from a normal distribution
  ci_params <- confint(ctf_model, level = 0.95, control = list(maxiter = 500))
  sigma <- as.list((ci_params[, 2] - ci_params[, 1]) / 2 / 1.96)
  # local function to plot the fit, uses the list fp (fitting parameters) from
  # the caller namespace (defined above)
  fitted_func <- function(x) {
    fp$a / (1 + exp(-fp$b * (log(x) - fp$c))) + fp$d
  }
  # local function to estimate tf from wgbs classifier score, uses fitting
  # parameters from the caller namespace
  impute_ctf <- function(wgbs_classifier_score) {
    exp(-1/fp$b * log(fp$a / (wgbs_classifier_score - fp$d) - 1) + fp$c)
  }
  # local function to draw tf from wgbs classifier score with fitting parameters
  # and estimated standard deviation of fitting parameters for a random draw
  impute_ctf_with_draw <- function(wgbs_classifier_score, fp, sigma) {
    a <- rnorm(1, mean = fp$a, sd = sigma$a)
    b <- rnorm(1, mean = fp$b, sd = sigma$b)
    c <- rnorm(1, mean = fp$c, sd = sigma$c)
    d <- rnorm(1, mean = fp$d, sd = sigma$d)
    exp(-1/b * log(a / (wgbs_classifier_score-d) - 1) + c)
  }
  impute_wgbs_classifier_score_min <- fp$d + score_frac_min * fp$a
  impute_wgbs_classifier_score_max <- fp$d + score_frac_max * fp$a

  if (draw_or_fixed == "fixed") {
    df_model$imputed_ctf <- sapply(df_model$wgbs_classifier_score,
                                    impute_ctf)
  } else {
    df_model$imputed_ctf <- sapply(df_model$wgbs_classifier_score,
                                    impute_ctf_with_draw,
                                    fp, sigma)
  }
  no_imputation <- !is.na(df_model$ctf) |
    df_model$wgbs_classifier_score < impute_wgbs_classifier_score_min |
    df_model$wgbs_classifier_score > impute_wgbs_classifier_score_max
  df_model$imputed_ctf[no_imputation] <- NA

  df_plot <- df_model %>%
    dplyr::mutate(ctf_source = factor(ifelse(is.na(ctf),
                                              "imputed", "measured"),
                                       levels = c("measured", "imputed")),
                  tf_plot = ifelse(is.na(ctf),
                                   imputed_ctf, ctf))

  hPlot <- ggplot2::ggplot(df_plot) +
    ggplot2::stat_function(ggplot2::aes(x = ctf),
                           fun = fitted_func) +
    ggplot2::geom_point(ggplot2::aes(x = tf_plot, y = wgbs_classifier_score,
                                     color = clinical_stage, shape = ctf_source)) +
    ggplot2::scale_x_log10()

  if(replace_with_imputed) {
    df_model <- df_model %>%
      dplyr::mutate(
        ctf_source = factor(ifelse(is.na(ctf),
                                    "imputed", "measured"),
                             levels = c("measured", "imputed")),
        ctf = ifelse(is.na(ctf) & !is.na(imputed_ctf),
                      imputed_ctf, ctf))
  }
  list(df_model = df_model, hPlot = hPlot)
}

#' Assigned standardized appearance and style to a figure and export it in
#' standard image file formats eps and svg.
#' @param hPlot: ggplot2 figure handle
#' @param title: Set title of figure to this string.
#' @param axis_x: Set x-axis label to this string.
#' @param axis_y: Set y-axis label to this string.
#' @param log_x: Transform x-axis into log space.
#' @param log_y: Transform y-axis into log space.
#' @param color_legend: Set the label for colors shown in the legend box.
#' @param shape_legend: Set the label for shapes shown in the legend box.
#' @param legend_position: Set position of the legend box (can be moved inside of
#'   the graph area). Must be a vector of two numerical values in [0, 1].
#' @param color_points: Set the palette for plots / points to a predefined palette.
#' @param x_axis_vertical: Display labels on x-Axis vertically.
#' @param export: Save figure with this name in eps and svg file formats.
#'
#' All parameters are optional and only affect the figure when selected.
#' @return ggplot2 handle to the updated figure to save or render in a markdown.
assign_figure_style <- function(hPlot, ...) {
  l_selections <- list(...)
  predefined_colors <- list(
    black = "#000000",
    red = "#ED1C24",
    mgreen = "#28b069",
    beige = "#B89466",
    mblue = "#1F6683",
    dpurple = "#612B5B",
    orange = "#ECA400",
    mgray = "#949494",
    dblue = "#3B1C52",
    dgray = "#404040",
    lgray = "#D9D9D9")
  # Set standard, basic figure properties and style.
  hPlot <- hPlot +
    ggplot2::theme_minimal() +
    ggplot2::theme(panel.grid.major = ggplot2::element_line(size = 0.3,
                                                            linetype = "solid",
                                                            colour = "light gray"),
                   axis.line = ggplot2::element_line(size = 0.5,
                                                     linetype = "solid",
                                                     colour = "black"),
                   legend.box.background = ggplot2::element_rect(color="white"))
  selections <- names(l_selections)
  if ("title" %in% selections) {
    hPlot <- hPlot +
      ggplot2::ggtitle(ggplot2::element_text(l_selections$title))
  }
  if("axis_x" %in% selections) {
    hPlot <- hPlot +
      ggplot2::xlab(l_selections$axis_x)
  }
  if("axis_y" %in% selections) {
    hPlot <- hPlot +
      ggplot2::ylab(l_selections$axis_y)
  }
  if("log_x" %in% selections) {
    hPlot <- hPlot +
      ggplot2::scale_x_log10()
  }
  if("log_y" %in% selections) {
    hPlot <- hPlot +
      ggplot2::scale_y_log10()
  }
  if("color_legend" %in% selections) {
    hPlot <- hPlot +
      ggplot2::labs(colour = l_selections$color_legend)
  }
  if("shape_legend" %in% selections) {
    hPlot <- hPlot +
      ggplot2::labs(shape = l_selections$shape_legend)
  }
  if("legend_position" %in% selections) {
    hPlot <- hPlot +
      ggplot2::theme(legend.position = l_selections$legend_position)
  }
  if("color_points" %in% selections) {
    hPlot <- hPlot +
      ggplot2::scale_color_manual(values=as.character(predefined_colors))
  }
  if("x_axis_vertical" %in% selections) {
    hPlot <- hPlot +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, hjust = 1))
  }
  if("export" %in% selections) {
    ggplot2::ggsave(filename = sprintf("../figures/%s.eps", l_selections$export),
                    plot = hPlot,
                    dpi = 600, device = "eps")
    ggplot2::ggsave(filename = sprintf("../figures/%s.svg", l_selections$export),
                    plot = hPlot,
                    dpi = 600, device = "svg")
  }
  hPlot
}

#' Derive total volume or surface area of primary lesions given (partially missing)
#' lesion size and focality information.
#' @param df_model: Input data frame with one row per case. Must have at least
#'   these columns:
#'   size_1, size_2, size_3, ..., size_6: Maximum diameter of a tumor focus in mm
#'   with up to 6 measured lesions for solid tumor foci.
#'   Sizes are expected sorted by size and can be NA.
#'   n_primary_lesions: Number of tumor lesions, taking focality and for breast
#'   cancer laterality into account. There can be less non-NA size measurements
#'   than primary lesions, and this function imputes missing sizes of smaller
#'   lesions when at least the size of an index lesion is provided.
#' @param l_imputation: List with imputation paramters. If this parameter is NULL,
#'   then imputation parameters are computed from the incoming data (training set).
#'   If not NULL, then imputation parameters are taken from this list (validation set).
#'   Imputation parameters are the size fraction of a 2nd and 3rd non-index lesion
#'   to the index lesion, and imputation then estimates missing size of possible
#'   2nd, 3rd and further non-index lesions based on the index lesion size
#' @param draw_or_fixed (optional, default "fixed"): If "fixed", then imputed
#'   sizes for 2nd, 3rd and further non-index lesions are set to index lesion size
#'   multiplied with the median of the ratio between index lesion size and 2nd or
#'   3rd observed lesions, respectively. If set to "draw", then the ratio is drawn
#'   randomly from the 2nd and 3rd quartile (non-extreme values) for each imputed case.
#' @param dimension (optional, default 3): If 3, then volume of all primary lesions
#'   is computed (used for breast and lung cancers).
#'   If 2, then surface area of all primary lesions is computed (used for colorectal
#'   cancers).
#' @return list with items:
#'   df_model: Incoming data frame with an additional column primary_volume or
#'     primary_area that contains the total primary tumor size, summed up over
#'     all lesions in multifocal or multilateral disease.
#'   l_size_imputation: Parameters for imputation passed from training to
#'     validation data
derive_primary_size <- function(
  df_model,
  l_imputation = NULL,
  draw_or_fixed = c("fixed", "draw"),
  dimension = 3) {
  draw_or_fixed <- match.arg(draw_or_fixed)
  # This function works with up to 6 lesions (for example bilateral, multifocal
  # breast cancers) that are not always specified in the input data. Missing
  # columns are added with no data
  v_missing_cols <- setdiff(c("size_1", "size_2", "size_3",
                              "size_4", "size_5", "size_6"),
                            colnames(df_model))
  df_model[, v_missing_cols] <- NA
  if (is.null(l_imputation)) {
    l_imputation <- list()
    vec_frac_lesion2 <- df_model$size_2 / df_model$size_1
    vec_frac_lesion2 <- vec_frac_lesion2[!is.na(vec_frac_lesion2)]
    if (length(vec_frac_lesion2) == 0) {
      vec_frac_lesion2 <- c(0.5)
    }
    med_frac_lesion2 <- median(vec_frac_lesion2)
    vec_frac_lesion3 <- df_model$size_3 / df_model$size_1
    vec_frac_lesion3 <- vec_frac_lesion3[!is.na(vec_frac_lesion3)]
    if (length(vec_frac_lesion3) == 0) {
      vec_frac_lesion3 <- c(0.25)
    }
    med_frac_lesion3 <- median(vec_frac_lesion3)
    l_imputation[["frac_lesion2"]] <- med_frac_lesion2
    l_imputation[["frac_lesion3"]] <- med_frac_lesion3
    if (length(vec_frac_lesion2) >= 5) {
      quartiles <- quantile(vec_frac_lesion2, c(0.25, 0.75))
      vec_frac_lesion2 <- vec_frac_lesion2[vec_frac_lesion2 >= quartiles[1] &
                                             vec_frac_lesion2 <= quartiles[2]]
    }
    if (length(vec_frac_lesion3) >= 5) {
      quartiles <- quantile(vec_frac_lesion3, c(0.25, 0.75))
      vec_frac_lesion3 <- vec_frac_lesion3[vec_frac_lesion3 >= quartiles[1] &
                                             vec_frac_lesion3 <= quartiles[2]]
    }
    l_imputation[["vec_frac_lesion2"]] <- vec_frac_lesion2
    l_imputation[["vec_frac_lesion3"]] <- vec_frac_lesion3
  } else {
    med_frac_lesion2 <- l_imputation[["frac_lesion2"]]
    med_frac_lesion3 <- l_imputation[["frac_lesion3"]]
    vec_frac_lesion2 <- l_imputation[["vec_frac_lesion2"]]
    vec_frac_lesion3 <- l_imputation[["vec_frac_lesion3"]]
  }

  #' Macro-like helper function to compute volume with or without random draw
  #' for imputation
  #' @function primary_volume
  #' @param n_primary_lesions: Number of primary tumor lesions, determined from
  #'   laterality and focality of primary tumor(s)
  #' @param size1, ..., size6: Known primary lesion diameter / largest extent
  #'   or NA if no such measurement exists
  #' @return Sum of volumes of all primary lesions
  #' This function uses the variables draw_or_fixed, med_frac_lesion2,
  #' med_frac_lesion3, vec_frac_lesion2, vec_frac_lesion3 from the caller's
  #' namespace.
  primary_size <- function(n_primary_lesions,
                           size_1, size_2, size_3, size_4, size_5, size_6) {
    v_sizes <- c(size_1, size_2, size_3, size_4, size_5, size_6)
    if (draw_or_fixed  == "fixed") {
      frac_lesion2 <- med_frac_lesion2
      frac_lesion3 <- med_frac_lesion3
    } else {
      # draw from central Q2 and Q3 quartiles
      frac_lesion2 <- vec_frac_lesion2[sample.int(length(vec_frac_lesion2), 1)]
      frac_lesion3 <- vec_frac_lesion3[sample.int(length(vec_frac_lesion3), 1)]
    }
    # create vector with ratio of largest (index), 2nd, 3rd and further primary
    # lesions to the index (largest) lesions. The same fraction is estimated
    # for non-primary lesions 3, 4, 5, 6.
    v_frac_lesions <- c(1.0, frac_lesion2, frac_lesion3,
                        frac_lesion3, frac_lesion3, frac_lesion3)
    # imputed sizes from all lesions are taken by scaling the index lesion
    v_imputed_sizes <- size_1 * v_frac_lesions
    # missing sizes are imputed by resp. fraction of index lesion
    v_replace <- which(is.na(v_sizes[1:n_primary_lesions]))
    v_sizes[v_replace] <- v_imputed_sizes[v_replace]
    if (dimension == 2) {
      scale <- 1/4 * pi
    }
    if (dimension == 3) {
      scale <- 1/6 * pi
    }
    sum(scale * v_sizes[1:n_primary_lesions] ^ dimension)
  }
  df_model$primary_size <- mapply(primary_size,
                                  df_model$n_primary_lesions,
                                  df_model$size_1, df_model$size_2, df_model$size_3,
                                  df_model$size_4, df_model$size_5, df_model$size_6)
  if (dimension == 2) {
    df_model$primary_area <- df_model$primary_size
  }
  if (dimension == 3) {
    df_model$primary_volume <- df_model$primary_size
  }
  list(df_model = df_model,
       l_imputation = l_imputation)
}

#' Create a figure showing violin plots separated by a category.
#' Typical example is cTF over clinical stage.
#' @param df_data: Data frame with at least a column name_x with a categorical
#'   variable and a column name_y with a numerical variable.
#' @param name_x: Categorical variable that defines data sets for which violin plots
#'   are created.
#' @param name_y: Numerical variable. Distribution of this variable is used to create
#'   the violin plots.
#' @param ...: Additional parameters passed on to ggplots::aes_string to define
#'   data mappings for this plot.
#' @return ggplot2 figure handle to the generated figure to save or render in a markdown
#'   Typically, this handle is passed on to assign_figure_style() to obtain consistent
#'   appearance for all figures.
basic_violin <- function(df_data, name_x, name_y, ...) {
  jitter_mapping <- ggplot2::aes_string(x = name_x, y = name_y, ...)
  ggplot2::ggplot(df_data,
                  ggplot2::aes_string(x = name_x, y = name_y)) +
    ggplot2::geom_violin(scale = "width") +
    ggplot2::geom_boxplot(width=0.1, outlier.shape = NA, fill = "light gray") +
    ggplot2::geom_jitter(jitter_mapping)
}

#' Derive and impute lymph node status and number of tumor-involved lymph nodes.
#' @param df_data Data frame with at least the columns ln_involved, n_stage, and
#'   clinical_stage
#' @return Incoming data frame with an additional column ln_status. This is
#'   0 if no tumor-involved lymph nodes present and 1 if tumor-involved lymph
#'     nodes are present. It is NA if information cannot be determined even by
#'     imputation from clinical stage alone.
#' The function first aggregates information of tumor-involved lymph nodes for
#' bilateral disease or multiple path reports that capture information for multiple
#' primary lesions.
#' Imputation imputes and adjudicates presence of tumor-involved lymph nodes
#' following this prioritized list of rules (later rules are only considered
#' if none of the previous rules triggered):
#' - Present if at least one path report lists at least one tumor-involved lymph
#'   node
#' - Present if clinical N-stage is N1, N2, N3
#' - Absent if clinical N-stage is N0
#' - Absent if path reports list the number of examined lymph nodes with no lymph
#'   node showing tumor involvement.
#' - Present if clinical stage is III or IV
#' - Absent if clinical stage is I or II
get_ln_status <- function(df_data) {
  derive_n_ln <- function(ln_spec) {
    ln_spec %>%
      stringr::str_remove_all("Other/missing,") %>%
      stringr::str_remove_all(",Other/missing") %>%
      stringr::str_split(",") %>%
      .[[1]] %>%
      as.integer() %>%
      sum
  }
  df_data %>%
    dplyr::rowwise() %>%
    dplyr::mutate(n_ln = derive_n_ln(ln_involved),
                  ln_status = dplyr::case_when(
                    n_ln > 0 ~ 1,
                    n_stage %in% c("N1", "N2", "N3") ~ 1,
                    n_stage == "N0" ~ 0,
                    n_ln == 0 ~ 0,
                    clinical_stage %in% c("III", "IV") ~ 1,
                    clinical_stage %in% c("I", "II") ~ 0,
                    TRUE ~ NA_real_))
}

#' Consecutively create and report characteristics of linear models after
#' selecting clinically significant covariates
#' @param df_model: Data frame with at least the columns specified in v_targets and
#'   v_variables
#' @param v_targets: Vector with two colum names. The first is the column name of
#'   the target variable for the analytical and predictive models. The is the second
#'   is the colum name of the target variable for a quantitative model.
#' @param v_variables: Vector of column names with the covariates / candidate inputs
#'   used for modeling.
#' @param v_constraints (optional): If provided, logical vector with the same length as
#'   v_variables. The vector selects (by setting an entry to TRUE) variables that
#'   have a real-world non-negativity constraint.
#'   Variables with a negative effect that are selected to be constrained are rejected
#'   as violating a non-negative constraint. The analysis model is updated without the
#'   constraint-violating parameters. This approach is chosen as compromise to enable
#'   physical non-negativity constraints while also obtaining interpretable p-values
#'   from a lienar model.
#' @return list with these items:
#'   analysis_model: Linear model with all variables that did not violate non-
#'     negativity contraints fitted to the first target variable.
#'   prediction_model: Linear model with only the variables with a p-value < 0.05
#'     in the analysis model fitted to the first target variable.
#'   quantitative_model: Linear model with only the variables with a p-value < 0.05
#'     in the analysis model fitted to the second target variable.
#'   df_analysis_model: Data frames that captures the output of summary(model)
#'     and an estimate or relative importance in a data frame to characterize the
#'     analysis model in reports and markdowns.
#'   violate_constraints: Vector of input variables that are not part of the
#'     analysis model because they violated a non-negativity constraint.
create_model <-function(df_model,
                        v_targets = c("ctf", "ctdna_amount_ge"),
                        v_variables,
                        v_constraints,
                        data_name = NULL) {
  lret <- list()
  # Create analysis model.
  variables <- paste(v_variables, collapse = " + ")
  formula <- sprintf("%s ~ %s", v_targets[1], variables)
  analysis_model <- lm(data = df_model, as.formula(formula))
  # Re-compute analysis model when non-negativity constraints are violated.
  if (!missing(v_constraints)) {
    v_coef <- coef(analysis_model)
    is_intercept <- (names(v_coef) == "(Intercept)")
    v_coef <- v_coef[!is_intercept]
    # The caller must ensure that number of covariates and number of entries in
    # v_constraint match.
    stopifnot(length(v_coef) == length(v_constraints))
    v_violations <- (v_coef < 0) & v_constraints
    if (sum(v_violations) > 0) {
      lret$violate_constraints <- names(v_coef[v_violations])
      v_variables <- v_variables[!v_violations]
      variables <- paste(v_variables, collapse = " + ")
      formula <- sprintf("%s ~ %s", v_targets[1], variables)
      analysis_model <- lm(data = df_model, as.formula(formula))
    } else {
      lret$violate_constraints <- c()
    }
  }
  df_analysis_model <- summarize_model(analysis_model)
  lret$analysis_model <- analysis_model
  lret$df_analysis_model <- df_analysis_model

  prediction_variables <- df_analysis_model %>%
    dplyr::filter(p_value <= 0.05) %>%
    dplyr::pull(variable)
  if (length(prediction_variables) > 0) {
    # pick variables and intercept for prediction model based on p-value of analytical model
    if ("(Intercept)" %in% prediction_variables) {
      prediction_variables <- prediction_variables[prediction_variables != "(Intercept)"]
      variables <- paste(prediction_variables, collapse = " + ")
    } else {
      variables <- paste(c("0", prediction_variables), collapse = " + ")
    }
    formula <- sprintf("%s ~ %s", v_targets[1], variables)
    prediction_model <- lm(data = df_model, as.formula(formula))
    # quantitative model uses the same variables, but a different target
    formula <- sprintf("%s ~ %s", v_targets[2], variables)
    quantitative_model <- lm(data = df_model, as.formula(formula))

    lret$prediction_model <- prediction_model
    lret$quantitative_model <- quantitative_model
  }
  lret
}

#' Create a "pretty print" table that summarizes a linear model
#' @param model: Fitted model as for example obtained from lm
#' @return Data frame with the summary of the model and an additional column with the relative
#'   importance of variables from the relaimpo package (only provided if the input model has an
#'   intercept estimate)
summarize_model <- function(model) {
  model_summary <- summary(model)
  # "pretty-print" data frame from the lm model summary
  df_summary <- as.data.frame(model_summary$coefficients) %>%
    dplyr::rename(estimate = Estimate,
                  std_error = `Std. Error`,
                  t_value = `t value`,
                  p_value = `Pr(>|t|)`) %>%
    dplyr::mutate(variable = rownames(model_summary$coefficients),
                  significance = dplyr::case_when(
                    p_value < 0.001 ~ "***",
                    p_value < 0.01 ~ "**",
                    p_value < 0.05 ~ "*",
                    p_value < 0.1 ~ ".",
                    TRUE ~ " "),
                  estimate = sprintf("%.3g", estimate),
                  print_p_value = sprintf("%.4g", p_value))
  # additional column from relaimpo package, which requires the model to have an intercept
  if ("(Intercept)" %in% names(model$coefficients)) {
    relative_importance <- relaimpo::calc.relimp.lm(model,
                                                    nk = TRUE,
                                                    rela = TRUE)
    df_summary <- df_summary %>%
      dplyr::left_join(data.frame(variable = names(relative_importance$lmg),
                                  rel_importance = relative_importance$lmg),
                       by = "variable") %>%
      dplyr::mutate(rel_importance = sprintf("%.3f", rel_importance)) %>%
      dplyr::select(variable, estimate, p_value, print_p_value,
                    significance, rel_importance)
  } else {
    df_summary <- df_summary %>%
      dplyr::select(variable, estimate, p_value, print_p_value, significance)
  }
  df_summary
}

#' Derive assay quantitation information.
#' @param df_model: Data frame with at least the columns weight_kg, height_m,
#'   sex, ctf, cfdna_conc_ng_ml
#' @return updated df_model data frame with additional columns:
#'   blood_volume_l: Estimated patient whole blood volume in liters
#'   plasma_volume_l: Estimated patient plasma volume in liters
#'   cfdna_conc_ge_ml: Concentration of all cfDNA in plasma in genome equivalents / ml
#'   ctdna_conc_ng_ml: Concentration of circulating tumor DNA in plasma in ng / ml
#'   ctdna_conc_ge_ml: Concentration of circulating tumor DNA in plasma in genome equivalents / ml
#'   ctdna_amount_ng: Total amount of ctDNA in patient circulation in ng
#'   ctdna_amount_ge: Total amount of ctDNA in patient circulation in genome equivalents
assay_quant <- function(df_model, df_assay) {
  # [Nad62] Nadler SB, Hidalgo JH, Bloch T. Prediction of blood volume in normal
  #   human adults. Surgery. 1962 Feb;51(2):224–32.
  # [Pio19] 1. Piovesan A, Pelleri MC, Antonaros F, Strippoli P, Caracausi M, Vitale L.
  #   On the length, weight and GC content of the human genome. BMC Res Notes
  #   [Internet]. 2019 Feb 27.
  df_model <- df_model %>%
    dplyr::mutate(blood_volume_l = ifelse(sex == "Female",  # [Nad62]
                                          0.3561 * height_m^3 + 0.03308 * weight_kg + 0.1833,
                                          0.3669 * height_m^3 + 0.03219 * weight_kg + 0.6041),
                  plasma_volume_l = blood_volume_l * 0.55,
                  ctdna_conc_ng_ml = cfdna_conc_ng_ml * ctf,
                  ctdna_amount_ng = ctdna_conc_ng_ml * plasma_volume_l * 1e3,
                  cfdna_conc_ge_ml = cfdna_conc_ng_ml / 6.5 * 1000,  # [Pio19]
                  ctdna_conc_ge_ml = ctdna_conc_ng_ml / 6.5 * 1000,
                  ctdna_amount_ge = ctdna_amount_ng / 6.5 * 1000)
}

#' Create a waterfall plot with one additional color-coded category
#' @param df_data: Data frame with at least the columns name_y, color_1, color_2
#' @param name_y name of column in df_data that holds a numeric value for height in the
#'   waterfall plot.
#' @param color_1 name of a column in df_data that holds a categorical value used to define
#'   the color of each bar in the waterfall plot.
#' @param color_2 name of a column in df_data that holds a categorical value used to define
#'   the color of an additional color bar below the waterfall plot.
#' @return ggplot2 figure handle with the created plot
waterfall_plot <-function(df_data, name_y, color_1, color_2) {
  v_colors_1 <- c("yellow", "blue", "red", "cyan", "magenta", "green", "black", "white")
  v_colors_2 <- c("black", "gray", "white", "dark grey", "light grey")
  n_classes_1 <- length(unique(df_data[[color_1]]))
  n_classes_2 <- length(unique(df_data[[color_2]]))

  used_colors <- c(v_colors_2[1:n_classes_2],
                   v_colors_1[1:n_classes_1])
  # For the display of the waterfall plot with tumor fraction on log scale
  # and boxes / markers colored by stage below each tumor fraction bar, the
  # display ranges of the boxes are individually computed.
  # Cases are sorted by tumor fraction. The left and right coordinates of the
  # boxes are rank and rank+1, resp. To order with highest bar on the left,
  # negative rank is used.
  df_data$xmin <- -rank(df_data[[name_y]], ties.method = "random")
  df_data$xmax <- df_data$xmin + 1
  # The top position of each bar is tumor fraction - the "height" of the bar
  df_data$ymax <- df_data[[name_y]]
  tf_top <- max(df_data$ymax, na.rm = TRUE)
  # The bottom of all bars is the lowest tumor fraction value scaled down by
  # 0.9 to give even the lowest bar a fixed height > 0.
  tf_bottom <- min(df_data$ymax, na.rm = TRUE) * 0.9
  df_data$ymin <- tf_bottom
  # The boxes colored by stage are located below the bars (their top coordinate
  # is the bottom coordinate of the bars), and their size is 10% of the highest
  # bar. The plot will be in log scale, so the height of the stage boxes has to
  # be computed to be 10% of tumor fraction bar height after log transform
  df_data$stage_bottom <- 10^(1.1 * log10(tf_bottom) - 0.1 * log(tf_top))
  df_data[[color_2]] <- paste0(" ", df_data[[color_2]])

  ggplot2::ggplot(data = df_data %>%
                    dplyr::filter(!is.na(xmin))) +
    ggplot2::geom_rect(mapping = ggplot2::aes_string(xmin = "xmin", xmax = "xmax",
                                                     ymin = "ymin", ymax = "ymax",
                                                     fill = color_1)) +
    ggplot2::geom_rect(mapping = ggplot2::aes_string(xmin = "xmin", xmax = "xmax",
                                                     ymin = "stage_bottom", ymax = "ymin",
                                                     fill = color_2)) +
    ggplot2::xlab("Case") +
    ggplot2::scale_y_log10() +
    ggplot2::theme_minimal() +
    ggplot2::scale_fill_manual(values = used_colors) +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_line(size = 0.3,
                                                              linetype = "solid",
                                                              colour = "light gray"),
                   panel.grid.major.x = ggplot2::element_blank(),
                   panel.grid.minor.x = ggplot2::element_blank())
}

#' Create an ROC plot, compute AUC and resp. CI, perform Wilcoxon test.
#' @param df_data: Data frame with at least the columns detected (aka test result)
#'   and feature.
#' @param detected: Name of a column in df_data with the test result as logical value.
#' @param feature: Name of a column in df_data with a feature that is input to the test.
#' @return list with
#'   \itemize {
#'    \item hPlot: ggplot2 plot handle with an ROI curve. The plot handle can get
#'      passed on to assign_figure_style for a standardized appearance.
#'    \item auc: Area under curve
#'    \item auc_ci: 95% confidence interval of area under curve
#'    \item wilcox_p: p-Value of a two sample one-sided Wilcoxon Rank Sum test to
#'      determine if detected cases have a higher feature value than non-detected cases.
#'  }
simple_roc <- function(df_data,
                       detected = "detected",
                       feature = "ctf_predicted") {
  controls <- df_data %>%
    dplyr::filter(!!rlang::sym(detected) == FALSE) %>%
    dplyr::pull(!!rlang::sym(feature))
  cases <- df_data %>%
    dplyr::filter(!!rlang::sym(detected) == TRUE) %>%
    dplyr::pull(!!rlang::sym(feature))
  object_roc <- pROC::roc(controls = controls,
                          cases = cases,
                          smooth = FALSE)
  df_roc <- data.frame(TPR = object_roc$sensitivities,
                       FPR = 1.0 - object_roc$specificities) %>%
    dplyr::arrange(TPR, FPR)
  auc <- pROC::auc(object_roc)
  auc_ci <- pROC::ci.auc(object_roc)
  hPlot <- ggplot2::ggplot(df_roc %>%
                             dplyr::arrange(TPR),
                           ggplot2::aes(x = FPR, y = TPR)) +
    ggplot2::geom_line() +
    ggplot2::geom_segment(ggplot2::aes(x = 0, xend = 1, y = 0, yend = 1),
                          color="grey", linetype="dashed") +
    ggplot2::xlim(0.0, 1.0) + ggplot2::ylim(0.0, 1.0)
  wilcox_test <- wilcox.test(cases, controls, alternative = "greater")
  return(list(hPlot = hPlot,
              auc = auc,
              auc_ci = auc_ci,
              wilcoxon_p = wilcox_test$p.value))
}

#' Create a ROC plot with additional break-down by a category, compute AUCs and
#' resp. CIs, perform Wilcoxon test.
#' @param df_data: Data frame with at least the columns detected (aka test result),
#'   feature, and breakdown.
#' @param detected: Name of a column in df_data with the test result as logical value.
#' @param feature: Name of a column in df_data with a feature that is input to the test.
#' @param breakdown: Name of a categorical variable (for example, clinical stage)
#'   to present additional ROC curves per category.
#' @return list with
#'   hPlot: ggplot2 plot handle with an ROI curve
#'   auc_xx: Area under curve
#'   auc_ci_xx: 95% confidence interval of area under curve
#'   wilcoxon_p_xx: p-Value of a two sample one-sided Wilcoxon Rank Sum test to
#'     determine if detected cases have a higher feature value than non-detected cases.
#'   xx is "all" or a category value from the breakdown column. The values are
#'   computed for all cases and then separately for each category breakdown, for
#'   example separately for each clinical stage.
combined_roc <- function(df_data,
                         detected = "detected",
                         feature = "tf_predicted",
                         breakdown = "clinical_stage") {
  scores <- df_data[[feature]]  # classifier score, model prediction, etc.
  labels <- df_data[[detected]]  # TRUE or FALSE for detected or not detected
  df_plot <- data.frame(TPR = numeric(),
                        FPR = numeric(),
                        breakdown = character())
  v_breakdown <- c("all", sort(unique(df_data[[breakdown]])))
  lret <- list()
  for (selected_breakdown in v_breakdown) {
    if (selected_breakdown != "all") {
      df_for_roc <- df_data %>%
        dplyr::filter(!!rlang::sym(breakdown) == selected_breakdown)
    } else {
      df_for_roc <- df_data
    }
    controls <- df_for_roc %>%
      dplyr::filter(!!rlang::sym(detected) == FALSE) %>%
      dplyr::pull(!!rlang::sym(feature))
    cases <- df_for_roc %>%
      dplyr::filter(!!rlang::sym(detected) == TRUE) %>%
      dplyr::pull(!!rlang::sym(feature))
    object_roc <- pROC::roc(controls = controls,
                            cases = cases,
                            smooth = FALSE)
    df_individual_roc <- data.frame(TPR = object_roc$sensitivities,
                                    FPR = 1.0 - object_roc$specificities) %>%
      dplyr::arrange(TPR)
    lret[[paste0("auc_", selected_breakdown)]] <- pROC::auc(object_roc)
    ci <- pROC::ci.auc(object_roc)
    lret[[paste0("auc_ci_", selected_breakdown)]] <- ci
    wilcox_test <- wilcox.test(cases, controls, alternative = "greater")
    lret[[paste0("wilcoxon_p_", selected_breakdown)]] <- wilcox_test$p.value
    df_plot <- rbind(df_plot,
                     df_individual_roc %>%
                       dplyr::transmute(TPR = TPR,
                                        FPR = FPR,
                                        breakdown = selected_breakdown))
  }
  hPlot <- ggplot2::ggplot() +
    ggplot2::geom_line(data = df_plot,
                       mapping = ggplot2::aes(x = FPR, y = TPR, color = breakdown)) +
    ggplot2::geom_segment(mapping = ggplot2::aes(x = 0, xend = 1, y = 0, yend = 1),
                          color="grey", linetype="dashed") +
    ggplot2::xlim(0.0, 1.0) + ggplot2::ylim(0.0, 1.0)
  lret[["hPlot"]] <- hPlot
  lret
}


#' Count the number of cases that require imputation of cTF, tumor size of a
#' non-index lesion in multifocal disease, or presence of tumor-involved lymph
#' nodes
#' @param df_data: Data frame with at least the columns ctf, ctf_source,
#'   size_1...size_6, n_primary_lesions, n_ln, n_stage, clinical_stage that are
#'   used in the functions impute_ctf, derive_primary_size, and get_ln_status,
#'   resp.
#' @return list with counts of cases with measured or imputed data:
#'   n_measured_ctf: Number of cases with measured cTF.
#'   n_imputed_ctf: Number of cases with imputed cTF.
#'   n_complete_size: Number of cases with complete size information (size for
#'     all lesions when multifocal disease is reported)
#'   n_imputed_size: Number of cases with multifocal disease where at least one
#'     size of a non-index lesion is imputed.
#'     n_ln_by_path: Number of cases where lymph node status was determined from
#'       number of tumor-involved lymph nodes reported on a pathology report.
#'     n_ln_by_n_stage: Number of cases where lymph node status was determined
#'       from N stage information.
#'     n_ln_by_clinical_stage: Number of cases where lymph node status was
#'       imputed from clinical stage.
report_imputation <- function(df_data) {
  l_ret <- list()
  l_ret$n_measured_ctf <-
    sum(!is.na(df_data$ctf) & df_data$ctf_source == "measured")
  l_ret$n_imputed_ctf <-
    sum(!is.na(df_data$ctf) & df_data$ctf_source == "imputed")
  v_all_size_columns <- c("size_1", "size_2", "size_3",
                          "size_4", "size_5", "size_6")
  df_data[, setdiff(v_all_size_columns, colnames(df_data))] <- NA
  df_data$n_size_measurements <- rowSums(!is.na(df_data[, v_all_size_columns]))
  l_ret$n_complete_size <-
    sum(df_data$n_size_measurements >= df_data$n_primary_lesions &
          !is.na(df_data$size_1))
  l_ret$n_imputed_size <-
    sum(df_data$n_size_measurements < df_data$n_primary_lesions &
          !is.na(df_data$size_1))
  df_data$ln_source <-
    dplyr::case_when(df_data$n_ln > 0 ~ "n LN",
                     df_data$n_stage %in% c("N1", "N2", "N3") ~ "N stage",
                     df_data$n_stage == "N0" ~ "N stage",
                     df_data$n_ln == 0 ~ "n LN",
                     df_data$clinical_stage %in% c("III", "IV") ~ "Clinical stage",
                     df_data$clinical_stage %in% c("I", "II") ~ "Clinical stage",
                     TRUE ~ "no information")
  l_ret$n_ln_by_path <- sum(df_data$ln_source == "n LN")
  l_ret$n_ln_by_n_stage <- sum(df_data$ln_source == "N stage")
  l_ret$n_ln_by_clinical_stage <- sum(df_data$ln_source == "Clinical stage")

  l_ret
}
