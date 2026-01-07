library(tidyverse)

## ============================================================
## 1. LOAD REP1 DATA, CLEAN, AND COMPUTE SUM/AREA
## ============================================================

data_dir  <- "/Users/nataliabarbosa/Frydman Lab Dropbox/Natalia Moreira Barbosa/Natalia B. with Judith/eIF5A and RQC/RQC_ANALYSIS_data/rep1"
csv_files <- list.files(path = data_dir, pattern = "\\.csv$", full.names = TRUE)

# Load one CSV from rep1
load_sample_data <- function(file_path) {
  file_name   <- basename(file_path)
  # Example filename: WT_1mMGC7_12h_A2-04 Region GFP...
  sample_name <- stringr::str_extract(file_name, "(?<=WT_1mMGC7_12h_).*?(?= Region GFP)")
  
  raw <- read.csv(file_path, check.names = FALSE)
  
  df <- raw %>%
    dplyr::select(
      `Area::Area!!R`,
      `IntensityMean_EGFP-T3::Intensity Mean Value of channel 'EGFP-T3'!!R`,
      `IntensityMean_mRF12-T2::Intensity Mean Value of channel 'mRF12-T2'!!R`,
      `IntensitySum1_EGFP-T3::Intensity Sum of channel 'EGFP-T3'!!R`,
      `IntensitySum1_mRF12-T2::Intensity Sum of channel 'mRF12-T2'!!R`
    ) %>%
    dplyr::mutate(
      Sample    = sample_name,  # acquisition ID (A2-04, A4-02, ...)
      Replicate = "rep1"
    )
  
  # Rename numeric columns
  colnames(df)[1:5] <- c("area", "EGFP_mean", "mRF12_mean",
                         "EGFP_sum", "mRF12_sum")
  df
}

# Bind all files
all_data <- purrr::map_dfr(csv_files, load_sample_data)

# Drop the units row: first row within each file (Sample × Replicate)
all_data <- all_data %>%
  dplyr::group_by(Sample, Replicate) %>%
  dplyr::slice(-1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    area       = as.numeric(area),
    EGFP_mean  = as.numeric(EGFP_mean),
    mRF12_mean = as.numeric(mRF12_mean),
    EGFP_sum   = as.numeric(EGFP_sum),
    mRF12_sum  = as.numeric(mRF12_sum),
    EGFP_sum_per_area  = EGFP_sum / area,
    mRF12_sum_per_area = mRF12_sum / area
  )

## ============================================================
## 2. MAP SAMPLE IDs → BIOLOGICAL LABELS, ADD CONDITION/REPORTER
## ============================================================

# Map well IDs to biological names
all_data$Sample <- forcats::fct_recode(
  all_data$Sample,
  "WT_RFP-GFP"                       = "A2-04",
  "WT_RFP-(AAA)12x-GFP"              = "A3-01",
  "WT_RFP-(AAG)12x-GFP"              = "A4-02", 
  "WT_RFP-(R)12x-GFP_exp"           = "A5_exp-01",
  "WT_RFP-GFP +1mMGC7"              = "B2-05",
  "WT_RFP-(AAA)12x-GFP +1mMGC7"     = "B3-02",
  "WT_RFP-(AAG)12x-GFP +1mMGC7"     = "B4-02", 
  "WT_RFP-(R)12x-GFP_exp +1mMGC7"   = "B5_exp-01"
)

# Condition and reporter identity
all_data <- all_data %>%
  dplyr::mutate(
    Condition = dplyr::if_else(
      stringr::str_detect(as.character(Sample), "1mMGC7"),
      "GC7", "noGC7"
    ),
    Reporter  = dplyr::case_when(
      stringr::str_detect(as.character(Sample), "AAA") ~ "(AAA)12x",
      stringr::str_detect(as.character(Sample), "AAG") ~ "(AAG)12x",
      stringr::str_detect(as.character(Sample), "\\(R\\)") ~ "(R)12x",
      TRUE ~ "Control"   # WT_RFP-GFP
    )
  )

## ============================================================
## 3. SUBSET TO CONTROL & AAG, PREP LONG TABLE FOR MEAN SIGNALS
## ============================================================

rep1_clean <- all_data %>%
  dplyr::filter(Reporter %in% c("Control", "(AAG)12x")) %>%
  dplyr::mutate(
    Reporter  = factor(Reporter,  levels = c("Control", "(AAG)12x")),
    Condition = factor(Condition, levels = c("noGC7", "GC7"))
  )

rep1_long <- rep1_clean %>%
  tidyr::pivot_longer(
    cols      = c(EGFP_mean, mRF12_mean),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Metric = factor(Metric, levels = c("EGFP_mean", "mRF12_mean")),
    Metric_color = dplyr::case_when(
      Metric == "EGFP_mean"  ~ "#4daf4a",  # GFP
      Metric == "mRF12_mean" ~ "#e41a1c"   # RFP
    )
  )

## ============================================================
## 4. FINAL PLOT — VIOLIN + BOX, GC7 DARKER
## ============================================================

ggplot(rep1_long,
       aes(x = interaction(Reporter, Condition), y = Value)) +
  
  # Violin layer (GC7 darker via alpha mapping)
  geom_violin(
    aes(fill = Metric_color, alpha = Condition),
    trim  = FALSE,
    color = NA
  ) +
  
  # Boxplot overlay
  geom_boxplot(
    aes(fill = Metric_color, alpha = Condition),
    width        = 0.18,
    color        = "black",
    outlier.size = 0.6
  ) +
  
  # Use hex colors directly; alpha is mapped to Condition
  scale_fill_identity() +
  scale_alpha_manual(values = c("noGC7" = 0.35, "GC7" = 0.8)) +
  
  # Order and label x-axis as Control/AAG × Condition
  scale_x_discrete(
    labels = c(
      "Control.noGC7"  = "Control\nnoGC7",
      "(AAG)12x.noGC7" = "AAG\nnoGC7",
      "Control.GC7"    = "Control\nGC7",
      "(AAG)12x.GC7"   = "AAG\nGC7"
    )
  ) +
  
  facet_wrap(
    ~ Metric,
    scales   = "free_y",
    labeller = as_labeller(c(
      EGFP_mean  = "GFP (EGFP_mean)",
      mRF12_mean = "RFP (mRF12_mean)"
    ))
  ) +
  
  labs(
    title = "Rep 1 — Mean Fluorescence Distribution (EGFP_mean / mRF12_mean)",
    x     = "Sample (Reporter × Condition)",
    y     = "Mean fluorescence (a.u.)",
    alpha = "Condition"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
    strip.text  = element_text(size = 13, face = "bold"),
    legend.title = element_blank()
  )












library(tidyverse)

## ============================================================
##  REP 2 – LOAD, CLEAN, ANNOTATE, AND PLOT FLUORESCENCE MEANS
## ============================================================

## ---------- 1. LOAD REP2 RAW DATA ---------------------------------

data_dir2  <- "/Users/nataliabarbosa/Frydman Lab Dropbox/Natalia Moreira Barbosa/Natalia B. with Judith/eIF5A and RQC/RQC_ANALYSIS_data/rep2"
csv_files2 <- list.files(path = data_dir2, pattern = "\\.csv$", full.names = TRUE)

# Function to load a single CSV from rep2
load_sample_rep2 <- function(file_path) {
  file_name  <- basename(file_path)
  # Example: WT_1mMGC7_12h_A2-04-Orthogonal Projection-108 Region Class 2
  # We want just the well ID: A2-04, A4-02, B2-05, B4-02, etc.
  sample_name <- stringr::str_extract(
    file_name,
    "(?<=WT_1mMGC7_12h_).*?(?=-Orthogonal Projection)"
  )
  
  raw <- read.csv(file_path, check.names = FALSE)
  
  df <- raw %>%
    dplyr::select(
      `Area::Area!!R`,
      `IntensityMean_EGFP-T3::Intensity Mean Value of channel 'EGFP-T3'!!R`,
      `IntensityMean_mRF12-T2::Intensity Mean Value of channel 'mRF12-T2'!!R`,
      `IntensitySum1_EGFP-T3::Intensity Sum of channel 'EGFP-T3'!!R`,
      `IntensitySum1_mRF12-T2::Intensity Sum of channel 'mRF12-T2'!!R`
    ) %>%
    dplyr::mutate(
      Sample    = sample_name,   # raw acquisition ID (A2-04, A4-02, ...)
      Replicate = "rep2"
    )
  
  colnames(df)[1:5] <- c("area", "EGFP_mean", "mRF12_mean",
                         "EGFP_sum", "mRF12_sum")
  df
}

# Bind all files for rep2
all_data_rep2 <- purrr::map_dfr(csv_files2, load_sample_rep2)

## ---------- 2. DROP UNIT ROWS + MAKE NUMERIC ----------------------

# Each file has a "units" row in the first line → drop first row per Sample × Replicate
all_data_rep2 <- all_data_rep2 %>%
  dplyr::group_by(Sample, Replicate) %>%
  dplyr::slice(-1) %>%
  dplyr::ungroup() %>%
  dplyr::mutate(
    area       = as.numeric(area),
    EGFP_mean  = as.numeric(EGFP_mean),
    mRF12_mean = as.numeric(mRF12_mean),
    EGFP_sum   = as.numeric(EGFP_sum),
    mRF12_sum  = as.numeric(mRF12_sum),
    EGFP_sum_per_area  = EGFP_sum / area,
    mRF12_sum_per_area = mRF12_sum / area
  )

## ---------- 3. MAP WELL IDs → BIOLOGICAL SAMPLE NAMES --------------

# Same mapping as rep1, now applied to rep2
all_data_rep2$Sample <- forcats::fct_recode(
  all_data_rep2$Sample,
  "WT_RFP-GFP"                       = "A2-04",
  "WT_RFP-(AAA)12x-GFP"              = "A3-01",
  "WT_RFP-(AAG)12x-GFP"              = "A4-02", 
  "WT_RFP-(R)12x-GFP_exp"           = "A5_exp-01",
  "WT_RFP-GFP +1mMGC7"              = "B2-05",
  "WT_RFP-(AAA)12x-GFP +1mMGC7"     = "B3-02",
  "WT_RFP-(AAG)12x-GFP +1mMGC7"     = "B4-02", 
  "WT_RFP-(R)12x-GFP_exp +1mMGC7"   = "B5_exp-01"
)

## ---------- 4. ADD CONDITION & REPORTER ----------------------------

all_data_rep2 <- all_data_rep2 %>%
  dplyr::mutate(
    Condition = dplyr::if_else(
      stringr::str_detect(as.character(Sample), "1mMGC7"),
      "GC7", "noGC7"
    ),
    Reporter = dplyr::case_when(
      stringr::str_detect(as.character(Sample), "AAA") ~ "(AAA)12x",
      stringr::str_detect(as.character(Sample), "AAG") ~ "(AAG)12x",
      stringr::str_detect(as.character(Sample), "\\(R\\)") ~ "(R)12x",
      TRUE ~ "Control"
    )
  )

## ---------- 5. SUBSET TO CONTROL & AAG -----------------------------

rep2_clean <- all_data_rep2 %>%
  dplyr::filter(Reporter %in% c("Control", "(AAG)12x")) %>%
  dplyr::mutate(
    Reporter  = factor(Reporter,  levels = c("Control", "(AAG)12x")),
    Condition = factor(Condition, levels = c("noGC7", "GC7"))
  )

## ---------- 6. LONG FORMAT FOR MEAN SIGNALS ------------------------

rep2_long <- rep2_clean %>%
  tidyr::pivot_longer(
    cols      = c(EGFP_mean, mRF12_mean),
    names_to  = "Metric",
    values_to = "Value"
  ) %>%
  dplyr::mutate(
    Metric = factor(Metric, levels = c("EGFP_mean", "mRF12_mean")),
    Metric_color = dplyr::case_when(
      Metric == "EGFP_mean"  ~ "#4daf4a",  # GFP
      Metric == "mRF12_mean" ~ "#e41a1c"   # RFP
    )
  )

## ---------- 7. FINAL PLOT – VIOLIN + BOX, GC7 DARKER ---------------

ggplot(rep2_long,
       aes(x = interaction(Reporter, Condition), y = Value)) +
  
  # Violin: color = GFP/RFP, GC7 darker via alpha
  geom_violin(
    aes(fill = Metric_color, alpha = Condition),
    trim  = FALSE,
    color = NA
  ) +
  
  # Boxplot overlay
  geom_boxplot(
    aes(fill = Metric_color, alpha = Condition),
    width        = 0.18,
    color        = "black",
    outlier.size = 0.6
  ) +
  
  scale_fill_identity() +
  scale_alpha_manual(values = c("noGC7" = 0.35, "GC7" = 0.8)) +
  
  scale_x_discrete(
    labels = c(
      "Control.noGC7"  = "Control\nnoGC7",
      "(AAG)12x.noGC7" = "AAG\nnoGC7",
      "Control.GC7"    = "Control\nGC7",
      "(AAG)12x.GC7"   = "AAG\nGC7"
    )
  ) +
  
  facet_wrap(
    ~ Metric,
    scales   = "free_y",
    labeller = as_labeller(c(
      EGFP_mean  = "GFP (EGFP_mean)",
      mRF12_mean = "RFP (mRF12_mean)"
    ))
  ) +
  
  labs(
    title = "Rep 2 — Mean Fluorescence Distribution (EGFP_mean / mRF12_mean)",
    x     = "Sample (Reporter × Condition)",
    y     = "Mean fluorescence (a.u.)",
    alpha = "Condition"
  ) +
  
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
    strip.text   = element_text(size = 13, face = "bold"),
    legend.title = element_blank()
  )


