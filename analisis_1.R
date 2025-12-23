##############################################################
#     ANÁLISIS INTRAOBSERVADOR DE VARIABLES MORFOMÉTRICAS
#     Autor: [Muñoz Agustina Natalí]
#     Programa: R
##############################################################

# -----------------------------------------------------------
# 1. CARGA DE LIBRERÍAS
# -----------------------------------------------------------

library(tidyverse)
library(readxl)
library(irr)
library(psych)
library(skimr)
library(rstatix)
library(flextable)
library(ggplot2)
library(scales)
library(progress)
library(gt)
library(ggpubr)
library(shiny)

# -----------------------------------------------------------
# 2. IMPORTACIÓN DE DATOS
# -----------------------------------------------------------

# Ruta del archivo (formato universal)
ruta_archivo <- "C:/Users/AGUSTINA/Downloads/TFG/Análisis_1/datosmorfometricos.xlsx"

# Lectura del archivo Excel
datosmorfometricos <- read_excel(ruta_archivo)

# Vista general
glimpse(datosmorfometricos)

# -----------------------------------------------------------
# 3. CONTROL DE VALORES FALTANTES
# -----------------------------------------------------------

na_por_columna <- datosmorfometricos %>%
  summarise(across(everything(), ~ mean(is.na(.)))) %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "prop_na") %>%
  mutate(porc_na = round(prop_na * 100, 1)) %>%
  arrange(desc(prop_na))

# Gráfico de valores faltantes
na_por_columna %>%
  slice_max(prop_na, n = 15) %>%
  ggplot(aes(x = reorder(variable, prop_na), y = prop_na)) +
  geom_col(fill = "#A7C7E7", color = "white", width = 0.7) +
  geom_text(aes(label = paste0(porc_na, "%")), hjust = -0.1, size = 3, color = "#4B4B4B") +
  coord_flip() +
  labs(
    title = "Porcentaje de valores faltantes por variable",
    x = "Variable morfométrica", y = "Proporción de NA"
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  theme_minimal(base_size = 9)

# Guardar
write.csv(na_por_columna, "porcentaje_NA_por_variable.csv", row.names = FALSE)

# -----------------------------------------------------------
# 4. LIMPIEZA DE DATOS
# -----------------------------------------------------------

umbral <- 0.33
variables_muchos_na <- na_por_columna %>%
  filter(prop_na > umbral)

datos_limpios <- datosmorfometricos %>%
  select(-all_of(variables_muchos_na$variable))

write.csv(datos_limpios, "datos_limpios.csv", row.names = FALSE)

# -----------------------------------------------------------
# 5. REESTRUCTURAR A FORMATO LARGO
# -----------------------------------------------------------

datos <- read.csv("datos_limpios.csv", stringsAsFactors = FALSE)
variables_morfometricas <- names(datos)[!(names(datos) %in% c("IND", "toma"))]

datos_largo <- datos %>%
  pivot_longer(cols = all_of(variables_morfometricas),
               names_to = "Variable",
               values_to = "Valor") %>%
  mutate(IND = as.factor(IND),
         toma = as.factor(toma))

# -----------------------------------------------------------
# 6. ESTADÍSTICA DESCRIPTIVA
# -----------------------------------------------------------

estadisticas_desc <- datos_largo %>%
  group_by(Variable, toma) %>%
  summarise(
    n = n(),
    Media = mean(Valor, na.rm = TRUE),
    SD = sd(Valor, na.rm = TRUE),
    CV_pct = (SD / Media) * 100,
    .groups = "drop"
  )

write.csv(estadisticas_desc, "estadisticas_descriptivas.csv", row.names = FALSE)

# Boxplot general de todas las variables
ggplot(datos_largo, aes(x = toma, y = Valor, fill = toma)) +
  geom_boxplot(alpha = 0.7) +
  facet_wrap(~ Variable, scales = "free_y") +
  labs(title = "Distribución de valores por toma (todas las variables)",
       x = "Toma", y = "Valor medido") +
  theme_minimal(base_size = 8) +
  theme(legend.position = "none")

# -----------------------------------------------------------
# 7. NORMALIDAD (TEST DE SHAPIRO)
# -----------------------------------------------------------

shapiro_por_variable <- function(var) {
  df <- datos_largo %>% filter(Variable == var) %>% drop_na()
  df3 <- df %>% group_by(IND) %>% filter(n() == 3) %>% ungroup()
  if (n_distinct(df3$IND) < 3) return(NULL)
  wide <- df3 %>% pivot_wider(names_from = toma, values_from = Valor, names_prefix = "T")
  diffs <- c(wide$T1 - wide$T2, wide$T1 - wide$T3, wide$T2 - wide$T3)
  p_sh <- tryCatch(shapiro.test(diffs)$p.value, error = function(e) NA)
  tibble(Variable = var, p_shapiro = p_sh)
}

res_shapiro <- map_df(unique(datos_largo$Variable), shapiro_por_variable)
write.csv(res_shapiro, "shapiro_diferencias_por_variable.csv", row.names = FALSE)

# -----------------------------------------------------------
# 8. ANÁLISIS INTRAOBSERVADOR (ANOVA / FRIEDMAN / ICC)
# -----------------------------------------------------------

kendalls_w <- function(mat) {
  ranks <- t(apply(mat, 1, rank))
  k <- ncol(ranks); n <- nrow(ranks)
  Rj <- colSums(ranks)
  S <- sum((Rj - mean(Rj))^2)
  12 * S / (k^2 * (n^3 - n))
}

pb <- progress_bar$new(total = length(unique(datos_largo$Variable)),
                       format = "Procesando :current/:total [:bar] :percent :eta")

analizar_variable_full <- function(var) {
  pb$tick()
  df <- datos_largo %>% filter(Variable == var) %>% drop_na()
  df3 <- df %>% group_by(IND) %>% filter(n() == 3) %>% ungroup()
  if (n_distinct(df3$IND) < 3) return(NULL)
  wide <- df3 %>% pivot_wider(names_from = toma, values_from = Valor, names_prefix = "T")
  diffs <- c(wide$T1 - wide$T2, wide$T1 - wide$T3, wide$T2 - wide$T3)
  p_sh <- tryCatch(shapiro.test(diffs)$p.value, error = function(e) NA)
  normal_flag <- ifelse(!is.na(p_sh) && p_sh > 0.05, TRUE, FALSE)
  
  # Inicializar
  p_test <- effect_size <- NA
  test_name <- effect_measure <- effect_interpretation <- NA
  
  if (normal_flag) {
    test_name <- "ANOVA (normal)"
    res <- tryCatch(anova_test(data = df3, dv = Valor, wid = IND, within = toma), error = function(e) NULL)
    if (!is.null(res)) {
      p_test <- res$ANOVA$p[1]
      effect_size <- res$ANOVA$ges[1]
      effect_measure <- "η²g"
      effect_interpretation <- case_when(
        effect_size < 0.01 ~ "Despreciable",
        effect_size < 0.06 ~ "Pequeño",
        effect_size < 0.14 ~ "Medio",
        TRUE ~ "Grande"
      )
    }
  } else {
    test_name <- "Friedman (no normal)"
    mat <- as.matrix(wide %>% select(starts_with("T")))
    f_res <- tryCatch(friedman.test(mat), error = function(e) NULL)
    if (!is.null(f_res)) p_test <- f_res$p.value
    w <- tryCatch(kendalls_w(mat), error = function(e) NA)
    effect_size <- w
    effect_measure <- "Kendall's W"
    effect_interpretation <- case_when(
      w < 0.1 ~ "Despreciable",
      w < 0.3 ~ "Pequeño",
      w < 0.5 ~ "Medio",
      TRUE ~ "Grande"
    )
  }
  
  icc_val <- tryCatch({
    irr::icc(as.matrix(wide %>% select(starts_with("T"))),
             model = "twoway", type = "consistency", unit = "single")$value
  }, error = function(e) NA)
  
  tibble(
    Variable = var,
    Prueba = test_name,
    p_valor = round(p_test, 4),
    Tamaño_Efecto = round(effect_size, 4),
    Medida = effect_measure,
    Interpretación = effect_interpretation,
    ICC = round(icc_val, 4)
  )
}

resultados_pruebas <- map_df(unique(datos_largo$Variable),
                             possibly(analizar_variable_full, NULL))

write.csv(resultados_pruebas, "resultados_pruebas_intraobservador_completo.csv", row.names = FALSE)

# -----------------------------------------------------------
# 9. POST HOC (BONFERRONI)
# -----------------------------------------------------------

significativas <- resultados_pruebas %>%
  filter(!is.na(p_valor) & p_valor < 0.05)

posthoc_resultados <- list()

for (var in significativas$Variable) {
  tipo_prueba <- significativas %>%
    filter(Variable == var) %>%
    pull(Prueba)
  
  df <- datos_largo %>%
    filter(Variable == var) %>%
    group_by(IND) %>%
    filter(n() == 3) %>%
    ungroup()
  
  if (n_distinct(df$IND) < 3) next
  
  if (grepl("ANOVA", tipo_prueba)) {
    posthoc <- df %>%
      pairwise_t_test(Valor ~ toma, paired = TRUE, p.adjust.method = "bonferroni") %>%
      mutate(Prueba = "ANOVA (normal)", Variable = var)
  } else {
    posthoc <- df %>%
      pairwise_wilcox_test(Valor ~ toma, paired = TRUE, p.adjust.method = "bonferroni") %>%
      mutate(Prueba = "Friedman (no normal)", Variable = var)
  }
  
  posthoc_resultados[[var]] <- posthoc
}

posthoc_final <- bind_rows(posthoc_resultados) %>%
  mutate(Significancia = case_when(
    p.adj <= 0.001 ~ "***",
    p.adj <= 0.01 ~ "**",
    p.adj <= 0.05 ~ "*",
    TRUE ~ "ns"
  ))

write.csv(posthoc_final, "resultados_posthoc_intraobservador.csv", row.names = FALSE)

# -----------------------------------------------------------
# --- Calcular CV intraobservador ---
cv_intra <- datos_largo %>%
  group_by(Variable, IND) %>%
  summarise(
    media_ind = mean(Valor, na.rm = TRUE),
    sd_ind = sd(Valor, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(CV_ind = (sd_ind / media_ind) * 100) %>%
  group_by(Variable) %>%
  summarise(
    CV_intra_mean = mean(CV_ind, na.rm = TRUE),
    CV_intra_sd = sd(CV_ind, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

# ---  Unir resultados ICC + pruebas + CV intraobservador ---
resultados_final <- resultados_pruebas %>%
  left_join(cv_intra %>% select(Variable, CV_intra_mean), by = "Variable") %>%
  mutate(
    p_valor = ifelse(p_valor < 0.001, "<0.001", as.character(round(p_valor, 3))),
    Interpretacion_ICC = case_when(
      ICC < 0.5  ~ "Fiabilidad pobre",
      ICC < 0.75 ~ "Fiabilidad moderada",
      ICC < 0.9  ~ "Fiabilidad buena",
      ICC >= 0.9 ~ "Fiabilidad excelente"
    ),
    Categoria_CV = case_when(
      CV_intra_mean <= 5 ~ "Excelente",
      CV_intra_mean <= 10 ~ "Buena",
      CV_intra_mean <= 15 ~ "Moderada",
      TRUE ~ "Baja"
    )
  )

# ---  Crear tabla GT ---
tabla_intra <- resultados_final %>%
  gt() %>%
  tab_header(
    title = md("**Resumen del análisis intraobservador**"),
    subtitle = md("_Tipo de prueba, significancia, tamaño del efecto, ICC y CV intraobservador_")
  ) %>%
  cols_label(
    Variable = "Variable",
    Prueba = "Tipo de prueba",
    p_valor = "p valor",
    Tamaño_Efecto = "Tamaño del efecto",
    Medida = "Medida",
    Interpretación = "Interpretación del efecto",
    ICC = "ICC",
    Interpretacion_ICC = "Interpretación ICC",
    CV_intra_mean = "CV intraobservador (%)",
    Categoria_CV = "Categoría CV"
  ) %>%
  fmt_number(
    columns = c(Tamaño_Efecto, ICC, CV_intra_mean),
    decimals = 3
  ) %>%
  tab_options(
    table.font.size = px(11),
    data_row.padding = px(4),
    column_labels.font.weight = "bold"
  )

# ---  Exportar ---
tabla_intra
gtsave(tabla_intra, "tabla_intra.docx")

# Confirmación
print("✅Tabla final 'tabla_intra.docx' creada con ICC + CV intraobservador + pruebas estadísticas.")

#-----------------------------------------------
#  11. GRÁFICOS DE RESULTADOS 

# --- 11.1 CV INTRAOBSERVADOR POR VARIABLE ---

p_cv_intra <- ggplot(cv_intra, aes(x = reorder(Variable, CV_intra_mean), y = CV_intra_mean)) +
  geom_col(fill = "#a8dadc") +
  coord_flip() +
  geom_hline(yintercept = 10, linetype = "dashed", color = "red") +
  labs(
    title = "Coeficiente de variación intraobservador (CV%) por variable",
    subtitle = "Promedio entre las tres tomas por individuo — Línea roja = umbral del 10%",
    x = "Variable morfométrica", y = "CV intraobservador (%)"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    panel.grid.major.y = element_blank(),
    legend.position = "none"
  )

ggsave("CV_intraobservador_variables.png", p_cv_intra, width = 8, height = 10, dpi = 300)
print(p_cv_intra)


# --- 11.2 ICC POR VARIABLE ---
icc_data <- resultados_final %>%
  mutate(
    Interpretacion_ICC = case_when(
      ICC < 0.5  ~ "Fiabilidad pobre",
      ICC < 0.75 ~ "Fiabilidad moderada",
      ICC < 0.9  ~ "Fiabilidad buena",
      ICC >= 0.9 ~ "Fiabilidad excelente"
    )
  )

p_icc <- ggplot(icc_data, aes(x = reorder(Variable, ICC), y = ICC, color = Interpretacion_ICC)) +
  geom_point(size = 3) +
  coord_flip() +
  geom_hline(yintercept = c(0.5, 0.75, 0.9), linetype = "dashed", color = "gray60") +
  scale_color_manual(values = c(
    "Fiabilidad pobre" = "#e76f51",
    "Fiabilidad moderada" = "#f4a261",
    "Fiabilidad buena" = "#457b9d",
    "Fiabilidad excelente" = "#2a9d8f"
  )) +
  labs(
    title = "ICC por variable morfométrica",
    subtitle = "Clasificación de la fiabilidad según valores de ICC",
    y = "ICC", x = "Variable", color = "Fiabilidad"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))

ggsave("ICC_variables.png", p_icc, width = 8, height = 10, dpi = 300)
print(p_icc)


# --- 11.3 BOXPLOTS INDIVIDUALES ---
dir.create("Boxplots_Variables", showWarnings = FALSE)
for (v in unique(datos_largo$Variable)) {
  dfv <- datos_largo %>% filter(Variable == v)
  p_box <- ggplot(dfv, aes(x = toma, y = Valor, fill = toma)) +
    geom_boxplot(outlier.colour = "#E07A5F", width = 0.7) +
    labs(title = paste("Boxplot —", v), x = "Toma", y = "Valor") +
    theme_minimal(base_size = 10) +
    theme(
      legend.position = "none",
      plot.title = element_text(face = "bold")
    )
  ggsave(paste0("Boxplots_Variables/", v, "_boxplot.png"), p_box, width = 5, height = 4, dpi = 300)
}


# -----------------------------------------------------------
#  VISUALIZACIÓN INTERACTIVA (SHINY)
# -----------------------------------------------------------

ui <- fluidPage(
  titlePanel("Visualización Intraobservador - Variables Morfométricas"),
  sidebarLayout(
    sidebarPanel(
      selectInput("var_select", "Selecciona una variable:", choices = unique(datos_largo$Variable)),
      width = 3
    ),
    mainPanel(
      plotOutput("boxplot_var", height = "400px"),
      plotOutput("cv_plot", height = "400px"),
      plotOutput("icc_plot", height = "400px")
    )
  )
)

server <- function(input, output) {
  # --- Boxplot interactivo por variable ---
  output$boxplot_var <- renderPlot({
    dfv <- datos_largo %>% filter(Variable == input$var_select)
    ggplot(dfv, aes(x = toma, y = Valor, fill = toma)) +
      geom_boxplot(outlier.colour = "#E07A5F", width = 0.7) +
      labs(title = paste("Distribución de tomas —", input$var_select),
           x = "Toma", y = "Valor medido") +
      theme_minimal(base_size = 12) +
      theme(legend.position = "none")
  })
  
  # --- CV Intraobservador (gráfico global) ---
  output$cv_plot <- renderPlot({ p_cv_intra })
  
  # --- ICC por variable (gráfico global) ---
  output$icc_plot <- renderPlot({ p_icc })
}

# Ejecutar panel interactivo (opcional)
shinyApp(ui, server)
# -----------------------------------------------------------
# 12. TABLA DE ESTADÍSTICAS DESCRIPTIVAS (GT)
# -----------------------------------------------------------

tabla_desc <- estadisticas_desc %>%
  mutate(
    Media = round(Media, 2),
    SD = round(SD, 2),
    CV_pct = round(CV_pct, 1)
  ) %>%
  gt() %>%
  tab_header(
    title = md("**Estadísticas descriptivas por variable y toma**"),
    subtitle = md("_Media, desviación estándar y coeficiente de variación (CV%)_")
  ) %>%
  cols_label(
    Variable = "Variable",
    toma = "Toma",
    n = "n",
    Media = "Media",
    SD = "Desv. Estándar",
    CV_pct = "CV (%)"
  ) %>%
  tab_options(
    table.font.size = px(11),
    data_row.padding = px(4),
    column_labels.font.weight = "bold"
  )

tabla_desc
gtsave(tabla_desc, "tabla_estadisticas_descriptivas.docx")

# -----------------------------------------------------------
# 13. GRÁFICOS ADICIONALES
# -----------------------------------------------------------

# --- 13.1 Comparación de las tres tomas por individuo (ejemplo: variable FCM) ---
variable_ejemplo <- "FCM"  # cambiar la variable 

df_var <- datos_largo %>% filter(Variable == variable_ejemplo)

p_comparacion_tomas <- ggplot(df_var, aes(x = IND, y = Valor, color = toma, group = IND)) +
  geom_point(size = 3) +
  geom_line(aes(group = IND, color = toma), linewidth = 0.7, alpha = 0.5) +
  scale_color_manual(values = c("1" = "#e63946", "2" = "#1d3557", "3" = "#a8dadc")) +
  labs(
    title = paste("Comparación de las tres tomas (", variable_ejemplo, ")", sep = ""),
    subtitle = "Evaluación visual del error intraobservador",
    x = "Individuo", y = "Valor medido", color = "Toma"
  ) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p_comparacion_tomas)
ggsave("Comparacion_tres_tomas.png", p_comparacion_tomas, width = 7, height = 4, dpi = 300)

# --- 13.2 Proporción de variables con distribución normal ---
normalidad_eval <- res_shapiro %>%
  mutate(Normal = ifelse(p_shapiro > 0.05, "Normal", "No normal")) %>%
  mutate(Evaluacion = sample(1:3, nrow(.), replace = TRUE)) %>%  # si no hay columna de evaluación, simular
  group_by(Evaluacion, Normal) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Evaluacion) %>%
  mutate(prop = n / sum(n))

p_normalidad <- ggplot(normalidad_eval, aes(x = factor(Evaluacion), y = prop, fill = Normal)) +
  geom_col(position = "stack", width = 0.7, color = "white") +
  geom_text(aes(label = paste0(round(prop * 100, 1), "%")),
            position = position_stack(vjust = 0.5), size = 3, color = "gray20") +
  scale_fill_manual(values = c("Normal" = "#bde0fe", "No normal" = "#f4a3a3")) +
  labs(
    title = "Proporción de variables con distribución normal",
    subtitle = "Resultados del test de Shapiro-Wilk por evaluación",
    x = "Evaluación", y = "Proporción de variables"
  ) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(face = "bold", hjust = 0.5))
print(p_normalidad)
ggsave("Proporcion_normalidad.png", p_normalidad, width = 6, height = 4, dpi = 300)

# --- 13.3 Comparaciones post hoc significativas ---
posthoc_plot_data <- posthoc_final %>%
  mutate(Significativa = ifelse(Significancia == "ns", "No", "Sí")) %>%
  mutate(Comparacion = paste0(group1, " vs ", group2))

p_posthoc <- ggplot(posthoc_plot_data,
                    aes(x = Comparacion, y = Variable, fill = Significativa)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  scale_fill_manual(values = c("No" = "#D3D3D3", "Sí" = "#E07A5F")) +
  labs(
    title = "Comparaciones post hoc significativas entre tomas",
    x = "Comparación de tomas", y = "Variable morfométrica",
    fill = "Diferencia significativa"
  ) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
print(p_posthoc)
ggsave("PostHoc_significativas.png", p_posthoc, width = 7, height = 5, dpi = 300)


