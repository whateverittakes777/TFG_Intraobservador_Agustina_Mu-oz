# 🦴 MORFOMETRÍA TRADICIONAL Y FOTOGRAMETRÍA DE CORTO ALCANCE  
### Validación metodológica para el análisis osteométrico en restos óseos humanos de Norpatagonia

**Autora:** Muñoz Agustina Natalí
**Directora:** Béguelin Marien
**Co Director:** Citton Paolo
**Carrera:** Licenciatura en Criminología y Ciencias Forenses  
**Universidad:** Universidad Nacional de Río Negro (UNRN)  
**Año:** 2025  

---

## 📘 Descripción del proyecto

Este repositorio acompaña el Trabajo Final de Grado titulado **“Morfometría tradicional y fotogrametría de corto alcance: validación metodológica para el análisis osteométrico en restos óseos humanos de Norpatagonia”**.

El objetivo principal del estudio es evaluar y comparan dos formas de registro de variables métricas, las medidas morfométricas tomadas de manera analógica y las mismas medidas tomadas de manera digital, a partir de su modelo tridimensional. Todo el relevamiento se desarrolla en una muestra esquelética contemporánea perteneciente a la Colección Osteológica de referencia Norpatagonica (CORN, Vázquez et al., 2025), de modo que paralelamente se genere un repositorio digital de la colección. Los análisis se realizaron bajo un enfoque estadístico reproducible implementado en **R**.

El análisis del error intraobservador en el relevamiento de las variables morfométricas incluye el análisis de consistencia de medidas, pruebas de normalidad, ANOVA/Friedman, tamaño del efecto, coeficiente de correlación intraclase (ICC) y visualizaciones gráficas.

---

## 📂 Estructura del repositorio

📁 Analisis_Morfometrico/
├── Analisis_Intraobservador.R           # Script principal en R
├── resultados_pruebas_intraobservador_completo.csv
├── resultados_posthoc_intraobservador.csv
├── estadisticas_descriptivas.csv
├── shapiro_diferencias_por_variable.csv
├── CV_promedio.png                      # Gráfico CV
├── ICC_por_variable.png                 # Gráfico ICC
├── Boxplots_Variables/                  # Boxplots individuales
└── README.md                            # Descripción del proyecto

## 🧠 Metodología estadística

El script `Analisis_Intraobservador.R` implementa los siguientes pasos:

1. **Carga y limpieza de datos** desde un archivo Excel (`datosmorfometricos.xlsx`).
2. **Control de valores faltantes** y filtrado de variables con más del 33% de NA.
3. **Reestructuración** a formato largo para análisis.
4. **Estadística descriptiva** (media, SD, CV%).
5. **Prueba de normalidad (Shapiro-Wilk)** sobre las diferencias de tomas.
6. **Análisis intraobservador**  
   - ANOVA de medidas repetidas (si distribución normal).  
   - Prueba de Friedman y W de Kendall (si no normal).  
   - Cálculo del ICC (modelo dos vías, consistencia, unidad individual).  
7. **Análisis Post Hoc** (t-test pareado o Wilcoxon con corrección Bonferroni).  
8. **Visualización de resultados** mediante gráficos y tablas interactivas.  

---

## 📊 Visualizaciones incluidas

- **Coeficiente de variación promedio (CV%) por variable**  

- **Fiabilidad intraobservador (ICC) por variable**  

- **Boxplots individuales por variable y toma**  
  (disponibles en la carpeta `Boxplots_Variables/`)

---

## ⚙️ Requisitos

Para ejecutar el script, se necesita tener instalado **R (≥ 4.2)** y las siguientes librerías:

```r
tidyverse, readxl, irr, psych, skimr, rstatix, flextable,
ggplot2, scales, progress, gt, ggpubr
Instalación rápida:

install.packages(c("tidyverse", "readxl", "irr", "psych", "skimr",
                   "rstatix", "flextable", "ggplot2", "scales",
                   "progress", "gt", "ggpubr"))
## 🚀 Ejecución

1. Clonar el repositorio:
   ```bash
   git clone https://github.com/whateverittakes777/TFG_Intraobservador_Agustina_Mu-oz.git


🧾 Licencia
Este trabajo se distribuye bajo la licencia Creative Commons Atribución-NoComercial 4.0 Internacional (CC BY-NC 4.0).
Esto significa que cualquier persona puede compartir y adaptar el material, siempre que se otorgue el crédito correspondiente y no se utilice con fines comerciales.
