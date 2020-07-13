# Análisis espacial del dengue en Costa Rica 

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

Este repositorio contiene toda la información referida al proyecto de investigación del curso SP-1626: Estadística Bayesiana, desarrollado por el estudiante de la Maestría Académica en Estadística: Pablo Vivas Corrales (<pablo.vivas@ucr.ac.cr>) para el primer semestre del 2020. 

## Archivos de repositorio

* Archivos
  * [`Análisis.R`](#análisis)
  * [`Escrito.pdf`](#escrito)
  * [LICENSE](#licencia)

* Carpetas
  * [Bibliografía](#bibliografía)
  * [Datos](#datos)
  * [Trabajo escrito](#trabajo-escrito)


### Archivos
---
### Análisis

El archivo `análisis.r` es un documento donde se presenta el código `R` utilizado para generar los resultados de este artículo. Comprende desde la lectura de los datos hasta la creación de los modelos y las figuras. Se utilizó el software R en su versión 3.6.3 (2020-02-29) y el IDE RStudio en su versión 1.2.5033.

### Escrito

Es una copia del archivo `Escrito.pdf` de la carpeta [Trabajo Escrito](#trabajo-escrito) donde se presenta el artículo. Se pone a disposición este archivo, para mayor facilidad en la búsqueda de la versión final este artículo. 

### Licencia

El código usado y presentado en este repositorio tiene una licencia [MIT](https://opensource.org/licenses/MIT)

---
### Carpetas
---
### Bibliografía

En esta carpeta se encuentran archivos `.pdf` que se utilizaron como principales referencias bibliográficas para el proyecto. Estas referencias fueron introducidas al software [Mendeley](https://www.mendeley.com/?interaction_required=true) para crear el archivo `Referencias.bib` de la carpeta [Trabajo Escrito](#trabajo-escrito) utilizado para hacer las citas bibliográficas correspondientes.

### Datos

Aquí se almacenan los datos utilizados en el análisis. Destacan principalmente los siguientes archivos:

* `Cantones`: Información geográfica de la división territorial de Costa Rica en cantones. Esta información se obtuvo de la página [http://daticos-geotec.opendata.arcgis.com](http://daticos-geotec.opendata.arcgis.com/datasets/249bc8711c33493a90b292b55ed3abad_0)
* `Datos cantonales.xlsx`: Compendio de indicadores cantonales que se utilizaron en el análisis, junto con la información espacial. Los indicadores utilizados son:
  * `pob`: Población del 2019. INEC. [LINK](http://services.inec.go.cr/proyeccionpoblacion/frmproyec.aspx)
  * `casos`: Casos de dengue. Ministerio de Salud. [LINK](https://www.ministeriodesalud.go.cr/index.php/vigilancia-de-la-salud/analisis-de-situacion-de-salud)
  * `dengue`: Tasa de Dengue (100.000 habitantes) del 2019. Ministerio de Salud. [LINK](https://www.ministeriodesalud.go.cr/index.php/vigilancia-de-la-salud/analisis-de-situacion-de-salud)
  * `tugurio`: Porcentaje de viviendas de tipo tugurio en el 2011. INEC. [LINK](http://sistemas.inec.cr:8080/bincri/RpWebEngine.exe/Portal?BASE=2011&lang=esp)
  * `densidad`: Densidad de la población del 2011. INEC. [LINK](http://sistemas.inec.cr:8080/bincri/RpWebEngine.exe/Portal?BASE=2011&lang=esp)
  * `residuos`: Porcentaje de viviendas que eliminan los residuos sólidos por camión recolector en el 2011. INEC. [LINK](http://sistemas.inec.cr:8080/bincri/RpWebEngine.exe/Portal?BASE=2011&lang=esp)
  * `acueducto`: Porcentaje de viviendas con acueducto en el 2011. INEC. [LINK](http://sistemas.inec.cr:8080/bincri/RpWebEngine.exe/Portal?BASE=2011&lang=esp)

### Trabajo Escrito

En esta carpeta se almacenan todos los archivos necesarios para la creación del documento `Escrito.pdf`. Estos archivos se crearon mediante el software [TexMaker](https://www.xm1math.net/texmaker/). Además contiene los gráficos creados en el análisis, estos últimos son de la forma `F#.pdf` y se describen a continuación:

* `F11.pdf`: Tasa de dengue (100.000 habitantes) por cantón, 2019.
* `F12.pdf`: Tasa de dengue en quintiles por cantón, 2019.
* `F21.pdf`: Estructuras de vecinos: Reina.
* `F22.pdf`: Estructuras de vecinos: Knn(4).
* `F31.pdf`: Tasa de dengue y tasa de dengue espacialmente rezagada.
* `F32.pdf`: Valores de influencia por cantón, 2019.
* `F4.pdf`: Residuales de los cuatro modelos ajustados
* `F5.pdf`: Riesgo relativo por cantón, 2019.
* `F6.pdf`: Intervalos de confianza del 95% para el riesgo relativo de dengue por cantón, 2019.
* `FA1.pdf`: Casos de dengue por cantón, 2019.
* `FA2.pdf`: Desviación estándar de la tasa de dengue por 100.000 habitantes por cantón, 2019.
* `FA3.pdf`: Porcentaje de viviendas de tipo tugurio por cantón, 2011.
* `FA4.pdf`: Densidad de la población por cantón, 2011.
* `FA5.pdf`: Porcentaje de viviendas que eliminan los residuos sólidos por camión recolector por cantón, 2011.
* `FA6.pdf`: Porcentaje de viviendas con acueducto por cantón, 2011.
* `FA7.pdf`: Probabilidad de cada supuesto en tasa de dengue por 100.000 habitantes
* `FA8.pdf`: Valores observados y valores esperados de la tasa de dengue por 100.000 habitantes
* `FA9.pdf`: Mapa de probabilidad de Chownoysky
* `FA10.pdf`: Estimador empírico global de Bayes (mm), Suavizado empírico de Bayes (ml) con parámetros $\nu = 0,255$ y $\alpha = 0,253$; estimador empírico local de Bayes (mm.local) y SMR