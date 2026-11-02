# microbiomaSuelo

Análisis de amplicones **16S** o **ITS** de suelo, desde los archivos
`*.extendedFrags.fastq.gz` (FLASH) hasta un reporte HTML.

## Instalación

El repositorio es privado: se necesita un token de GitHub con acceso de lectura.

```r
install.packages(c("remotes", "BiocManager"))
Sys.setenv(GITHUB_PAT = "SU_TOKEN")
remotes::install_github("braddmg/microbiomaSuelo",
                        repos = BiocManager::repositories())
```

## Ejemplo

```r
library(microbiomaSuelo)

analizar_microbioma(
  tipo_datos        = "16S",                # o "ITS"
  carpeta_entrada   = "Data 16s",
  archivo_metadata  = "mdata.csv",
  variable_analisis = "Tratamiento",
  ruta_taxonomia    = "Databases/silva_nr99_v138.2_toGenus_trainset.fa.gz",
  ruta_especies     = "Databases/silva_v138.2_assignSpecies.fa.gz"
)
```

- **Metadata:** CSV con la columna `SampleID` (el texto antes de
  `.extendedFrags.fastq.gz`) y la variable a comparar.
- **ITS:** use `ruta_taxonomia = "Databases/sh_general_release_19.02.2025.tgz"`,
  sin `ruta_especies`.
- **Taxones de interés (opcional):** agregue `archivo_taxones = "lista.csv"`,
  con las columnas `Grupo`, `Genero`, `Especie`, `Clasificacion` y `Funcion`.

El resultado queda en `Resultados_16S/Reporte_final.html`. Para compartirlo,
envíe la carpeta completa: los Krona se abren desde `05_Krona/`.

## Licencia

MIT. Incluye Krona (Ondov et al. 2011, *BMC Bioinformatics* 12:385) bajo su
propia licencia, en `inst/krona/LICENSE_Krona.txt`.
