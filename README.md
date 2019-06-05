# sargazo-plot
Script en Python para graficar datos en mallas no estructuradas.
Utiliza multiprocesamiento (varios CPU's para realizar la graficación)
Toma datos de elevación del mar de un archivo netCDF y posición de partículas (sargazo) de un archivo .mat y los grafica.
Además agrega información tipo *shape* para delimitar los mapas (estados, países, costas, etc.)
Agrega una barra de color personalizada (my_color.py)

## Bibliotecas necesarias 
* Cartopy
* netCDF4
* scipy
* numpy
* matplotlib

## Probando
Se deben indicar los siguientes parámetros:
variable: zeta, temp
título: texto de título
ruta: Ruta del archivo .nc
output: nombre base de los archivos de salida (sin extensión)
n_cpus: número de CPU que se usará para el procesamiento

## Para contribuir
Para información sobre cómo ayudar o publicar errores vea CONTRIBUTING.md.

## Autores
https://github.com/ma-robles

## Agradecimientos
