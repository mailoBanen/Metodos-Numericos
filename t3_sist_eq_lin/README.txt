# Instrucciones para ejecutar los programas

La carpeta incluye un codigo por cada problema, ademas de dos archivos para la libreria
que contiene las funciones que solucionan los sistemas.

Los algoritmos están en operations_mat.c 
Tambien se incluyen los txt para cada sistema.


Todos los programas de la tarea se corren de la misma manera.
Desde el mismo directorio donde están los archivos, compilar cada archivo por separado con: 

" gcc main.c libreria.c -o main -lm "

gcc p1_t3.c operations_mat.c -o main -lm
gcc p2_t3.c operations_mat.c -o main -lm
gcc p3_t3.c operations_mat.c -o main -lm
gcc p5_t3.c operations_mat.c -o main -lm
gcc pextra_t3.c operations_mat.c -o main -lm

Una vez compilado, ejecutar con: 

./main
