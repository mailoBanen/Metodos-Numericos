#include <stdio.h>
#include <stdlib.h>


int main(){

    size_t n=3;

    FILE *file = fopen("archivo.txt", "r");

    if(file == NULL){
        return 1;
    }

    double **arr = malloc(n * sizeof(double *));

    for(int i=0; i<n; i++){
        //Creo memoria para cada renglon y la almaceno en arr
        arr[i] = malloc(n * sizeof(double));

        for(int j=0; j<n; j++){

            if(fscanf(file, "%lf", &arr[i][j]) != 1){ //Seria lo mismo si pongo *(arr+i)+j

                printf("Error\n");
                break;
            }
        }
    }

    for(int i=0; i<n; i++){
        for(int j=0; j<n; j++){

                printf("%lf \n", arr[i][j]);
            
        }
    }

    fclose(file);

    for(int i=0; i<n; i++){

        free(arr[i]);
    }
    free(arr);

    return 0;
}