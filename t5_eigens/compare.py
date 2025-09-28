import numpy as np
from scipy.linalg import eigh

def leer_y_calcular_eigenvalores(archivo, k, mayores=True):
    """
    Versión simplificada para un archivo específico
    """
    # Leer matriz
    with open(archivo, 'r') as f:
        n = int(f.readline().strip().split()[0])
        matriz = []
        for i in range(n):
            linea = f.readline().strip()
            fila = [float(x) for x in linea.split()]
            matriz.append(fila)
    
    matriz = np.array(matriz)
    
    # Calcular todos los eigenvalores
    eigenvalores = eigh(matriz, eigvals_only=True)
    
    # Seleccionar k más grandes o más pequeños
    if mayores:
        return np.sort(eigenvalores)[-k:][::-1]  # k más grandes (orden descendente)
    else:
        return np.sort(eigenvalores)[:k]  # k más pequeños (orden ascendente)


archivo = "archivos_tarea5/Eigen_50x50.txt"
k = 7

# Calcular los 5 eigenvalores más grandes
eigenvalores_grandes = leer_y_calcular_eigenvalores(archivo, k, mayores=True)
print(f"Los {k} eigenvalores más grandes: {eigenvalores_grandes}")

# Calcular los 3 eigenvalores más pequeños
eigenvalores_pequenos = leer_y_calcular_eigenvalores(archivo, k, mayores=False)
print(f"Los {k} eigenvalores más pequeños: {eigenvalores_pequenos}")