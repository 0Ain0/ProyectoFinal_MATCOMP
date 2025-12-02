def print_matrix(name, M):
    print(f"{name} =")
    for fila in M:
        print(" ", fila)
    print()


def lu_decomposition(A):
    """
    Descomposición LU sin pivoteo.
    Imprime el proceso paso por paso con índices desde 1 en TODO.
    """
    n = len(A)

    U = [fila[:] for fila in A]  
    L = [[0.0] * n for _ in range(n)]
    for i in range(n):
        L[i][i] = 1.0

    print("\n=== MATRIZ INICIAL A ===")
    print_matrix("A", A)

    
    for k in range(n - 1):
        pivot = U[k][k]
        if pivot == 0:
            raise ValueError("Pivote cero: no se puede hacer LU sin pivoteo.")

        print(f"\n========== PASO k = {k+1} ==========")
        print(f"Pivote = U[{k+1}][{k+1}] = {pivot}\n")

        for i in range(k + 1, n):
            print(f"\n--- Eliminando entrada U[{i+1}][{k+1}] ---")
            print(f"Valor actual U[{i+1}][{k+1}] = {U[i][k]}")

            l_ik = U[i][k] / pivot
            L[i][k] = l_ik

            print(f"Factor l[{i+1}][{k+1}] = {l_ik}")

            print(f"Fila U[{i+1}] antes: {U[i][:]}")

            for j in range(k, n):
                U[i][j] -= l_ik * U[k][j]

            print(f"Fila U[{i+1}] después: {U[i][:]}")
            print_matrix("L", L)
            print_matrix("U", U)

    return L, U


# # SUGERENCIA PARA PODER RESOLVER EL SISTEMA, NO SOLO DESCOMPONER

def solve_lu(L, U, b):
    """
    Resuelve el sistema Ax = b usando las matrices L y U obtenidas.
    1. Ly = b (Sustitución hacia adelante)
    2. Ux = x (Sustitución hacia atrás)
    """
    n = len(L)
    
    # 1. Resolver Ly = b
    y = [0.0] * n
    for i in range(n):
        suma = sum(L[i][j] * y[j] for j in range(i))
        y[i] = b[i] - suma
    
    # 2. Resolver Ux = y
    x = [0.0] * n
    for i in range(n - 1, -1, -1):
        suma = sum(U[i][j] * x[j] for j in range(i + 1, n))
        x[i] = (y[i] - suma) / U[i][i]
        
    return x

# #######################################################################


# if __name__ == "__main__":
#     A = [
#         [2, 1, 1],
#         [4, -6, 0],
#         [-2, 7, 2]
#     ]

#     L, U = lu_decomposition(A)

#     print("\n===== MATRIZ FINAL L =====")
#     print_matrix("L", L)

#     print("===== MATRIZ FINAL U =====")
#     print_matrix("U", U)


# NUEVO MAIN PARA RESOLVER LOS DOS PROBLEMAS DE APLICACIÓN

if __name__ == "__main__":
    # === PROBLEMA 1 ===
    print("\n\n################ PROBLEMA 1 ################")
    A1 = [
        [1, 1, 1],
        [1, 2, 2],
        [1, 2, 3]
    ]
    b1 = [6, 9, 10]
    
    L1, U1 = lu_decomposition(A1)
    x1 = solve_lu(L1, U1, b1)
    
    print("\n>>> SOLUCIÓN PROBLEMA 1 (x):", x1)
    # Debería dar [1.0, 1.0, 3.0] o similar

    # === PROBLEMA 2 (CIRCUITOS) ===
    print("\n\n################ PROBLEMA 2 (CIRCUITOS) ################")
    A2 = [
        [10, -5, 0],
        [-5, 15, -5],
        [0, -5, 10]
    ]
    b2 = [10, 0, 0] # 10V en la primera malla, 0 en las otras
    
    L2, U2 = lu_decomposition(A2)
    x2 = solve_lu(L2, U2, b2)
    
    print("\n>>> SOLUCIÓN PROBLEMA 2 (Corrientes):", x2)