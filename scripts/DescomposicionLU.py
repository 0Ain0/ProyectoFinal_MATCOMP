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

if __name__ == "__main__":
    A = [
        [2, 1, 1],
        [4, -6, 0],
        [-2, 7, 2]
    ]

    L, U = lu_decomposition(A)

    print("\n===== MATRIZ FINAL L =====")
    print_matrix("L", L)

    print("===== MATRIZ FINAL U =====")
    print_matrix("U", U)
