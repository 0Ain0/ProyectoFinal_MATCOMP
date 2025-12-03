def sor(A, b, x0, omega, tol=1e-6, max_iter=1000):
    """
    Parámetros:
        A        : matriz (lista de listas) n x n, diagonalmente dominante
        b        : vector (lista) de tamaño n
        x0       : vector de aproximación inicial (lista) de tamaño n
        omega    : factor de relajación (0 < omega < 2)
        tol      : tolerancia para criterio de paro
        max_iter : máximo número de iteraciones

    Regresa:
        x        : aproximación de la solución
        iteraciones : número de iteraciones realizadas
        error    : error final (máxima diferencia entre iteraciones)
    """
    n = len(A)
    x = x0[:] 

    iteraciones = 0
    error = tol + 1.0

    print("\n=== INICIANDO MÉTODO SOR ===")
    print(f"omega = {omega}, tol = {tol}, max_iter = {max_iter}")
    print(f"x0 = {x0}\n")

    while error > tol and iteraciones < max_iter:
        print(f"\n===== Iteración {iteraciones} =====")
        error = 0.0

        for i in range(n):
            suma = 0.0
            xprev = x[i]

            for j in range(n):
                if j != i:
                    suma += A[i][j] * x[j]

            x[i] = (1 - omega) * x[i] + omega * (b[i] - suma) / A[i][i]

            diff = abs(x[i] - xprev)
            if diff > error:
                error = diff

            print(f"x[{i+1}] previo: {xprev:.6f}  →  nuevo: {x[i]:.6f}   (cambio = {diff:.6e})")

        print(f"Error de la iteración: {error:.6e}")
        print(f"Vector x = {x}")

        iteraciones += 1

    print("\n=== FIN DEL MÉTODO SOR ===")
    return x, iteraciones, error

if __name__ == "__main__":
    A = [
        [ 4.0,  1.0,  0.0],
        [ 1.0, 4.0,  1.0],
        [0.0,  1.0,  4.0]
    ]

    b = [5.0, 6.0, 5.0]

    x0 = [0.0, 0.0, 0.0]

    omega = 1.2

    x, iters, err = sor(A, b, x0, omega, tol=0.07, max_iter=1000)

    print("Solución aproximada x:", x)
    print("Iteraciones:", iters)
    print("Error final:", err)

