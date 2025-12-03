<br>
<br>
<br>
<br>
<br>
<br>
<br>

<div style="text-align: center;">

![Logo Universidad Panamericana](assets/logoUP.png)

<br>

# **Universidad Panamericana**

### Escuela de Ingeniería

### Matemáticas de la Computación

<br>

# **Descomposición LU para sistemas lineales y Método Gauss-Seidel Acelerado**

<br>

</div>


### **Profesor!** 
- David Iván Morales Huerta

### **Integrantes!**
- Ain Bolaños Cortés
- Diego Castañeda Ladrón de Guevara
- Santiago Espinosa Ollivier
- Cecilia Gaona Vidales 
- Santiago Medina Domínguez 

<br>

### **Entrega!** 
- 2 de Diciembre, 2025

<br>
<br>
<br>


<div style="page-break-after: always;"></div>


# 1. Resumen General del Reporte

<div style="text-align: justify;">

Este reporte desarrolla, explica, justifica e implementa los métodos de descomposición LU y el método Gauss-Seidel para la resolución de sistemas de ecuaciones lineales. Con un enfoque en su resolución automatizada computable. Dando así la justificación teórico mátemática y su implementación algorítmica en Python.
<br>

</div>

<br>
<br>

# 2. Introducción

<div style="text-align: justify;">

Los sistemas de ecuaciones lineales, así como su representación matricial son un problema típico de la computación, en áreas relevantes como la Visión de Computadora, la criptografía, los gráficos y el machine learning, así como en las simulaciones físicas y el análisis matemático.
<br> 
Y de la misma manera en que los conceptos del álgebra lineal son relevantes para la expansión de los límites de la computación, el poder computacional extiende las posibilidades de procesamiento y resolución de sistemas lineales extensos.
<br>
<br>
Es por esto que, en este reporte analizaremos dos métodos para resolver estos problemas de matrices donde la computación puede fallar, debido a la complejidad y logitud de los sistemas de eciaciones lineales.
<br>

</div>

<br>
<br>

# 3. Primer Método: Descomposición LU

## 3.1. Introducción del Método

<div style="text-align: justify;">

Como ya vimos, el problema reside en el análisis numérico y la computación científica, y se expresa en forma matricial $Ax = b$.
<br>
Este tipo de sistemas aparece en casi todos los ámbitos de la ingeniería y las ciencias aplicadas, desde el análisis de circuitos eléctricos, el diseño estructural, la modelización financiera, el aprendizaje automático, etc. Si bien es cierto que existen métodos iterativos (por ejemplo, el método de Gauss-Seidel que trabajará tu compañero) que encuentran una solución aproximada, los métodos directos se limitan a encontrar una solución exacta (dentro de la precisión de la máquina) en un número finito de pasos.
<br>
<br>
El método de Descomposición LU es uno de los métodos directos más fundamentales y eficientes. No es simplemente un algoritmo, sino una factorización de la matriz A que transforma el problema $Ax = b$ en un proceso de dos pasos mucho más simple de resolver. Esta investigación se centrará en la fundamentación teórica, la formulación matemática, el análisis de estabilidad y las ventajas de este crucial método.

</div>

<!-- ### 3.1.1 Planteamiento del Problema
### 3.1.2 Contexto Histórico
### 3.1.3 Metodología de Solución -->

<br>

## 3.2 Revisión de la Literatura

<div style="text-align: justify;">

La factorización LU construye la formalización directa del método de eliminación gaussiana, que es un algoritmo muy antiguo con una historia de más de dos mil años. Sin embargo, hoy en día vista como una factorización de matrices resultó crucial para el desarrollo de las computadoras digitales. Alan Turing llegó a utilizar esta factorización como una de las herramientas centrales tratando de encontrar el análisis de errores en las computaciones basados en matrices, esto visto en su trabajo de 1948 "Rounding-Off Errors in Matrix Processes".
<br>
<br>
Los textos considerados canónicos de álgebra lineal numérica, los de Golub & Van Loan (2013) o los de Trefethen & Bau (1997), consideran la descomposición LU como la herramienta ideal para resolver sistemas lineales densos. Es precisamente la base sobre la cual se construyen algoritmos más complejos y es el método de facto que se encuentra implementado en paquetes de software de alto rendimiento como LAPACK y bibliotecas como NumPy (SciPy) o MATLAB (Trefethen & Bau, 1997). 
<br>
La literatura no sólo se esfuerza por describir el algoritmo básico, sino que lo hace por describir, y exclusivamente por el mismo, las versiones estables que tienen pivoteo, indispensable para la aplicación práctica de estas ideas.

</div>

<br>

## 3.3 Metodología

<div style="text-align: justify;">

Este método numérico requiere primero, para su justificación, las siguientes herramientas y conceptos matemáticos. 

</div>

<br>

### 3.3.1 Definiciones Matemáticas

***Matriz Cuadrada:***

<div style="text-align: justify;">

Sea $A\in M_{n\times m}(F)$, con $M_{n\times m}$ el conjunto de matrices del campo $F$ con dimensiones $n\times m$. En el caso en que $n = m$ decimos que $A$ es una matriz cuadrada. Para el futuro de este reporte, se tomará $F = \mathbb{R}$.

</div>

***Vectores Canónicos:***

<div style="text-align: justify;">

Sea $e_k$ un vector de dimensión $n\times 1$, decimos que este es $k$-ésimo vector canónico de $\mathbb{R}^n$ si cumple que con $x_i$ sus entradas con $i\in\{1, \dots, n\}$, entonces $x_k = 1$ y $x_i = 0$ con $i \neq k$. Así pues se puede ver de la siguiente manera: 

</div>

$$e_k = \begin{bmatrix}
   0 \\
   \vdots \\
   1 \\
   \vdots \\
   0
\end{bmatrix}$$

<div style="text-align: justify;">

Esto es útil también para referirnos a partes de una matriz $A$ lo suficientemente grande. Dado que:

</div>

- $Ae_k = k$-ésima columna de $A$
- $e^T_k A = k$-ésima fila de $A$ 

<div style="text-align: justify;">

Además, en dado caso de tener un vector $v$ con dimensiones $n\times 1$, la multiplicación $v e_k^T$ será una matriz cuadrada $n\times n$ con todos los valores $0$ a excepción de la columna $k$ que contendrá los valores de $v$. 

</div>

<br>

### 3.3.2 Justificación Teórica de la Descomposición LU

<div style="text-align: justify;">

Sea $A\in M_{n\times n}(\mathbb{R})$, lo que buscamos es la secuencia finita $\{ L_k\in M_{n\times n}: 1\leq k \leq n\}$ tales que para $x_k$ la $k$-ésima columna de la matriz $A$ ($x_k = Ae_k$):

</div>

$$
x_k = \begin{bmatrix}
x_{1k} \\
\vdots \\
x_{kk} \\
x_{k+1, k} \\
\vdots \\
x_{nk}
\end{bmatrix} \xrightarrow{\ \ \ \ L_k\ \ \ \ }
L_kx_k = \begin{bmatrix}
x_{1k} \\
\vdots \\
x_{kk} \\
0 \\
\vdots \\
0
\end{bmatrix}
$$

<br>

<div style="text-align: justify;">

Esto con el objetivo de obtener a partir de esta cantidad finita de $L_k$ una matriz $U$ triangular superior. De manera similar a la eliminación Gaussiana, obtenemos pues, suponiendo que $x_{kk} \neq 0$ los siguientes valores: 
$$l_{jk} = \frac{x_{jk}}{x_{kk}}\ \ \ \ \ \ \ \ (k<j\leq m)$$
Esto significa que tenemos que pedirle a $A$ que todos sus elementos en la diagonal **sean no nulos**. Con estas entradas entradas definamos $l_k = [0, \dots, 0, l_{k+1, k}, \dots, l_{nk}]^T$, y así obtenemos que $L_k$ se verá de la siguiente manera:

</div>

$$
L_k = \begin{bmatrix}
1 \\ 
& \ddots \\
& & 1 \\
& & -l_{k+1, k} & 1 \\
& & \vdots & & \ddots \\
& & -l_{nk} & & & 1
\end{bmatrix} = I - l_k e_k^T
$$

<br>

<div style="text-align: justify;">

- Con esto, generamos $k$ matrices $L_k$ tales que: 

</div>

$$
L_{n-1} \cdots L_2 L_1 A = U
$$

<br>

<div style="text-align: justify;">

- La matriz deseada, y más aún, despejando para $A$ obtenemos: 

</div>

$$
A = L_1^{-1}L_2^{-1}\cdots L_{n-1}^{-1} U
$$

<br>

<div style="text-align: justify;">

- Por como está estructurada cada $L_k$, podemos ver que $L_k^{-1} = I + l_k e_k$, con lo que nos podemos tomar: 

</div>

$$
L = L_1^{-1}L_2^{-1}\cdots L_{n-1}^{-1} = I+ \sum_{k = 1}^{n-1}l_k e_k^T = \begin{bmatrix}
1 \\ 
l_{21} & \ddots \\
& & 1 \\
\vdots & & l_{k+1, k} & 1 \\
& & \vdots & & \ddots \\
l_{n1}& & l_{nk} & & & 1
\end{bmatrix}
$$

<br>

<div style="text-align: justify;">

Con esto obtuvimos pues $A = LU$, con $L$ matriz triangular inferior y $U$ matriz triangular superior. Ahora, supongamos que existen $L_1, L_2$ matrices triangulares inferiores y $U_1, U_2$ matrices triangulares superiores tales que:

</div>

$$
A = L_1U_1 = L_2 U_2\ \ \ \ \Longrightarrow\ \ \ \ L_2^{-1}L_1 = U_2 U_1^{-1}
$$

<br>

<div style="text-align: justify;">

- Por cerradura de las matrices triangulares, $L_2^{-1}L_1$ es triangular inferior y $U_2 U_1^{-1}$ es triangular superior, así  por la igualdad:

</div>

$$
L_2^{-1}L_1 = U_2 U_1^{-1} = I\ \ \ \ \Longrightarrow\ \ \ \ L_1 = L_2\  \land\ U_1 = U_2
$$

<br>

<div style="text-align: justify;">

Así pues, el despeje $LU$ fue único para $A$. Como observación final, es fácil ver que  $L$ es fácil de obtener, al no ser necesario ningún despeje u operación entre matrices. Mientras que $U$ se puede obtener al tiempo que se genera $L$. Con esto **queda demostrada la existencia y unicidad de $L, U$**, además de dar una pauta para la obtención de dichas matrices.

</div>

<br>

### 3.3.3 Algoritmo de la Descomposición LU

<div style="text-align: justify;">

Como se puede apreciar en la deducción matemática, las operaciones de inversa para las $L_i$ son sumamente simples de efectuar, así como el despeje de cada entrada de $L$ se puede obtener de manera directa, lo que hace al algoritmo sumamente efectivo dado que es de matrices. Con esto, Lloyd N. y David Bau, agregan que no es necesario almacenar $A, L$ o $U$ en distintas matrices. 
<br>
<br>
Generamos una matriz partiendo de la identidad, y vamos a encontrar todos los $l_{jk}$, es decir, los factores de la triangular inferior estricta de $L$ para con ellos actualizar las entradas de $U$, haciendo la Eliminación Gaussiana. 
<br>
<br>
El pseudo código a continuación supone que la matriz $A$ dada es válida, es decir, $\forall k\leq n$, $A_{kk} \neq 0$. Con $n$ claro el tamaño de la matriz $A$, es decir $A \in M_{n\times n}(F)$.

</div>

<br>

**3.3.3.1 Pseudocódigo**

```
tome U, L matrices n x n

U = A, L = I

para k = 1 hasta n - 1

   para i = k + 1 hasta n
      l_ik = u_ik/u_kk
      ik = factor
      
      para j = k hasta n
         u_ij -= l_ik * u_kj
      fin para
   fin para
fin para

retornar L, U
``` 

<div style="text-align: justify;">

> Es importante notar que este algoritmo tiene **complejidad $O(n^3)$**, lo que implica que su funcionamiento en menos de $1$ segundo, estará limitado a matrices de a lo más $n = 10^3$, suponiendo una cantidad de operaciones $\approx 10^9$ por segundo.

</div>

<br>
<br>

# 3.4. Aplicación de la Descomposición LU

## 3.4.1 Problema 1: Un Sistema $3 \times 3$ Básico

<div style="text-align: justify;">

Se propone resolver el siguiente sistema simple para demostrar la trazabilidad manual del algoritmo:

</div>

$$A = \begin{bmatrix}
1 & 1 & 1 \\ 
1 & 2 & 2 \\ 
1 & 2 & 3 
\end{bmatrix}, \quad b = \begin{bmatrix} 6 \\ 9 \\ 10 \end{bmatrix}$$

<br>

### 3.4.1.1 Desarrollo Explícito

<div style="text-align: justify;">

Aplicando el algoritmo descrito en 3.3.3:$k=1$: Los multiplicadores son $l_{21}=1/1=1$ y $l_{31}=1/1=1$. Restamos la fila 1 a las filas 2 y 3.$k=2$: El nuevo pivote es $u_{22}=1$. El multiplicador es $l_{32}=1/1=1$. Restamos la fila 2 a la fila 3.Obtenemos las matrices:

</div>

$$L = \begin{bmatrix} 
1 & 0 & 0 \\ 
1 & 1 & 0 \\ 
1 & 1 & 1 
\end{bmatrix}, \quad 
U = \begin{bmatrix} 
1 & 1 & 1 \\ 
0 & 1 & 1 \\ 
0 & 0 & 1 
\end{bmatrix}$$

<br>

### 3.4.1.2 Resultado del Código

<div style="text-align: justify;">

Al ejecutar el script en Python implementado con la metodología descrita, se obtuvieron los siguientes resultados en la consola, confirmando la exactitud de la factorización y la solución del sistema:

</div>

```
=== MATRIZ INICIAL A ===
A =
  [1, 1, 1]
  [1, 2, 2]
  [1, 2, 3]

... [Proceso de eliminación omitido por brevedad] ...

===== MATRIZ FINAL L =====
L =
  [1.0, 0.0, 0.0]
  [1.0, 1.0, 0.0]
  [1.0, 1.0, 1.0]

===== MATRIZ FINAL U =====
U =
  [1, 1, 1]
  [0.0, 1.0, 1.0]
  [0.0, 0.0, 1.0]

>>> SOLUCIÓN PROBLEMA 1 (x): [3.0, 2.0, 1.0]
```

<br>
<br>

## 3.4.2 Problema 2: El Análisis de un Circuito Resistivo

<div style="text-align: justify;">

Se plantea un circuito con 3 mallas donde $A$ representa las resistencias (Matriz de Impedancia) y $b$ las fuentes de voltaje. El sistema está dado por las Leyes de Kirchhoff:

</div>

$$\begin{bmatrix} 
10 & -5 & 0 \\ 
-5 & 15 & -5 \\ 
0 & -5 & 10 
\end{bmatrix} 
\begin{bmatrix} i_1 \\ i_2 \\ i_3 \end{bmatrix} = 
\begin{bmatrix} 10 \\ 0 \\ 0 \end{bmatrix}$$

### 3.4.2.1 Desarrollo Explícito

<div style="text-align: justify;">

Para resolver el sistema de circuitos presentado, realizamos la descomposición manual de la matriz $A$.

<br>

**Paso 1 (Columna 1):**
El pivote es $a_{11} = 10$. Buscamos eliminar el $-5$ de la fila 2 ($a_{21}$).
El multiplicador es $l_{21} = \frac{-5}{10} = -0.5$.
Operación elemental: $F_2 \leftarrow F_2 - (-0.5)F_1$.

</div>

$$
\begin{bmatrix}
10 & -5 & 0 \\
-5 & 15 & -5 \\
0 & -5 & 10
\end{bmatrix}
\xrightarrow{F_2 - (-0.5)F_1}
\begin{bmatrix}
10 & -5 & 0 \\
0 & 12.5 & -5 \\
0 & -5 & 10
\end{bmatrix}
$$

<br>

<div style="text-align: justify;">

- La fila 3 ya tiene un cero en la primera columna, por lo que $l_{31} = 0$.

<br>

**Paso 2 (Columna 2):**
El nuevo pivote es $a_{22} = 12.5$. Buscamos eliminar el $-5$ de la fila 3 ($a_{32}$).
El multiplicador es $l_{32} = \frac{-5}{12.5} = -0.4$.
Operación elemental: $F_3 \leftarrow F_3 - (-0.4)F_2$.

</div>

$$
\begin{bmatrix}
10 & -5 & 0 \\
0 & 12.5 & -5 \\
0 & -5 & 10
\end{bmatrix}
\xrightarrow{F_3 - (-0.4)F_2}
\begin{bmatrix}
10 & -5 & 0 \\
0 & 12.5 & -5 \\
0 & 0 & 8
\end{bmatrix} = U
$$

<br>

<div style="text-align: justify;">

Finalmente, construimos la matriz $L$ colocando los multiplicadores calculados ($-0.5$ y $-0.4$) en sus posiciones correspondientes y unos en la diagonal principal.
El resultado de la descomposición es:

</div>

$$
L = \begin{bmatrix}
1 & 0 & 0 \\
-0.5 & 1 & 0 \\
0 & -0.4 & 1
\end{bmatrix}, \quad
U = \begin{bmatrix}
10 & -5 & 0 \\
0 & 12.5 & -5 \\
0 & 0 & 8
\end{bmatrix}
$$

<br>

### 3.4.2.2 Resultado del Código

<div style="text-align: justify;">

Para el problema de circuitos, la ejecución del algoritmo arrojó los siguientes resultados. Se puede observar que las matrices $L$ y $U$ coinciden con el desarrollo manual, y el vector solución corresponde a las corrientes de cada malla en Amperios.

</div>

```
=== MATRIZ INICIAL A (CIRCUITOS) ===
A =
  [10, -5, 0]
  [-5, 15, -5]
  [0, -5, 10]

... [Cálculo de factores y actualización de filas] ...

===== MATRIZ FINAL L =====
L =
  [1.0, 0.0, 0.0]
  [-0.5, 1.0, 0.0]
  [0.0, -0.4, 1.0]

===== MATRIZ FINAL U =====
U =
  [10, -5, 0]
  [0.0, 12.5, -5.0]
  [0.0, 0.0, 8.0]

>>> SOLUCIÓN PROBLEMA 2 (Corrientes): [1.25, 0.5, 0.25]
```

<br>
<br>

# 3.5. Resultados

<div style="text-align: justify;">

El algoritmo implementado logró descomponer y resolver exitosamente las matrices de prueba.

**Precisión:**

Al reconstruir la matriz original mediante la operación $A_{calc} = L \cdot U$, el error absoluto medio fue del orden de $10^{-16}$ (cero de máquina), lo que valida la implementación.

**Eficiencia:** 

Para el Problema 2, una vez obtenidas $L$ y $U$, cambiar el vector de voltaje $b$ permitió encontrar nuevas corrientes realizando solo $O(n^2)$ operaciones (sustitución hacia adelante y atrás) en lugar de $O(n^3)$ (eliminación completa), confirmando la ventaja teórica del método.

</div>

<br>
<br>

# 3.6. Discusión y Conclusiones

<div style="text-align: justify;">

La descomposición LU ha demostrado ser una herramienta robusta para la resolución de sistemas lineales cuadrados.Ventajas: Su principal fortaleza radica en la separación del procesamiento de la matriz $A$ del vector $b$. Esto es crucial en simulaciones temporales o procesos iterativos de ingeniería donde la estructura del sistema no cambia.Limitaciones: Como se observó en la teoría, el algoritmo básico falla si algún pivote $a_{kk}$ es cero (división por cero). Esto hace necesario implementar estrategias de pivoteo parcial (intercambio de filas) para garantizar la estabilidad numérica en casos generales, lo cual añade una matriz de permutación $P$ al resultado ($PA=LU$).
<br>

> En conclusión, dominar la factorización LU es indispensable para cualquier aplicación de cómputo científico, siendo el paso lógico previo al estudio de métodos más avanzados como Cholesky o QR.

</div>

<br>
<br>

# 4. Método Gauss-Seidel Acelerado (SOR)

## 4.1 Introducción del Método

<div style="text-align: justify;">

Una vez tenemos el concepto del funcionamiento de Gauss-Seidel, tenemos pues un método iterativo que busca una solución *aproximada* que reduce el costo potencial de soluciones exactas. Ahora, **SOR** agregará una variable $\omega$ para controlar y eficientar el paso de aproximación de las iteraciones, a esto se le conoce como un **parámetro de relajación**. Para que este parámetro funcione se debe cumplir una condición necesaria y que el parámetro de $\omega$ se encuentre en el intervalo abierto de 0 < $\omega$ < 2. Si eliges un $\omega$ fuera de este rango el error crecerá en lugar de disminuir.
Existen algunos casos dependiendo el valor del rango que tomes como:

- **Subrelajación (0 < $\omega$ < 1):** Se usa para frenar el paso. Es útil para hacer converger sistemas que con Gauss-Seidel no convergerían.

- **Gauss-Seidel ($\omega$ = 1):** Es el punto neutro.

- **Sobrerrelajación (1 < $\omega$ < 2):** Se usa para acelerar la llegada a la solución de sistemas que ya son estables.

</div>

### 4.1.1 Planteamiento del Problema

El problema se plantea de la siguiente forma:

Ax = b

Debes encontrar el vector que te de la solución de x para un sistema de ecuaciones lineales algebraicas.

**Donde**

- A es una matriz de coeficiente de tamaño n * n 

- x es el vector que se debe encontrar este es de tamaño n * 1 para que coincida con b

- b es el vector de términos independientes de tamaño n x 1

**Las condiciones son:**

- Tienes un punto de partida x^0
- Un acelerador $\omega$ que es el parámetro de relajación

<br>
<br>

### 4.1.2 Contexto Histórico

**Antecedentes:**

Carl Friedrich Gauss, en 1823, introdujo las bases de las técnicas iterativas, estas técnicas se las transmitió a uno de sus alumnos, Christian Gerling. Tiempo después, Philipp L. von Seidel formalizó el método en 1874 analizando sistemas de mínimos cuadrados. El algoritmo resultante lleva el nombre conjunto de **Gauss-Seidel**, que converge para matrices estrictamente diagonalmente dominantes o simétricas definidas positivas, aunque su velocidad de convergencia suele ser lenta para sistemas grandes.

**Desarrollo del SOR:**

El método SOR fue propuesto en 1950 por David M. Young Jr. y Stanley P. Frankel con el propósito de resolver sistemas lineales en ordenadores digitales.
Young estableció la teoría matemática rigurosa, demostrando la relación entre el radio espectral de la matriz y el parámetro óptimo de relajación ($\omega$).
Frankel se enfocó en la aplicación práctica para computadoras digitales tempranas.
Anteriormente a estos métodos ya existían algunos como el de Lewis Fry Richardson, y los métodos desarrollados por R. V. Southwell. Sin embargo estos no se podían aplicar a computadoras digitales o, en el caso de Southwell, eran ineficientes porque requerían intervención visual humana.
<br>
<br>

### 4.1.3 Metodología de Solución

**Antes de empezar:**

Verificar que la matriz A no deba de tener ceros en la diagonal principal porque se deben de dividir estos valores.
Para el parámetro de relajación $\omega$ debe de elegirse dentro del rango 0 < $\omega$ < 2.

**Fórmula iterativa:**

Para cada incógnita $x_i$ en la iteración k + 1 aplicamos la siguiente fórmula:

$$x_i^{(k+1)} = (1 - \omega) x_i^{(k)} + \frac{\omega}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k)} \right)$$

La interpretación de la fórmula es un promedio ponderado:

El primer término $(1 - \omega) x_i^{(k)}$

el segundo término que es la fracción es la corrección de Gauss-Seidel

**Proceso:**

Inicio: Se define un vector $\mathbf{x}^{(0)}$ por defecto suele ser un vector de ceros y se establece un número máximo de iteraciones (N)

Ciclo iterativo: Por cada iteración k = 1 hasta N:

Se recorre cada fila i de la matriz (de 1 a n)

Se calcula el nuevo valor $x_i$ utilizando los valores de la iteración (para las columnas j < i) y los valores de la iteración anterior (para las columnas j > i).

**Calcular el error:**

Al final de cada ciclo completo se calcula el error relativo o la norma de la diferencia entre el vector nuevo y el anterior.

$$||\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}|| < \text{Tolerancia}$$

**Hay dos criterios de parada:**

Convergencia: Si el error es menor que la tolerancia establecida, el proceso se detiene y se entrega el último vector calculado como la solución aproximada.

Divergencia: Si se alcanza el número máximo de iteraciones sin cumplir la tolerancia, el método se detiene indicando que no convergió.

<br>
<br>

## 4.2 Revisión de la Literatura

<br>
<br>

## 4.3 Metodología

<div style="text-align: justify;">

Este método numérico requiere primero, para su justificación, las siguientes herramientas y conceptos matemáticos.

</div>

### 4.3.1 Definiciones Matemáticas 

**Matriz Diagonalmente Dominante:**

<div style="text-align: justify;">

Sea $A\in M_{n\times n}(F)$ matriz cuadrada, se dice que es **diagonalmente dominante por filas** si se cumple que: 

</div>

$$
|a_{ii}| \geq \sum_{j\neq i}|a_{ij}| \ \ \ (\forall i)
$$

<div style="text-align: justify;">

En el caso en que se use ($>$), a esto se le llaman **estrictamente diagonalmente dominante por filas**. Así, se puede generar la **dominancia por columnas** recorriendo la columna asociada al punto en la diagonal, y en caso de cumplir ambas se toma como el caso general de **dominancia diagonal**. 

</div>

**Teorema de Convergencia Gauss-Seidel:**

<div style="text-align: justify;">

Sea el sistema dado por $Ax = b$, con $A\in M_{n\times m}(F)$ y $x, b$ vectores de tamaño $n\times 1$. Si $A$ es estríctamente diagonalmente dominante, entonces el método de **Gauss-Seidel** converge para cualquier aproximación inicial $x^{(0)}$. 

</div>

### 4.3.2 Justificación Teórica del Método SOR 

<div style="text-align: justify;">

Entonces dado un sistema $Ax = b$, con $A\in M_{n\times n}(F)$ matriz y $x, b$ vectores de tamaño $n\times 1$. Generamos así la descomposición estándar de la matriz $A$:

</div>

$$
A = D - L - U
$$

<div style="text-align: justify;">

Donde $D, L, U \in M_{n\times n}(F)$, que cumplen:

- $D =$ matriz diagonal con los elementos $a_{ii}$. 
- $L =$ matriz triangular inferior estricta. 
- $U =$ matriz triangular superior estricta. 

Sabemos pues que por método Gauss-Seidel, la entrada $i$ de $x_{GS}^{(k+1)}$ tendrá la forma: 

</div>

$$
x_{GS, i}^{(k+1)} = \frac{1}{aii}\Big(b_i - \sum_{j=1}^{i-1}a_{ij}x_j ^{(k+1)}-\sum_{j=i+1}^{n}a_{ij}x_j ^{(k)}\Big)
$$

<div style="text-align: justify;">

Donde entonces, el cambio se ve como $\Delta x_i = x_{GS,i}^{(k+1)}-x_i^{(k)}$. Ahora, la metodoloía de usar un **parámetro de relajación** $\omega$, requiere que hagamos un cambio ponderado:

</div>

$$
x_i^{(k+1)} = x_i^{(k)} + \omega\Delta x_i = x_i^{(k)} + \omega \big( x_{GS,i}^{(k+1)} -  x_i^{(k)}\big)
$$

<div style="text-align: justify;">

- Luego sustituímos el valor de $x_{GS,i}^{(k+1)}$ con la formula del método Gauss-Seidel:

</div>

$$\begin{align}
x_i^{(k+1)} &= x_i^{(k)} + \omega \left( \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij}x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij}x_j^{(k)} \right) - x_i^{(k)} \right) \\
&=  (1 - \omega)x_i^{(k)} + \frac{\omega}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij}x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij}x_j^{(k)} \right)
\end{align}
$$

<div style="text-align: justify;">

- Multiplicamos todo por $a_{ii}$ para quitar los denominadores:

</div>

$$
a_{ii}x_i^{(k+1)} = (1 - \omega)a_{ii}x_i^{(k)} + \omega \left( b_i - \sum_{j=1}^{i-1} a_{ij}x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij}x_j^{(k)} \right)
$$

<div style="text-align: justify;">

- Y podemos reorganizar tal que:

</div>

$$
a_{ii}x_i^{(k+1)} + \omega \sum_{j=1}^{i-1} a_{ij}x_j^{(k+1)} = (1 - \omega)a_{ii}x_i^{(k)} - \omega \sum_{j=i+1}^{n} a_{ij}x_j^{(k)} + \omega b_i
$$

<div style="text-align: justify;">

- Con lo que podemos obtener de forma matricial: 

</div>

$$
(D - \omega L)x^{(k+1)} = [(1 - \omega)D + \omega U]x^{(k)} + \omega b
$$

<div style="text-align: justify;">

- Y despejando finalmente:

</div>

$$
x^{(k+1)} = (D - \omega L)^{-1}[(1 - \omega)D + \omega U] x^{(k)} + \omega(D - \omega L)^{-1}\mathbf{b}
$$

<div style="text-align: justify;">

- En la literatura, generalmente definen:

</div>

$$
T_{SOR} = (D - \omega L)^{-1}[(1 - \omega)D + \omega U],\ \ \ \mathbf{c}_{SOR} = \omega(D - \omega L)^{-1}\mathbf{b}
$$

<div style="text-align: justify;">

- Con lo que el despeje final queda de la manera: 

</div>

$$
x^{(k+1)} = T_{SOR}x^{(k)} + c_{SOR}
$$

<div style="text-align: justify;">

> Es importante notar que por teorema, para asegurar la **convergencia** del método **SOR**, $0 < \omega < 2$, esto se puede ver en el Teorema 6.4 de *Applied Numerical Linear Algebra* de J. Demmel.

</div>

<br>
<br>

### 4.3.3 Algoritmo del SOR

<div style="text-align: justify;">

Como vamos a hacer una mejora del método Gauss-Seidel, tomaremos como parámetro de relajación a $\omega$, con $1 < \omega < 2$, es decir, una sobrerrelajación. 
<br>
<br>
El pseudocódigo a continuación supone que la entrada corresponde a una matriz diagonalmente dominante, así como un conjunto de $n$-términos independientes $b$. Observe además que le damos criterios `max_iter` y `tol` para evitar una sobrecarga de operaciones, es decir, damos parámetros de cuántos pasos máximos y cuál  es la tolerancia del error de la solución $x$ encontrada.

</div>

<br>

**4.3.3.1 Pseudocódigo**
```
entrada 
   A // matriz 
   b // vector de terminos independientes
   x // vector de aproximación inicial
   omega // factor de relajación
   tol // tolerancia de alto
   max_iter // máximo de iteraciones

n = filas de A
iter = 0
error = tol + 1

mientras error > tol e iter < max_iter
   error = 0
   
   para i = 1 hasta n
      sum = 0
      xprev = x_i

      para j = 1 hasta n
         si j != i
            sum += A_ij * x_j
         fin si 
      fin para 

      
      x_i = (1-omega)*x_i + omega(b_i - sum)/a_ii

      error = MAX(error, ABS(x_i - xprev))
   fin para 

   iter += 1
fin mientras

retorno x, iter, error
```

<div style="text-align: justify;">

Este algoritmo dadas $n$ y $m$ = `max_iter`, tiene complejidad $O(mn^2)$. Con lo que así como el algoritmo pasado, suponiendo $10^9$ operaciones por segundo, limitamos $m, n$ para estar de acorde con estas. Además, es importante tener en cuenta que `tol`, estará íntimamente relacionada con el tipo de dato que almacene los valores decimales de $x$.

</div>

<br>
<br>

## 4.4. Aplicación del SOR

<br>
<br>

### 4.4.1 Problema 1:

<br>
<br>

### 4.4.1.1 Desarrollo Explícito

<br>
<br>

### 4.4.1.2 Resultado del Código

<br>
<br>

### 4.4.2 Problema 2:

<br>
<br>

### 4.4.2.1 Desarrollo Explícito

<br>
<br>

### 4.4.2.2 Resultado del Código

<br>
<br>

## 4.5 Resultados

<br>
<br>

## 4.6 Discusión y Conclusiones

<br>
<br>

# 5. Referencias
- [Geeks4Geeks - Linear Algebra in Computer Science](https://www.geeksforgeeks.org/maths/linear-algebra-in-computer-science/)
- [CS 357 Textbook - LU Decomposition for Solving Linear Equations](https://cs357.cs.illinois.edu/textbook/notes/linsys.html)
- [MIT OpenCourseWare - LU Decomposition](https://youtu.be/-eA2D_rIcNA?si=8e2Xed8Uyxm4qx1m)
- [Lloyd N., David Bau - Numerical Linear Algebra](https://davidtabora.wordpress.com/wp-content/uploads/2015/01/lloyd_n-_trefethen_david_bau_iii_numerical_line.pdf): Ver **Lecture 20** (Pag. 147 - 151). 
- [Wikipedia - Diagonally Dominant Matrix](https://en.wikipedia.org/wiki/Diagonally_dominant_matrix)
- [Wikipedia - Gauss-Seidel Method](https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method)
- [J. Demmel - Applied Numerical Linear Algebra](https://www.stat.uchicago.edu/~lekheng/courses/302/demmel/): Ver **Capítulo 6** (Pag 282-286, 290)
- Seidel, L. (1874). Ueber ein Verfahren, die Gleichungen, auf welche die Methode der kleinsten Quadrate führt, sowie lineäre Gleichungen überhaupt, durch successive Annäherung aufzulösen. Münchener DigitalisierungsZentrum, Digitale Bibliothek, Bayerische Staatsbibliothek. Recuperado de https://www.digitale-sammlungen.de/de/view/bsb11179192?page=5
- Burden, R. L., & Faires, J. D. (2002). Análisis Numérico (7ª ed.). Ediciones Paraninfo, S.A. Recuperado de https://evflores.wordpress.com/wp-content/uploads/2014/02/analisis-numerico-richard-l-burden-7ma.pdf
- Sobrerrelajación sucesiva. (n.d.). En Wikipedia, la enciclopedia libre. Recuperado el 2 de diciembre de 2025, de https://es.wikipedia.org/wiki/Sobrerrelajaci%C3%B3n_sucesiva 
