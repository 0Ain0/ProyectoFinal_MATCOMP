<br>
<br>
<br>
<br>
<br>
<br>
<br>
<div align="center">

![Logo Universidad Panamericana](assets/logoUP.png)

<br>

# **Universidad Panamericana**
## **Escuela de Ingeniería**
## Matemáticas de la Computación

<br>
<br>
<br>

# **Descomposición LU para sistemas lineales y Método Gauss-Seidel Acelerado**

<br>
<br>
<br>
<br>
</div>

### **Profesor:** David Iván Morales Huerta

### **Equipo:**
- Ain Bolaños Cortés
- Casta
- Santiago Espinosa Ollivier
- Cecilia Gaona Vidales 
- Santiago Medina Domínguez 


<br>

**Entrega:** 3 de Diciembre del 2025

<div style="page-break-after: always;"></div>

# 1. Resumen
Este reporte desarrolla, explica, justifica e implementa los métodos de descomposición LU y el método Gauss-Seidel para la resolución de sistemas de ecuaciones lineales. Con un enfoque en su resolución automatizada computable. Dando así la justificación teórico mátemática y su implementación algorítmica en Python. 

<br>
<br>


# 2. Introducción
Los sistemas de ecuaciones lineales, así como su representación matricial son un problema típico de la computación, en áreas relevantes como la Visión de Computadora, la criptografía, los gráficos y el machine learning, así como en las simulaciones físicas y el análisis matemático. Y de la misma manera en que los conceptos del álgebra lineal son relevantes para la expansión de los límites de la computación, el poder computacional extiende las posibilidades de procesamiento y resolución de sistemas lineales extensos.   

<br>
<br>

# 3. Descomposición LU

# 3.1. Introducción del Método
## 3.1.1 Planteamiento del Problema
## 3.1.2 Contexto Histórico
## 3.1.3 Metodología de Solución

<br>
<br>

# 3.2 Revisión de la Literatura

<br>
<br>

# 3.3 Metodología
Este método numérico requiere primero, para su justificación, las siguientes herramientas y conceptos matemáticos. 

## 3.3.1 Definiciones Matemáticas
### Matriz Cuadrada 
Sea $A\in M_{n\times m}(F)$, con $M_{n\times m}$ el conjunto de matrices del campo $F$ con dimensiones $n\times m$. En el caso en que $n = m$ decimos que $A$ es una matriz cuadrada. Para el futuro de este reporte, se tomará $F = \mathbb{R}$. 

### Vectores Canónicos
Sea $e_k$ un vector de dimensión $n\times 1$, decimos que este es $k$-ésimo vector canónico de $\mathbb{R}^n$ si cumple que con $x_i$ sus entradas con $i\in\{1, \dots, n\}$, entonces $x_k = 1$ y $x_i = 0$ con $i \neq k$. Así pues se puede ver de la siguiente manera: 
$$e_k = \begin{bmatrix}
   0 \\
   \vdots \\
   1 \\
   \vdots \\
   0
\end{bmatrix}$$
Esto es útil tambien para referirnos a partes de una matriz $A$ lo suficientemente grande. Dado que: 
- $Ae_k = k$-ésima columna de $A$
- $e^T_k A = k$-ésima fila de $A$ 

Además, en dado caso de tener un vector $v$ con dimensiones $n\times 1$, la multiplicación $v e_k^T$ será una matriz cuadrada $n\times n$ con todos los valores $0$ a excepción de la columna $k$ que contendrá los valores de $v$. 

## 3.3.2. Justificación Teórica de la Descomposición LU
Sea $A\in M_{n\times n}(\mathbb{R})$, lo que buscamos es la secuencia finita $\{ L_k\in M_{n\times n}: 1\leq k \leq n\}$ tales que para $x_k$ la $k$-ésima columna de la matríz $A$ ($x_k = Ae_k$): 
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
Esto con el objetivo de obtener a partir de esta cantidad finita de $L_k$ una matriz $U$ triangular superior. De manera similar a la eliminación Gaussiana, obtenemos pues, suponiendo que $x_{kk} \neq 0$ los siguientes valores: $$l_{jk} = \frac{x_{jk}}{x_{kk}}\ \ \ \ \ \ \ \ (k<j\leq m)$$
Esto significa que tenemos que pedirle a $A$ que todos sus elementos en la diagonal **sean no nulos**. Con estas entradas entradas definamos $l_k = [0, \dots, 0, l_{k+1, k}, \dots, l_{nk}]^T$, y así obtenemos que $L_k$ se verá de la siguiente manera: 
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
Con esto, generamos $k$ matrices $L_k$ tales que: 
$$
L_{n-1} \cdots L_2 L_1 A = U
$$
La matriz deseada, y más aún, despejando para $A$ obtenemos: 
$$
A = L_1^{-1}L_2^{-1}\cdots L_{n-1}^{-1} U
$$
Por como está estructurada cada $L_k$, podemos ver que $L_k^{-1} = I + l_k e_k$, con lo que nos podemos tomar: 
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
Con esto obtuvimos pues $A = LU$, con $L$ matriz triangular inferior y $U$ matriz triangular superior. Ahora, supongamos que existen $L_1, L_2$ matrices triangulares inferiores y $U_1, U_2$ matrices triangulares superiores tales que: 
$$
A = L_1U_1 = L_2 U_2\ \ \ \ \Longrightarrow\ \ \ \ L_2^{-1}L_1 = U_2 U_1^{-1}
$$
Por cerradura de las matrices triangulares, $L_2^{-1}L_1$ es triangular inferior y $U_2 U_1^{-1}$ es triangular superior, así  por la igualdad: 
$$
L_2^{-1}L_1 = U_2 U_1^{-1} = I\ \ \ \ \Longrightarrow\ \ \ \ L_1 = L_2\  \land\ U_1 = U_2
$$
Así pues, el despeje $LU$ fue único para $A$. Como observación final, es facil ver que  $L$ es facil de obtener, al no ser necesario ningún despeje u operación entre matrices. Mientras que $U$ se puede obtener al tiempo que se genera $L$. Con esto **queda demostrada la existencia y unicidad de $L, U$**, además de dar una pauta para la obtención de dichas matrices.

## 3.3.3. Algoritmo de la Descomposición LU
Como se puede apreciar en la deducción matemática, las opraciones de inversa para las $L_i$ son sumamente simples de efectuar, así como el despeje de cada entrada de $L$ se puede obtener de manera directa, lo que hace al algoritmo sumamente efectivo dado que es de matrices. Con esto, Lloyd N. y David Bau, agregan que no es necesario almacenar $A, L$ o $U$ en distintas matrices. 

Generamos una matriz partiendo de la identidad, y vamos a encontrar todos los $l_{jk}$, es decir, los factores de la triangular inferior estricta de $L$ para con ellos actualizar las entradas de $U$, haciendo la Eliminación Gaussiana. 

El pseudo código a continuación supone que la matríz $A$ dada es válida, es decir, $\forall k\leq n$, $A_{kk} \neq 0$. Con $n$ claro el tamaño de la matriz $A$, es decir $A \in M_{n\times n}(F)$.  

## 3.3.3.1. Pseudocódigo
```
tome U, L matrices n x n

U = A, L = I

para k = 1 hasta n - 1

   para i = k + 1 hasta n
      l_ik = u_ik/u_kk
      
      para j = k hasta n
         u_ij -= l_ik * u_kj
      fin para
   fin para
fin para

retornar L, U
``` 
Es importante notar que este algoritmo tiene **complejidad $O(n^3)$**, lo que implica que su funcionamiento en menos de $1$ segundo, estará limitado a matrices de a lo más $n = 10^3$, suponiendo una cantidad de operaciones $\approx 10^9$ por segundo. 
<br>
<br>

# 3.4. Aplicación de la Descomposición LU
## 3.4.1 Problema 1:
## 3.4.1.1 Desarrollo Explícito
## 3.4.1.2 Resultado del Código

<br>
<br>

## 3.4.1 Problema 2:
## 3.4.1.1 Desarrollo Explícito
## 3.4.1.2 Resultado del Código

<br>
<br>

# 3.5. Resultados

<br>
<br>

# 3.6. Discusión y Conclusiones
<br>
<br>

---

<br>

# 4. Método Gauss-Seidel Acelerado (SOR)

# 4.1. Introducción del Método
Una vez tenemos el concepto del funcionamiento de Gauss-Seidel, tenemos pues un método iterativo que busca una solución *aproximada* que reduce el costo potencial de soluciones exactas. Ahora, **SOR** agregará una variable $\omega$ para controlar y eficienciar el paso de aproximamiento de las iteraciones, a esto se le conoce como un **parámetro de relajación**. Para que este parámetro funcione se debe de cumplir una condición necesaria y que el parámetro de $\omega$ se encuentre en el intervalo abierto de 0 < $\omega$ < 2. Si eliges un $\omega$ fuera de este rango el error crecerá en lugar de disminuir.
Existen algunos casos dependiendo el valor del rango que tomes como:

- **Subrelajación (0 < $\omega$ < 1):** Se usa para frenar el paso. Es útil para hacer converger sistemas que con Gauss-Seidel no convergerían.

- **Gauss-Seidel ($\omega$ = 1):** Es el punto neutro.

- **Sobrerrelajación (1 < $\omega$ < 2):** Se usa para acelerar la llegada a la solución de sistemas que ya son estables.

## 4.1.1 Planteamiento del Problema
El problema se plantea de la siguiente forma:

Ax = b

Debes de encontrar el vector que te de la solución de x para un sistema de ecuaciones lineales algebraicas.

**Donde**

- A es una matriz de coeficiente de tamaño n * n 

- x es el vector que se debe de encontrar este es de tamaño n * 1 para que coincida con b

- b es el vector de términos independientes de tamaño n x 1

- **Las condiciones son:**
- Tienes un punto de partida x^0
- Un acelerador $\omega$ que es el parámetro de relajación


## 4.1.2 Contexto Histórico

- **Antecedentes**

Carl Friedrich Gauss, en 1823, introdujo las bases de las técnicas iterativas, estas técnicas se las transmitió a uno de sus alumnos, Christian Gerling. Tiempo después, Philipp L. von Seidel formalizó el método en 1874 analizando sistemas de mínimos cuadrados. El algoritmo resultante lleva el nombre conjunto de **Gauss-Seidel**, que converge para matrices estrictamente diagonalmente dominantes o simétricas definidas positivas, aunque su velocidad de convergencia suele ser lenta para sistemas grandes.

- **Desarrollo del SOR**

El método SOR fue propuesto en 1950 por David M. Young Jr. y Stanley P. Frankel con el propósito de resolver sistemas lineales en ordenadores digitales.
Young estableció la teoría matemática rigurosa, demostrando la relación entre el radio espectral de la matriz y el parámetro óptimo de relajación ($\omega$).
Frankel se enfocó en la aplicación práctica para computadoras digitales tempranas.
Anteriormente a estos métodos ya existían algunos como el de Lewis Fry Richardson, y los métodos desarrollados por R. V. Southwell. Sin embargo estos no se podían aplicar a computadoras digitales o, en el caso de Southwell, eran ineficientes porque requerían intervención visual humana.


## 4.1.3 Metodología de Solución   

- **Antes de empezar**
Verificar que la matriz A no deba de tener ceros en la diagonal principal porque se deben de dividir estos valores.
Para el parámetro de relajación $\omega$ debe de elegirse dentro del rango 0 < $\omega$ < 2.

- **Fórmula iterativa**
Para cada incógnita $x_i$ en la iteración k + 1 aplicamos la siguiente fórmula:

$$x_i^{(k+1)} = (1 - \omega) x_i^{(k)} + \frac{\omega}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k)} \right)$$

La interpretación de la fórmula es un promedio ponderado:

El primer término $(1 - \omega) x_i^{(k)}$

el segundo término que es la fracción es la correción de Gauss-Seidel

**Proceso**
Inicio: Se define un vector $\mathbf{x}^{(0)}$ por defecto suele ser un vector de ceros y se establece un número máximo de iteraciones (N)

Ciclo iterativo: Por cada iteración k = 1 hasta N:

Se recorre cada fila i de la matriz (de 1 a n)

Se calcula el nuevo valor $x_i$ utilizando los valores de la iteración (para las columnas j < i) y los valores de la iteración anterior (para las columnas j > i).

**Calcular el error**

Al final de cada ciclo completo se calcula el error relativo o la norma de la diferencia entre el vector nuevo y el anterior.

$$||\mathbf{x}^{(k+1)} - \mathbf{x}^{(k)}|| < \text{Tolerancia}$$

**Hay dos criterios de parada**

Convergencia: Si el error es menor que la tolerancia establecida, el proceso se detiene y se entrega el último vector calculado como la solución aproximada.

Divergencia: Si se alcanza el número máximo de iteraciones sin cumplir la tolerancia, el método se detiene indicando que no convergió.

<br>
<br>

# 4.2 Revisión de la Literatura

# 4.3 Metodología
Este método numérico requiere primero, para su justificación, las siguientes herramientas y conceptos matemáticos.

## 4.3.1 Definiciones Matemáticas 
### Matriz Diagonalmente Dominante
Sea $A\in M_{n\times n}(F)$ matríz cuadrada, se dice que es **diagonalmente dominante por filas** si se cumple que: 
$$
|a_{ii}| \geq \sum_{j\neq i}|a_{ij}| \ \ \ (\forall i)
$$
En el caso en que se use ($>$), a esto se le llaman **estrictamente diagonalmente dominante por filas**. Así, se puede generar la **dominancia por columnas** recorriendo la columna asociada al punto en la diagonal, y en caso de cumplir ambas se toma como el caso general de **dominancia diagonal**. 

### Teorema de Convergencia Gauss-Seidel
Sea el sistema dado por $Ax = b$, con $A\in M_{n\times m}(F)$ y $x, b$ vectores de tamaño $n\times 1$. Si $A$ es estríctamente diagonalmente dominante, entonces el método de **Gauss-Seidel** converge para cualquier aproximación inicial $x^{(0)}$. 

## 4.3.2. Justificación Teórica del Método SOR 
Entonces dado un sistema $Ax = b$, con $A\in M_{n\times n}(F)$ matriz y $x, b$ vectores de tamaño $n\times 1$. Generamos así la descomposición estándar de la matriz $A$: 
$$
A = D - L - U
$$
Donde $D, L, U \in M_{n\times n}(F)$, que cumplen: 
- $D =$ matriz diagonal con los elementos $a_{ii}$. 
- $L =$ matriz triangular inferior estricta. 
- $U =$ matriz triangular superior estricta. 

Sabemos pues que por método Gauss-Seidel, la entrada $i$ de $x_{GS}^{(k+1)}$ tendrá la forma: 
$$
x_{GS, i}^{(k+1)} = \frac{1}{aii}\Big(b_i - \sum_{j=1}^{i-1}a_{ij}x_j ^{(k+1)}-\sum_{j=i+1}^{n}a_{ij}x_j ^{(k)}\Big)
$$
Donde entonces, el cambio se ve como $\Delta x_i = x_{GS,i}^{(k+1)}-x_i^{(k)}$. Ahora, la metodoloía de usar un **parámetro de relajación** $\omega$, requiere que hagamos un cambio ponderado: 
$$
x_i^{(k+1)} = x_i^{(k)} + \omega\Delta x_i = x_i^{(k)} + \omega \big( x_{GS,i}^{(k+1)} -  x_i^{(k)}\big)
$$
Luego sustituímos el valor de $x_{GS,i}^{(k+1)}$ con la formula del método Gauss-Seidel: 
$$\begin{align}
x_i^{(k+1)} &= x_i^{(k)} + \omega \left( \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij}x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij}x_j^{(k)} \right) - x_i^{(k)} \right) \\
&=  (1 - \omega)x_i^{(k)} + \frac{\omega}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij}x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij}x_j^{(k)} \right)
\end{align}
$$
Multiplicamos todo por $a_{ii}$ para quitar los denominadores: 
$$
a_{ii}x_i^{(k+1)} = (1 - \omega)a_{ii}x_i^{(k)} + \omega \left( b_i - \sum_{j=1}^{i-1} a_{ij}x_j^{(k+1)} - \sum_{j=i+1}^{n} a_{ij}x_j^{(k)} \right)
$$
Y podemos reorganizar tal que: 
$$
a_{ii}x_i^{(k+1)} + \omega \sum_{j=1}^{i-1} a_{ij}x_j^{(k+1)} = (1 - \omega)a_{ii}x_i^{(k)} - \omega \sum_{j=i+1}^{n} a_{ij}x_j^{(k)} + \omega b_i
$$
Con lo que podemos obtener de forma matricial: 
$$
(D - \omega L)x^{(k+1)} = [(1 - \omega)D + \omega U]x^{(k)} + \omega b
$$
Y despejando finalmente: 
$$
x^{(k+1)} = (D - \omega L)^{-1}[(1 - \omega)D + \omega U] x^{(k)} + \omega(D - \omega L)^{-1}\mathbf{b}
$$
En la literatura, generalmente definen:
$$
T_{SOR} = (D - \omega L)^{-1}[(1 - \omega)D + \omega U],\ \ \ \mathbf{c}_{SOR} = \omega(D - \omega L)^{-1}\mathbf{b}
$$
Con lo que el despeje final queda de la manera: 
$$
x^{(k+1)} = T_{SOR}x^{(k)} + c_{SOR}
$$

> Es importante notar que por teorema, para asegurar la **convergencia** del método **SOR**, $0 < \omega < 2$, esto se puede ver en el Teorema 6.4 de *Applied Numerical Linear Algebra* de J. Demmel. 

<br>
<br>

## 4.3.3. Algoritmo del SOR
Como vamos a hacer una mejora del método Gauss-Seidel, tomaremos como parámetro de relajación a $\omega$, con $1 < \omega < 2$, es decir, una sobrerrelajación. 

El pseudocódigo a continuación supone que la entrada corresponde a una matríz diagonalmente dominante, así como un conjunto de $n$-términos independientes $b$. Observe además que le damos criterios `max_iter` y `tol` para evitar una sobrecarga de operaciones, es decir, damos parámetros de cuántos pasos máximos y cuál  es la tolerancia del error de la solución $x$ encontrada.

## 4.3.3.1. Pseudocódigo
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
Este algoritmo dadas $n$ y $m$ = `max_iter`, tiene complejidad $O(mn^2)$. Con lo que así como el algoritmo pasado, suponiendo $10^9$ operaciones por segundo, limitamos $m, n$ para estar de acorde con estas. Además, es importante tener en cuenta que `tol`, estará íntimamente relacionada con el tipo de dato que almacene los valores decimales de $x$. 

<br>
<br>

# 4.4. Aplicación del SOR


## 4.4.1. Problema 1:
## 4.4.1.1. Desarrollo Explícito
## 4.4.1.2. Resultado del Código

<br>
<br>

## 4.4.1. Problema 2:
## 4.4.1.1. Desarrollo Explícito
## 4.4.1.2. Resultado del Código

<br>
<br>

# 4.5. Resultados

<br>
<br>

# 4.6. Discusión y Conclusiones

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
