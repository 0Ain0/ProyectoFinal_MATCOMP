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
- Cecilia 
- Casta
- Santiago Medina Domínguez  
- Santiago Espinosa Ollivier
- Ain Bolaños Cortés

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

El problema de resolver sistemas de ecuaciones lineales reside, en el núcleo del análisis numérico y la computación científica, y se expresa en forma matricial Ax = b. Este tipo de sistemas aparece en casi todos los ámbitos de la ingeniería y las ciencias aplicadas, desde el análisis de circuitos eléctricos, el diseño estructural, la modelización financiera, el aprendizaje automático, etc. Si bien es cierto que existen métodos iterativos (por ejemplo, el método de Gauss-Seidel que trabajará tu compañero) que encuentran una solución aproximada, los métodos directos se limitan a encontrar una solución exacta (dentro de la precisión de la máquina) en un número finito de pasos.
El método de Descomposición LU es uno de los métodos directos más fundamentales y eficientes. No es simplemente un algoritmo, sino una factorización de la matriz A que transforma el problema Ax = b en un proceso de dos pasos mucho más simple de resolver. Esta investigación se centrará en la fundamentación teórica, la formulación matemática, el análisis de estabilidad y las ventajas de este crucial método.

## 3.1.2 Contexto Histórico

La factorización LU es la formalización moderna de la eliminación gaussiana, un algoritmo con raíces que se remontan a más de dos milenios (como en los "Nueve Capítulos sobre el Arte Matemático" chinos). Sin embargo, su concepción actual como una factorización matricial explícita fue fundamental para el nacimiento de la computación digital. Alan Turing, en su trabajo seminal de 1948 "Rounding-Off Errors in Matrix Processes", utilizó esta factorización como herramienta central para analizar los errores de redondeo en procesos matriciales automatizados. Desde entonces, ha evolucionado para convertirse en el estándar de facto en bibliotecas de software de alto rendimiento como LAPACK.


## 3.1.3 Metodología de Solución

La estrategia general consiste en descomponer la matriz $A$ en el producto de dos matrices triangulares: L (Lower, triangular inferior) y U (Upper, triangular superior), tal que A=LU.
Esta factorización permite desacoplar el sistema $Ax=b$ en dos etapas sucesivas: primero se resuelve un sistema para un vector intermedio utilizando L (sustitución hacia adelante) y posteriormente se obtiene la solución final utilizando $U$ (sustitución hacia atrás). Para garantizar la viabilidad en computadoras, se incorpora una matriz de permutación P (pivoteo), resultando en la forma PA=LU.

<br>
<br>

# 3.2 Revisión de la Literatura

La literatura canónica en álgebra lineal numérica posiciona a la descomposición LU como la herramienta fundamental para sistemas densos. Autores como Golub y Van Loan (2013) y Trefethen y Bau (1997) describen este método como la "navaja suiza" del cálculo matricial.La investigación teórica destaca dos aspectos críticos:Eficiencia: Una vez calculada la factorización (con costo cúbico $O(n^3)$), resolver para múltiples vectores $b$ es cuadrático $O(n^2)$, lo cual es una ventaja decisiva sobre la inversión de matrices explícita.Estabilidad: Los textos clásicos (Burden et al., 2017) enfatizan que la factorización pura es inestable sin intercambio de filas. La literatura moderna se centra en variantes con pivoteo parcial para acotar el crecimiento de los errores de redondeo, asegurando que los multiplicadores en la eliminación no excedan la unidad en magnitud.

La factorización LU constituye la formalización directa del método de eliminación gaussiana, que es un algoritmo muy antiguo con una historia de más de dos mil años. No obstante, su concepción moderna como una factorización de matrices resultó crucial para el desarrollo de las computadoras digitales. Alan Turing llegó a utilizar esta factorización como una de las herramientas centrales en su trabajo de 1948 "Rounding-Off Errors in Matrix Processes", tratando de encontrar el análisis de errores en las computaciones basados en matrices (Turing, 1948).

Los textos considerados canónicos de álgebra lineal numérica, los de Golub & Van Loan (2013) o los de Trefethen & Bau (1997), consideran la descomposición LU como la "navaja suiza" para la resolución de sistemas lineales densos. Es precisamente la base sobre la cual se construyen algoritmos más complejos y es el método de facto que se encuentra implementado en paquetes de software de alto rendimiento como LAPACK y bibliotecas como NumPy (SciPy) o MATLAB (Trefethen & Bau, 1997). La literatura no sólo se esfuerza por describir el algoritmo básico, sino que lo hace por describir, y exclusivamente por el mismo, las versiones estables que tienen pivoteo, indispensable para la aplicación práctica de estas ideas.




# 3.3 Metodología

Este método numérico requiere primero, para su justificación, las siguientes herramientas y conceptos matemáticos. 

## 3.3.1 Definiciones Matemáticas

El método se implementa mediante la factorización con pivoteo parcial para garantizar la estabilidad numérica:$$PA = LU$$Donde $P$ es una matriz de permutación, L es triangular inferior con unos en la diagonal, y U es triangular superior.Al sustituir en Ax=b, obtenemos LUx = Pb. Definimos b' = Pb$ y un vector auxiliar z = Ux. El problema se resuelve en dos pasos:Paso 1: Sustitución hacia adelante ($Lz = b'$).Se despeja$z comenzando desde la primera fila:z_i = b'_i - \sum_{j=1}^{i-1} L_{ij}z_j, \quad \text{para } i=2, \dots, n \quad (\text{con } z_1 = b'_1)Paso 2: Sustitución hacia atrás (Ux = z).Se despeja x comenzando desde la última fila:$$   x_i = \frac{1}{U_{ii}} \left( z_i - \sum_{j=i+1}^{n} U_{ij}x_j \right), \quad \text{para } i=n-1, \dots, 1 \quad (\text{con } x_n = \frac{z_n}{U_{nn}})$$

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
Así pues, el despeje $LU$ fue único para $A$. Como observación final, es facil ver que  $L$ es facil de obtener, al no ser necesario ningún despeje u operación entre matrices. Mientras que $U$ se puede obtener al tiempo que se genera $L$. Con esto queda demostrada la existencia y unicidad de $L, U$, además de dar una pauta para la obtención de dichas matrices.

<br>
<br>

## 3.3.3. Algoritmo de la Descomposición LU
## 3.3.3.1. Pseudocódigo
```

``` 

## 3.3.3.2. Script en Python
```python
print("Buenas buenas")
```

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

La descomposición LU, reforzada con pivoteo parcial, constituye el pilar de los métodos directos. Su elegancia teórica transforma un problema denso en dos sistemas triangulares triviales. Su verdadera fortaleza radica en la eficiencia para resolver sistemas con múltiples cargas ($b$) constantes y en su estabilidad controlada, lo que la convierte en la herramienta predilecta para simulaciones ingenieriles donde la precisión exacta es prioritaria sobre la velocidad de iteración aproximada.

---

<br>

# 4. Método Gauss-Seidel Acelerado (SOR)

# 4.1. Introducción del Método
Una vez tenemos el concepto del funcionamiento de Gauss-Seidel, tenemos pues un método iterativo que busca una solución *aproximada* que reduce el costo potencial de soluciones exactas. Ahora, **SOR** agregará una variable $\omega$ para controlar y eficienciar el paso de aproximamiento de las iteraciones, a esto se le conoce como un **parámetro de relajación**. 

## 4.1.1 Planteamiento del Problema
## 4.1.2 Contexto Histórico
## 4.1.3 Metodología de Solución

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

# 4.3.2. Justificación Teórica del Método SOR 
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
x^{(k+1)} = T_{SOR}\mathbf{x}^{(k)} + c_{SOR}
$$

> Es importante notar que por teorema, para asegurar la **convergencia** del método **SOR**, $0 < \omega < 2$, esto se puede ver en el Teorema 6.4 de *Applied Numerical Linear Algebra* de J. Demmel. 

<br>
<br>

## 4.3.3. Algoritmo de la Descomposición LU
## 4.3.3.1. Pseudocódigo
```

``` 

## 4.3.3.2. Script en Python
```python
print("Buenas buenas")
```

# 4.4. Aplicación de la Descomposición LU
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
- Burden, R. L., Faires, J. D., & Burden, A. M. (2017). Análisis numérico (10ª ed.). Cengage Learning.
- Golub, G. H., & Van Loan, C. F. (2013). Matrix computations (4ª ed.). Johns Hopkins University Press.
- Trefethen, L. N., & Bau, D. (1997). Numerical linear algebra. SIAM.
- Turing, A. M. (1948). Rounding-Off Errors in Matrix Processes. The Quarterly Journal of Mechanics and Applied Mathematics, 1(1), 287–308.