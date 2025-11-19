<br>
<br>
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

### **Profesor:** [Nombre del Profesor]

### **Equipo:**
- Cecilia Casta
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


# 2. Objetivos

## 2.1. Objetivo General
[Descripción del objetivo general]

## 2.2. Objetivos Específicos
- [Objetivo específico 1]
- [Objetivo específico 2]
- [Objetivo específico 3]

<br>
<br>

# 3. Introducción
Los sistemas de ecuaciones lineales, así como su representación matricial son un problema típico de la computación, en áreas relevantes como la Visión de Computadora, la criptografía, los gráficos y el machine learning, así como en las simulaciones físicas y el análisis matemático. Y de la misma manera en que los conceptos del álgebra lineal son relevantes para la expansión de los límites de la computación, el poder computacional extiende las posibilidades de procesamiento y resolución de sistemas lineales extensos.   

<br>
<br>

# 4. Descomposición LU

# 4.1. Planteamiento del Problema: 
[Descripción del problema para Descomposición LU]

<br>
<br>

# 4.2. Marco Teórico: 
[Fundamentos teóricos de la Descomposición LU]

## 4.2.1 Definiciones Matemáticas

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


<br>
<br>

# 4.3. Desarrollo Teórico: 
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
Así pues, el despeje $LU$ fue único para $A$. Como observación final, es facil ver que  $L$ es facil de obtener, al no ser necesario ningún despeje u operación entre matrices. Mientras que $U$ se puede obtener al tiempo que se genera $L$. 

<br>
<br>

# 4.4. Método Práctico: 
[Aplicación práctica del método]

<br>
<br>

# 4.5. Algoritmo: 
## 4.5.1. Pseudocódigo
```

``` 

## 4.5.2. Script en Python
```python
print("Buenas buenas")
```

## 4.5.3. Resultados

<br>
<br>

---

<br>

# 5. Método Gauss-Seidel Acelerado (SOR)
# 5.1. Planteamiento del Problema: 
[Descripción del problema para Descomposición LU]

<br>
<br>

# 5.2. Marco Teórico: 
[Fundamentos teóricos de la Descomposición LU]

<br>
<br>

# 5.3. Desarrollo Teórico: 
[Desarrollo matemático y teórico]

<br>
<br>

# 5.4. Método Práctico: 
[Aplicación práctica del método]

<br>
<br>

# 5.5. Algoritmo: 
## 5.5.1. Pseudocódigo
```

``` 

## 5.5.2. Script en Python
```python
print("Buenas buenas")
```

## 5.5.3. Resultados

# 6. Concluciones 

# 7. Referencias
- [Geeks4Geeks - Linear Algebra in Computer Science](https://www.geeksforgeeks.org/maths/linear-algebra-in-computer-science/)
- [CS 357 Textbook - LU Decomposition for Solving Linear Equations](https://cs357.cs.illinois.edu/textbook/notes/linsys.html)
- [MIT OpenCourseWare - LU Decomposition](https://youtu.be/-eA2D_rIcNA?si=8e2Xed8Uyxm4qx1m)
- [Lloyd N., David Bau - Numerical Linear Algebra](https://davidtabora.wordpress.com/wp-content/uploads/2015/01/lloyd_n-_trefethen_david_bau_iii_numerical_line.pdf): Ver **Lecture 20** (Pag. 147 - 151). 