# QR-ColPiv-LeastSquares
C/C++ code to solve the least-squares problem via column-pivoted QR decomposition

1. [Mathematical Overview](#mathematical-overview)
2. [Core Algorithm](#core-algorithm)
3. [Key Equations](#key-equations)
4. [Numerical Stability](#numerical-stability)
5. [Validation Metrics](#validation-metrics)
6. [Comparison with SVD](#comparison-with-svd)


## Mathematical Overview  

### Linear Least-Squares Problem  

&zwnj;**Solves the linear least-squares problem:**&zwnj;

$$\min_{{x}} {||}{A}{x} - {b}{||}_2^2$$


&zwnj;**Where:**&zwnj;
- ${A} \in \mathbb{R}^{m \times n}$ (Full column rank)  
- ${b} \in \mathbb{R}^{m}$  
- ${x} \in \mathbb{R}^{n}$  

### Decomposition  
$${A}{\Pi} = {Q}{R}$$  

| Symbol                | Description                                   |
|-----------------------|-----------------------------------------------|
| ${\Pi}$            | Column permutation matrix                     |
| ${Q}$              | Orthogonal matrix (${Q}^\top{Q} = {I}$)       |  
| ${R}$              | Upper triangular matrix                       |  


## Core Algorithm

### Step 1: Column Pivoting  
&zwnj;**At iteration  k :**&zwnj;
1. Compute remaining column norms:  

$${\gamma}_j = \|{R}(k:m, j)\|_2^2 \quad \text{for } j = k:n$$  

2. Swap column $p$ with $k$-th column:  

$$p = \arg\max_j \gamma_j$$

### Step 2: Householder Transformation  
&zwnj;**For column $k$:**&zwnj;
1. Compute Householder vector ${v}$:  

$${v} = {x} \pm {||}{x}\{||}_2 {e}_k, \quad {x} = {R}(k:m, k)$$
     
2. Apply reflection to submatrix:  

$${R}(k:m, k:n) \leftarrow {R}(k:m, k:n) - {v} \left( \frac{2}{{v}^\top{v}} {v}^\top {R}(k:m, k:n) \right)$$
   
3. Update orthogonal matrix:  

$${Q} \leftarrow {Q} - {Q}{v} \left( \frac{2}{{v}^\top{v}} {v}^\top \right)$$  


## Key Equations

### Triangular System Solving & Validation  

#### After Decomposition  
&zwnj;**Triangular System Solving:**&zwnj;

$${R}_{11} {z} = {Q}_1^\top {b}$$  

- ${R}_{11}$: Leading n $\times$ n submatrix of ${R}$  

### Solution Recovery  

$${x} = {\Pi} {z}$$  

### Residual Orthogonality (Validation)  

$$({A}{x})^\top ({b} - {A}{x}) \approx 0 \quad \text{(Machine precision $\sim 10^{-16}$)}$$  


## Numerical Stability

### Key Mechanisms  

1. Column Pivoting  

$$\|r_{kk}\| \geq \|r_{ij}\| \quad \text{for } i \geq k, \, j \geq k$$  

2. Householder Reflections  
&zwnj;**Backward stable**&zwnj; with error bound:  

$${{||}{\Delta}A{||}_F \over {||}A{||}_F} {\leq} O (mn {\epsilon}\text{(machine)})$$  

3. Explicit Orthogonality  

$$\|{Q}^\top{Q} - {I}\|_F \leq 10^{-14} \quad \text{(Typical in 2025 hardware)}$$  


## Validation Metrics

| Metrics           | Formula                                                                  | Acceptable Threshold|
|:------------------|:------------------------------------------------------------------------:|------------------: |
| QR Accuracy       | $\frac{\|{A}{\Pi} - {Q}{R}\|_F}{\|{A}\|_F} < 10^{-12}$                  | $<10^{-12}$    |
| Orthogonality     | $\|{Q}^\top{Q} - {I}\|_F < 10^{-14}$                            | $<10^{-14}$    |
| Solution Error    | $\|{A}{x} - {b}\|_2$                                            | Problem-dependent  |


## Comparison with SVD

| &zwnj;**Property**&zwnj;            | &zwnj;**QR (Column-Pivoted)**&zwnj;              | &zwnj;**SVD**&zwnj;                    |
|-------------------------|---------------------------------------|-----------------------------|
| &zwnj;**Flops**&zwnj;               | $2mn^2 - \frac{2}{3}n^3$          | $4m^2n + 22n^3$         |
| &zwnj;**Rank Detection**&zwnj;      | Threshold-based                       | Exact                       |
| &zwnj;**Hardware Utilization**&zwnj;| Optimized for multi-core CPUs         | GPU-friendly                |
| &zwnj;**Industry Adoption**&zwnj;   | 85% of real-time systems              | 95% of medical imaging      |


