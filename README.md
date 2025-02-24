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
$$
\min_{\bm{x}} \|\bm{A}\bm{x} - \bm{b}\|_2^2
$$  

&zwnj;**Where:**&zwnj;
- $\bm{A} \in \mathbb{R}^{m \times n}$ (Full column rank)  
- $\bm{b} \in \mathbb{R}^{m}$  
- $\bm{x} \in \mathbb{R}^{n}$  

### Decomposition  
$$
\bm{A}\bm{\Pi} = \bm{Q}\bm{R}
$$  

| Symbol                | Description                                   |
|-----------------------|-----------------------------------------------|
| $\bm{\Pi}$            | Column permutation matrix                     |
| $\bm{Q}$              | Orthogonal matrix ($\bm{Q}^\top\bm{Q} = \bm{I}$) |  
| $\bm{R}$              | Upper triangular matrix                      |  


## Core Algorithm

### Step 1: Column Pivoting  
&zwnj;**At iteration $ k $:**&zwnj;
1. Compute remaining column norms:  
   $$
   \gamma_j = \|\bm{R}(k:m, j)\|_2^2 \quad \text{for } j = k:n
   $$  
2. Swap column $ p $ with $ k $-th column:  
   $$
   p = \arg\max_j \gamma_j
   $$

### Step 2: Householder Transformation  
&zwnj;**For column $ k $:**&zwnj;
1. Compute Householder vector $ \bm{v} $:  
   $$
   \bm{v} = \bm{x} \pm \|\bm{x}\|_2 \bm{e}_k, \quad \bm{x} = \bm{R}(k:m, k)
   $$  
2. Apply reflection to submatrix:  
   $$
   \bm{R}(k:m, k:n) \leftarrow \bm{R}(k:m, k:n) - \bm{v} \left( \frac{2}{\bm{v}^\top\bm{v}} \bm{v}^\top \bm{R}(k:m, k:n) \right)
   $$  
3. Update orthogonal matrix:  
   $$
   \bm{Q} \leftarrow \bm{Q} - \bm{Q}\bm{v} \left( \frac{2}{\bm{v}^\top\bm{v}} \bm{v}^\top \right)
   $$  


## Key Equations

### Triangular System Solving & Validation  

#### After Decomposition  
&zwnj;**Triangular System Solving:**&zwnj;
$$
\bm{R}_{11} \bm{z} = \bm{Q}_1^\top \bm{b}
$$  
- $\bm{R}_{11}$: Leading $n \times n$ submatrix of $\bm{R}$  

### Solution Recovery  
$$
\bm{x} = \bm{\Pi} \bm{z}
$$  

### Residual Orthogonality (Validation)  
$$
(\bm{A}\bm{x})^\top (\bm{b} - \bm{A}\bm{x}) \approx 0 \quad \text{(Machine precision $\sim 10^{-16}$)}
$$  


## Numerical Stability

### Key Mechanisms  

1. Column Pivoting  
$$
\|r_{kk}\| \geq \|r_{ij}\| \quad \text{for } i \geq k, \, j \geq k
$$  

2. Householder Reflections  
&zwnj;**Backward stable**&zwnj; with error bound:  
$$
\frac{\|\Delta\bm{A}\|_F}{\|\bm{A}\|_F} \leq O(mn \epsilon_{\text{machine}})
$$  

3. Explicit Orthogonality  
$$
\|\bm{Q}^\top\bm{Q} - \bm{I}\|_F \leq 10^{-14} \quad \text{(Typical in 2025 hardware)}
$$  


## Validation Metrics

| Metrics           | Formula                                                                  | Acceptable Threshold|
|:------------------|:------------------------------------------------------------------------:|------------------: |
| QR Accuracy       | $  \frac{\|\bm{A}\bm{\Pi} - \bm{Q}\bm{R}\|_F}{\|\bm{A}\|_F} < 10^{-12} $ | $ <10^{-12} $    |
| Orthogonality     | $\|\bm{Q}^\top\bm{Q} - \bm{I}\|_F < 10^{-14}$                            | $ <10^{-14} $    |
| Solution Error    | $\|\bm{A}\bm{x} - \bm{b}\|_2$                                            | Problem-dependent  |


## Comparison with SVD

| &zwnj;**Property**&zwnj;            | &zwnj;**QR (Column-Pivoted)**&zwnj;              | &zwnj;**SVD**&zwnj;                    |
|-------------------------|---------------------------------------|-----------------------------|
| &zwnj;**Flops**&zwnj;               | $ 2mn^2 - \frac{2}{3}n^3 $          | $ 4m^2n + 22n^3 $         |
| &zwnj;**Rank Detection**&zwnj;      | Threshold-based                       | Exact                       |
| &zwnj;**Hardware Utilization**&zwnj;| Optimized for multi-core CPUs         | GPU-friendly                |
| &zwnj;**Industry Adoption**&zwnj;   | 85% of real-time systems              | 95% of medical imaging      |


