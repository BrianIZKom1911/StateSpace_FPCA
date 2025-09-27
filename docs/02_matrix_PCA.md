Reference: Jolliffe, I.T. 2002. *Principal Component Analysis*, 2nd Ed. Springer, New York, NY.

# Background
## Derivation of Principal Components

Let $X=(X_1, X_2, \ldots, X_m)^T$ be the m-dimensional *random vector*. To do standard PCA, we are seeking the linear combinations called **principal components** (PC) satisfying: 

(1) $w_k^T X=(w_{k1}, \ldots, w_{km})(X_1, \ldots, X_m)^T$ has maximum variance; 

(2) $w_k^T X$ is uncorrelated with all preceding PCs, for $k=1, 2, \ldots, m$. 

- Up to $m$ PCs could be found, but it is hoped, in general, that most of the variation in $X$ will be accounted for by $p$ PCs, $p<<m$.

- For now, presume $\Sigma=var(X)$ is known, then the k-th PC turns out to be $Z_k=a_k^T X$ where $a_k$ is an eigenvector of $\Sigma$ corresponding to its k-th largest eigenvalue $\lambda_k$. Furthermore, if we choose $a_k$ to have unit length, $a_k^T a_k=1$, then $var(Z_k)=\lambda_k$.

**Derivation.**
(The maximization is subject to a constant norm of the coefficient vector; otherwise, the inner product can be as large as one wishes. For simplicity, just set it to unit length.)

> $$\underset{w_1\in \mathbb{R}}{\max} var(w_1^T X)=w_1^T\Sigma w_1 \text{ subject to } w_1^T w_1=1$$
>
> Using Lagrange multipliers, 
>
> $$\underset{w_1\in \mathbb{R}}{\max} var(w_1^T X)=w_1^T\Sigma w_1 - \lambda (w_1^Tw_1-1)$$
> 
> Differentiation w.r.t. $w_1$ gives $\Sigma w_1-\lambda w_1=0$ or $(\Sigma -\lambda I_m) w_1=0$. Thus, $\lambda$ is an eigenvalue of $\Sigma$ and $w_1$ is the corresponding eigenvector. Then, the quantity to be maximized is $w_1^T\Sigma w_1=w_1^T\lambda w_1=\lambda$, so $\lambda$ must be as large as possible, which is the largest eigenvalue of $\Sigma$, i.e., $\lambda =\lambda_1$. From now on, let $a_j$ denote the unit-length eigenvector associated with the j-th largest eigenvalue $\lambda_j$.
>
> Now, for $k≥2$, suppose the first $k$ PCs are given by $Z_j=Xa_j$; $a_j^Ta_j=1$; $cov(a_j^T X, a_i^T X)=0$ for any $i<j$, $j=k-1, \ldots, 1$. (For $k=1$ we do not need to assume any conditions beyond the first step and prove directly the nature of $w_2$.) We want to prove the coefficient defining the (k+1)th PC is $a_{k+1}$. That is, $a_{k+1}$ maximizes $w_{k+1}^T\Sigma w_{k+1}$ subject to $w_{k+1}^T w_{k+1}=1$ and $cov(w_{k+1}^T X, a_j^T X)=0$. Note that
>
> $$cov(w_{k+1}^T X,a_j^T X)=w_{k+1}^T\Sigma a_j=w_{k+1}^T\lambda_j a_j=\lambda_j w_{k+1}^T a_j$$
>
> Thus the zero correlation constraints are equivalent to $w_{k+1}^T a_j=0$ for any $j\leq k$.
>
> $$\underset{w_{k+1}}{\max} w_{k+1}^T\Sigma w_{k+1}-\lambda (w_{k+1}^T w_{k+1}-1)-\sum_{j\leq k}\delta_j (w_{k+1}^T a_j)$$
>
> Differentiation w.r.t. $w_{k+1}$ gives $\Sigma w_{k+1}-\lambda w_{k+1}-\sum \delta_j a_j=0$. Left multiplying $a_1^T$ on both sides gives
>
> $$\begin{align}
> a_1^T \Sigma w_{k+1}-\lambda a_1^T w_{k+1}-\sum \delta_j a_1^T a_j &=0\\
> \sum \delta_j a_1^T a_j &=0
> \end{align}$$
>
> But $a_1^T a_j=0$ for all $j\neq 1$ and $a_1^T a_1=1$. So $\delta_1=0$.
> Doing so with $a_2^T,\ldots,a_k^T$, it must be that $\delta_j=0$ for any $j\leq k$. Therefore, $\Sigma w_{k+1}-\lambda w_{k+1}=0$, or $(\Sigma -\lambda I_m) w_{k+1}=0$. $\lambda$ is an eigenvalue of $\Sigma$. Again, $w_{k+1}^T \Sigma w_{k+1}=\lambda$ must be as large as possible. Assuming $\Sigma$ does not have repeated eigenvalues*, $\lambda$ cannot equal any of $\lambda_1, \ldots, \lambda_k$. If it did, it would follow that $w_{k+1}=a_k$ violating $w_{k+1}^T a_k=0$. So, $\lambda =\lambda_{k+1}, w_{k+1}=a_{k+1}$.
> Hence, the proof follows from mathematical induction. ||

- $\Sigma$ might have repeated eigenvalues. That is, some eigenvalue's (algebraic or geometric) multiplicity $=q>1$. From [Spectral Theorem II](01_matrix_factorization/#spectral-theorem), since $\Sigma$ is symmetric (real and Hermitian), all its eigenvalues are real and the algebraic and geometric multiplicities agree.
In this case, the $q$ eigenvectors corresponding to the $q$ equal eigenvalues span a certain unique q-dimensional space, but, within this space, they are arbitrary.
Therefore, these $q$ PCs are not uniquely determined---Any set of orthogonal vectors that span this eigenspace can be chosen as the direction of each of the PCs. But the total variance explained by the components in this subspace is fixed. What we lose is the interpretability.
More discussions can be found in Jolliffe (2002), Sections 2.4, 3.7. 

## Properties (based on a known population variance)
Let $Z^q=(Z_1, Z_2, \ldots, Z_q)^T (q\leq m)$. $Z_k$ is the k-th PC with k-th largest variance $var(Z_k)=\lambda_k$ which equals the k-th (largest) eigenvalues. Then $Z^q={A^q}^T X$, where $A^q$ is the first $q$ columns of $A$, the $m\times m$ orthogonal matrix $(a_1 \ldots a_m)$ whose k-th column is the k-th eigenvectors of $\Sigma$. That is, the PCs are defined by an orthonormal linear transformation of $X$. We have

$$\Sigma A^q=A^q \Lambda_q$$

where $\Lambda_q=diag(\lambda_1,\ldots,\lambda_q)$. But only when $q=m$ do we have these equations:

$$\begin{align}
A^T\Sigma A &=\Lambda \\
\Sigma &=A\Lambda A^T
\end{align}$$

(Trivially, $A^1\equiv a_1; A^m\equiv A$)

> The derivation can be found in [Spectral Decomposition](01_matrix_factorization/#spectral-theorem).

**Property A1.** For any integer $q$, $1\leq q\leq m$, consider the orthonormal linear transformation $\xi=B^T X$ of $X$, where $\xi$ is a q-dimensional vector and $B$ is an $(m×q)$ matrix. Let $var(\xi)≔\Sigma_\xi=B^T \Sigma B$. Then $\tr(\Sigma_\xi)$ is maximized by taking $B=A^q$.

**Property A2.**
Let $\xi=B^T X$, then $tr(\Sigma_\xi)$ is minimized by taking $B=A_*^q$ where $A_*^q$ consists of the last $q$ columns of $A$.

**Property A3.** (spectral decomposition)

$$\Sigma =\lambda_1 a_1 a_1^T+\lambda_2 a_2 a_2^T+\ldots+\lambda_m a_m a_m^T$$

> This is a direct implication from $\Sigma =A\Lambda A^T$ by expanding the RHS matrix product.

- Looking at diagonal elements, we see $var(X_j)=\sum_{k=1}^m \lambda_k a_kj^2$. Not only can we decompose the variances of all the elements of $X$ into decreasing contributions $(\lambda_1,\ldots,\lambda_m)$ due to each PC, we can also decompose the whole covariance matrix $\Sigma$ into contributions $\lambda_k a_k a_k^T$ from each PC.
- Property A1 emphasizes that the PCs explain, successively, as much as possible of $tr(\Sigma)$. The current property shows, intuitively, that they also do a good job of explaining the off-diagonal elements of $\Sigma$.

**Corollary 1.** The linear combination of $X$ that has maximum variance conditional on the first $q$ PCs, $Z^q$, is precisely the (q+1)th PC. That is,

$$a_{q+1}= \underset{w}{\text{argmax }} var(Xw|Z^q)$$

Notice: This proposition is true ONLY IF $E(X|Z^q)$ is linear in $Z^q$ and $var(X|Z^q)$ is constant, which is a necessary condition but not pointed out in the textbook. Otherwise, the RHS is just $var(X-BZ^q)$, the variance of the best linear predictor of $X$ from $Z^q$, $B:=\Sigma_{xz}\Sigma_{zz}^{-1}$.

> *Sketched proof.* Under the above condition, we have a useful result that $var(X|Z^q)$ can be decomposed as $\Sigma -\Sigma_{xz}\Sigma_{zz}^{-1} \Sigma_{zx}$, where $\Sigma_{zz}=var(Z^q)$, $\Sigma_{xz}=cov(X,Z^q)$, $\Sigma_{zx}=\Sigma_{xz}^T$.
> The k-th column of $\Sigma_{xz}$ is $\lambda_k a_k$. The matrix $\Sigma_{zz}^{-1}=diag(\lambda_1^{-1},\ldots,\lambda_q^{-1}$), so
> 
> $$\Sigma_{xz}\Sigma_{zz}^{-1} \Sigma_{zx}=\sum_{k=1}^q\lambda_k a_k a_k^T$$
>
> and
>
> $$\Sigma -\Sigma_{xz}\Sigma_{zz}^{-1} \Sigma_{zx}=\sum_{k=q+1}^m\lambda_k a_k a_k^T$$
>
> Finding a linear function of $X$ having maximum conditional variance reduces to finding the eigenvalues and eigenvectors of the conditional covariance matrix. It's easy to verify that these are simply $(\lambda_{q+1},a_{q+1}),(\lambda_{q+2},a_{q+2}),\ldots,(\lambda_m,a_m)$. The eigenvector associated with the largest
of these eigenvalues is $a_{q+1}$. ||
