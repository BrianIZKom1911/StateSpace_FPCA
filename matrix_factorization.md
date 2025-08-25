# Matrix Factorization
## Spectrum and Multiplicities
**Definition.** The **spectrum** of a square matrix is the set of its eigenvalues.
An square matrix for which all eigenvalues are distinct is said to have a *simple* spectrum. 

**Definitions.**
Algebraic multiplicity $\mu(\lambda)$ is the number of roots of the characteristic polynomial $det⁡(A-xI)$ equal to $\lambda$.
Geometric multiplicity $m(\lambda):=\text{dim}⁡(\ker(A-\lambda I))$ is the dimension of the eigenspace associated with eigenvalue $\lambda$. Some easy properties:
- $m(\lambda)\leq\mu(\lambda)$. When $m(\lambda)<\mu(\lambda)$, the eigenvalue $\lambda$ is said to be defective.
- If $\lambda$ is not repeated, then $m(\lambda)=\mu(\lambda)=1$
> *Proof.* Suppose the geometric multiplicity $m(\lambda)=g$, so there are g LID eigenvectors, say $v_1,\ldots,v_g$, associated with $\lambda$.
> We can extend this set to form a basis for the entire vector space $\mathbb{C}^n$. Denote the basis by ${v_1,\ldots,v_g,w_{g+1},\ldots,w_n}$ and define a change-of-basis matrix
> $$P=(v_1 \ldots v_g  w_{g+1} \ldots w_n)$$
> 
> Since the columns of $P$ form a basis, $P$ is invertible. Consider the matrix $A'$ (not a transpose) defined by $A'=P^{-1}AP$. $A'$ is similar to $A$, so it has the same charac-polys.
> Hence, we can examine the algebraic multiplicity of $\lambda$ for $A'$. 
> The columns of AP are:
> 
> $$Av_1=\lambda_1 v_1,\ldots,Av_g=\lambda_g v_g,$$
> $$Aw_j=\sum_{i=1}^g c_{ij}v_i+\sum_{i=g+1}^n c_{ij}w_i, \forall j>g$$
> 
> When we transform $AP$ into the new basis by multiplying by $P^{-1}$, the first $g$ columns of $AP$ becomes $P^{-1}Av_j=\lambda P^{-1}v_j$ where $P^{-1}v_j=e_j$
> (the column vector with 1 in the $j$th element). So,  
>
> $$
> A'=
> \begin{pmatrix}
> \lambda I_g & B\\
> 0 & C
> \end{pmatrix}
> $$
>  
> $$
> \begin{align}
> p_{A'}(x)&:=\det⁡(A'-xI)=\det⁡\begin{pmatrix}
> (\lambda-x)I_g & B\\
> 0 & C-xI_{n-g}
> \end{pmatrix}\\
> &=\det⁡((\lambda-x)I_g)\det⁡(C-xI_{n-g})\\
> &=(\lambda-x)^g\det⁡(C-xI_{n-g})
> \end{align}
> $$
> 
> So, $p_{A'}(x)$ has at least $g$ roots equal to $\lambda$. In other words, $\mu(\lambda)geqg$. ||

**Theorem** (*Schur decomposition*). Square matrix $A\in \mathbb{C}^{n\times n}$ can be expressed as

$$A=QUQ^{*}\quad\text{(called \textit{Schur form})}$$

for some unitary matrix $Q$, and some upper triangular matrix $U$.
- Unitary means the inverse is also the conjugate transpose, $Q^{-1}=Q^{*}$ (also denoted by $\bar{Q}^T$).
- Since $A$ and $U$ are similar, they have the same **spectrum**. And since $U$ is triangular, its eigenvalues are the diagonal entries.
> *Proof.* Suppose $\lambda_1, \lambda_2, \ldots, \lambda_n\in\mathbb{C}$ are the $n$ eigenvalues of $A$. Let $q_1$ be an eigenvector of norm 1 associated with $\lambda_1$.
> Pick any $n-1$ vectors that are of norm 1 and orthogonal to $q_1$. They form an orthonormal basis $Q_1=(q_1 q_2 \ldots q_n)$ of $\mathbb{C}^n$, $Q_1^{*}Q_1=I$.
> 
> $$
> Q_1^{*}AQ_1=\begin{pmatrix}
> \lambda_1 & A_{12}\\
> 0 & A_{22}
> \end{pmatrix}
> $$
> 
> Clearly $A_{22}$ has eigenvalues $\lambda_2, \ldots, \lambda_n$ because $\det⁡(A-\lambda I)=\det⁡(Q_1^{*}AQ_1-\lambda I)=(\lambda-\lambda_1)\det⁡(A_{22}-\lambda I)$.
> Now we prove the theorem by **induction**: 
>
> The theorem is trivial for $n=1$. Assume the theorem holds for $n=k$. For $n=k+1$, we proceed as above and apply the theorem to $A_{22}$ which is $k\times k$.
> So, $\exists Q_2$ that is unitary and $U_2$ upper triangular s.t. $A_{22}=Q_2 U_2 Q_2^{*}$. Define
> 
> $$
> Q=Q_1\begin{pmatrix}
> 1 & 0\\
> 0 & Q_2
> \end{pmatrix}
> $$
> 
> Then
> 
> $$\begin{align}
> AQ&=AQ_1\begin{pmatrix}
> 1 & 0\\
> 0 & Q_2
> \end{pmatrix}=Q_1\begin{pmatrix}
> \lambda_1 & A_{12}\\
> 0 & A_{22}
> \end{pmatrix}\begin{pmatrix}
> 1 & 0\\
> 0 & Q_2\\
> \end{pmatrix}\\
> &=Q_1\begin{pmatrix}
> \lambda_1 & A_{12}Q_2\\
> 0 & A_{22}Q_2
> \end{pmatrix}=Q_1\begin{pmatrix}
> \lambda_1 & A_{12}Q_2\\
> 0 & Q_2U_2
> \end{pmatrix}\\
> &=Q_1\begin{pmatrix}
> 1 & 0\\
> 0 & Q_2
> \end{pmatrix}\begin{pmatrix}
> \lambda_1 & A_{12}Q_2\\
> 0 & U_2
> \end{pmatrix}=Q\begin{pmatrix}
> \lambda_1 & A_{12}Q_2\\
> 0 & U_2
> \end{pmatrix}
> \end{align}$$
>
> Let $U$ be the second matrix, which is upper triangular because $U_2$ is. Hence, $AQ=QU \Leftrightarrow A=QUQ^{*}$. ||

## Diagonalization

