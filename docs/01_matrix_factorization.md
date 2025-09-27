# Matrix Factorization
## Spectrum and Multiplicities
**Definition.** The **spectrum** of a square matrix is the set of its eigenvalues.
An square matrix for which all eigenvalues are distinct is said to have a *simple* spectrum. 

**Definitions.**

Algebraic multiplicity $\mu(\lambda)$ is the number of roots of the characteristic polynomial $det⁡(A-xI)$ equal to $\lambda$.

Geometric multiplicity $m(\lambda):=\text{dim}⁡(\ker(A-\lambda I))$ is the dimension of the eigenspace associated with eigenvalue $\lambda$. Some easy properties:

- $m(\lambda)\leq\mu(\lambda)$. When $m(\lambda)<\mu(\lambda)$, the eigenvalue $\lambda$ is said to be *defective*.

- If $\lambda$ is not repeated, then $m(\lambda)=\mu(\lambda)=1$

> *Proof.* Suppose the geometric multiplicity $m(\lambda)=g$, so there are g LID eigenvectors, say $v_1, \ldots, v_g$, associated with $\lambda$.
> We can extend this set to form a basis for the entire vector space $\mathbb{C}^n$. Denote the basis by ${v_1, \ldots, v_g, w_{g+1}, \ldots, w_n}$ and define a change-of-basis matrix
> $$P=(v_1 \ldots v_g, w_{g+1} \ldots w_n)$$
> 
> Since the columns of $P$ form a basis, $P$ is invertible. Consider the matrix $A'$ (not a transpose) defined by $A'=P^{-1}AP$. $A'$ is similar to $A$, so it has the same charac-polys.
> Hence, we can examine the algebraic multiplicity of $\lambda$ for $A'$. 
> The columns of AP are:
> 
> $$Av_1=\lambda_1 v_1, \ldots, Av_g=\lambda_g v_g,$$
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
> So, $p_{A'}(x)$ has at least $g$ roots equal to $\lambda$. In other words, $\mu(\lambda)\geq g$. ||


**Theorem** (*Schur decomposition*). Square matrix $A\in \mathbb{C}^{n\times n}$ can be expressed as

$$A=QUQ^{*}\quad\text{(called \textit{Schur form})}$$

for some unitary matrix $Q$, and some upper triangular matrix $U$.
- Unitary means the inverse is also the conjugate transpose, $Q^{-1}=Q^{\*}$ (also denoted by $\bar{Q}^T$).
- Since $A$ and $U$ are similar, they have the same **spectrum**. And since $U$ is triangular, its eigenvalues are the diagonal entries.
> *Proof.* Suppose $\lambda_1, \lambda_2, \ldots, \lambda_n\in\mathbb{C}$ are the $n$ eigenvalues of $A$. Let $q_1$ be an eigenvector of norm 1 associated with $\lambda_1$.
> Pick any $n-1$ vectors that are of norm 1 and orthogonal to $q_1$. They form an orthonormal basis $Q_1=(q_1, q_2, \ldots, q_n)$ of $\mathbb{C}^n$, $Q_1^{\*}Q_1=I$.
> 
> $$
> Q_1^{\*}AQ_1=\begin{pmatrix}
> \lambda_1 & A_{12}\\
> 0 & A_{22}
> \end{pmatrix}
> $$
> 
> Clearly $A_{22}$ has eigenvalues $\lambda_2, \ldots, \lambda_n$ because $\det⁡(A-\lambda I)=\det⁡(Q_1^{\*}AQ_1-\lambda I)=(\lambda-\lambda_1)\det⁡(A_{22}-\lambda I)$.
> Now we prove the theorem by **induction**: 
>
> The theorem is trivial for $n=1$. Assume the theorem holds for $n=k$. For $n=k+1$, we proceed as above and apply the theorem to $A_{22}$ which is $k\times k$.
> So, $\exists Q_2$ that is unitary and $U_2$ upper triangular s.t. $A_{22}=Q_2 U_2 Q_2^{\*}$. Define
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
> Let $U$ be the second matrix, which is upper triangular because $U_2$ is. Hence, $AQ=QU \Leftrightarrow A=QUQ^{\*}$. ||

## Diagonalization
**Definition.** Square matrix $A\in\mathbb{C}^{n\times n}$ is **diagonalizable** if there exists a diagonal matrix $\Lambda$ such that

$$A=P\Lambda P^{-1}$$

Easy properties:
- If $A$ can be diagonalized, then $A^k=P\Lambda^k P^{-1}, k\in \mathbb{R}$.
- An $n\times n$ matrix can be diagonalized iff there exists $n$ linearly independent eigenvectors.
- Special case 1: An $n\times n$ matrix has $n$ distinct eigenvalues $\Rightarrow$ diagonalizable.

> *Proof.* Suppose $n\times n$ matrix $A$ has $n$ distinct eigenvalues. We aim to prove $A$ has $n$ linearly independent (LID) eigenvectors. Proof by **induction**: 
> 
> Trivially, $A$ has at least 1 LID eigenvector $v_1$.
> Now, assume $A$ has $k$ LID eigenvectors, $1\leq k\leq n$, and we show that it implies $A$ has $k+1$ LID eigenvectors.
> By assumption, there exists a set ${v_1,\ldots,v_{k}}$ of LID eigenvectors associated with eigenvalues $\lambda_1,\ldots,\lambda_k$.
> For an eigenvector $v_{k+1}$ according to a different eigenvalue $\lambda_{k+1}$, consider the equation (1)
>
> $$c_1 v_1+c_2 v_2+\ldots+c_{k+1} v_{k+1}=\mathbf{0}\quad (1)$$
>
> to see its solution.
> 
> Multiplying both sides by $A$ gives
>
> $$\begin{align}
> c_1 Av_1+c_2 Av_2+\ldots+c_{k+1}Av_{k+1}&=A\mathbf{0}\\
> c_1\lambda_1 v_1+c_2\lambda_2 v_2+\ldots+c_{k+1}\lambda_{k+1}v_{k+1}&=\mathbf{0}
> \end{align}$$
>
> Multiplying both sides by $\lambda_{k+1}$ gives
> 
> $$c_1\lambda_{k+1}v_1+c_2\lambda_{k+1}v_2+\ldots+c_{k+1}\lambda_{k+1}v_{k+1}=\mathbf{0}$$
>
> Subtracting the two equations above, we get
> 
> $$c_1(\lambda_1-\lambda_{k+1})v_1+c_2(\lambda_2-\lambda_{k+1})v_2+\ldots+c_k (\lambda_k-\lambda_{k+1})v_k=\mathbf{0}$$
>
> All coefficients must be zero since $v_1,v_2,\ldots,v_k$ are LID. But $\lambda_1\neq \lambda_{k+1},\ldots,\lambda_k\neq \lambda_{k+1}$ since they are distinct.
> So, $c_1=c_2=\ldots=c_k=0$. Thus equation (1) becomes $c_{k+1}v_{k+1}=\mathbf{0}$ which implies $c_{k+1}=\mathbf{0}$.
> Hence, (1) only has trivial solution; $v_1, \ldots, v_k, v_{k+1}$ are LID. ||

- Special case 2: An Hermitian [in the real case: symmetric] matrix is diagonalizable. Actually, it can be unitarily [orthogonally] diagonalized,
$A=P\Lambda P^{\*}\ [A=P\Lambda P^T]$ --> See Spectral Theorem I.


**Eigendecomposition**

Let's look at the factorization. Let $A$ be **diagonalizable** $n\times n$ matrix, then it has $n$ linearly independent eigenvectors $q_1,\ldots,q_n$. 
Then define $Q=(q_1 \ldots q_n)$, $\Lambda =\text{diag}(\lambda_1,\ldots,\lambda_n)$, we find

$$\begin{align}
AQ=A(q_1 \ldots q_n)=(\lambda_1 q_1 \ldots \lambda_n q_n)&=Q\Lambda \\
\Leftrightarrow A&=Q\Lambda Q^{-1}
\end{align}$$

($Q$ is invertible because its columns are LID.)

The LID eigenvectors ${q_i: i=1,\ldots,n}$ with nonzero eigenvalues form a **basis** (not necessarily orthonormal) for all possible products $Aw, \forall w\in\mathbb{C}^n$.

## Spectral Theorem
**Theorem** (spectral theorem I). Let $A$ be an $n\times n$ Hermitian matrix. There exists an unitary matrix $Q$ and a diagonal matrix $D$ such that $A=QDQ^{\*}$.

> *Proof.* By **Schur Theorem**, there exists an unitary matrix $Q$ s.t. $D=Q^{\*}AQ$ is an **upper triangular** matrix. D is also Hermitian:
>
> $$D^{\*}=(Q^{\*}AQ)^{\*}=Q^{\*}A^{\*}Q=Q^{\*}AQ=D$$
>
> An triangular and Hermitian matrix must be **diagonal**. ||

- An direct observation is that $Q^{\*}=Q^{-1}$, so $A$ is in fact unitarily diagonalized by $Q$, so we have

**Corollary.** Hermitian matrices are (unitarily) diagonalizable.


**Theorem** (spectral theorem II). For an $n\times n$ Hermitian matrix $A$:

(1) All its $n$ eigenvalues are real. (Though they can be repeated.) 

(2) Eigenspaces for distinct eigenvalues are orthogonal.

(3) Their algebraic multiplicities and geometric multiplicities agree. 

> *Proof.* (1) Suppose $\lambda$ is a complex eigenvalue, so $Av=\lambda v$ for some $v\neq 0$, then
> 
> $$\bar{v}^T Av=\lambda\bar{v}^T v=\lambda\lVert v\rVert^2$$
> 
> Taking conjugate transpose of both sides gives
>
> $$\bar{v}^T\bar{A}^T v=\bar{\lambda}\lVert v\rVert^2$$
>
> Since $A$ is Hermitian, $\bar{A}^T=A$, so $\bar{\lambda}\lVert v\rVert^2=\lambda\lVert v\rVert^2 \Leftrightarrow \lambda =\bar{\lambda}$, which proves $\lambda$ is a real number.
>
> (2) Suppose $Av=\lambda v$ and $Aw=\mu w$ with real numbers $\lambda\neq \mu$, then $\bar{v}^T\bar{A}^T=\bar{\lambda}\bar{v}^T=\lambda\bar{v}^T$. And
>
> $$\lambda \bar{v}^T w=\bar{v}^T\bar{A}^T w=\bar{v}^T Aw=\mu\bar{v}^T w$$
>
> Since $\lambda\neq \mu$, it must be that $\bar{v}^T w=0$. 
>
> (3) Again, use Schur Theorem. Consider the same $D=Q^{\*}AQ$. Since $Q^{\*}=Q^{-1}$, $D=Q^{-1}AQ$ is similar to $A$, and we can analyze the multiplicities for $D$ in place of $A$.
> 
> Since $D$ is triangular (here, diagonal), its diagonal entries must be the $n$ (real) eigenvalues of $A$.
> The algebraic multiplicity $\mu(\lambda)$ is the number of times it appears on the diagonal of $D$. 
> The kernel space (or null space) $\ker(D-\lambda I)$ is straightforward to find. The dimension of $\ker(D-\lambda I)$ is just the number of zeros on its diagonal,
> which is equal to the times $\lambda$ appears in $D$. Thus, $m(\lambda)=\mu(\lambda)$. ||


**Theorem** (spectral theorem III). Let $A$ be an $n\times n$ Hermitian matrix, then there exists orthogonal projections $P_1, \ldots, P_r$ 
and real numbers $\lambda_1, \ldots, \lambda_r$ such that

$$\begin{align}
(1)& \quad P_i P_j=0 \text{ if } i\neq j\\
(2)& \quad I_n=P_1+\ldots+P_r\\
(3)& \quad AP_i=\lambda_i P_i
\end{align}$$

This is the decomposition of the identity into **eigenspace projections**. In these notations, we have

$$A=\lambda_1 P_1+\ldots+\lambda_r P_r$$

> *Proof.* Let $\lambda_1, \ldots, \lambda_r$ be the distinct eigenvalues of $A$ which are real (by Spectral Theorem II). 
> Let $P_i$ be the orthogonal projection onto the eigenspace $\ker⁡(A-\lambda_i I)$, then we have $AP_i=\lambda_i P_i$. 
> Moreover, by Theorem II, the different eigenspaces are mutually orthogonal, so $P_i P_j=0$. Finally, $I_n=P_1+\ldots+P_r$ is equivalent to saying there is a basis of eigenvectors.

**Spectral Decomposition**

Now, suppose $(n\times n)$ Hermitian matrix $A$ has $n$ distinct eigenvalues $\lambda_1, \ldots, \lambda_n$, then there is such a decomposition $A=Q\Lambda Q^*$ where $Λ=\text{diag}(\lambda_1, \ldots, \lambda_n)$, $Q=(q_1 \ldots q_n)$ where each column is the corresponding eigenvectors of unit length. 
Expanding the RHS yields the spectral decomposition of $A$:

$$A=\sum_{i=1}^n \lambda_i q_i\bar{q}_i^T$$

> *Proof.* Since $A$ is Hermitian and eigenvalues $\lambda_1,\ldots,\lambda_n$ are distinct, the $n$ eigenvectors are orthogonal (and hence LID).
> Thus, we have the eigendecomposition $A=Q\Lambda Q^{-1}$ by the stated definitions of $\Lambda$ and $Q$.
> Now, we find
> 
> $$
> Q^{\*}Q=(\bar{q}_1^T \ldots \bar{q}_n^T)\begin{pmatrix}
> q_1\\
> \vdots\\
> q_n
> \end{pmatrix}=\begin{pmatrix}
> \bar{q}_1^T\bar{q}_1 & \ldots & \bar{q}_1^Tq_n\\
> \vdots & & \vdots\\
> \bar{q}_n^Tq_1 & \ldots & \bar{q}_n^Tq_n
> \end{pmatrix}=I_n
> $$
>
> since $\bar{q}_i^Tq_j=0$ for any $i\leq j$ (orthogonality) and $\bar{q}_j^Tq_j=1$ for all $j=1,\ldots,n$. Thus, $Q^{\*}=Q^{-1}$. ||

**Normal Matrix**

**Definition.** Complex square matrix $A$ is said **normal** if it commutes with its conjugate transpose $A^{\*}$:

$$A^{\*}A=AA^{\*}$$

Easy properties:

- $A$ is normal *iff* it is **unitarily similar** to a diagonal matrix.

> *Proof.*
>
> (Necessity $\Leftarrow$) It is obvious. If $\exists$ unitary $Q$ and diagonal $D$ s.t. $A=QDQ^{\*}$, then $A^{\*} A=(QD^{\*} Q^{\*})QDQ^{\*}=QD^{\*}DQ^{\*}$,
> while $AA^{\*}=QDQ^{\*}(QD^{\*}Q^{\*})=QDD^{\*}Q^{\*}$. Then conclusion follows from $D^{\*}D=\text{diag}(|d_{11}|^2,\ldots,|d_{nn}|^2 )=DD^{\*}$. 
>
> (Sufficiency $\Rightarrow$) We can use **Spectral Theorem I** or not.
> *Using it enables an easier proof, but a proof without calling the stricter version of the statement makes their logic clear.* So we will do it the hard way. 
>
> By Schur decomposition, $A=QUQ^{\*}$ which means $A$ is unitarily similar to an upper triangular $U$. It's easy to show $U$ is also normal. 
> We now show that being normal and upper triangular implies **diagonality**.
> 
> Proof by **induction**:
> 
> $n=1$ is a trivial case. When $n=2$,
> 
> $$U=\begin{pmatrix}
> u_{11} & u_{12}\\
> 0 & u_{22}
> \end{pmatrix}, \quad \text{Then } 
> U^{\*}=\begin{pmatrix}
> \bar{u}\_{11} & 0\\
> \bar{u}\_{12} & \bar{u}\_{22}
> \end{pmatrix}$$
> 
> $$\begin{align}
> 0=UU^{\*}-U^{\*}U &= \begin{pmatrix}
> |u_{11}|^2+|u_{12}|^2 & u_{12}\bar{u}\_{22}\\
> u_{22}\bar{u}\_{12} & |u_{22}|^2 
> \end{pmatrix} - \begin{pmatrix}
> |u_{11}|^2 & \bar{u}\_{11}u_{12}\\
> \bar{u}\_{12}u_{11} & |u_{12}|^2+|u_{22}|^2
> \end{pmatrix}\\
> &= \begin{pmatrix}
> |u_{12}|^2 & u_{12}\bar{u}\_{22} - \bar{u}\_{11}u_{12}\\
> u_{22}\bar{u}\_{12}-\bar{u}\_{12}u_{11} & -|u_{12}|^2
> \end{pmatrix}
> \end{align}$$
> 
> Therefore $|u_{12}|^2=0 \Leftrightarrow u_{12}=0$; $U$ is diagonal.
> 
> Assume it holds for $n=k-1$. Now we show that it holds for $n=k$. Write
>
> $$
> U=\begin{pmatrix}
> u_{11} & U_{12}\\
> 0 & U_{22}
> \end{pmatrix}
> $$
>
> where $U_{12}$ is $1\times (k-1)$, and $U_{22}$ is $(k-1)\times (k-1)$ upper triangular matrix. Similarly,
>
> $$
> 0=UU^{\*}-U^{\*} U=\begin{pmatrix}
> U_{12}U_{12}^{\*} & U_{12}U_{22}^{\*}-\bar{u}\_{11}U_{12}\\
> U_{22}U_{12}^{\*}-U_{12}^{\*}u_{11} & U_{22}U_{22}^{\*}-U_{12}U_{12}^{\*}-U_{22}^{\*}U_{22}\end{pmatrix}
> $$
>
> Thus $U_{12}U_{12}^{\*}=0 \Leftrightarrow U_{12}=0$. So $U_{22}U_{22}^{\*}-U_{22}^{\*}U_{22}=0$, or $U_{22}$ is normal.
> By the assumption of induction, $U_{22}$ is diagonal. Finally, $U$ is diagonal. ||

- That is to say, normal matrix is (unitarily) diagonalizable.

Obviously, Hermitian matrix is normal. Therefore, the spectral theorem I is just a stricter version of this property of normal matrix, whose proof would be redundant if we showed this first.

## Singular Value Decomposition (SVD)
The Schur decomposition and the spectral decomposition are special cases of the **singular value decomposition** for square matrices and Hermitian matrices, resp. It will be extended to matrices of any shape.

**Definition.** Suppose matrix $A$ represent the linear transformation from vector spaces $V^m$ to $V^n$. 
An non-negative real number $\sigma$⁠ is called a **singular value** if there exist unit vectors $u\in V^m$ and $v\in V^n$ s.t. 

$$Av=\sigma u,A^{\*}u=\sigma v$$

The vectors $u$ and $v$ are called **left-singular** and **right-singular vectors** for $\sigma$, resp.

**Theorem** (SVD). 
For an $m\times n$ matrix $A$, $\exists U, V$ that are $m\times m$ and $n\times n$ unitary matrix and $R$ that is $m\times n$ rectangular diagonal with non-negative real entries, s.t. 

$$A=URV^*$$

- Rectangular diagonal means all entries but those at $(i, i)$ are zero.

> *Proof of existence.*
> 
> WLoG, suppose $m\leq n$. Since $A^{\*}A$ is positive semi-definite and Hermitian, by spectral theorem  there exists an ⁠n\times n unitary matrix $V$ such that
>
> $$V^{\*}A^{\*}AV=D=\begin{pmatrix} 
> D_{11} & 0\\ 
> 0 & 0 
> \end{pmatrix}$$
>
> where $D_{11}$ is diagonal and positive definite, of dimension $\ell\times\ell$, with $\ell$ the number of non-zero eigenvalues of $A^{\*}A$. (Obviously $\ell\leq m$). 
> Here $V$ is by definition a matrix whose ith column is the $i$-th eigenvector of $A^{\*}A$, corresponding to the eigenvalue $d_{ii}=\lambda_i$.
> Moreover, the $k$-th column of $V$ for $k\geq\ell+1$ is an eigenvector of $A^{\*}A$ with eigenvalue $d_{kk}=0$.
> Thus we can write $V=(V_1\quad V_2)$  where the columns of $V_1$ and $V_2$ contain the eigenvectors of $A^{\*}A$
> corresponding to non-zero and zero eigenvalues, resp. $V$ being unitary implies $V_1^{\*}V_1=I_1, V_2^{\*}V_2=I_2$ and $V_1 V_1^{\*}+V_2 V_2^{\*}=I_n$
> where subscripts in $I_1$ and $I_2$ denote the dimension equals $\ell$ and $n-\ell$, resp. 
>
> Using this partition of $V$ the equation becomes
>
> $$\begin{pmatrix}
> V_1^* \\
> V_2^{\*}
> \end{pmatrix}A^{\*}A\begin{pmatrix}
> V_1 & V_2
> \end{pmatrix}=\begin{pmatrix} 
> V_1^{\*}A^{\*}AV_1 & V_1^{\*}A^{\*}AV_2\\
> V_2^{\*}A^{\*}AV_1 & V_2^{\*}A^{\*}AV_2
> \end{pmatrix}=\begin{pmatrix} D_{11} & 0\\
> 0 & 0
> \end{pmatrix}$$
>
> which implies that
>
> $$V_1^{\*}A^{\*}AV_1=D_{11}, V_2^{\*}A^{\*}AV_2=0$$
> 
> The second equation implies $AV_2=0$. [It's easy to show that $A^{\*}A=0 \Leftrightarrow A=0$]
>
> Let's now define $U_1=AV_1 D_{11}^{-1/2}$. Then 
> 
> $$
> U_1 D_{11}^{1/2}V_1^{\*}=AV_1 D_{11}^{-1/2}D_{11}^{1/2}V_1^{\*}=AV_1 V_1^{\*}=A(I_n-V_2 V_2^{\*})=A-(AV_2)V_2^{\*}=A
> $$
>
> ($U_1$ and $V_1$ are in general not unitary, since they might not be square.)
> However, we know that $U_1$ is $m\times\ell$ with $m\geq\ell$. Also, since 
>
> $$
> U_1^{\*}U_1=D_{11}^{-1/2}V_1^{\*}A^{\*}AV_1 D_{11}^{-1/2}=D_{11}^{-1/2}D_{11}D_{11}^{-1/2}=I_1
> $$
> 
> the columns in $U_1$ are orthonormal and can be extended to an orthonormal basis. This means we can choose $U_2$ s.t. $U=(U_1\quad U_2)$ is unitary. 
> Now, define 
>
> $$
> R=\begin{pmatrix}
> D_{11}^{1/2} & 0_{l\times(n-\ell)}\\
> 0_{(m-l)\times\ell} & 0_{(m-l)\times(n-l)}
> \end{pmatrix}  
> $$
> 
> that has $(n-m)$ zero rows removed from $D$ and thus has $(m-l)$ zero rows. Hence
> 
> $$
> \begin{pmatrix}
> U_1 & U_2
> \end{pmatrix}\begin{pmatrix}
> D_{11}^{1/2} & 0 \\
> 0 & 0
> \end{pmatrix}\begin{pmatrix}
> V_1^{\*} \\
> V_2^{\*}
> \end{pmatrix}=\begin{pmatrix}
> U_1 & U_2
> \end{pmatrix}\begin{pmatrix}
> D_{11}^{1/2} V_1^{\*} \\
> 0
> \end{pmatrix}=U_1 D_{11}^{1/2} V_1^{\*}=A
> $$
