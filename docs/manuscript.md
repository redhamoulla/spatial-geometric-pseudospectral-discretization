# Geometric pseudo-spectral discretization: theory and simulations

This manuscript brings together a reading of the paper by Moulla, Lefèvre and
Maschke (2012), a reconstruction of its spatial method, its coupling to
projected-effort time integration, and the repository's reproducible simulations.

The spatial construction combines mixed polynomial spaces, an effort quotient
and boundary traces. The time discretization targets the balance

```math
H_{n+1}-H_n=\text{boundary work}-\text{dissipation},
```

under the quadrature and solver assumptions stated in the proofs.
The temporal ingredients are related to AVF, HBVM and Petrov–Galerkin methods.
The results distinguish reproduction of the paper, trajectory accuracy,
energy preservation and symplecticity.

The theory, proofs, figures and bibliography form a single continuous document.
Links to the code and data make each experiment verifiable.
Installation and execution commands are provided in the
[repository README](../README.md).

## Contents

1. [The 2012 paper and conventions](#1-the-2012-paper-and-conventions)
2. [Spatial construction and proofs](#2-spatial-construction-and-proofs)
3. [Time discretization](#3-time-discretization)
4. [Space-time assembly and interconnection](#4-space-time-assembly-and-interconnection)
5. [Numerical protocol and results](#5-numerical-protocol-and-results)
6. [References and relation to prior work](#6-references-and-relation-to-prior-work)

## 1. The 2012 paper and conventions

### 1.1. Source examined

R. Moulla, L. Lefèvre and B. Maschke, **Pseudo-spectral methods for the spatial
symplectic reduction of open systems of conservation laws**, *Journal of
Computational Physics* 231(4), 1272–1292 (2012).
[DOI](https://doi.org/10.1016/j.jcp.2011.10.008) ·
[HAL manuscript](https://hal.science/hal-01625008).

The reading is based on the author's copy of the HAL manuscript:
38 pages including the HAL cover, with the text dated December 16, 2010,
and deposited in HAL in 2017. Equation, proposition and table numbers below
refer to that copy. The publisher's typeset version has not been audited
separately. This repository contains an original reconstruction and commentary;
the source article file is not redistributed.

### 1.2. What the paper constructs

| Component | Location in the manuscript | Implementation in this repository |
|---|---|---|
| Differential-form spaces and open system | §2.3, Definitions 11–14 | Explicit sign convention and ports |
| Efforts in $`\mathbb P_N`$, flows and states in $`\mathbb P_{N-1}`$ | (31)–(36) | LagrangeBasis, MixedDiscretization |
| Exact differentiation and boundary traces | (36)–(38) | D, trace |
| Rectangular power pairing | (39)–(40) | M |
| Isotropy before reduction | Proposition 18 | Stokes identity test |
| Effort quotient by $`\ker M^\top`$ | (44), Proposition 20 | Reduced efforts and trace completions |
| Input–output realization | Proposition 21, (50) | Direct assembly of J, B, C, S |
| Constitutive closure | §4.1 and §4.2 | Quadratic, quartic and Saint-Venant energies |
| Ideal transmission-line spectrum | §6, Tables 1, 3, 4 | paper_spectrum.csv and an independent test |

The central property concerns the **spaces**, the power pairing and the
quotient; it does not require a particular family of Gauss nodes.
The geometry must survive an invertible change of polynomial coordinates.
The implementation checks Gauss, Lobatto, interior uniform and perturbed nodes.
Including boundary nodes through a Lobatto grid is a natural extension of the
polynomial argument; the manuscript's initial definitions use interior nodes.

The paper's term “spatial symplecticity” refers to preserving this structure
after reduction. By itself, it does not establish symplecticity of the time-step
map $`x_n\mapsto x_{n+1}`$.

### 1.3. Checking the conventions

We fix

```math
\partial_t q=-\partial_z e_p,\qquad
\partial_t p=-\partial_z e_q,\qquad
\dot H=e_q(a)e_p(a)-e_q(b)e_p(b).
```

With the definition $`D_{ki}=\varphi_i'(z_k)`$, integration by parts necessarily gives

```math
MD+D^\top M^\top=t_b^\top t_b-t_a^\top t_a.
```

Equation (43) in the HAL copy displays the opposite sign for the boundary terms
when combined literally with that definition of $`D`$.
The manuscript also alternates between some conventions for flows and state
derivatives. The code uses the convention above throughout, rather than
transplanting an isolated sign. The tests independently evaluate boundary traces,
power and integration by parts.

Additional points to keep in mind:

- The boundary factors in (67) are reversed relative to the stated conditions;
  (70) gives compatible factors.
- The final implication in the proof of Theorem 24 must retain the quadrature
  metric. Weighted skew-symmetry does not generally imply Euclidean
  skew-symmetry by simply dropping the weights.
- Remark 25 explicitly states a conjecture for classical Chebyshev collocation.
  The repository does not turn that conjecture into a theorem.
- The rank argument used here relies on injectivity of differentiation together
  with a trace, followed by the dimension of the reduced graph.

These observations concern the examined copy and the implementation conventions;
they do not undermine the central geometric construction.

### 1.4. Precise numerical reproduction

The published values correspond to a length $`\ell=2`$, capacitance
$`C=3`$, inductance per unit length $`L_{\rm ind}=2`$, eight states per field,
and homogeneous inputs. The length is identified here from the frequencies
in Table 1:

```math
\omega_j=\frac{(j+\frac12)\pi}{\ell\sqrt{C L_{\rm ind}}},\qquad j=0,\ldots,7.
```

The test compares **all eight** computed frequencies with the printed digits
in Table 4, including the high modes that a model of order 16 approximates
poorly. Reproducing the table does not mean that all eight modes accurately
approximate the continuum: the CSV records the reproduction error separately
from the error relative to the analytical frequencies.

The repository does not claim to reproduce the ill-conditioned classical
collocation variants in Tables 2 and 5. The nonlinear, open and interconnected
cases are identified as additional experiments.

## 2. Spatial construction and proofs

### 2.1. Spaces, differentiation and power

On $`[a,b]`$, choose two sets of distinct nodes. The sets need not coincide or
be nested. Their associated Lagrange bases are

```math
\varphi_i\in\mathbb P_N,\quad i=0,\ldots,N,\qquad
\psi_k\in\mathbb P_{N-1},\quad k=0,\ldots,N-1.
```

The states $`q,p`$ and their time derivatives use $`\psi`$; spatial efforts use
$`\varphi`$. Define

```math
D_{ki}=\varphi_i'(z_k),\qquad
M_{ik}=\int_a^b\varphi_i\psi_k\,dz,\qquad
t_a=(\varphi_i(a))_i,\quad t_b=(\varphi_i(b))_i.
```

Here $`D\in\mathbb R^{N\times(N+1)}`$ and $`M\in\mathbb R^{(N+1)\times N}`$.
The derivative of an element of $`\mathbb P_N`$ belongs exactly to
$`\mathbb P_{N-1}`$, so no derivative projection is required.

**Proposition S1 — Discrete Stokes identity.**

```math
MD+D^\top M^\top=t_b^\top t_b-t_a^\top t_a.
```

**Proof.** Entry $`(i,j)`$ of the first term is
$`\int\varphi_i\varphi_j'`$, since interpolation of $`\varphi_j'`$ in $`\psi`$
is exact. Add the transposed term and integrate
$`(\varphi_i\varphi_j)'`$. The identity is independent of the nodes. $`\square`$

In the implementation, polynomials are represented in a Legendre basis.
Gauss quadrature with $`N+1`$ points exactly integrates $`M`$, the mass matrices
and the required products. The **quadrature** nodes do not replace the user's
**representation** nodes.

### 2.2. Why the quotient is necessary

The internal power is

```math
(e_q^{\rm full})^\top M\dot q+(e_p^{\rm full})^\top M\dot p.
```

This pairing is degenerate on the full effort space. The efforts conjugate
to the state coefficients are

```math
g_q=M^\top e_q^{\rm full},\qquad
g_p=M^\top e_p^{\rm full}.
```

**Proposition S2 — Rank and completion.** $`M`$ has rank $`N`$.
The kernel of $`M^\top`$ is one-dimensional and represents the degree-$`N`$
Legendre polynomial mapped to $`[a,b]`$. The matrices

```math
C_a=\begin{bmatrix}M^\top\\t_a\end{bmatrix},\qquad
C_b=\begin{bmatrix}M^\top\\t_b\end{bmatrix}
```

are invertible.

**Proof.** If $`Mv=0`$, the polynomial $`\psi v`$ is orthogonal to all of
$`\mathbb P_N`$, hence to itself, so it vanishes. The transpose kernel
therefore represents the orthogonal complement of $`\mathbb P_{N-1}`$ in
$`\mathbb P_N`$. Its generator has nonzero boundary values proportional to
$`(-1)^N`$ and $`1`$. One trace removes this remaining degree of freedom. $`\square`$

This rank property explains why the two reduced effort vectors and two boundary
inputs suffice to reconstruct the full efforts.

### 2.3. Ports and the Dirac graph

Choose

```math
u=\begin{bmatrix}e_q(a)\\e_p(b)\end{bmatrix},\qquad
y=\begin{bmatrix}e_p(a)\\-e_q(b)\end{bmatrix}.
```

At each evaluation,

```math
e_q^{\rm full}=C_a^{-1}\begin{bmatrix}g_q\\u_1\end{bmatrix},\qquad
e_p^{\rm full}=C_b^{-1}\begin{bmatrix}g_p\\u_2\end{bmatrix},
```

```math
\dot q=-De_p^{\rm full},\qquad \dot p=-De_q^{\rm full}.
```

The outputs are obtained from the same traces. This elimination yields

```math
\dot x=Jg+Bu,\qquad y=Cg+Su,\qquad x=(q,p).
```

**Proposition S3 — Reduced structure.**

```math
J^\top=-J,\qquad C=B^\top,\qquad S^\top=-S.
```

**Proof.** Applying S1 to the reconstructed efforts gives
$`g^\top\dot x=u^\top y`$ for every $`(g,u)`$. By polarization, the matrix

```math
\mathcal J=\begin{bmatrix}J&B\\-C&-S\end{bmatrix}
```

is skew-symmetric. Its graph is isotropic and has dimension $`2N+2`$ in a bond
space of dimension $`4N+4`$, so it is maximally isotropic.
It defines a Dirac structure for the canonical pairing. $`\square`$

The term $`S`$ is generally nonzero. Removing it before interconnecting
subsystems changes the boundary relations. The code retains this term.
It measures the raw skew-symmetry residual; no operation $`(J-J^\top)/2`$
is used to hide an assembly error.

### 2.4. Constitutive closure

Restrict the **energy** to the state space, then differentiate it:

```math
H_N(x)=\int_a^b\mathcal H(\psi q,\psi p,z)\,dz,\qquad g=\nabla H_N(x).
```

For the ideal transmission line,

```math
H_N=\tfrac12 q^\top Q_Cq+\tfrac12p^\top Q_Lp,\qquad
(Q_C)_{ij}=\int_a^b\frac{\psi_i\psi_j}{C(z)}\,dz.
```

For Saint-Venant flow in a normalized rectangular channel with $`\rho=B=1`$,
the variables used are $`(h,v)`$, rather than $`(h,hv)`$:

```math
H_N(h,v)=\int_a^b\left(\tfrac g2 h_N^2+\tfrac12h_Nv_N^2\right)\,dz.
```

Thus $`e_h=gh+v^2/2`$ and $`e_v=hv`$. The energy is cubic in the coefficients.
Spatial quadrature is chosen to integrate this polynomial, its gradient and
its Hessian exactly. Underintegrating a constitutive energy defines a different
discrete energy; preserving that energy does not establish accuracy relative
to the continuum energy.

A dissipative closure adds $`-Rg`$, with $`R=R^\top\succeq0`$:

```math
\dot x=(J-R)\nabla H_N+Bu,\qquad
\dot H_N=u^\top y-g^\top Rg.
```

For the constant-parameter line, resistance on the second field is
$`R=\mathrm{diag}(0,rM_\psi^{-1})`$, where
$`M_\psi=\int\psi^\top\psi`$. This realizes the current term $`-r p/L_{\rm ind}`$.

### 2.5. Coordinate invariance and conditioning

A change of state coordinates $`x'=Tx`$ transforms the objects according to

```math
J'=TJT^\top,\quad B'=TB,\quad S'=S,\quad
H'_N(x')=H_N(T^{-1}x').
```

For a quadratic energy, $`Q'=T^{-\top}QT^{-1}`$. The dynamics are similar,
and power and the spectrum are identical. This conclusion assumes the same
constitutive functional, or quadratures yielding the same restriction.

The invariance is algebraic; it does not provide immunity to roundoff.
Uniform and perturbed nodes become poorly conditioned as $`N`$ increases.
The results report condition numbers and the corresponding residuals.
Gauss and Lobatto are the repository's usual computational choices.

## 3. Time discretization

### 3.1. Problem and relation to existing methods

After spatial reduction,

```math
\dot x=(J-R)\nabla H_N(x)+Bu,\qquad
y=B^\top\nabla H_N(x)+Su,
```

with $`J^\top=-J`$, $`R^\top=R\succeq0`$ and $`S^\top=-S`$.
The implementation treats these matrices as constant. $`H_N`$ may be nonlinear
and $`u(t)`$ time-dependent. An extension to state-dependent matrices is stated
in §3.7; it is not implemented.

The time discretization uses a **continuous polynomial state and a projected
effort**. For a constant conservative structure, it belongs to the AVF/HBVM
family. Its Petrov–Galerkin interpretation and its extension to ports have
explicit precedents, discussed in the
[annotated references](#6-references-and-relation-to-prior-work).
The contribution of this repository is a complete, verifiable and reproducible
construction coupled to the 2012 spatial method; no priority is claimed
for these time integrators.

### 3.2. Two temporal spaces

On a step $`[t_n,t_n+h]`$, write $`c=(t-t_n)/h`$ and choose

```math
X(c)\in\mathbb P_s([0,1];\mathbb R^d),\quad X(0)=x_n,\qquad
V(c)=\frac1hX'(c)\in\mathbb P_{s-1}.
```

Let $`P_0,\ldots,P_{s-1}`$ be the orthonormal shifted Legendre polynomials,
with $`P_0=1`$, and let $`I_j(c)=\int_0^cP_j(\tau)\,d\tau`$. Parameterize

```math
X(c)=x_n+h\sum_{j=0}^{s-1}I_j(c)a_j,\qquad
V(c)=\sum_{j=0}^{s-1}P_j(c)a_j.
```

The initial trace removes the constant kernel of differentiation.
The unknowns are the $`sd`$ coefficients $`a_j`$, and
$`x_{n+1}=X(1)=x_n+ha_0`$.

Define the temporal effort by orthogonal projection:

```math
\gamma_j=\int_0^1P_j(c)\nabla H_N(X(c))\,dc,\qquad
E(c)=\sum_{j=0}^{s-1}P_j(c)\gamma_j.
```

In general, $`E`$ is not pointwise equal to $`\nabla H_N(X)`$.
It is the effort **conjugate to the discrete flow** $`V`$ in the temporal
power pairing. This distinction is necessary for a nonquadratic energy.

### 3.3. Computable scheme

A positive quadrature rule $`(c_\ell,b_\ell)_{\ell=1}^k`$, exact on
$`\mathbb P_{2s-2}`$, defines the weak dynamics:

**(T1) Stage equations.**

```math
a_i=\sum_{\ell=1}^k b_\ell P_i(c_\ell)
\big[(J-R)E(c_\ell)+Bu(t_n+hc_\ell)\big].
```

The outputs used in the boundary work are

**(T2) Port outputs.**

```math
Y_\ell=B^\top E(c_\ell)+Su_\ell.
```

The implementation also uses the same Gauss quadrature to compute $`\gamma_i`$:

**(T3) Quadrature projection.**

```math
\gamma_i^Q=\sum_\ell b_\ell P_i(c_\ell)\nabla H_N(X(c_\ell)).
```

The theoretical identity with exact projection must be distinguished from
its implementation with a quadrature approximation.

**Proposition T1 — Energy balance over one step.** Suppose the coefficients
$`\gamma_i`$ are computed exactly and T1 is solved exactly. Then

**(T4) One-step energy balance.**

```math
\boxed{
H_N(x_{n+1})-H_N(x_n)
=h\sum_\ell b_\ell
\left[u_\ell^\top Y_\ell-E_\ell^\top R E_\ell\right].
}
```

**Proof.** Applying the continuous chain rule along the polynomial path gives

```math
\Delta H_N=\int_0^1\nabla H_N(X(c))^\top X'(c)\,dc
=h\sum_i\gamma_i^\top a_i.
```

Insert T1 and interchange the sums:

```math
\Delta H_N=h\sum_\ell b_\ell E_\ell^\top
\big[(J-R)E_\ell+Bu_\ell\big].
```

The terms $`E_\ell^\top JE_\ell`$ and $`u_\ell^\top Su_\ell`$ vanish.
Use T2. $`\square`$

Without input, energy cannot increase; without input or dissipation, it is
preserved at the **step endpoints**. No claim is made that it remains constant
at every point along the interior polynomial. To interpret the balance as
a passivity inequality with nonnegative storage, $`H_N`$ must be bounded below
on the relevant domain.

### 3.4. Quadrature requirement and special cases

If $`H_N`$ is a polynomial of total degree $`r`$, then
$`P_i(c)\nabla H_N(X(c))`$ has degree at most $`rs-1`$.
Gauss quadrature with $`k`$ points is therefore sufficient when

**(T5) Sufficient quadrature conditions.**

```math
k\ge s,\qquad 2k\ge rs.
```

| Energy | Degree $`r`$ | Temporal state degree $`s=1`$ | Temporal state degree $`s=2`$ |
|---|---:|---:|---:|
| Ideal transmission line | 2 | $`k=1`$ | $`k=2`$ |
| Rectangular-channel Saint-Venant | 3 | $`k=2`$ | $`k=3`$ |
| Line with quartic storage | 4 | $`k=2`$ | $`k=4`$ |

With $`s=1`$, the effort is the average gradient along the segment joining the
two states, giving AVF. For a quadratic Hamiltonian, this average equals
the gradient at the midpoint, recovering the implicit midpoint rule.

With $`k=s`$, the projected effort and the gradient coincide at the Gauss nodes:
the scheme is classical Gauss collocation. With $`k>s`$, it is the
$`\mathrm{HBVM}(k,s)`$ realization for a constant conservative structure.
For smooth autonomous problems in this class, the endpoint order is $`2s`$.
The simulations measure orders 2 and 4 against independent trajectory references.

The proof of T4 is not a general proof of order $`2s`$ for nonautonomous,
variable-matrix or constrained extensions. The repository's order claims
are limited to the cited classes and the measured cases.

For a nonpolynomial Hamiltonian, no finite $`k`$ generally guarantees exactness
of T3. Quadrature must be refined, or an exact discrete-gradient closure adopted.
The chain-rule defect is then a quantity to monitor, rather than set
artificially to zero.

### 3.5. Numerical certificate and three error sources

Let $`\rho_i`$ be the residual of T1 with sign convention
$`\rho_i=a_i-\mathrm{RHS}_i`$, and define

```math
\varepsilon_{\rm ch}=\Delta H_N-h\sum_i(\gamma_i^Q)^\top a_i.
```

For an exactly skew-symmetric structure, the following a posteriori identity holds:

**(T6) Residual balance identity.**

```math
\Delta H_N-W_n+D_n
=\varepsilon_{\rm ch}+h\sum_i(\gamma_i^Q)^\top\rho_i.
```

In floating-point arithmetic, the spatial structure defect contributes an
additional term. The program independently computes endpoint energy, boundary
work, dissipation, the chain-rule defect and the Newton residual.
It does not obtain energy by accumulating the very balance being checked.

The solver uses Newton iterations with an analytical Hessian and a line search.
Failure to converge raises an exception. Excellent energy preservation alone
does not guarantee accurate phase or trajectories, which is why the convergence
tables and comparison with the continuum are essential.

If $`H_N`$ is $`C^2`$ and the input is sufficiently regular, the derivative of the
residual with respect to $`a`$ is the identity at $`h=0`$. The implicit function
theorem gives a locally unique branch for sufficiently small steps.
This establishes neither global existence of an implicit solution for every
step size nor Newton convergence from every initial guess.

### 3.6. Symplecticity is a separate property

For a closed system with invertible $`J`$, a map $`\Phi_h`$ is symplectic if

```math
(D\Phi_h)^\top J^{-1}D\Phi_h=J^{-1}.
```

Gauss schemes satisfy this property but do not generally preserve a nonquadratic
Hamiltonian. The AVF/HBVM variants used here preserve polynomial energy with
adequate quadrature but are not generally symplectic. The repository measures
both defects on a quartic oscillator, using a tangent matrix obtained by
implicit differentiation and checked against finite differences.

Energy preservation alone also does not establish preservation of every Casimir,
positive water depth, absence of oscillations or stability of a noncoercive system.

### 3.7. Algebraic extension to variable matrices

In T1, replace $`J,R,B,S`$ with their values at $`X(c_\ell)`$.
If $`J(X)^\top=-J(X)`$, $`R(X)\succeq0`$ and $`S(X)^\top=-S(X)`$ at each stage,
the same proof gives T4 with the stage matrices: it only uses pointwise
skew-symmetry and positivity.

This still requires projecting the effort in **every** power term and exactly
integrating that projection whenever an exact balance is claimed.
Regularity, order and solvability of the resulting equations require analysis
for the particular model. The current solver, whose Jacobian exploits constant
matrices, does not implement this extension.

## 4. Space-time assembly and interconnection

### 4.1. Why a tensor product is not sufficient

Exact differentiation and a discrete Stokes identity address the **bilinear**
geometry. They do not automatically supply the chain rule for a nonquadratic
Hamiltonian. A complete construction therefore includes:

1. a compatible polynomial complex and its traces;
2. an effort quotient that makes the power pairing nondegenerate;
3. a constitutive closure derived from the spatial energy;
4. a temporal closure using a projected effort;
5. ports evaluated using that same effort.

The integrator implements items 3–5 on each time slab.
The tensor identities below check items 1–2 at the geometric level.
The tensorization tests are distinct from the dynamical simulations.

### 4.2. Polynomial complex on a rectangle

Write $`E_t=\mathbb P_s`$, $`F_t=\mathbb P_{s-1}\,dt`$,
$`E_z=\mathbb P_N`$ and $`F_z=\mathbb P_{N-1}\,dz`$. On one slab,

```math
V^0=E_t\otimes E_z,\qquad
V^1=(F_t\otimes E_z)\oplus(E_t\otimes F_z),\qquad
V^2=F_t\otimes F_z.
```

With orientations $`dt,dz`$, the matrices of the complex are

```math
d_0=
\begin{bmatrix}
D_t\otimes I_{E_z}\\ I_{E_t}\otimes D_z
\end{bmatrix},\qquad
d_1=
\begin{bmatrix}
-I_{F_t}\otimes D_z&D_t\otimes I_{F_z}
\end{bmatrix}.
```

**Proposition ST1.** $`d_1d_0=0`$.

**Proof.** The two terms in the product are
$`-D_t\otimes D_z`$ and $`D_t\otimes D_z`$. They cancel exactly.
The property does not depend on a particular node family. $`\square`$

### 4.3. Boundary identity on a slab

For $`\alpha\in\{t,z\}`$, define

```math
S_\alpha=M_\alpha D_\alpha,\qquad
G_\alpha=\int\varphi_\alpha^\top\varphi_\alpha,\qquad
B_\alpha=t_{\alpha,+}^\top t_{\alpha,+}
-t_{\alpha,-}^\top t_{\alpha,-}.
```

Here $`S_\alpha`$ denotes an integration-by-parts matrix, unrelated to
the port feedthrough matrix $`S`$.
The one-dimensional matrices satisfy
$`S_\alpha+S_\alpha^\top=B_\alpha`$. It follows that

```math
\mathcal K=S_t\otimes G_z+G_t\otimes S_z,
```

```math
\boxed{\mathcal K+\mathcal K^\top
=B_t\otimes G_z+G_t\otimes B_z.}
```

**Proof.** Transpose the Kronecker products and use the symmetry of
$`G_t,G_z`$. $`\square`$

The right-hand side contains the initial and final faces and the two spatial
boundaries. Internal faces of adjacent slabs cancel when they use the same
traces with opposite orientations.

This is an exact identity on polynomial spaces, not a proof of multisymplectic
conservation for every nonlinear closure. A multisymplectic law would require
studying the variations and two-forms associated with the discrete dynamics.
No such claim is made here.

### 4.4. Balance theorem for the complete simulation

**Proposition ST2.** Under the assumptions of Propositions S1–S3 and T1,
for step sizes $`h_n>0`$, the simulation satisfies

```math
H_N(x_M)-H_N(x_0)
=\sum_{n=0}^{M-1}h_n\sum_\ell b_\ell
\left[u_{n\ell}^\top Y_{n\ell}
-E_{n\ell}^\top R E_{n\ell}\right].
```

**Proof.** Spatial reduction provides the matrices from S3 and the closure
$`g=\nabla H_N`$. Apply T4 on each slab and sum: energies at the internal
endpoints telescope. $`\square`$

The balance remains valid for variable step sizes under the same integration
and solver assumptions on each step. It concerns $`H_N`$, the energy of the
spatially discrete model. The error relative to the continuum is evaluated
separately.

### 4.5. Composition of open systems

Assemble several systems in block-diagonal form, then connect the ports through

```math
u=Ky+v,\qquad K^\top=-K.
```

At each stage, $`u^\top y=v^\top y`$: the interconnection injects no power.
Since $`y=B^\top E+Su`$,

```math
(I-KS)u=KB^\top E+v.
```

If $`I-KS`$ is invertible, set $`F=(I-KS)^{-1}`$. The composed system is

```math
J_{\rm c}=J+BFKB^\top,\qquad
B_{\rm c}=BF,\qquad S_{\rm c}=SF.
```

Retain $`R`$ and add the Hamiltonians.
The identities $`FK=KF^\top`$ and $`SF=F^\top S`$ imply
$`J_{\rm c}^\top=-J_{\rm c}`$ and $`S_{\rm c}^\top=-S_{\rm c}`$;
the output is $`B_{\rm c}^\top E+S_{\rm c}v`$.

**Proposition ST3.** The sum of the component balances T4 is the balance of the
composed system when the components share the same time stages and the same
power quadrature.

**Proof.** The internal powers cancel through $`y^\top Ky=0`$ at every stage,
before summing over time. $`\square`$

The experiment connects the right boundary of line A to the left boundary
of line B: current and voltage are continuous, with the orientations defined
in §2.3. The remaining two ports have zero inputs. Energy actually transfers
from one line to the other.

If $`I-KS`$ is singular, its inverse must not be replaced with an arbitrary
pseudoinverse: constraints may arise and require a differential-algebraic
port-Hamiltonian formulation. The code reports this case.
Interconnections with different clocks are not supported.

### 4.6. Open research questions

| Topic | Status in this version |
|---|---|
| One-dimensional spatial geometry with arbitrary nodes | Proof and implementation |
| Nonlinear temporal energy balance | Proof, implementation and simulations |
| Tensor-product complex and Stokes identity | Proved and numerically checked |
| Synchronous port interconnection | Proof and simulation |
| PH systems with variable matrices | Algebraic balance proved; solver not implemented |
| PH systems with algebraic constraints | Not solved in this version |
| Two- and three-dimensional spatial geometries | Not implemented |
| Shocks and guaranteed positivity for Saint-Venant | Not established |
| Multisymplecticity of the complete scheme | Not claimed |
| Publishable mathematical novelty | Not established by this construction alone |

A concrete extension would preserve energetic composition when subdomains use
different time steps. Power cancellation would then need to be transferred
between distinct temporal spaces, rather than assumed at shared stages.

## 5. Numerical protocol and results

The results were computed on October 10, 2026, by
[experiments/reproduce.py](../experiments/reproduce.py).
The [manifest](../results/manifest.json) records the Python and NumPy versions,
platform and SHA-256 hashes of the computational sources.

### 5.1. Reproducing the paper

Ideal transmission line: $`\ell=2`$, $`C=3`$, $`L_{\rm ind}=2`$, $`N=8`$,
with homogeneous ports. All eight frequencies in Table 4 are recovered for
four node families. The maximum discrepancy from the published digits is
$`1.26\,10^{-12}`$.

Transformations from Gauss coordinates to perturbed-node coordinates give
relative defects of $`2.07\,10^{-13}`$ for $`J`$ and $`2.47\,10^{-15}`$ for $`Q`$.
The constitutive quadratures are exact.

[Frequencies and errors relative to the continuum](../results/tables/paper_spectrum.csv).

### 5.2. Geometry and conditioning

| Check | Result |
|---|---:|
| Maximum relative Dirac defect, four node families, $`N\le12`$ | $`4.59\,10^{-11}`$ |
| Same measure at $`N=16`$ | $`1.91\,10^{-8}`$ |
| Tensor-product Stokes identity | $`1.02\,10^{-14}`$ |
| Identity $`d_1d_0=0`$, normalized residual | $`7.78\,10^{-18}`$ |

The deterioration at $`N=16`$ comes from uniform or perturbed coordinates:
the condition number of a completion reaches approximately $`4.1\,10^8`$.
This does not contradict invariance in exact arithmetic.
The structural validation threshold of $`10^{-8}`$ is explicitly applied to
$`N\le12`$; the results at $`N=16`$ are still published.
The time-domain simulations use Gauss coordinates.

[Structure and condition numbers](../results/tables/spatial_structure.csv) ·
[Space-time complex](../results/tables/space_time_geometry.csv).

### 5.3. Trajectory convergence

Linear temporal convergence is measured at fixed spatial resolution against
the exponential of the semidiscrete model, computed independently through
diagonalization in energy coordinates.

Nonlinear temporal convergence uses a fine explicit RK4 reference independent
of the implicit solver. Increasing the reference from 6,000 to 12,000 steps
changes the final state by $`2.86\,10^{-14}`$ in the Euclidean norm.

| Model | Temporal degree | Order measured from the two finest step sizes |
|---|---:|---:|
| Ideal transmission line | 1 | 1.999 |
| Ideal transmission line | 2 | 3.986 |
| Quartic line | 1 | 2.000 |
| Quartic line | 2 | 3.987 |

![Temporal convergence](../results/figures/time_convergence.svg)

An analytical standing wave independently measures the spatial and space-time
errors. At $`N=12`$, the spatial $`L^2`$ error is $`1.46\,10^{-8}`$;
with temporal degree 2 and $`h=0.025`$, the total error is also approximately
$`1.46\,10^{-8}`$. This is a smooth analytical test case, not a uniform estimate
for every port-Hamiltonian equation.

![Spatial convergence](../results/figures/space_convergence.svg)

[Temporal data](../results/tables/time_convergence.csv) ·
[Nonlinear data](../results/tables/nonlinear_convergence.csv) ·
[Spatial data](../results/tables/space_convergence.csv).

### 5.4. Quartic Hamiltonian and positive and negative controls

Normalized closed model:

```math
H_N=\int_0^1\left[\tfrac12(q_N^2+p_N^2)+\tfrac14q_N^4\right]dz.
```

$`N=6`$, $`h=0.03`$, 1,000 steps, $`T=30`$.
This experiment evaluates integration of the semidiscrete model;
it is not a PDE convergence study over that time horizon.

| Method | $`(s,k)`$ | Maximum relative energy error |
|---|---:|---:|
| Implicit midpoint | (1,1) | $`2.71\,10^{-3}`$ |
| AVF with exact quadrature | (1,2) | $`9.23\,10^{-15}`$ |
| Two-stage Gauss | (2,2) | $`4.92\,10^{-6}`$ |
| Degree-2 projected effort | (2,4) | $`2.42\,10^{-14}`$ |

![Nonlinear energy preservation](../results/figures/nonlinear_energy.svg)

The quartic oscillator distinguishes energy preservation from symplecticity.
For one step of size $`0.2`$, implicit midpoint has a symplectic defect of
$`2.86\,10^{-17}`$ but an energy error of $`2.88\,10^{-4}`$.
AVF preserves energy to roundoff but has a symplectic defect of $`2.58\,10^{-3}`$.

[Oscillator measurements](../results/tables/symplectic_vs_energy.csv).

### 5.5. Active ports and dissipation

The same quartic storage is used with resistance $`r=0.15`$ on the second field,
two nonzero time-dependent inputs, $`h=0.02`$, $`T=10`$, and method $`(s,k)=(2,4)`$.
Energy is computed directly from the state; work and dissipation are integrated
independently at the stages.

```math
\max_n|H_n-H_0-W_{0:n}+D_{0:n}|=1.23\,10^{-15}.
```

Total dissipation is approximately $`0.10906`$, and all its increments are positive.
The maximum single-step balance defect is $`5.42\,10^{-17}`$.

![Open-system energy balance](../results/figures/forced_balance.svg)

[Step-by-step balance](../results/tables/forced_damped.csv).

### 5.6. Nonlinear Saint-Venant equations

Normalized rectangular channel, $`g=1`$, $`N=12`$:

```math
h(z,0)=1+0.02\sin(\pi z/2),\qquad v(z,0)=0.
```

The left input imposes the Bernoulli effort $`e_h=1`$; the right input imposes
the discharge $`e_v=0`$. The balance includes the actual work at the left boundary.
Method $`(s,k)=(2,3)`$, time step $`h_t=0.002`$, $`T=0.6`$.

- Cumulative balance defect: $`4.24\,10^{-16}`$.
- Vector-field residual at the lake-at-rest state: $`8.15\,10^{-13}`$.
- Minimum depth at 401 sampled points at every time-step endpoint: $`0.999866`$.
- Final RMS depth difference against $`N=18`$ and $`h_t=0.001`$: $`3.30\,10^{-7}`$.

![Saint-Venant dynamics](../results/figures/shallow_water.svg)

Positive **sampled** depth does not prove positivity throughout the spatial
interval or between time steps. The refinement changes space and time
simultaneously, so it does not measure their orders separately.
The experiment remains smooth; shocks, dry beds and hydraulic jumps are
not treated.

[Balance](../results/tables/shallow_water_ledger.csv) ·
[Profiles](../results/tables/shallow_water_profiles.csv).

### 5.7. Interconnection

Two lines with $`N_A=6`$, $`N_B=5`$ and different constitutive coefficients.
Initially, energy is present only in A. Port feedthrough is retained during
the elimination.

- Maximum total balance defect: $`3.47\,10^{-16}`$.
- Maximum internal power cancellation defect: $`8.33\,10^{-17}`$.
- Maximum energy reached in B: $`0.03839`$.

![Energy exchange](../results/figures/interconnection.svg)

[Data](../results/tables/interconnection.csv).

### 5.8. Scope of the validation

The unit tests examine the identities and several failure modes:
underintegration, incorrect ports, confusion between energy and symplecticity,
and Newton failure. The 22 global criteria are defined in validate() and recorded
in [summary.json](../results/summary.json).

The thresholds are explicit and allow for platform-dependent roundoff.
They validate the stated experiments and identities under the assumptions
of Sections 2–4; they do not provide a general certification of
port-Hamiltonian models.

## 6. References and relation to prior work

References were consulted on October 10, 2026.
This list focuses on the construction in this repository; it does not claim
to establish priority through an exhaustive review of the literature.

### 6.1. Spatial method

**Moulla, Lefèvre, Maschke (2012).** *Pseudo-spectral methods for the spatial
symplectic reduction of open systems of conservation laws*. JCP 231(4),
1272–1292. [DOI](https://doi.org/10.1016/j.jcp.2011.10.008).

The direct source for mixed spaces, the rectangular pairing, the quotient
and ports. The complete HAL manuscript was read. The audit and mapping to the
code are given in [Section 1](#1-the-2012-paper-and-conventions).

### 6.2. Time integration

**Hairer (2010).** *Energy-preserving variant of collocation methods*.
[Author's text](https://www.unige.ch/~hairer/preprints/crete.pdf).

Develops collocation variants extending AVF to high order.
It distinguishes energy preservation from symplectic properties.
The repository carries this distinction into its negative controls.

**Brugnano, Iavernaro, Trigiante.** *Analysis of Hamiltonian Boundary Value
Methods (HBVMs): A class of energy-preserving Runge–Kutta methods for the
numerical solution of polynomial Hamiltonian systems*.
[2009 preprint](https://arxiv.org/abs/0909.5659) ·
[2015 publication](https://doi.org/10.1016/j.cnsns.2014.05.030).

The conservative constant-matrix cases of our integrator are $`\mathrm{HBVM}(k,s)`$.
Path degree and quadrature size play different roles: increasing quadrature
can make the balance exact without increasing the number of polynomial unknowns.

**Celledoni, Høiseth (2017, preprint).** *Energy-Preserving and
Passivity-Consistent Numerical Discretization of Port-Hamiltonian Systems*.
[arXiv](https://arxiv.org/abs/1706.08621).

Presents discrete-gradient and splitting discretizations with discrete passivity.
This is a direct precedent for any general claim of energy-preserving
port-Hamiltonian integration.

**Kotyczka, Lefèvre (2018, preprint).** *Discrete-time port-Hamiltonian
systems: A definition based on symplectic integration*.
[arXiv](https://arxiv.org/abs/1811.07852).

Constructs a discrete-time Dirac structure through collocation and symplectic
integration, distinguishing discrete balance from energy error.
This framework does not justify calling the balance of an original nonquadratic
Hamiltonian “exact” under insufficient quadrature.

**Giesselmann, Karsai, Tscherpel (2025).** *Energy-consistent Petrov–Galerkin
time discretization of port-Hamiltonian systems*. SMAI-JCM 11, 335–367.
[Article](https://smai-jcm.centre-mersenne.org/articles/10.5802/smai-jcm.127/) ·
[Full text](https://smai-jcm.centre-mersenne.org/item/10.5802/smai-jcm.127.pdf).

Scheme 2.2 and Proposition 2.5 project the effort to obtain a balance at temporal
endpoints. Remark 2.4 addresses state-dependent operators.
The text distinguishes exact projection from numerical quadrature.
This relationship is explicit: temporal effort projection is not a new idea
introduced by this repository. Their framework includes nonlinear operators
beyond the constant-matrix implementation provided here.

### 6.3. Concrete contribution of this version

- A traceable reconstruction of the 2012 spatial method.
- An explicit connection between its quotients, ports and energy and a
  compatible time integrator.
- A space-time balance with proved quadrature conditions.
- Simulations distinguishing reproduction, accuracy, structure, energy,
  passivity and symplecticity.

A new publication would require an additional distinguishing result, such as
composition across different time grids, treatment of port-Hamiltonian
constraints, or uniform convergence in space and time.
None of these results is presented as already achieved.
