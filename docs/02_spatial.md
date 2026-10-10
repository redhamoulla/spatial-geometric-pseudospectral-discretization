# Construction spatiale et preuves

## 1. Espaces, dérivée et puissance

Sur $[a,b]$, choisir deux ensembles de points distincts, sans imposer qu'ils
soient les mêmes ou emboîtés. Les bases de Lagrange associées sont

$$
\varphi_i\in\mathbb P_N,\quad i=0,\ldots,N,\qquad
\psi_k\in\mathbb P_{N-1},\quad k=0,\ldots,N-1.
$$

Les états $q,p$ et leurs dérivées temporelles utilisent $\psi$ ; les efforts
spatiaux utilisent $\varphi$. Définir

$$
D_{ki}=\varphi_i'(z_k),\qquad
M_{ik}=\int_a^b\varphi_i\psi_k\,dz,\qquad
t_a=(\varphi_i(a))_i,\quad t_b=(\varphi_i(b))_i.
$$

$D\in\mathbb R^{N\times(N+1)}$ et $M\in\mathbb R^{(N+1)\times N}$.
La dérivée d'un élément de $\mathbb P_N$ appartient exactement à
$\mathbb P_{N-1}$ ; aucune projection de la dérivée n'est nécessaire.

**Proposition S1 — Stokes discret.**

$$
MD+D^\top M^\top=t_b^\top t_b-t_a^\top t_a.
$$

**Preuve.** L'entrée $(i,j)$ du premier terme est
$\int\varphi_i\varphi_j'$ puisque l'interpolation de $\varphi_j'$ dans $\psi$
est exacte. Ajouter le terme transposé puis intégrer
$(\varphi_i\varphi_j)'$. L'identité est indépendante des points. $\square$

Pour le calcul, les polynômes sont exprimés dans une base de Legendre ; une
quadrature de Gauss à $N+1$ points intègre exactement $M$, les matrices de masse
et les produits nécessaires. Les nœuds de **quadrature** ne remplacent pas les
nœuds de **représentation** de l'utilisateur.

## 2. Le quotient n'est pas facultatif

La puissance interne est

$$
(e_q^{\rm full})^\top M\dot q+(e_p^{\rm full})^\top M\dot p.
$$

Elle est dégénérée sur les efforts complets. Les efforts conjugués aux
coefficients des états sont

$$
g_q=M^\top e_q^{\rm full},\qquad
g_p=M^\top e_p^{\rm full}.
$$

**Proposition S2 — Rang et complétion.** $M$ a rang $N$.
Le noyau de $M^\top$ est unidimensionnel et représente le polynôme de Legendre
de degré $N$ transporté sur $[a,b]$. Les matrices

$$
C_a=\begin{bmatrix}M^\top\\t_a\end{bmatrix},\qquad
C_b=\begin{bmatrix}M^\top\\t_b\end{bmatrix}
$$

sont inversibles.

**Preuve.** Si $Mv=0$, le polynôme $\psi v$ est orthogonal à tout
$\mathbb P_N$, donc à lui-même : il est nul. Le noyau transposé représente
donc l'orthogonal de $\mathbb P_{N-1}$ dans $\mathbb P_N$. Son générateur a
des valeurs de bord non nulles, proportionnelles à $(-1)^N$ et $1$.
Une trace élimine ce dernier degré de liberté. $\square$

Cette propriété de rang explique pourquoi deux efforts réduits et deux
entrées de bord suffisent à reconstruire les efforts complets.

## 3. Ports et graphe de Dirac

Choisir

$$
u=\begin{bmatrix}e_q(a)\\e_p(b)\end{bmatrix},\qquad
y=\begin{bmatrix}e_p(a)\\-e_q(b)\end{bmatrix}.
$$

À chaque évaluation :

$$
e_q^{\rm full}=C_a^{-1}\begin{bmatrix}g_q\\u_1\end{bmatrix},\qquad
e_p^{\rm full}=C_b^{-1}\begin{bmatrix}g_p\\u_2\end{bmatrix},
$$

$$
\dot q=-De_p^{\rm full},\qquad \dot p=-De_q^{\rm full}.
$$

Les sorties proviennent des mêmes traces. Cette élimination fournit

$$
\dot x=Jg+Bu,\qquad y=Cg+Su,\qquad x=(q,p).
$$

**Proposition S3 — Structure réduite.**

$$
J^\top=-J,\qquad C=B^\top,\qquad S^\top=-S.
$$

**Preuve.** S1 appliquée aux efforts reconstruits donne
$g^\top\dot x=u^\top y$ pour tout $(g,u)$. Par polarisation, la matrice

$$
\mathcal J=\begin{bmatrix}J&B\\-C&-S\end{bmatrix}
$$

est antisymétrique. Son graphe est isotrope et de dimension $2N+2$ dans
un espace de liaison de dimension $4N+4$ ; il est donc maximal isotrope.
Il définit une structure de Dirac pour le produit canonique. $\square$

Le terme $S$ n'est généralement pas nul. Le supprimer avant d'interconnecter
des sous-systèmes change les relations de bord. Le code conserve ce terme.
Il mesure l'antisymétrie brute ; aucune opération $(J-J^\top)/2$ ne vient
masquer une erreur d'assemblage.

## 4. Fermeture constitutive

Restreindre **l'énergie** à l'espace des états, puis la différencier :

$$
H_N(x)=\int_a^b\mathcal H(\psi q,\psi p,z)\,dz,\qquad g=\nabla H_N(x).
$$

Pour la ligne idéale :

$$
H_N=\tfrac12 q^\top Q_Cq+\tfrac12p^\top Q_Lp,\qquad
(Q_C)_{ij}=\int_a^b\frac{\psi_i\psi_j}{C(z)}\,dz.
$$

Pour Saint-Venant, dans un canal rectangulaire normalisé avec $\rho=B=1$,
les variables utilisées sont $(h,v)$, et non $(h,hv)$ :

$$
H_N(h,v)=\int_a^b\left(\tfrac g2 h_N^2+\tfrac12h_Nv_N^2\right)\,dz.
$$

Ainsi $e_h=gh+v^2/2$ et $e_v=hv$. L'énergie est cubique dans les coefficients.
La quadrature spatiale est choisie pour intégrer exactement ce polynôme,
ses gradients et sa Hessienne. Une énergie constitutive sous-intégrée est
une autre énergie discrète ; la conservation de celle-ci ne prouve pas la
précision vis-à-vis de l'énergie du continuum.

Une fermeture dissipative ajoute $-Rg$, avec $R=R^\top\succeq0$ :

$$
\dot x=(J-R)\nabla H_N+Bu,\qquad
\dot H_N=u^\top y-g^\top Rg.
$$

Pour la ligne à paramètres constants, la résistance du second champ est
$R=\operatorname{diag}(0,rM_\psi^{-1})$, où
$M_\psi=\int\psi^\top\psi$. Elle réalise le terme de courant $-r p/L_{\rm ind}$.

## 5. Invariance et conditionnement

Un changement de coordonnées d'état $x'=Tx$ transporte les objets selon

$$
J'=TJT^\top,\quad B'=TB,\quad S'=S,\quad
H'_N(x')=H_N(T^{-1}x').
$$

Pour une énergie quadratique,
$Q'=T^{-\top}QT^{-1}$. Les dynamiques sont similaires, les puissances et le
spectre sont identiques. Cette conclusion suppose la même fonctionnelle
constitutive, ou des quadratures qui donnent la même restriction.

Cette invariance est algébrique, pas une immunité à l'arrondi. Les nœuds
uniformes et perturbés se conditionnent mal lorsque $N$ augmente.
Les résultats publient les conditionnements et les résidus correspondants.
Gauss et Lobatto sont les choix de calcul usuels du dépôt.

Suite : [discrétisation temporelle](03_temporal.md).
