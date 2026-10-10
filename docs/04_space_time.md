# Assemblage espace-temps et interconnexion

## 1. Pourquoi un produit tensoriel ne suffit pas

La dérivation exacte et une identité de Stokes discrète règlent la géométrie
**bilinéaire**. Elles ne créent pas automatiquement la règle de chaîne d'un
Hamiltonien non quadratique. Une construction complète comporte donc :

1. un complexe polynomial compatible et ses traces ;
2. un quotient des efforts donnant un produit de puissance non dégénéré ;
3. une fermeture constitutive issue de l'énergie spatiale ;
4. une fermeture temporelle par effort projeté ;
5. des ports évalués avec ce même effort.

L'intégrateur réalise les points 3–5 sur chaque tranche temporelle.
Les identités tensorielles ci-dessous vérifient les points 1–2 au niveau
géométrique. Les tests de tensorisation sont distincts des simulations dynamiques.

## 2. Complexe polynomial du rectangle

Notons $E_t=\mathbb P_s$, $F_t=\mathbb P_{s-1}\,dt$,
$E_z=\mathbb P_N$, $F_z=\mathbb P_{N-1}\,dz$. Sur une tranche :

$$
V^0=E_t\otimes E_z,\qquad
V^1=(F_t\otimes E_z)\oplus(E_t\otimes F_z),\qquad
V^2=F_t\otimes F_z.
$$

Avec les orientations $dt,dz$, les matrices du complexe sont

$$
d_0=
\begin{bmatrix}
D_t\otimes I_{E_z}\\ I_{E_t}\otimes D_z
\end{bmatrix},\qquad
d_1=
\begin{bmatrix}
-I_{F_t}\otimes D_z&D_t\otimes I_{F_z}
\end{bmatrix}.
$$

**Proposition ST1.** $d_1d_0=0$.

**Preuve.** Les deux termes du produit sont
$-D_t\otimes D_z$ et $D_t\otimes D_z$. Ils s'annulent exactement.
La propriété n'utilise pas une famille particulière de nœuds. $\square$

## 3. Identité de bord sur la tranche

Pour $\alpha\in\{t,z\}$, poser

$$
S_\alpha=M_\alpha D_\alpha,\qquad
G_\alpha=\int\varphi_\alpha^\top\varphi_\alpha,\qquad
B_\alpha=t_{\alpha,+}^\top t_{\alpha,+}
-t_{\alpha,-}^\top t_{\alpha,-}.
$$

Le symbole $S_\alpha$ désigne ici une matrice d'intégration par parties,
sans rapport avec le terme de transmission directe $S$ des ports.
Les matrices unidimensionnelles satisfont
$S_\alpha+S_\alpha^\top=B_\alpha$. Il en découle

$$
\mathcal K=S_t\otimes G_z+G_t\otimes S_z,
$$

$$
\boxed{\mathcal K+\mathcal K^\top
=B_t\otimes G_z+G_t\otimes B_z.}
$$

**Preuve.** Transposer les produits de Kronecker et utiliser la symétrie
de $G_t,G_z$. $\square$

Le membre droit contient les faces initiale et finale ainsi que les deux
frontières spatiales. Les faces internes de deux tranches accolées
s'annulent si elles utilisent les mêmes traces avec des orientations opposées.

Il s'agit d'une identité exacte sur les espaces polynomiaux, pas d'une preuve
de conservation multisymplectique pour toutes les fermetures non linéaires.
Une loi multisymplectique requerrait l'étude des variations et des deux-formes
associées à la dynamique discrète. Elle n'est pas revendiquée ici.

## 4. Théorème de bilan pour la simulation complète

**Proposition ST2.** Sous les hypothèses des propositions S1–S3 et T1, pour
des pas $h_n>0$, la simulation satisfait

$$
H_N(x_M)-H_N(x_0)
=\sum_{n=0}^{M-1}h_n\sum_\ell b_\ell
\left[u_{n\ell}^\top Y_{n\ell}
-E_{n\ell}^\top R E_{n\ell}\right].
$$

**Preuve.** La réduction spatiale donne les matrices de S3 et la fermeture
$g=\nabla H_N$. Appliquer T4 sur chaque tranche et sommer : les énergies
aux extrémités internes se télescopent. $\square$

Le bilan reste valable avec des pas variables, sous les mêmes conditions
d'intégration et de résolution à chaque pas. Il concerne $H_N$, donc
l'énergie du modèle spatial discret ; l'erreur vis-à-vis du continuum est
évaluée séparément.

## 5. Composition de systèmes ouverts

Assembler plusieurs systèmes par blocs diagonaux, puis connecter les ports
avec

$$
u=Ky+v,\qquad K^\top=-K.
$$

Au niveau de chaque étage, $u^\top y=v^\top y$ : l'interconnexion n'injecte
aucune puissance. Comme $y=B^\top E+Su$,

$$
(I-KS)u=KB^\top E+v.
$$

Si $I-KS$ est inversible, poser $F=(I-KS)^{-1}$. Le système composé est

$$
J_{\rm c}=J+BFKB^\top,\qquad
B_{\rm c}=BF,\qquad S_{\rm c}=SF.
$$

On conserve $R$ et on additionne les Hamiltoniens.
Les identités $FK=KF^\top$ et $SF=F^\top S$ donnent
$J_{\rm c}^\top=-J_{\rm c}$, $S_{\rm c}^\top=-S_{\rm c}$ ; la sortie
s'exprime comme $B_{\rm c}^\top E+S_{\rm c}v$.

**Proposition ST3.** La somme des bilans T4 des composants est le bilan
du système composé si les composants partagent les mêmes étages temporels
et la même quadrature de puissance.

**Preuve.** Les puissances internes sont annulées par $y^\top Ky=0$
à chaque étage, avant la sommation en temps. $\square$

L'expérience connecte la frontière droite de la ligne A à la frontière gauche
de la ligne B : continuité du courant et de la tension, avec les orientations
définies dans la note spatiale. Les deux autres ports ont des entrées nulles.
L'énergie passe réellement d'une ligne à l'autre.

Si $I-KS$ est singulière, il ne faut pas remplacer l'inverse par une
pseudo-inverse arbitraire : des contraintes peuvent apparaître et requérir
une formulation PH différentielle-algébrique. Le programme signale ce cas.
Les interconnexions avec des horloges différentes ne sont pas prises en charge.

## 6. Ce qui reste un travail de recherche

| Sujet | Statut de cette version |
|---|---|
| Géométrie spatiale 1D à points arbitraires | Démonstration et code |
| Bilan d'énergie temporel non linéaire | Démonstration, code et simulations |
| Complexe et Stokes tensoriels | Identités démontrées et vérifiées |
| Interconnexion synchrone avec ports | Démonstration et simulation |
| PH à matrices variables | Bilan algébrique démontré, solveur non implémenté |
| PH avec contraintes algébriques | Non résolu dans cette version |
| Géométries spatiales 2D/3D | Non implémentées |
| Chocs et positivité garantie de Saint-Venant | Non établis |
| Multisymplecticité du schéma complet | Non revendiquée |
| Nouveauté mathématique publiable | Non établie par cette seule construction |

Le prolongement le plus concret consiste à préserver la composition énergétique
avec des pas différents entre sous-domaines : la cancellation de puissance
doit alors être transportée entre espaces temporels distincts, plutôt que
supposée à des étages communs.

Suite : [protocole et résultats](05_validation.md).
