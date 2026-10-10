# Extension temporelle : efforts projetés et bilan exact

## 1. Problème et statut de la construction

Après réduction spatiale :

$$
\dot x=(J-R)\nabla H_N(x)+Bu,\qquad
y=B^\top\nabla H_N(x)+Su,
$$

avec $J^\top=-J$, $R^\top=R\succeq0$, $S^\top=-S$.
Le code traite ces matrices constantes. $H_N$ peut être non linéaire et
$u(t)$ variable. La généralisation à des matrices dépendant de l'état est
énoncée en §7, sans prétendre qu'elle est déjà implémentée.

Le choix temporel est une méthode à **état polynomial continu et effort
projeté**. Pour une structure constante conservative, elle appartient à la
famille AVF/HBVM. Son interprétation Petrov–Galerkin et son extension aux ports
ont des antécédents explicites, présentés dans la
[bibliographie commentée](06_literature.md). L'apport de ce dépôt est la
construction complète, vérifiable et reproductible avec la méthode spatiale
de 2012 ; aucune priorité sur ces intégrateurs temporels n'est revendiquée.

## 2. Deux espaces temporels

Sur un pas $[t_n,t_n+h]$, écrire $c=(t-t_n)/h$ et choisir

$$
X(c)\in\mathbb P_s([0,1];\mathbb R^d),\quad X(0)=x_n,\qquad
V(c)=\frac1hX'(c)\in\mathbb P_{s-1}.
$$

Soient $P_0,\ldots,P_{s-1}$ les polynômes de Legendre décalés orthonormaux,
$P_0=1$, et $I_j(c)=\int_0^cP_j(\tau)\,d\tau$. Paramétrer

$$
X(c)=x_n+h\sum_{j=0}^{s-1}I_j(c)a_j,\qquad
V(c)=\sum_{j=0}^{s-1}P_j(c)a_j.
$$

La trace initiale élimine le noyau constant de la dérivée. Les inconnues
sont les $sd$ coefficients $a_j$, et
$x_{n+1}=X(1)=x_n+ha_0$.

Définir l'effort temporel par projection orthogonale :

$$
\gamma_j=\int_0^1P_j(c)\nabla H_N(X(c))\,dc,\qquad
E(c)=\sum_{j=0}^{s-1}P_j(c)\gamma_j.
$$

$E$ n'est pas en général $\nabla H_N(X)$ point par point. Il est l'effort
**conjugué au flux discret** $V$ dans le produit de puissance temporel.
Cette différence est nécessaire pour une énergie non quadratique.

## 3. Schéma calculable

Une quadrature positive $(c_\ell,b_\ell)_{\ell=1}^k$, exacte pour
$\mathbb P_{2s-2}$, fournit la dynamique faible

$$
a_i=\sum_{\ell=1}^k b_\ell P_i(c_\ell)
\big[(J-R)E(c_\ell)+Bu(t_n+hc_\ell)\big].
\tag{T1}
$$

Les sorties utilisées dans le travail sont

$$
Y_\ell=B^\top E(c_\ell)+Su_\ell.
\tag{T2}
$$

Dans le code, la même quadrature de Gauss calcule aussi $\gamma_i$ :

$$
\gamma_i^Q=\sum_\ell b_\ell P_i(c_\ell)\nabla H_N(X(c_\ell)).
\tag{T3}
$$

Il faut distinguer l'identité théorique avec projection exacte et
l'implémentation avec projection quadraturée.

**Proposition T1 — Bilan d'énergie d'un pas.** Supposons que les coefficients
$\gamma_i$ soient calculés exactement et que T1 soit résolu exactement. Alors

$$
\boxed{
H_N(x_{n+1})-H_N(x_n)
=h\sum_\ell b_\ell
\left[u_\ell^\top Y_\ell-E_\ell^\top R E_\ell\right].
}
\tag{T4}
$$

**Preuve.** La règle de chaîne continue appliquée au chemin polynomial donne

$$
\Delta H_N=\int_0^1\nabla H_N(X(c))^\top X'(c)\,dc
=h\sum_i\gamma_i^\top a_i.
$$

Insérer T1 et intervertir les sommes :

$$
\Delta H_N=h\sum_\ell b_\ell E_\ell^\top
\big[(J-R)E_\ell+Bu_\ell\big].
$$

Les termes $E_\ell^\top JE_\ell$ et $u_\ell^\top Su_\ell$ sont nuls.
Utiliser T2. $\square$

Sans entrée, l'énergie ne croît pas ; sans entrée ni dissipation, elle est
conservée aux **extrémités** des pas. Il n'est pas affirmé qu'elle reste
constante à tout point du polynôme intérieur. Pour interpréter le bilan comme
une inégalité de passivité avec stockage positif, $H_N$ doit être borné
inférieurement sur le domaine pertinent.

## 4. Règle de quadrature et cas particuliers

Si $H_N$ est un polynôme de degré total $r$, alors
$P_i(c)\nabla H_N(X(c))$ a degré au plus $rs-1$.
La quadrature de Gauss à $k$ points est donc suffisante si

$$
k\ge s,\qquad 2k\ge rs.
\tag{T5}
$$

| Énergie | Degré $r$ | État temporel $s=1$ | État temporel $s=2$ |
|---|---:|---:|---:|
| Ligne idéale | 2 | $k=1$ | $k=2$ |
| Saint-Venant rectangulaire | 3 | $k=2$ | $k=3$ |
| Ligne avec stockage quartique | 4 | $k=2$ | $k=4$ |

Avec $s=1$, l'effort est la moyenne du gradient sur le segment joignant les
deux états. On obtient AVF. Pour un Hamiltonien quadratique, cette moyenne
est le gradient au milieu ; on retrouve le point milieu implicite.

Avec $k=s$, les valeurs de l'effort projeté aux nœuds de Gauss coïncident
avec celles du gradient : le schéma est la collocation de Gauss classique.
Avec $k>s$, on obtient la réalisation HBVM$(k,s)$ pour une structure constante
conservative. Pour les problèmes autonomes lisses de cette classe, l'ordre
aux extrémités est $2s$. Les simulations mesurent effectivement les ordres 2 et 4,
avec des références de trajectoire indépendantes.

La preuve de T4 n'est pas une preuve générale d'ordre $2s$ pour toutes les
extensions non autonomes, à matrices variables ou avec contraintes. Les
affirmations d'ordre du dépôt se limitent aux classes citées et aux cas mesurés.

Pour un Hamiltonien non polynomial, aucune valeur finie de $k$ ne garantit
en général l'exactitude de T3. Il faut augmenter la quadrature ou adopter
une fermeture par gradient discret exact. Le défaut de chaîne est alors
une quantité à contrôler, pas à mettre artificiellement à zéro.

## 5. Certificat numérique : trois sources d'erreur

Soit $\rho_i$ le résidu de T1, avec le signe
$\rho_i=a_i-\mathrm{RHS}_i$, et définissons

$$
\varepsilon_{\rm ch}=\Delta H_N-h\sum_i(\gamma_i^Q)^\top a_i.
$$

À structure antisymétrique exacte, on a l'identité a posteriori

$$
\Delta H_N-W_n+D_n
=\varepsilon_{\rm ch}+h\sum_i(\gamma_i^Q)^\top\rho_i.
\tag{T6}
$$

En arithmétique flottante s'ajoute le défaut de la structure spatiale.
Le programme calcule séparément l'énergie aux extrémités, le travail des ports,
la dissipation, le défaut de chaîne et le résidu de Newton. Il n'obtient pas
l'énergie par cumul du bilan qu'il est censé vérifier.

Le solveur est Newton avec Hessienne analytique et recherche linéaire. Un
échec de convergence lève une exception. Une énergie très bien conservée
ne garantit pas une bonne phase ou une bonne trajectoire : les tableaux de
convergence et la comparaison au continuum sont donc indispensables.

Si $H_N$ est $C^2$ et l'entrée suffisamment régulière, la fonction résiduelle
en $a$ a une dérivée égale à l'identité pour $h=0$. Le théorème des fonctions
implicites donne une branche locale unique pour un pas assez petit. Ceci
n'établit ni l'existence globale de la solution implicite pour tout pas,
ni la convergence de Newton pour toute initialisation.

## 6. Symplecticité : propriété distincte

Pour un système fermé avec $J$ inversible, une application $\Phi_h$ est
symplectique si

$$
(D\Phi_h)^\top J^{-1}D\Phi_h=J^{-1}.
$$

Les schémas de Gauss satisfont cette propriété, mais ne conservent pas
généralement un Hamiltonien non quadratique. Les variantes AVF/HBVM retenues
conservent l'énergie polynomiale avec la quadrature adéquate, mais ne sont
pas généralement symplectiques. Le dépôt mesure ces deux défauts sur un
oscillateur quartique avec une matrice tangente obtenue par différentiation
implicite, elle-même vérifiée par différences finies.

La conservation énergétique ne suffit pas non plus à établir la préservation
de tous les Casimirs, la positivité de la hauteur d'eau, l'absence d'oscillations,
ou la stabilité d'un système non coercif.

## 7. Extension algébrique aux matrices variables

Remplacer, dans T1, $J,R,B,S$ par leurs valeurs à $X(c_\ell)$.
Si $J(X)^\top=-J(X)$, $R(X)\succeq0$ et $S(X)^\top=-S(X)$ à chaque étage,
la même preuve donne T4 avec les matrices d'étage : elle n'a utilisé
l'antisymétrie et la positivité que point par point.

Cela exige toujours la projection de l'effort dans **tous** les termes de
puissance et une quadrature exacte de cette projection lorsque l'on revendique
un bilan exact. La régularité, l'ordre et la résolution des équations résultantes
restent à analyser pour le modèle considéré. Le solveur actuel, dont la
Jacobienne exploite des matrices constantes, ne met pas cette extension en œuvre.

Suite : [assemblage espace-temps et interconnexion](04_space_time.md).
