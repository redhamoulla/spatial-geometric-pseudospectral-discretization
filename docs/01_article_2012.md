# Lecture critique et traçabilité de l'article de 2012

## Source effectivement étudiée

R. Moulla, L. Lefèvre et B. Maschke, **Pseudo-spectral methods for the spatial
symplectic reduction of open systems of conservation laws**, *Journal of
Computational Physics* 231(4), 1272–1292 (2012).
[DOI](https://doi.org/10.1016/j.jcp.2011.10.008) ·
[Manuscrit HAL](https://hal.science/hal-01625008).

La lecture porte sur la copie du manuscrit HAL fournie par l'auteur :
38 pages avec couverture HAL, texte daté du 16 décembre 2010, dépôt HAL de 2017.
Les numéros d'équations, propositions et tableaux ci-dessous sont ceux de cette
copie. La version mise en page par l'éditeur n'a pas été auditée séparément.
Le dépôt contient une reconstruction et des commentaires originaux ; le fichier
source de l'article n'est pas redistribué.

## Ce que construit l'article

| Étape | Localisation dans le manuscrit | Réalisation dans ce dépôt |
|---|---|---|
| Espaces de formes et système ouvert | §2.3, définitions 11–14 | Convention de signe et ports explicites |
| Efforts dans $\mathbb P_N$, flux et états dans $\mathbb P_{N-1}$ | (31)–(36) | LagrangeBasis, MixedDiscretization |
| Dérivation exacte et traces de bord | (36)–(38) | D, trace |
| Produit de puissance rectangulaire | (39)–(40) | M |
| Isotropie avant réduction | Proposition 18 | Test de Stokes |
| Quotient des efforts par $\ker M^\top$ | (44), proposition 20 | Efforts réduits et complétions par traces |
| Réalisation avec entrées et sorties | Proposition 21, (50) | J, B, C, S assemblés directement |
| Fermeture constitutive | §4.1 et §4.2 | Énergies quadratique, quartique, Saint-Venant |
| Spectre de la ligne idéale | §6, tableaux 1, 3, 4 | paper_spectrum.csv et test indépendant |

La propriété centrale concerne les **espaces**, le produit de puissance et le
quotient, et non une famille obligatoire de points de Gauss. La géométrie doit
survivre à un changement de coordonnées polynomiales inversible. Le programme
vérifie Gauss, Lobatto, points uniformes intérieurs et points perturbés. Les
points de Lobatto aux frontières sont une extension naturelle de l'argument
polynomial ; les nœuds sont intérieurs dans les définitions initiales du manuscrit.

Le terme « symplecticité spatiale » de l'article désigne la conservation de cette
structure après réduction. Il ne suffit pas, à lui seul, à établir la
symplecticité de l'application temporelle $x_n\mapsto x_{n+1}$.

## Conventions vérifiées plutôt que recopiées

Nous fixons

$$
\partial_t q=-\partial_z e_p,\qquad
\partial_t p=-\partial_z e_q,\qquad
\dot H=e_q(a)e_p(a)-e_q(b)e_p(b).
$$

Avec la définition $D_{ki}=\varphi_i'(z_k)$, l'intégration par parties donne
nécessairement

$$
MD+D^\top M^\top=t_b^\top t_b-t_a^\top t_a.
$$

L'équation (43) de la copie HAL affiche le signe opposé des termes de bord
lorsqu'on la combine littéralement avec cette définition de $D$.
Le manuscrit alterne aussi certaines conventions de flux et de dérivée d'état.
Le code adopte la convention ci-dessus de bout en bout ; il ne transplante pas
isolément un signe d'une formule. Les tests évaluent indépendamment les traces,
la puissance et l'intégration par parties.

Autres précautions de lecture :

- Les facteurs de bord de (67) sont inversés par rapport aux conditions
  annoncées ; (70) donne les facteurs compatibles.
- La dernière implication de la preuve du théorème 24 doit conserver
  soigneusement la métrique de quadrature. Une antisymétrie pondérée ne donne
  pas en général une antisymétrie euclidienne par simple suppression des poids.
- La remarque 25 est explicitement une conjecture pour la collocation classique
  de Chebyshev. Le dépôt ne transforme pas cette conjecture en théorème.
- L'argument de rang utilisé ici s'appuie sur l'injectivité de la dérivée
  accompagnée d'une trace, puis sur la dimension du graphe réduit.

Ces observations concernent la copie étudiée et les conventions d'implémentation ;
elles ne remettent pas en cause la construction géométrique centrale.

## Reproduction numérique précise

Les valeurs publiées correspondent à une longueur $\ell=2$, une capacité
$C=3$, une inductance linéique $L_{\rm ind}=2$, huit états par champ et des
entrées homogènes. La longueur est ici identifiée par les fréquences du tableau 1 :

$$
\omega_j=\frac{(j+\frac12)\pi}{\ell\sqrt{C L_{\rm ind}}},\qquad j=0,\ldots,7.
$$

Le test compare les **huit** fréquences calculées aux décimales du tableau 4,
y compris les modes élevés mal approximés par un modèle d'ordre 16.
Reproduire le tableau ne signifie pas que ses huit modes approchent tous bien
le continuum : le CSV conserve séparément l'erreur de reproduction et l'erreur
par rapport aux fréquences analytiques.

Le dépôt ne prétend pas reproduire les variantes de collocation classique
mal conditionnées des tableaux 2 et 5. Les nouveaux cas non linéaires, ouverts et
interconnectés sont identifiés comme des expériences supplémentaires.

Suite : [construction spatiale](02_spatial.md).
