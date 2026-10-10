# Positionnement bibliographique

Consultation effectuée le 10 octobre 2026. Cette liste est ciblée sur la
construction du dépôt ; elle ne prétend pas établir une priorité par une
revue exhaustive de toute la littérature.

## Source spatiale

**Moulla, Lefèvre, Maschke (2012).** *Pseudo-spectral methods for the spatial
symplectic reduction of open systems of conservation laws*. JCP 231(4),
1272–1292. [DOI](https://doi.org/10.1016/j.jcp.2011.10.008).

Source directe : espaces mixtes, produit rectangulaire, quotient et ports.
Le manuscrit HAL complet a été lu. L'audit et la correspondance avec le code
sont dans [01_article_2012.md](01_article_2012.md).

## Intégration temporelle

**Hairer (2010).** *Energy-preserving variant of collocation methods*.
[Texte de l'auteur](https://www.unige.ch/~hairer/preprints/crete.pdf).

Développe des variantes de collocation prolongeant AVF à ordre élevé.
Ce travail distingue conservation de l'énergie et propriétés symplectiques.
Le dépôt reprend cette distinction dans ses contrôles négatifs.

**Brugnano, Iavernaro, Trigiante.** *Analysis of Hamiltonian Boundary Value
Methods (HBVMs): A class of energy-preserving Runge–Kutta methods for the
numerical solution of polynomial Hamiltonian systems*.
[Prépublication 2009](https://arxiv.org/abs/0909.5659) ·
[Publication 2015](https://doi.org/10.1016/j.cnsns.2014.05.030).

Les cas conservatifs à matrice constante de notre intégrateur relèvent des
HBVM$(k,s)$. Le degré du chemin et le nombre de points de quadrature jouent
des rôles différents : augmenter la quadrature peut rendre le bilan exact
sans augmenter le nombre d'inconnues polynomiales.

**Celledoni, Høiseth (2017, prépublication).** *Energy-Preserving and
Passivity-Consistent Numerical Discretization of Port-Hamiltonian Systems*.
[arXiv](https://arxiv.org/abs/1706.08621).

Présente des discrétisations par gradients discrets et par séparation
d'opérateurs, avec passivité discrète. C'est un antécédent direct à toute
revendication générale d'intégration PH conservant l'énergie.

**Kotyczka, Lefèvre (2018, prépublication).** *Discrete-time port-Hamiltonian
systems: A definition based on symplectic integration*.
[arXiv](https://arxiv.org/abs/1811.07852).

Construit une structure de Dirac en temps discret par collocation et
intégration symplectique ; distingue bilan discret et erreur énergétique.
Ce cadre ne justifie pas d'appeler « exact » le bilan de l'énergie originale
d'un Hamiltonien non quadratique sous une quadrature insuffisante.

**Giesselmann, Karsai, Tscherpel (2025).** *Energy-consistent Petrov–Galerkin
time discretization of port-Hamiltonian systems*. SMAI-JCM 11, 335–367.
[Article](https://smai-jcm.centre-mersenne.org/articles/10.5802/smai-jcm.127/) ·
[Texte intégral](https://smai-jcm.centre-mersenne.org/item/10.5802/smai-jcm.127.pdf).

Le schéma 2.2 et la proposition 2.5 projettent l'effort pour obtenir un bilan
aux extrémités temporelles. La remarque 2.4 traite des opérateurs dépendant
de l'état. Le texte distingue projection exacte et quadrature numérique.
Cette proximité est explicite : la projection temporelle n'est pas une idée
nouvelle du présent dépôt. Leur cadre inclut des opérateurs non linéaires
au-delà de l'implémentation matricielle constante livrée ici.

## Contribution concrète de cette version

- Reconstruction traçable de la méthode spatiale de 2012.
- Raccordement explicite de ses quotients, ports et énergie à un intégrateur
  temporel compatible.
- Bilan espace-temps avec conditions de quadrature démontrées.
- Simulations distinguant reproduction, précision, structure, énergie,
  passivité et symplecticité.

Une publication nouvelle demanderait un résultat supplémentaire différenciant :
composition avec grilles temporelles différentes, traitement de contraintes
PH, ou convergence uniforme en espace et en temps. Aucun de ces résultats
n'est annoncé comme déjà obtenu.
