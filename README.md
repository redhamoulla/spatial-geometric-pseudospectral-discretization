# Spatial geometric pseudospectral discretization

**Méthode de Moulla–Lefèvre–Maschke (2012), extension en temps et simulations
port-hamiltoniennes reproductibles.**

Ce dépôt reconstruit la discrétisation spatiale géométrique de l'article
[JCP 2012](https://doi.org/10.1016/j.jcp.2011.10.008), puis la complète avec
une intégration temporelle à efforts projetés. Les modèles ouverts conservent
leur bilan discret :

$$
H_{n+1}-H_n
=\text{travail des ports}-\text{dissipation}.
$$

La construction utilise des espaces polynomiaux mixtes, un quotient des efforts
et des traces de bord. Elle conserve les Hamiltoniens polynomiaux non linéaires
si la quadrature est suffisante et la résolution implicite convergée.
Les ingrédients temporels sont rattachés à AVF, HBVM et Petrov–Galerkin :
leur nouveauté n'est pas revendiquée.

## Résultats exécutés

| Vérification | Résultat |
|---|---:|
| Huit fréquences du tableau 4 de 2012 | erreur max. $1.3\,10^{-12}$ |
| Ordres temporels, linéaire et non linéaire | environ 2 et 4 |
| Énergie quartique, degré temporel 2, 1 000 pas | erreur relative max. $2.5\,10^{-14}$ |
| Bilan avec deux ports actifs et résistance | défaut cumulé $1.3\,10^{-15}$ |
| Saint-Venant lisse non linéaire | défaut de bilan $4.3\,10^{-16}$ |
| Deux lignes interconnectées | défaut total $3.5\,10^{-16}$ |

**15 tests unitaires et 22 critères de simulation passent** sur l'exécution
de référence. Les seuils, données et limites figurent dans
[le rapport de validation](docs/05_validation.md).

![Énergie non linéaire](results/figures/nonlinear_energy.svg)

## Lire la théorie

1. [Article de 2012 : lecture, conventions et reproduction](docs/01_article_2012.md)
2. [Construction spatiale : Stokes, quotient, ports, invariance](docs/02_spatial.md)
3. [Extension temporelle : schéma, preuves et quadrature](docs/03_temporal.md)
4. [Espace-temps : complexe polynomial et interconnexion](docs/04_space_time.md)
5. [Protocole, résultats et limites](docs/05_validation.md)
6. [Positionnement dans la littérature](docs/06_literature.md)

Les démonstrations sont en Markdown avec des équations rendues sur GitHub.
Les graphiques sont en SVG et les données en CSV/JSON.

## Reproduire

Python 3.10 ou plus récent. Les calculs dépendent de NumPy ; ReportLab
produit les figures vectorielles.

~~~bash
git clone https://github.com/redhamoulla/spatial-geometric-pseudospectral-discretization.git
cd spatial-geometric-pseudospectral-discretization
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[figures]"
python -m unittest discover -s tests -v
python experiments/reproduce.py --check
~~~

Sans les figures :

~~~bash
python -m pip install -e .
python experiments/reproduce.py --check --no-figures --output /tmp/gpps-results
~~~

Chaque exécution génère tableaux, validations et manifeste des sources.
Les différences d'arrondi entre plateformes sont attendues ; les critères
utilisent des tolérances explicites.

## Organisation

| Répertoire | Contenu |
|---|---|
| docs/ | Démonstrations, lecture du papier, références et résultats commentés |
| src/gpps/ | Polynômes, réduction spatiale, Hamiltoniens, ports et intégrateur |
| experiments/ | Campagne complète et production des figures |
| tests/ | Géométrie, énergie, convergence, interconnexion et contrôles négatifs |
| results/tables/ | Données des simulations exécutées |
| results/figures/ | Six figures vectorielles |
| results/summary.json | Mesures et critères de validation |
| results/manifest.json | Environnement numérique et empreintes des sources |
| .github/workflows/ | Vérification automatique des tests et simulations |

## Utilisation minimale

~~~python
import numpy as np
from gpps import MixedDiscretization, transmission_line, ProjectedGauss, integrate

mesh = MixedDiscretization(8, interval=(0, 2), family="gauss")
system = transmission_line(mesh, capacitance=3, inductance=2)
x0 = np.r_[np.sin(np.pi * mesh.z / 4), np.zeros(mesh.n)]
path, ledger = integrate(system, x0, 0.01, 200, ProjectedGauss(2, 2))
print(max(abs(row["balance"]) for row in ledger))
~~~

## Périmètre

Version de recherche 0.1 : systèmes 1D lisses, structure matricielle constante,
énergies quadratique/cubique/quartique, ports et interconnexions synchrones.
Conservation énergétique ne signifie pas automatiquement symplecticité,
positivité de Saint-Venant, précision des modes élevés ou traitement des chocs.
Les PH à contraintes et les géométries 2D/3D restent à développer.

La méthode spatiale est attribuée aux trois auteurs de 2012. Les
[références temporelles](docs/06_literature.md) doivent également être citées ;
voir [CITATION.cff](CITATION.cff).
