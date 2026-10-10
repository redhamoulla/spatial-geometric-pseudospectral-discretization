# Spatial geometric pseudospectral discretization

Théorie et simulations reproductibles de la méthode spatiale de
[Moulla–Lefèvre–Maschke (2012)](https://doi.org/10.1016/j.jcp.2011.10.008),
avec intégration temporelle à efforts projetés.

**[Lire le manuscrit complet](docs/manuscrit.md)** : étude de l'article,
construction spatiale, preuves, discrétisation temporelle, interconnexion,
simulations et bibliographie, dans un seul document avec sommaire.
Les équations sont rendues sur GitHub et les figures sont en SVG.

## Reproduire les simulations

Python 3.10 ou plus récent. NumPy assure les calculs ; ReportLab produit
les figures vectorielles.

~~~bash
git clone https://github.com/redhamoulla/spatial-geometric-pseudospectral-discretization.git
cd spatial-geometric-pseudospectral-discretization
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e ".[figures]"
python -m unittest discover -s tests -v
python experiments/reproduce.py --check
~~~

Sans générer les figures :

~~~bash
python -m pip install -e .
python experiments/reproduce.py --check --no-figures --output /tmp/gpps-results
~~~

L'exécution de référence passe **15 tests unitaires et 22 critères de
simulation**. Les huit fréquences publiées en 2012 sont reproduites avec un
écart maximal de $1.3\,10^{-12}$ ; les ordres temporels mesurés sont 2 et 4.
Le [chapitre des résultats](docs/manuscrit.md#5-protocole-et-résultats-numériques)
présente les paramètres, figures, erreurs et limites de chaque expérience.
Les seuils tolèrent les différences d'arrondi entre plateformes.

## Organisation

| Emplacement | Contenu |
|---|---|
| [docs/manuscrit.md](docs/manuscrit.md) | Manuscrit scientifique complet |
| src/gpps/ | Polynômes, réduction spatiale, Hamiltoniens, ports et intégrateur |
| experiments/ | Campagne de simulations et production des figures |
| tests/ | Vérification de la géométrie, de l'énergie et de la convergence |
| results/tables/ | Données numériques en CSV |
| results/figures/ | Six figures vectorielles en SVG |
| results/summary.json | Mesures et critères de validation |
| results/manifest.json | Environnement numérique et empreintes des sources |
| .github/workflows/ | Vérification automatique des tests et simulations |

## Périmètre et citation

Version de recherche 0.1 : systèmes 1D lisses, matrices de structure constantes,
énergies quadratique, cubique ou quartique et interconnexions synchrones.
La conservation énergétique ne garantit pas la symplecticité, la positivité
de Saint-Venant ou le traitement des chocs.

La méthode spatiale est attribuée aux auteurs de 2012. Les méthodes temporelles
AVF, HBVM et Petrov–Galerkin ont leurs propres antécédents, détaillés dans la
[bibliographie](docs/manuscrit.md#6-bibliographie-et-positionnement).
Voir [CITATION.cff](CITATION.cff) et [CONTRIBUTING.md](CONTRIBUTING.md).
