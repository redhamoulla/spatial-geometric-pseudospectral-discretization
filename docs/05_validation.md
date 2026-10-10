# Protocole et résultats numériques

Résultats exécutés le 10 octobre 2026 par
[experiments/reproduce.py](../experiments/reproduce.py). Le
[manifeste](../results/manifest.json) conserve les versions de Python et NumPy,
la plateforme et les empreintes SHA-256 des sources de calcul.

## 1. Reproduction de l'article

Ligne idéale : $\ell=2$, $C=3$, $L_{\rm ind}=2$, $N=8$, ports homogènes.
Les huit fréquences du tableau 4 sont retrouvées pour quatre familles de
points. L'écart maximal aux décimales publiées est $1.26\,10^{-12}$.

Les transformations entre les coordonnées de Gauss et les coordonnées
perturbées donnent un défaut relatif de $2.07\,10^{-13}$ pour $J$ et de
$2.47\,10^{-15}$ pour $Q$. Les quadratures constitutives sont exactes.

[Fréquences et erreurs par rapport au continuum](../results/tables/paper_spectrum.csv).

## 2. Géométrie et conditionnement

| Vérification | Résultat |
|---|---:|
| Défaut relatif de Dirac maximal, quatre familles, $N\le12$ | $4.59\,10^{-11}$ |
| Même mesure à $N=16$ | $1.91\,10^{-8}$ |
| Identité de Stokes tensorielle | $1.02\,10^{-14}$ |
| Identité $d_1d_0=0$, résidu normalisé | $7.78\,10^{-18}$ |

La dégradation à $N=16$ vient des coordonnées uniformes ou perturbées :
le conditionnement d'une complétion atteint environ $4.1\,10^8$.
Cela ne contredit pas l'invariance en arithmétique exacte. Le seuil
$10^{-8}$ de la validation structurelle est appliqué explicitement à
$N\le12$ ; les résultats à $N=16$ restent publiés.
Les simulations temporelles utilisent les coordonnées de Gauss.

[Structure et conditionnements](../results/tables/spatial_structure.csv) ·
[Complexe espace-temps](../results/tables/space_time_geometry.csv).

## 3. Convergence de trajectoire

La convergence temporelle linéaire est mesurée à espace fixé contre
l'exponentielle du modèle semi-discret, calculée indépendamment par
diagonalisation dans les coordonnées énergétiques.

La convergence temporelle non linéaire utilise un RK4 explicite fin, indépendant
du solveur implicite. Le passage de 6 000 à 12 000 pas change l'état de
référence de $2.86\,10^{-14}$ en norme euclidienne.

| Modèle | Degré temporel | Ordre mesuré sur les deux derniers pas |
|---|---:|---:|
| Ligne idéale | 1 | 1.999 |
| Ligne idéale | 2 | 3.986 |
| Ligne quartique | 1 | 2.000 |
| Ligne quartique | 2 | 3.987 |

![Convergence temporelle](../results/figures/time_convergence.svg)

Une onde stationnaire analytique mesure séparément l'erreur spatiale
et l'erreur espace-temps. À $N=12$, l'erreur spatiale $L^2$ est
$1.46\,10^{-8}$ ; avec le degré temporel 2 et $h=0.025$,
l'erreur totale est également de l'ordre de $1.46\,10^{-8}$.
Il s'agit d'un cas analytique lisse, pas d'une estimation uniforme pour
toute équation PH.

![Convergence spatiale](../results/figures/space_convergence.svg)

[Temps](../results/tables/time_convergence.csv) ·
[Non linéaire](../results/tables/nonlinear_convergence.csv) ·
[Espace](../results/tables/space_convergence.csv).

## 4. Hamiltonien quartique : contrôles positifs et négatifs

Modèle fermé normalisé :

$$
H_N=\int_0^1\left[\tfrac12(q_N^2+p_N^2)+\tfrac14q_N^4\right]dz.
$$

$N=6$, $h=0.03$, 1 000 pas, $T=30$.
Cette expérience évalue l'intégration du modèle semi-discret ; elle ne
constitue pas une étude de convergence PDE à cet horizon.

| Méthode | $(s,k)$ | Maximum de l'erreur relative d'énergie |
|---|---:|---:|
| Point milieu | (1,1) | $2.71\,10^{-3}$ |
| AVF quadraturé exactement | (1,2) | $9.23\,10^{-15}$ |
| Gauss à deux étages | (2,2) | $4.92\,10^{-6}$ |
| Effort projeté, degré 2 | (2,4) | $2.42\,10^{-14}$ |

![Énergie non linéaire](../results/figures/nonlinear_energy.svg)

L'oscillateur quartique distingue énergie et symplecticité. Sur un pas de
taille $0.2$, le point milieu a un défaut symplectique $2.86\,10^{-17}$ mais
une erreur d'énergie $2.88\,10^{-4}$. AVF conserve l'énergie à l'arrondi près
mais a un défaut symplectique $2.58\,10^{-3}$.

[Mesures sur l'oscillateur](../results/tables/symplectic_vs_energy.csv).

## 5. Ports actifs et dissipation

Même stockage quartique, résistance $r=0.15$ sur le second champ,
deux entrées temporelles non nulles, $h=0.02$, $T=10$, méthode $(s,k)=(2,4)$.
L'énergie est calculée directement à partir de l'état ; le travail et la
dissipation sont intégrés indépendamment aux étages.

$$
\max_n|H_n-H_0-W_{0:n}+D_{0:n}|=1.23\,10^{-15}.
$$

La dissipation totale vaut environ $0.10906$ et tous ses incréments sont
positifs. Le défaut maximal du bilan sur un pas est $5.42\,10^{-17}$.

![Bilan ouvert](../results/figures/forced_balance.svg)

[Bilan pas à pas](../results/tables/forced_damped.csv).

## 6. Saint-Venant non linéaire

Canal rectangulaire normalisé, $g=1$, $N=12$ :

$$
h(z,0)=1+0.02\sin(\pi z/2),\qquad v(z,0)=0.
$$

L'entrée gauche impose l'effort de Bernoulli $e_h=1$ ; la droite impose
le débit $e_v=0$. Le bilan tient compte du travail effectif de la frontière
gauche. Méthode $(s,k)=(2,3)$, $h_t=0.002$, $T=0.6$.

- Défaut cumulé du bilan : $4.24\,10^{-16}$.
- Résidu du champ de vecteurs à l'eau au repos : $8.15\,10^{-13}$.
- Hauteur minimale sur 401 points à chaque extrémité temporelle : $0.999866$.
- Écart RMS final des hauteurs avec $N=18$ et $h_t=0.001$ : $3.30\,10^{-7}$.

![Saint-Venant](../results/figures/shallow_water.svg)

La hauteur positive **échantillonnée** n'est pas une preuve de positivité
sur tout l'intervalle ni entre les pas. Le raffinement compare simultanément
l'espace et le temps : il ne mesure pas séparément leurs ordres.
L'expérience reste lisse ; aucun choc, fond sec ou ressaut n'est traité.

[Bilan](../results/tables/shallow_water_ledger.csv) ·
[Profils](../results/tables/shallow_water_profiles.csv).

## 7. Interconnexion

Deux lignes, $N_A=6$, $N_B=5$, coefficients constitutifs différents.
L'énergie initiale est uniquement dans A. Le terme de transmission directe
des ports est conservé lors de l'élimination.

- Défaut maximal du bilan total : $3.47\,10^{-16}$.
- Défaut maximal de cancellation des puissances internes : $8.33\,10^{-17}$.
- Énergie maximale atteinte dans B : $0.03839$.

![Échange d'énergie](../results/figures/interconnection.svg)

[Données](../results/tables/interconnection.csv).

## 8. Portée de la validation

Les tests unitaires examinent les identités et plusieurs erreurs possibles :
sous-intégration, mauvais ports, confusion énergie/symplecticité et échec de
Newton. Les 22 critères globaux sont définis dans validate() et consignés
dans [summary.json](../results/summary.json).

Les seuils sont explicites et tolèrent les différences d'arrondi entre
plateformes. Ils valident les expériences et les identités indiquées,
avec les hypothèses des notes théoriques, et ne constituent pas une
certification générale des modèles PH.
