# Contribuer

Toute modification de la géométrie doit conserver des conventions explicites
de signes, dimensions et produit de puissance. Ne pas réparer silencieusement
une matrice en l'antisymétrisant après assemblage.

Pour un nouveau Hamiltonien, fournir valeur, gradient, Hessienne et
quadrature constitutive. Séparer les erreurs de structure, de temps et d'espace.
Une bonne conservation de l'énergie n'est pas, seule, un test de trajectoire.

~~~bash
python -m unittest discover -s tests -v
python experiments/reproduce.py --check --no-figures --output /tmp/gpps-review
~~~

Documenter les résultats positifs et négatifs. Citer les méthodes existantes ;
ne pas qualifier de démontrée une propriété seulement observée.
