# Projet génomique

Le but est de retrouver la position de reads sur le génome de référence. Le programme utilise les k-mers et la BWT, puis compare les résultats au fichier BAM.

## Pour lancer

Ouvrir ce dossier dans VS Code, puis `ProjetMappingNotebook.ipynb`.
Les trois fichiers FASTA, FASTQ et BAM peuvent rester à côté du notebook, comme sur GitLab, ou être regroupés dans `donnees/`, sans changer leurs noms.

À la première utilisation, exécuter la cellule d'installation :

```python
%pip install -r requirements.txt
```

Redémarrer le noyau, puis lancer les cellules dans l'ordre.

- `NOMBRE_READS = 100` pour tester rapidement.
- `NOMBRE_READS = 0` pour tout le fichier. Cela peut prendre plusieurs heures.

Le notebook garde automatiquement 30 reads pour les illustrations et les petits essais sur k et le nettoyage. Le bilan principal utilise tous les reads demandés. Relancer la cellule de mapping recommence le calcul ; pour revoir un graphique dans la même session, relancer seulement sa cellule.

## Les fichiers

- `GenomToolBox.py` contient notamment DC3 et la BWT.
- `Workflow.py` fait le mapping.
- `evaluation.py` compare les placements au BAM.

## Comment le code est organisé

Le notebook lance les étapes dans cet ordre : lecture du génome, construction des index, lecture des reads, mapping, puis comparaison au BAM. Le BAM sert à vérifier les réponses, pas à trouver les positions.

Une **classe** regroupe des données et des fonctions. Un **objet** est une utilisation concrète de cette classe. Ici, chaque chromosome a son objet `BWT`, avec sa séquence et les tables pour chercher dedans. Le dictionnaire `indexes` associe le nom du chromosome à cet objet. Une fonction d'une classe s'appelle une **méthode**.

Dans `GenomToolBox.py` :

- `Preparation` vérifie la séquence et permet de faire les k-mers ou le complémentaire inversé.
- `DC3` construit le tableau des suffixes, utilisé par la BWT.
- `BWT` sert aux recherches : `compter` donne le nombre d'occurrences et `rechercher_adn` donne les positions sur les deux brins.
- `SeqMetrics` contient des outils de distance et d'alignement. Pour le mapping actuel, la vérification des candidats se fait avec `aligner_debut` dans `Workflow.py`.

Les fonctions principales du workflow sont :

- `lire_reference` et `lire_fastq` : lire les fichiers. Un read est un dictionnaire avec son identifiant, sa séquence et ses qualités.
- `construire_index` : créer ou recharger les objets BWT des chromosomes.
- `mapper_read` : chercher les k-mers d'un read, faire voter les positions candidates et vérifier les meilleurs placements.
- `mapper_flux` : traiter les reads un par un et écrire les résultats au fur et à mesure.
- `nettoyer_read` : couper les bases de mauvaise qualité à la fin, pour le petit essai de nettoyage.

Dans `evaluation.py`, `evaluer_fichier` compare les placements au BAM. `ajuster_verite_nettoyage` adapte les positions attendues quand un read a été raccourci.

## Ce qu'on peut modifier

Les réglages sont dans la cellule des paramètres du notebook :

- `K` : taille des k-mers. Plus petit, on trouve davantage de répétitions ; plus grand, les erreurs peuvent gêner la recherche.
- `ERREURS` : nombre maximal de substitutions, insertions et délétions acceptées au total.
- `MAX_OCCURRENCES` et `MAX_CANDIDATS` : limites pour éviter trop de recherches. Les augmenter peut retrouver plus de placements, mais prend plus de temps.
- `GRAINE` : changer le tirage de l'échantillon. Garder la même valeur permet de refaire le même essai.

Le seuil du nettoyage se change avec `SEUIL_QUALITE`, dans sa cellule. Pour utiliser d'autres données, modifier les noms `FASTA`, `FASTQ` et `BAM` au début du notebook, en gardant une référence et un BAM correspondant aux reads.

Les graphiques s'affichent dans le notebook. Les résultats sont enregistrés dans un dossier `sorties_...` : `placements.csv`, `placements.jsonl` et `bilan.json`.
