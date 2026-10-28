"""Mapping de reads par k-mers, puis vérification des régions candidates."""

import argparse
import csv
import gzip
import json
import random
from collections import Counter
from pathlib import Path
from time import perf_counter

from Bio import SeqIO
from GenomToolBox import BWT, Preparation


def ouvrir_texte(chemin):
    if str(chemin).lower().endswith(".gz"):
        return gzip.open(chemin, "rt")
    return open(chemin, "r")


def lire_reference(chemin):
    references = {}
    with ouvrir_texte(chemin) as fichier:
        for record in SeqIO.parse(fichier, "fasta"):
            if record.id in references:
                raise ValueError("Identifiant FASTA répété : " + record.id)
            references[record.id] = Preparation(str(record.seq)).seq
    if not references:
        raise ValueError("Le fichier FASTA est vide.")
    return references


def lire_fastq(chemin):
    """Lit un read à la fois pour limiter la mémoire."""
    with ouvrir_texte(chemin) as fichier:
        for record in SeqIO.parse(fichier, "fastq"):
            yield {"id": record.id, "sequence": str(record.seq).upper(),
                   "qualite": record.letter_annotations["phred_quality"]}


def nettoyer_read(read, seuil=20):
    """Retire seulement les bases de mauvaise qualité à l'extrémité 3'."""
    if type(seuil) is not int or seuil < 0:
        raise ValueError("Le seuil de qualité doit être un entier positif ou nul.")
    sequence, qualites = read["sequence"], read["qualite"]
    if len(sequence) != len(qualites):
        raise ValueError("La séquence et ses qualités doivent avoir la même longueur.")
    fin = len(sequence)
    while fin > 0 and qualites[fin - 1] < seuil:
        fin -= 1
    return {**read, "sequence": sequence[:fin], "qualite": list(qualites[:fin]),
            "bases_retirees": len(sequence) - fin}


def echantillonner_fastq(chemin, nombre=300, graine=42):
    """Échantillon réservoir : chaque read a la même chance d'être retenu."""
    if type(nombre) is not int or nombre <= 0:
        raise ValueError("Le nombre de reads doit être un entier positif.")
    hasard = random.Random(graine)
    reads = []
    total = 0
    for total, read in enumerate(lire_fastq(chemin), 1):
        if total <= nombre:
            reads.append(read)
        else:
            position = hasard.randrange(total)
            if position < nombre:
                reads[position] = read
    return reads, total


def construire_index(references, dossier_cache="index_bwt"):
    """Un index par chromosome : aucun placement ne traverse une frontière."""
    if not references:
        raise ValueError("Il faut au moins une séquence de référence.")
    dossier = Path(dossier_cache) if dossier_cache is not None else None
    if dossier is not None:
        dossier.mkdir(parents=True, exist_ok=True)
    indexes = {}
    for numero, (nom, sequence) in enumerate(references.items(), 1):
        sequence = Preparation(sequence).seq
        debut = perf_counter()
        chemin = dossier / f"chromosome_{numero}.npz" if dossier else None
        index = None
        if chemin is not None and chemin.exists():
            index = BWT.charger(chemin)
            # Un autre FASTA peut avoir été utilisé auparavant.
            if index.texte != sequence + "$":
                index = None
        action = "chargé"
        if index is None:
            index = BWT(sequence)
            if chemin is not None:
                index.sauvegarder(chemin)
            action = "construit"
        indexes[nom] = index
        print(f"{nom} : {len(sequence):,} bases, index {action} "
              f"en {perf_counter() - debut:.2f} s", flush=True)
    return indexes


def aligner_debut(read, reference, position, erreurs):
    """Aligne le read entier à partir de position, avec au plus erreurs écarts.

    La fin est libre entre L-d et L+d bases de référence. Une insertion
    désigne une base du read en plus par rapport à la référence.
    """
    longueur_reference = len(reference) - int(reference.endswith("$"))
    if position < 0 or position >= longueur_reference:
        return None
    n = len(read)
    fenetre = reference[position:min(position + n + erreurs, longueur_reference)]
    p = len(fenetre)
    if p < n - erreurs:
        return None

    # On stocke seulement la bande |i-j| <= erreurs, avec une case de marge.
    # La colonne j de la ligne i est à l'indice j-i+erreurs+1 dans cette bande.
    M = [[erreurs + 1] * (2 * erreurs + 3) for _ in range(n + 1)]
    for i in range(min(n, erreurs) + 1):
        M[i][erreurs + 1 - i] = i
    for j in range(min(p, erreurs) + 1):
        M[0][j + erreurs + 1] = j
    for i in range(1, n + 1):
        precedente, ligne = M[i - 1], M[i]
        for j in range(max(1, i - erreurs), min(p, i + erreurs) + 1):
            b = j - i + erreurs + 1
            ligne[b] = min(precedente[b + 1] + 1,
                           ligne[b - 1] + 1,
                           precedente[b] + (read[i - 1] != fenetre[j - 1]))
        if min(ligne) > erreurs:
            return None

    fins = range(max(1, n - erreurs), p + 1)
    j = min(fins, key=lambda fin: (M[n][fin - n + erreurs + 1], abs(fin - n), fin))
    distance = M[n][j - n + erreurs + 1]
    if distance > erreurs:
        return None
    fin = position + j
    i = n
    substitutions = insertions = deletions = 0
    while i > 0 or j > 0:
        b = j - i + erreurs + 1
        if (i > 0 and j > 0 and
                M[i][b] == M[i - 1][b] + (read[i - 1] != fenetre[j - 1])):
            substitutions += int(read[i - 1] != fenetre[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and M[i][b] == M[i - 1][b + 1] + 1:
            insertions += 1
            i -= 1
        else:
            deletions += 1
            j -= 1
    return {"position": position, "fin": fin, "distance": distance,
            "substitutions": substitutions, "insertions": insertions,
            "deletions": deletions}


def mapper_read(read_id, sequence, indexes, k=21, erreurs=3,
                max_occurrences=200, max_candidats=30):
    """Renvoie les meilleurs placements parmi les régions examinées.

    Positions à partir de 0, fin exclue, y compris sur le brin inverse.
    Les lettres IUPAC sont comparées littéralement, sans effet de joker.
    """
    for nom, valeur in [("k", k), ("max_occurrences", max_occurrences),
                         ("max_candidats", max_candidats)]:
        if type(valeur) is not int or valeur <= 0:
            raise ValueError(nom + " doit être un entier positif.")
    if type(erreurs) is not int or erreurs < 0:
        raise ValueError("erreurs doit être un entier positif ou nul.")
    if not indexes:
        raise ValueError("Aucun index de référence.")
    sequence = Preparation(sequence).seq
    resultat = {"id": read_id, "statut": "non_mappe", "alignements": [],
                "recherche_limitee": False, "kmers_ignores": 0}
    if len(sequence) < k:
        resultat["raison"] = "read_plus_court_que_k"
        return resultat
    if erreurs >= len(sequence):
        raise ValueError("Le nombre d'erreurs doit être inférieur à la longueur du read.")

    inverse = BWT.complementaire_inverse(sequence)
    orientations = [(sequence, "+"), (inverse, "-")]
    if sequence == inverse:
        orientations = [(sequence, "+/-")]

    # Si le read entier est exact, aucune distance plus petite n'est possible.
    exacts = []
    total_exacts = 0
    for read, brin in orientations:
        for chromosome, index in indexes.items():
            gauche, droite = index._intervalle(read)
            total_exacts += droite - gauche
            reste = max(0, max_occurrences - len(exacts))
            for position in index.tableau_suffixes[gauche:min(droite, gauche + reste)]:
                position = int(position)
                exacts.append({"chromosome": chromosome, "position": position,
                               "fin": position + len(read), "brin": brin,
                               "distance": 0, "substitutions": 0, "insertions": 0,
                               "deletions": 0, "support": len(read) - k + 1})
    if exacts:
        resultat["alignements"] = sorted(exacts, key=lambda a: (a["chromosome"], a["position"], a["brin"]))
        resultat["recherche_limitee"] = total_exacts > len(exacts)
        resultat["statut"] = "ambigu" if total_exacts > 1 or sequence == inverse else "mappe"
        return resultat

    # Les k-mers votent pour un début : position trouvée - décalage dans le read.
    votes = {}
    for read, brin in orientations:
        positions_kmers = {}
        for decalage, mot in enumerate(Preparation(read).kmer(k)):
            positions_kmers.setdefault(mot, []).append(decalage)
        for mot, decalages in positions_kmers.items():
            intervalles = []
            total = 0
            for chromosome, index in indexes.items():
                gauche, droite = index._intervalle(mot)
                total += droite - gauche
                intervalles.append((chromosome, gauche, droite))
            if total > max_occurrences:
                resultat["kmers_ignores"] += len(decalages)
                resultat["recherche_limitee"] = True
                continue
            for chromosome, gauche, droite in intervalles:
                index = indexes[chromosome]
                for position in index.tableau_suffixes[gauche:droite]:
                    for decalage in decalages:
                        debut = int(position) - decalage
                        if -erreurs <= debut < len(index.texte) - 1:
                            cle = (chromosome, brin, debut)
                            votes.setdefault(cle, set()).add(decalage)

    candidats = sorted(votes, key=lambda cle: (-len(votes[cle]), cle))
    if len(candidats) > max_candidats:
        resultat["recherche_limitee"] = True
    candidats = candidats[:max_candidats]
    sequences = {brin: read for read, brin in orientations}
    visites = set()
    alignements = []
    meilleure_distance = erreurs
    for chromosome, brin, debut in candidats:
        reference = indexes[chromosome].texte
        # Un indel peut décaler le vote de quelques bases.
        for position in range(max(0, debut - erreurs), debut + erreurs + 1):
            cle = (chromosome, brin, position)
            if cle in visites:
                continue
            visites.add(cle)
            # Une réponse plus mauvaise que la meilleure connue sera écartée.
            alignement = aligner_debut(sequences[brin], reference, position, meilleure_distance)
            if alignement is not None:
                meilleure_distance = alignement["distance"]
                soutiens = set()
                for origine in range(position - erreurs, position + erreurs + 1):
                    soutiens.update(votes.get((chromosome, brin, origine), set()))
                alignement.update(chromosome=chromosome, brin=brin, support=len(soutiens))
                alignements.append(alignement)
    if alignements:
        meilleur = min(a["distance"] for a in alignements)
        alignements = [a for a in alignements if a["distance"] == meilleur]
        resultat["alignements"] = sorted(alignements, key=lambda a: (a["chromosome"], a["position"], a["brin"]))
        resultat["statut"] = "ambigu" if len(alignements) > 1 or sequence == inverse else "mappe"
    else:
        resultat["raison"] = "aucun_candidat_valide"
    return resultat


def mapper_reads(reads, indexes, **parametres):
    resultats = []
    for read in reads:
        if not read["sequence"] and read.get("bases_retirees", 0) > 0:
            # Un read entièrement nettoyé reste compté parmi les non-mappés.
            resultat = {"id": read["id"], "statut": "non_mappe", "alignements": [],
                        "recherche_limitee": False, "kmers_ignores": 0,
                        "raison": "read_plus_court_que_k"}
        else:
            resultat = mapper_read(read["id"], read["sequence"], indexes, **parametres)
        resultats.append(resultat)
    return resultats


def sauvegarder_resultats(resultats, chemin):
    """Une ligne par placement ; une ligne vide de position pour un read non mappé."""
    colonnes = ["id", "statut", "chromosome", "position", "fin", "brin", "distance",
                "substitutions", "insertions", "deletions", "support",
                "recherche_limitee", "kmers_ignores", "raison"]
    with open(chemin, "w", newline="", encoding="utf-8") as fichier:
        sortie = csv.DictWriter(fichier, fieldnames=colonnes)
        sortie.writeheader()
        for numero, resultat in enumerate(resultats, 1):
            infos = {cle: resultat.get(cle, "") for cle in
                     ["id", "statut", "recherche_limitee", "kmers_ignores", "raison"]}
            for alignement in resultat["alignements"] or [{}]:
                sortie.writerow({**infos, **alignement})
            if numero % 1000 == 0:
                fichier.flush()


def mapper_flux(reads, indexes, chemin_csv, chemin_jsonl, apercu=300,
                graine=42, afficher_tous=1000, **parametres):
    """Traite tous les reads reçus, sans conserver tous les résultats en mémoire.

    Le CSV contient les placements, le JSONL un résultat complet par ligne.
    On ne garde qu'un petit aperçu pour les figures du notebook.
    """
    if type(apercu) is not int or apercu < 0:
        raise ValueError("apercu doit être un entier positif ou nul.")
    if type(afficher_tous) is not int or afficher_tous < 0:
        raise ValueError("afficher_tous doit être un entier positif ou nul.")
    if Path(chemin_csv).resolve() == Path(chemin_jsonl).resolve():
        raise ValueError("Il faut deux fichiers différents pour CSV et JSONL.")
    valeurs = {"k": 21, "erreurs": 3, "max_occurrences": 200, "max_candidats": 30}
    for cle, valeur in parametres.items():
        if cle not in valeurs:
            raise TypeError("Paramètre inconnu : " + cle)
        valeurs[cle] = valeur
    for cle, valeur in valeurs.items():
        minimum = 0 if cle == "erreurs" else 1
        if type(valeur) is not int or valeur < minimum:
            raise ValueError("Valeur invalide pour " + cle)
    if not indexes:
        raise ValueError("Aucun index de référence.")

    resume = {"total": 0, "statuts": Counter(), "recherches_limitees": 0,
              "longueurs": Counter(), "qualites": [],
              "reads_apercu": [], "resultats_apercu": []}
    hasard = random.Random(graine)
    debut = perf_counter()

    def resultats_progressifs():
        with open(chemin_jsonl, "w", encoding="utf-8") as fichier:
            for numero, read in enumerate(reads, 1):
                resultat = mapper_read(read["id"], read["sequence"], indexes, **parametres)
                fichier.write(json.dumps(resultat, ensure_ascii=False) + "\n")
                resume["total"] = numero
                resume["statuts"][resultat["statut"]] += 1
                resume["recherches_limitees"] += int(resultat["recherche_limitee"])
                resume["longueurs"][len(read["sequence"])] += 1

                # Un histogramme par position remplace la matrice de tous les reads.
                for position, qualite in enumerate(read["qualite"]):
                    if position == len(resume["qualites"]):
                        resume["qualites"].append(Counter())
                    resume["qualites"][position][int(qualite)] += 1

                # L'aperçu est tiré au hasard ; le mapping porte bien sur tous les reads.
                place = numero - 1 if numero <= apercu else hasard.randrange(numero)
                if place < apercu:
                    if numero <= apercu:
                        resume["reads_apercu"].append(read)
                        resume["resultats_apercu"].append(resultat)
                    else:
                        resume["reads_apercu"][place] = read
                        resume["resultats_apercu"][place] = resultat
                if afficher_tous and numero % afficher_tous == 0:
                    fichier.flush()
                    print(f"{numero:,} reads traités en {perf_counter() - debut:.1f} s", flush=True)
                yield resultat

    sauvegarder_resultats(resultats_progressifs(), chemin_csv)
    resume["temps_secondes"] = perf_counter() - debut
    resume["statuts"] = dict(resume["statuts"])
    resume["longueurs"] = dict(resume["longueurs"])
    resume["qualites"] = [dict(histogramme) for histogramme in resume["qualites"]]
    print(f"Terminé : {resume['total']:,} reads, {resume['temps_secondes']:.2f} s", flush=True)
    return resume


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", required=True)
    parser.add_argument("--fastq", required=True)
    parser.add_argument("--sortie", default="resultats.csv")
    parser.add_argument("--nombre", type=int, default=0, help="0 : tous les reads (défaut) ; sinon échantillon aléatoire")
    parser.add_argument("--graine", type=int, default=42)
    parser.add_argument("--k", type=int, default=21)
    parser.add_argument("--erreurs", type=int, default=3)
    parser.add_argument("--max-occurrences", type=int, default=200)
    parser.add_argument("--max-candidats", type=int, default=30)
    parser.add_argument("--cache", default="index_bwt")
    args = parser.parse_args()
    if args.nombre < 0:
        parser.error("--nombre doit être positif ou nul")
    if args.k <= 0 or args.erreurs < 0 or args.max_occurrences <= 0 or args.max_candidats <= 0:
        parser.error("k et les plafonds doivent être positifs ; erreurs doit être positif ou nul")
    debut = perf_counter()
    indexes = construire_index(lire_reference(args.reference), args.cache)
    print(f"Préparation des index : {perf_counter() - debut:.2f} s")
    if args.nombre == 0:
        reads = lire_fastq(args.fastq)
    else:
        reads, total = echantillonner_fastq(args.fastq, args.nombre, args.graine)
        print(f"Échantillon : {len(reads)} reads parmi {total}")
    chemin_jsonl = Path(args.sortie).with_suffix(".jsonl")
    resume = mapper_flux(reads, indexes, args.sortie, chemin_jsonl, graine=args.graine,
                         k=args.k, erreurs=args.erreurs,
                         max_occurrences=args.max_occurrences, max_candidats=args.max_candidats)
    with open(Path(args.sortie).with_suffix(".bilan.json"), "w", encoding="utf-8") as fichier:
        json.dump(resume, fichier, ensure_ascii=False, indent=2)
    print("Placements :", args.sortie)


if __name__ == "__main__":
    main()
