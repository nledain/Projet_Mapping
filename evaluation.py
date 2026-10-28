"""Comparaison de nos placements avec les positions du BAM fourni."""

import json
import re
import sqlite3
import tempfile
import warnings
from copy import deepcopy
from pathlib import Path

import bamnostic


def lire_verite_bam(chemin, identifiants):
    """Lit le BAM entier et garde uniquement les reads demandés.

    Les positions commencent à 0 et la fin est exclue, comme dans Python.
    Cette fonction concerne les reads simples du projet, pas les reads appariés.
    Le parcours est linéaire en taille du BAM, sans fichier d'index .bai.
    """
    identifiants = set(identifiants)
    verite = {}
    if not identifiants:
        return verite

    # L'index .bai n'est pas nécessaire pour un parcours séquentiel.
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="No supplied index file was not found.*",
                                category=UserWarning)
        fichier = bamnostic.AlignmentFile(str(chemin), "rb", require_index=False)
    with fichier as bam:
        for read in bam:
            if read.query_name not in identifiants:
                continue
            if read.is_secondary or read.is_supplementary:
                continue
            if read.query_name in verite:
                raise ValueError("Plusieurs alignements primaires pour " + read.query_name)

            if read.is_unmapped:
                verite[read.query_name] = {
                    "chromosome": None, "position": None, "fin": None,
                    "brin": None, "mappe": False, "cigar": None,
                }
            else:
                # M, D, N, = et X consomment des bases de la référence.
                longueur = sum(n for operation, n in read.cigartuples
                               if operation in (0, 2, 3, 7, 8))
                verite[read.query_name] = {
                    "chromosome": read.reference_name,
                    "position": read.pos,
                    "fin": read.pos + longueur,
                    "brin": "-" if read.is_reverse else "+",
                    "mappe": True,
                    "cigar": read.cigarstring,
                }

    manquants = identifiants - verite.keys()
    if manquants:
        warnings.warn(str(len(manquants)) + " identifiant(s) absent(s) du BAM.")
    return verite


def evaluer_mapping(resultats, verite, tolerance=0):
    """Compare le chromosome, la position gauche et le brin au BAM.

    L'origine est retrouvée si elle figure parmi les meilleurs placements.
    Cela ne signifie pas qu'un read ambigu est placé de manière certaine.
    La tolérance porte seulement sur la position, jamais sur le chromosome/brin.
    """
    if type(tolerance) is not int or tolerance < 0:
        raise ValueError("La tolérance doit être un entier positif ou nul.")

    bilan = {
        "tolerance": tolerance,
        "total": 0,
        "sans_verite": 0,
        "verite_non_mappee": 0,
        "evaluables": 0,
        "non_mappes": 0,
        "uniques": 0,
        "ambigus": 0,
        "recherches_limitees": 0,
        "uniques_evaluables": 0,
        "uniques_corrects": 0,
        "origine_retrouvee": 0,
        "mappes_sans_origine": 0,
        "details": [],
    }
    deja_vus = set()

    for resultat in resultats:
        identifiant = resultat["id"]
        if identifiant in deja_vus:
            raise ValueError("Résultat en double pour " + identifiant)
        deja_vus.add(identifiant)
        bilan["total"] += 1

        alignements = resultat["alignements"]
        # Seuls les ex aequo au meilleur score sont des réponses retenues.
        if alignements:
            distance = min(a["distance"] for a in alignements)
            meilleurs = [a for a in alignements if a["distance"] == distance]
        else:
            meilleurs = []

        # Deux copies de la même réponse ne rendent pas un read ambigu.
        positions = {(a["chromosome"], a["position"], a["fin"], a["brin"])
                     for a in meilleurs}
        unique = (len(positions) == 1
                  and next(iter(positions))[3] != "+/-"
                  and resultat.get("statut") != "ambigu")
        if not positions:
            statut = "non_mappe"
            bilan["non_mappes"] += 1
        elif unique:
            statut = "unique"
            bilan["uniques"] += 1
        else:
            statut = "ambigu"
            bilan["ambigus"] += 1

        limite = bool(resultat.get("recherche_limitee", False))
        bilan["recherches_limitees"] += int(limite)
        attendu = verite.get(identifiant)
        detail = {"id": identifiant, "statut": statut,
                  "recherche_limitee": limite, "origine_retrouvee": None}

        if attendu is None:
            bilan["sans_verite"] += 1
        elif not attendu.get("mappe", True):
            bilan["verite_non_mappee"] += 1
        else:
            bilan["evaluables"] += 1
            correct = attendu.get("origine_definie", True) and any(
                a["chromosome"] == attendu["chromosome"]
                and abs(a["position"] - attendu["position"]) <= tolerance
                and a["brin"] in (attendu["brin"], "+/-")
                for a in meilleurs
            )
            detail["origine_retrouvee"] = correct
            bilan["origine_retrouvee"] += int(correct)
            bilan["mappes_sans_origine"] += int(bool(positions) and not correct)
            if unique:
                bilan["uniques_evaluables"] += 1
                bilan["uniques_corrects"] += int(correct)
        bilan["details"].append(detail)

    # None signifie que le dénominateur est nul, pas que le score vaut zéro.
    n = bilan["total"]
    e = bilan["evaluables"]
    u = bilan["uniques_evaluables"]
    bilan["taux_mapping"] = (n - bilan["non_mappes"]) / n if n else None
    bilan["taux_origine_retrouvee"] = bilan["origine_retrouvee"] / e if e else None
    bilan["precision_uniques"] = bilan["uniques_corrects"] / u if u else None
    return bilan


def ajuster_verite_nettoyage(verite, reads_originaux, reads_nettoyes):
    """Adapte les positions du BAM après retrait de l'extrémité 3' des reads.

    Sur le brin inverse, cette coupe enlève le début de l'alignement sur la
    référence. Le CIGAR permet de compter séparément les bases du read et
    celles de la référence, notamment en présence d'insertions ou délétions.
    Coût linéaire dans les séquences vérifiées et les opérations du CIGAR.
    Les reads supprimés restent évaluables, mais sans origine retrouvable.
    """
    originaux = {}
    for read in reads_originaux:
        if read["id"] in originaux:
            raise ValueError("Read original en double : " + read["id"])
        originaux[read["id"]] = read
    resultat = {}
    vus = set()
    for read in reads_nettoyes:
        identifiant = read["id"]
        if identifiant not in originaux or identifiant in vus:
            raise ValueError("Identifiant nettoyé absent ou répété : " + identifiant)
        vus.add(identifiant)
        sequence = originaux[identifiant]["sequence"]
        longueur = len(sequence)
        conserve = len(read["sequence"])
        retire = longueur - conserve
        if retire < 0 or read["sequence"] != sequence[:conserve]:
            raise ValueError("Le nettoyage doit seulement retirer la fin du read : " + identifiant)
        if read.get("bases_retirees", retire) != retire:
            raise ValueError("Nombre de bases retirées incohérent : " + identifiant)
        if identifiant not in verite:
            continue
        attendu = deepcopy(verite[identifiant])
        resultat[identifiant] = attendu
        if not retire or not attendu.get("mappe", True):
            continue
        if not conserve:
            attendu.update(origine_definie=False, cigar=None)
            continue

        cigar = attendu.get("cigar") or ""
        operations = re.findall(r"(\d+)([MIDNSHP=X])", cigar)
        if (not operations or "".join(n + op for n, op in operations) != cigar
                or any(int(n) <= 0 for n, _ in operations)):
            raise ValueError("CIGAR absent ou invalide pour ajuster " + identifiant)
        longueur_cigar = sum(int(n) for n, op in operations if op in "MIS=X")
        if longueur_cigar != longueur:
            raise ValueError("Longueur FASTQ différente du CIGAR pour " + identifiant)
        if attendu["brin"] not in ("+", "-"):
            raise ValueError("Brin BAM invalide pour " + identifiant)

        # Intervalle conservé dans le read orienté comme la référence.
        debut_read, fin_read = (0, conserve) if attendu["brin"] == "+" else (retire, longueur)
        q, r = 0, attendu["position"]
        debut_ref = fin_ref = None
        for nombre, operation in operations:
            nombre = int(nombre)
            if operation in "M=X":
                gauche, droite = max(q, debut_read), min(q + nombre, fin_read)
                if gauche < droite:
                    if debut_ref is None:
                        debut_ref = r + gauche - q
                    fin_ref = r + droite - q
            if operation in "MIS=X":
                q += nombre
            if operation in "MDN=X":
                r += nombre
        if attendu.get("fin") is not None and r != attendu["fin"]:
            raise ValueError("Fin BAM différente du CIGAR pour " + identifiant)
        # Le CIGAR d'origine ne décrit plus le read raccourci.
        attendu["cigar"] = None
        attendu["origine_definie"] = debut_ref is not None
        if debut_ref is not None:
            attendu.update(position=debut_ref, fin=fin_ref)
    if vus != originaux.keys():
        raise ValueError("Les mêmes identifiants sont nécessaires avant et après nettoyage.")
    return resultat


def _ajouter_au_bilan(bilan, resultat, attendu):
    """Ajoute un read aux compteurs, sans conserver tous les détails."""
    verite = {resultat["id"]: attendu} if attendu is not None else {}
    partie = evaluer_mapping([resultat], verite, bilan["tolerance"])
    # On additionne les comptes ; les taux seront recalculés à la fin.
    for cle, valeur in partie.items():
        if cle != "tolerance" and type(valeur) is int:
            bilan[cle] += valeur

    detail = partie["details"][0]
    if detail["origine_retrouvee"] is None:
        categorie = "Sans vérité mappée"
    elif detail["statut"] == "non_mappe":
        categorie = "Non mappé"
    elif not detail["origine_retrouvee"]:
        categorie = "Placé ailleurs"
    elif detail["statut"] == "unique":
        categorie = "Unique, correct"
    else:
        categorie = "Ambigu, origine présente"
    bilan["categories"][categorie] = bilan["categories"].get(categorie, 0) + 1


def _calculer_taux(bilan):
    n, e, u = bilan["total"], bilan["evaluables"], bilan["uniques_evaluables"]
    bilan["taux_mapping"] = (n - bilan["non_mappes"]) / n if n else None
    bilan["taux_origine_retrouvee"] = bilan["origine_retrouvee"] / e if e else None
    bilan["precision_uniques"] = bilan["uniques_corrects"] / u if u else None


def evaluer_fichier(chemin_jsonl, chemin_bam, dossier_travail=None,
                    tolerance=3, progression=100000, identifiants_apercu=()):
    """Compare les résultats au BAM sans tout charger en mémoire.

    SQLite conserve les résultats sur disque pour les retrouver par identifiant,
    même si l'ordre du BAM est différent de celui du FASTQ. Le fichier temporaire
    est supprimé à la fin. On renvoie les bilans strict et avec tolérance.
    On peut garder la vérité de 300 identifiants au plus pour les petits essais.
    """
    if type(tolerance) is not int or tolerance < 0:
        raise ValueError("La tolérance doit être un entier positif ou nul.")
    if type(progression) is not int or progression < 0:
        raise ValueError("La progression doit être un entier positif ou nul.")
    apercu = set()
    for identifiant in identifiants_apercu:
        if not isinstance(identifiant, str) or not identifiant:
            raise ValueError("Un identifiant de l'aperçu est invalide.")
        apercu.add(identifiant)
        if len(apercu) > 300:
            raise ValueError("L'aperçu est limité à 300 identifiants.")
    bilan = evaluer_mapping([], {}, tolerance=0)
    bilan_tolere = evaluer_mapping([], {}, tolerance=tolerance) if tolerance else None
    bilans = (bilan, bilan_tolere) if tolerance else (bilan,)
    for mesure in bilans:
        mesure["categories"] = {}
        mesure["repartition_reference"] = {}
        mesure["lectures_bam"] = 0
        mesure["verite_apercu"] = {}

    dossier = Path(dossier_travail) if dossier_travail is not None else Path(chemin_jsonl).parent
    dossier = dossier.resolve()
    dossier.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="evaluation_", dir=dossier) as nom:
        base = sqlite3.connect(str(Path(nom) / "resultats.sqlite"))
        try:
            # Environ 8 Mio de cache ; les résultats restent sur disque.
            base.execute("PRAGMA cache_size = -8192")
            base.execute("CREATE TABLE reads (id TEXT PRIMARY KEY, resultat TEXT, vu INTEGER DEFAULT 0)")
            with open(chemin_jsonl, encoding="utf-8") as fichier:
                for numero, ligne in enumerate(fichier, 1):
                    resultat = json.loads(ligne)
                    identifiant = resultat["id"]
                    if not isinstance(identifiant, str) or not identifiant:
                        raise ValueError("Identifiant de read invalide à la ligne " + str(numero))
                    try:
                        base.execute("INSERT INTO reads (id, resultat) VALUES (?, ?)",
                                     (identifiant, ligne))
                    except sqlite3.IntegrityError as erreur:
                        raise ValueError("Résultat en double pour " + identifiant) from erreur
                    if numero % 10000 == 0:
                        base.commit()
            base.commit()

            with warnings.catch_warnings():
                warnings.filterwarnings("ignore", message="No supplied index file was not found.*",
                                        category=UserWarning)
                fichier_bam = bamnostic.AlignmentFile(str(chemin_bam), "rb", require_index=False)
            repartition = {}
            total_bam = 0
            avec_verite = 0
            with fichier_bam as bam:
                for total_bam, read in enumerate(bam, 1):
                    if progression and total_bam % progression == 0:
                        print(f"Évaluation : {total_bam:,} enregistrements BAM parcourus", flush=True)
                    if read.is_secondary or read.is_supplementary:
                        continue
                    ligne = base.execute("SELECT resultat, vu FROM reads WHERE id = ?",
                                         (read.query_name,)).fetchone()
                    if ligne is None and read.query_name not in apercu:
                        continue
                    if ((ligne is not None and ligne[1])
                            or read.query_name in bilan["verite_apercu"]):
                        raise ValueError("Plusieurs alignements primaires pour " + read.query_name)
                    if read.is_unmapped:
                        attendu = {"chromosome": None, "position": None, "fin": None,
                                   "brin": None, "mappe": False, "cigar": None}
                    else:
                        longueur = sum(n for operation, n in read.cigartuples
                                       if operation in (0, 2, 3, 7, 8))
                        attendu = {"chromosome": read.reference_name, "position": read.pos,
                                   "fin": read.pos + longueur,
                                   "brin": "-" if read.is_reverse else "+", "mappe": True,
                                   "cigar": read.cigarstring}
                    if read.query_name in apercu:
                        for mesure in bilans:
                            mesure["verite_apercu"][read.query_name] = dict(attendu)
                    if ligne is None:
                        continue
                    resultat = json.loads(ligne[0])
                    if attendu["mappe"]:
                        chromosome = attendu["chromosome"]
                        repartition[chromosome] = repartition.get(chromosome, 0) + 1
                    for mesure in bilans:
                        _ajouter_au_bilan(mesure, resultat, attendu)
                    base.execute("UPDATE reads SET vu = 1 WHERE id = ?", (read.query_name,))
                    avec_verite += 1
                    if avec_verite % 10000 == 0:
                        base.commit()
            base.commit()

            # Les reads absents du BAM restent dans le dénominateur total.
            for (ligne,) in base.execute("SELECT resultat FROM reads WHERE vu = 0"):
                resultat = json.loads(ligne)
                for mesure in bilans:
                    _ajouter_au_bilan(mesure, resultat, None)
            for mesure in bilans:
                mesure["lectures_bam"] = total_bam
                mesure["repartition_reference"] = dict(repartition)
                _calculer_taux(mesure)
            if bilan["sans_verite"]:
                warnings.warn(str(bilan["sans_verite"]) + " identifiant(s) absent(s) du BAM.")
        finally:
            base.close()
    return bilan, bilan_tolere if tolerance else deepcopy(bilan)
