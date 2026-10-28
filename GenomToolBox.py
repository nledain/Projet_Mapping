import numpy as np
from array import array
from collections import Counter


def counting_sort(items, key, k):
    """Tri par comptage STABLE de `items` selon key(item) dans [0, k[.

    Complexité O(len(items) + k).
    """
    compte = [0] * (k + 1)
    for x in items:
        compte[key(x) + 1] += 1
    for i in range(k):
        compte[i + 1] += compte[i]          # compte[c] = première place de la clé c
    sortie = [None] * len(items)
    for x in items:                         # parcours dans l'ordre -> stabilité
        c = key(x)
        sortie[compte[c]] = x
        compte[c] += 1
    return sortie


def radix_sort(items, keys, largeur, k):
    """Radix sort LSD : trie `items` selon le tuple keys(item) (longueur `largeur`,
    chaque composante dans [0, k[). Une passe de tri par comptage par composante,
    de la moins significative à la plus significative.
    """
    for d in range(largeur - 1, -1, -1):
        items = counting_sort(items, lambda x, d=d: keys(x)[d], k)
    return items


class SeqMetrics :
    
    def __init__(self,S,T):
        self.S = S.upper()
        self.T = T.upper()
        self.n = len(S)
        self.p = len(T)
         
    def Distance(self, Mode, afficher=True):
        """Renvoie la distance ; afficher=False permet de l'utiliser dans un calcul."""
        
        
        
        if Mode.upper() == "HAMMING" : 
            
            if  self.n != self.p :
                if afficher:
                    print("Two strings arn't the same length, please use Levenshtein Distance instead")
            
            else : 
                d = 0 
                for i in range(self.n):
                    if self.S[i] != self.T[i] :
                        d+=1
                if afficher:
                    print("Distance:",d)
                return d
        
        elif Mode.upper() == "LEVENSHTEIN" : 
            
            M = np.full((self.n+ 1, self.p+1), None, dtype=object)
            for i in range(self.n+1):
                for j in range(self.p+1):
                    if min((i,j)) == 0 :
                        M[i,j] = np.max((i,j))
                    elif self.S[i-1] == self.T[j-1] : 
                        M[i,j] = np.min((M[i-1,j]+ 1,M[i,j-1]+ 1,M[i-1,j-1]))
                    elif self.S[i-1] != self.T[j-1] :
                        M[i,j] = np.min((M[i-1,j]+1,M[i,j-1]+1,M[i-1,j-1] + 1))
            d = int(M[self.n,self.p])
            if afficher:
                print("Distance:", d)
            return d
        
        else :
            if afficher:
                print("Doesn't understand your request, please try Mode = 'LEVENSHTEIN' or Mode = 'HAMMING'")


    def Alignment(self, Mode):
        M = np.full((self.n + 1, self.p + 1), None, dtype=object)
        
        if Mode.upper() == "LOCAL": 
            
            for i in range(self.n + 1):
                for j in range(self.p + 1):
                    if i == 0 or j == 0:
                        M[i, j] = 0
                    elif self.S[i - 1] == self.T[j - 1]: 
                        M[i, j] = np.max((0, M[i - 1, j] - 1, M[i, j - 1] - 1, M[i - 1, j - 1] + 1))
                    elif self.S[i - 1] != self.T[j - 1]:
                        M[i, j] = np.max((0, M[i - 1, j] - 1, M[i, j - 1] - 1, M[i - 1, j - 1] - 1))
            
            
            max_val = -1
            max_i, max_j = 0, 0
            for i in range(self.n + 1):
                for j in range(self.p + 1):
                    if M[i, j] > max_val:
                        max_val = M[i, j]
                        max_i, max_j = i, j
            
            i, j = max_i, max_j
            aligned_S, aligned_T = [], []
            while i > 0 and j > 0 and M[i, j] > 0:
                score_diff = 1 if self.S[i - 1] == self.T[j - 1] else -1
                if M[i, j] == M[i - 1, j - 1] + score_diff:
                    aligned_S.append(self.S[i - 1])
                    aligned_T.append(self.T[j - 1])
                    i -= 1
                    j -= 1
                elif M[i, j] == M[i - 1, j] - 1:
                    aligned_S.append(self.S[i - 1])
                    aligned_T.append('-')
                    i -= 1
                elif M[i, j] == M[i, j - 1] - 1:
                    aligned_S.append('-')
                    aligned_T.append(self.T[j - 1])
                    j -= 1
                else:
                    break
            
            return "".join(reversed(aligned_S)), "".join(reversed(aligned_T)), int(max_val)
        
        elif Mode.upper() == "GLOBAL": 
            
            for i in range(self.n + 1):
                for j in range(self.p + 1):
                    if i == 0 or j == 0:
                        M[i, j] = -np.max((i, j))
                    elif self.S[i - 1] == self.T[j - 1]: 
                        M[i, j] = np.max((M[i - 1, j] - 1, M[i, j - 1] - 1, M[i - 1, j - 1] + 1))
                    elif self.S[i - 1] != self.T[j - 1]:
                        M[i, j] = np.max((M[i - 1, j] - 1, M[i, j - 1] - 1, M[i - 1, j - 1] - 1))
            
            
            i, j = self.n, self.p
            aligned_S, aligned_T = [], []
            while i > 0 or j > 0:
                if i > 0 and j > 0 and M[i, j] == M[i - 1, j - 1] + (1 if self.S[i - 1] == self.T[j - 1] else -1):
                    aligned_S.append(self.S[i - 1])
                    aligned_T.append(self.T[j - 1])
                    i -= 1
                    j -= 1
                elif i > 0 and M[i, j] == M[i - 1, j] - 1:
                    aligned_S.append(self.S[i - 1])
                    aligned_T.append('-')
                    i -= 1
                elif j > 0 and M[i, j] == M[i, j - 1] - 1:
                    aligned_S.append('-')
                    aligned_T.append(self.T[j - 1])
                    j -= 1
                else:
                    if i > 0:
                        aligned_S.append(self.S[i - 1])
                        aligned_T.append('-')
                        i -= 1
                    elif j > 0:
                        aligned_S.append('-')
                        aligned_T.append(self.T[j - 1])
                        j -= 1
            
            return "".join(reversed(aligned_S)), "".join(reversed(aligned_T)), int(M[self.n, self.p])
        
        else: 
            print("Doesn't understand your request, please try mode = 'GLOBAL' or mode = 'LOCAL'")


class Preparation : 

    def __init__(self, seq):
        if not isinstance(seq, str):
            raise TypeError("La séquence doit être une chaîne de caractères.")
        self.seq = seq.upper()
        if not self.seq or set(self.seq) - set("ACGTRYSWKMBDHVN"):
            raise ValueError("La séquence contient des caractères ADN invalides.")

    def Complementaire(self):
        correspondance = str.maketrans("ACGTRYSWKMBDHVN", "TGCAYRSWMKVHDBN")
        return self.seq.translate(correspondance)

    def ComplementaireInverse(self):
        return self.Complementaire()[::-1]

    def kmer(self, k):
        if type(k) is not int or k <= 0:
            raise ValueError("k doit être un entier strictement positif.")
        kmers = []
        for i in range(len(self.seq) - k + 1):
            kmers.append(self.seq[i:i + k])
        return kmers



class DC3:
    """Construction du tableau des suffixes par l'algorithme DC3 (Kärkkäinen-Sanders).

    Usage : DC3(sequence).DC3() renvoie la liste des positions des suffixes
    triés. `sequence` est une chaîne ASCII ou une liste d'entiers > 0
    (cas de la récursion).
    """

    def __init__(self, seq):
        self.seq = seq
        self._reset()

    def _reset(self):
        """Remet à zéro les structures de travail (permet de rappeler DC3())."""
        self.Clean_seq = []
        self.P0 = []
        self.P1 = []
        self.P2 = []
        self.P1_2 = []
        self.R1_2 = []
        self.R1_2_sorted = []
        self.Index = []
        self.Order = []
        self.T_prime = []
        self.Duplicates = None
        self.k = 0          # taille de l'alphabet (valeurs dans [0, k[)

    # Préparation des données
    def Clean(self):
        "Clean les datas"
        if isinstance(self.seq, str):
            # Chaque caractère devient son code ASCII (A=65, C=67, G=71, T=84 ...).
            # Tout alphabet ASCII est donc accepté, pas seulement ACGT.
            self.Clean_seq = [ord(c) for c in self.seq.upper()]
            if any(c <= 0 or c > 127 for c in self.Clean_seq):
                raise ValueError("La chaîne doit contenir des caractères ASCII non nuls.")
        else:
            self.Clean_seq = list(self.seq)

        "Ajout des sentinelles"
        self.Clean_seq.extend([0, 0, 0])
        self.k = max(self.Clean_seq) + 1

    def Build_Triplet(self, SEQ):
        n = len(SEQ)
        for i in range(n - 2):
            if i % 3 == 0:
                self.P0.append(i)
            elif i % 3 == 1:
                self.P1.append(i)
            else:
                self.P2.append(i)
        self.P1_2 = self.P1 + self.P2

        "Radix sort des triplets T[i],T[i+1],T[i+2] pour i dans P1_2"
        k = max(SEQ) + 1
        self.Index = radix_sort(
            self.P1_2,
            lambda i: (SEQ[i], SEQ[i + 1], SEQ[i + 2]),
            3, k)

        "Creation de R1_2 sorted + Index"
        self.R1_2_sorted = [(SEQ[i], SEQ[i + 1], SEQ[i + 2]) for i in self.Index]
        self.R1_2 = list(zip(self.R1_2_sorted, self.Index))

        "Creation de la liste ordre (rang 1, 2, ... ; égal si triplets égaux)"
        self.Order = []
        if self.R1_2_sorted:
            r = 1
            self.Order.append(r)
            for j in range(len(self.R1_2_sorted) - 1):
                if self.R1_2_sorted[j] != self.R1_2_sorted[j + 1]:
                    r += 1
                self.Order.append(r)

        "Modification du flag Duplicates"
        self.Duplicates = len(set(self.Order)) < len(self.Order)

    def Build_T_prime(self):
        rank_dict = dict(zip(self.Index, self.Order))
        self.T_prime = [rank_dict[i] for i in self.P1] + [rank_dict[i] for i in self.P2]

    def Backward(self, Index_prime):
        self.Index = [self.P1_2[i] for i in Index_prime]

    def Sort_P0(self, T, rank):
        """Trie les positions de P0 par (T[p], rang de p+1) avec un radix sort.

        rank : dict position -> rang dans le suffix array de P1_2.
        Le rang d'une position absente (au-delà de la fin) vaut -1, décalé ici
        à 0 pour avoir une clé positive : rang+1.
        """
        k1 = self.k
        k2 = len(rank) + 1
        return radix_sort(
            self.P0,
            lambda p: (T[p], rank.get(p + 1, -1) + 1),
            2, max(k1, k2))

    def DC3(self):

        self._reset()

        "Nettoyage"
        self.Clean()
        T = self.Clean_seq
        n = len(self.seq)

        "Trie des P"
        self.Build_Triplet(T)
        self.Build_T_prime()

        "Récursion"
        if self.Duplicates:
            sub_inst = DC3(self.T_prime)
            Index_prime = sub_inst.DC3()

            self.Backward(Index_prime)
            Index1_2 = self.Index.copy()

        else:
            Index1_2 = self.Index.copy()

        rank = {
            pos: r
            for r, pos in enumerate(Index1_2)
        }

        "Construction de P0 (tri par radix sort)"
        Index0 = self.Sort_P0(T, rank)

        "Fusion selon les règles"
        Index012 = []
        i = 0
        j = 0

        while i < len(Index0) and j < len(Index1_2):

            a = Index0[i]
            b = Index1_2[j]

            if a == n:
                Index012.append(a)
                i += 1

            elif b == n:
                Index012.append(b)
                j += 1

            elif b % 3 == 1:

                left = (T[a], rank.get(a + 1, -1))
                right = (T[b], rank.get(b + 1, -1))

                if left <= right:
                    Index012.append(a)
                    i += 1
                else:
                    Index012.append(b)
                    j += 1
            else:

                left = (T[a], T[a + 1], rank.get(a + 2, -1))
                right = (T[b], T[b + 1], rank.get(b + 2, -1))

                if left <= right:
                    Index012.append(a)
                    i += 1
                else:
                    Index012.append(b)
                    j += 1

        "ajout des dernières positions"
        Index012.extend(Index0[i:])
        Index012.extend(Index1_2[j:])

        "Retire les sentinelles"
        Index012 = [p for p in Index012 if p != n]

        return Index012


class BWT:
    """Transformée de Burrows-Wheeler et recherche exacte de motifs.

    Les positions commencent à 0.
    Le symbole $ est réservé à la fin du texte.
    Les insertions, délétions et substitutions ne sont pas traitées ici.
    """

    def __init__(self, texte):
        if not isinstance(texte, str):
            raise TypeError("Le texte doit être une chaîne de caractères.")

        if "$" not in texte:
            texte += "$"

        if texte.count("$") != 1 or not texte.endswith("$"):
            raise ValueError("Le symbole $ doit apparaître une seule fois, à la fin.")

        if any(ord(lettre) <= 36 or ord(lettre) > 127 for lettre in texte[:-1]):
            raise ValueError("Le texte doit contenir des lettres ASCII sans espaces.")

        self.texte = texte.upper()
        self.tableau_suffixes = self._construire_suffixes()
        self._construire_tables()

    def _construire_tables(self):
        """Construit bwt, C et les compteurs à partir de texte et tableau_suffixes."""
        # On prend la lettre qui précède chaque suffixe trié.
        # Pour la position 0, l'indice -1 correspond au $ final.
        self.bwt = "".join(
            self.texte[int(position) - 1]
            for position in self.tableau_suffixes
        )

        # C contient le nombre de caractères plus petits que chaque lettre.
        self.C = {}
        cumul = 0

        for lettre, nombre in sorted(Counter(self.bwt).items()):
            self.C[lettre] = cumul
            cumul += nombre

        # On conserve un compteur tous les 128 caractères pour limiter la mémoire.
        self.pas = 128
        self.compteurs = {}
        # Au plus 4**6 suffixes ADN : la mémoire ne grandit pas avec les reads.
        self.intervalles_courts = {}

        for lettre in self.C:
            compteurs = array("I", [0])
            cumul = 0

            for debut in range(0, len(self.bwt), self.pas):
                cumul += self.bwt.count(lettre, debut, debut + self.pas)
                compteurs.append(cumul)

            self.compteurs[lettre] = compteurs

    def _construire_suffixes(self):
        """Tableau des suffixes de texte$ construit avec l'algorithme DC3."""
        corps = self.texte[:-1]                 # texte sans le '$'
        sa = DC3(corps).DC3()                   # DC3 ne gère pas le '$' lui-même
        # Le suffixe '$' (position len(corps)) est le plus petit : il passe en tête.
        return array("I", [len(corps)] + sa)

    def occurrence(self, lettre, position):
        """Compte la lettre dans la BWT avant position, borne exclue."""
        if not isinstance(lettre, str) or len(lettre) != 1:
            raise ValueError("Il faut donner une seule lettre.")

        if type(position) is not int or not 0 <= position <= len(self.bwt):
            raise ValueError("La position est hors du texte.")

        if lettre not in self.compteurs:
            return 0

        bloc = position // self.pas
        debut = bloc * self.pas

        return (
            self.compteurs[lettre][bloc]
            + self.bwt.count(lettre, debut, position)
        )

    def _intervalle(self, motif):
        """Trouve les lignes correspondant au motif par recherche arrière."""
        if not isinstance(motif, str):
            raise TypeError("Le motif doit être une chaîne de caractères.")

        if not motif or "$" in motif:
            raise ValueError("Le motif doit être non vide et ne pas contenir $.")

        motif = motif.upper()
        suffixe = motif[-6:]
        if len(suffixe) == 6 and set(suffixe) <= set("ACGT"):
            intervalle = self.intervalles_courts.get(suffixe)
            if intervalle is None:
                intervalle = self._etendre_intervalle(suffixe, 0, len(self.bwt))
                self.intervalles_courts[suffixe] = intervalle
            # La recherche arrière commence souvent par les mêmes six lettres.
            return self._etendre_intervalle(motif[:-6], *intervalle)
        return self._etendre_intervalle(motif, 0, len(self.bwt))

    def _etendre_intervalle(self, motif, gauche, droite):
        """Continue une recherche arrière dont les bornes sont déjà valides."""
        if gauche >= droite:
            return 0, 0
        pas = self.pas
        compter = self.bwt.count

        for lettre in reversed(motif):
            if lettre not in self.C:
                return 0, 0
            compteurs = self.compteurs[lettre]
            bloc_gauche, bloc_droite = gauche // pas, droite // pas
            # Même calcul que occurrence ; les bornes proviennent de la BWT.
            gauche = self.C[lettre] + compteurs[bloc_gauche] + compter(lettre, bloc_gauche * pas, gauche)
            droite = self.C[lettre] + compteurs[bloc_droite] + compter(lettre, bloc_droite * pas, droite)

            if gauche >= droite:
                return 0, 0

        return gauche, droite

    def compter(self, motif):
        """Renvoie le nombre d'occurrences exactes du motif."""
        gauche, droite = self._intervalle(motif)
        return droite - gauche

    def rechercher(self, motif):
        """Renvoie les positions des occurrences exactes, dans l'ordre croissant."""
        gauche, droite = self._intervalle(motif)

        return sorted(
            int(position)
            for position in self.tableau_suffixes[gauche:droite]
        )

    def mot_a_position(self, position, k):
        """Renvoie le mot de longueur k commençant à la position demandée."""
        if type(position) is not int or type(k) is not int:
            raise ValueError("La position et k doivent être des entiers.")

        if position < 0 or k <= 0 or position + k > len(self.texte) - 1:
            raise ValueError("La position ou la longueur du mot est invalide.")

        return self.texte[position:position + k]

    @staticmethod
    def complementaire_inverse(sequence):
        """Renvoie le complémentaire inversé d'une séquence ADN."""
        if not isinstance(sequence, str):
            raise TypeError("La séquence doit être une chaîne de caractères.")

        sequence = sequence.upper()
        alphabet = set("ACGTRYSWKMBDHVN")

        if not sequence or set(sequence) - alphabet:
            raise ValueError("La séquence contient des caractères ADN invalides.")

        correspondance = str.maketrans(
            "ACGTRYSWKMBDHVN",
            "TGCAYRSWMKVHDBN"
        )

        return sequence.translate(correspondance)[::-1]

    def rechercher_adn(self, motif):
        """Recherche le motif sur les deux brins de la référence.

        La référence doit être en majuscules.
        Les lettres ambiguës sont recherchées littéralement, sans rôle de joker.
        La position reste celle de gauche sur la référence directe.
        Pour un palindrome, '+/-' indique les deux orientations possibles.
        """
        alphabet = set("ACGTRYSWKMBDHVN")

        if set(self.C) - alphabet - {"$"}:
            raise ValueError(
                "La référence doit être une séquence ADN en majuscules. "
                "Construis l'index avec BWT(sequence.upper())."
            )

        inverse = self.complementaire_inverse(motif)
        motif = motif.upper()

        brin = "+/-" if motif == inverse else "+"

        resultats = [
            {"position": position, "brin": brin}
            for position in self.rechercher(motif)
        ]

        if motif != inverse:
            resultats.extend(
                {"position": position, "brin": "-"}
                for position in self.rechercher(inverse)
            )

        return sorted(resultats, key=lambda resultat: resultat["position"])

    def inverser(self):
        """Reconstruit le texte depuis sa BWT, avec le $ final."""
        deja_vus = Counter()
        liens = array("I")

        # La première occurrence d'une lettre dans la dernière colonne
        # correspond à sa première occurrence dans la première colonne.
        for lettre in self.bwt:
            liens.append(self.C[lettre] + deja_vus[lettre])
            deja_vus[lettre] += 1

        ligne = self.bwt.index("$")
        resultat = []

        for _ in range(len(self.bwt)):
            resultat.append(self.bwt[ligne])
            ligne = liens[ligne]

        return "".join(reversed(resultat))


    def sauvegarder(self, chemin_fichier):
        """Enregistre le texte et le tableau des suffixes."""
        np.savez(
            chemin_fichier,
            texte=np.frombuffer(self.texte.encode("ascii"), dtype=np.uint8),
            suffixes=np.asarray(self.tableau_suffixes, dtype=np.uint32),
        )

    @classmethod
    def charger(cls, chemin_fichier):
        """Recharge un index sauvegardé, sans recalculer DC3."""
        with np.load(chemin_fichier, allow_pickle=False) as fichier:
            if not {"texte", "suffixes"}.issubset(fichier.files):
                raise ValueError("Le fichier ne contient pas un index BWT.")
            texte_code = fichier["texte"]
            suffixes = fichier["suffixes"]

        if texte_code.ndim != 1 or texte_code.dtype != np.uint8:
            raise ValueError("Le texte sauvegardé est invalide.")
        texte = texte_code.tobytes().decode("ascii")
        n = len(texte)
        if (texte.count("$") != 1 or not texte.endswith("$")
                or texte != texte.upper()
                or any(ord(c) <= 36 or ord(c) > 127 for c in texte[:-1])):
            raise ValueError("Le texte sauvegardé est invalide.")
        if (suffixes.ndim != 1 or suffixes.dtype.kind not in "iu"
                or len(suffixes) != n
                or np.any(suffixes < 0) or np.any(suffixes >= n)):
            raise ValueError("Le tableau des suffixes est invalide.")

        # Chaque position doit apparaître une fois dans le tableau.
        presents = np.zeros(n, dtype=bool)
        presents[suffixes] = True
        if not presents.all() or suffixes[0] != n - 1:
            raise ValueError("Le tableau des suffixes est invalide.")

        index = cls.__new__(cls)
        index.texte = texte
        index.tableau_suffixes = array("I", suffixes)
        index._construire_tables()
        return index

