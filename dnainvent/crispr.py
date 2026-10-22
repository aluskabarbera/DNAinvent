"""Edición genética con CRISPR/Cas9 (Cas9 de Streptococcus pyogenes).

Cómo funciona en el laboratorio, y cómo lo modela este módulo:

1) BUSCAR DIANAS  (find_targets)
   Cas9 solo puede cortar junto a una PAM: una secuencia NGG (cualquier base
   seguida de dos G). Los 20 nucleótidos justo antes de la PAM son el
   "protoespaciador", y el ARN guía se diseña con esa misma secuencia.

       5'- ...[ protoespaciador de 20 nt ][N G G]... -3'
                                    ^
                          corte a 3 pb de la PAM

   La PAM puede estar en cualquiera de las dos hebras. En la hebra de abajo
   se ve, leída desde la de arriba, como CCN (la complementaria inversa de NGG).

2) CORTAR
   Cas9 rompe las dos hebras 3 pares de bases antes de la PAM.

3) REPARAR  (repair_nhej / repair_hdr)
   La célula no deja el ADN roto. Tiene dos vías:
   - NHEJ (unión de extremos no homólogos): pega los extremos a toda prisa y
     suele equivocarse, insertando o borrando unas pocas bases. Si el cambio no
     es múltiplo de 3, desplaza el marco de lectura y suele inactivar el gen.
   - HDR (reparación dirigida por homología): si se aporta un ADN molde con
     brazos iguales a los lados del corte, la célula lo copia. Así se inserta
     la secuencia que se quiera.
"""

import random
from dataclasses import dataclass

from .sequences import reverse_complement
from .validation import InvalidSequenceError, validate_dna

PAM_LENGTH = 3
GUIDE_LENGTH = 20
CUT_OFFSET = 3          # Cas9 corta 3 pb antes de la PAM
HOMOLOGY_ARM = 30       # longitud de cada brazo del molde para HDR
MAX_NHEJ_DELETION = 10


def _is_pam(triplet: str) -> bool:
    """PAM de SpCas9: NGG."""
    return len(triplet) == PAM_LENGTH and triplet[1:] == "GG"


@dataclass(frozen=True)
class Target:
    """Una diana posible para Cas9.

    Todas las posiciones están en coordenadas de la hebra superior (la introducida),
    desde 0, y los intervalos son [inicio, fin) como en los slices de Python.

    strand         '+' si la PAM está en la hebra superior, '-' si está en la inferior.
    protospacer    Secuencia de la diana leída 5'->3' en su propia hebra.
    pam            PAM leída 5'->3' en su propia hebra (siempre NGG).
    start, end     Dónde está el protoespaciador en la hebra superior.
    pam_start, pam_end  Dónde está la PAM en la hebra superior.
    cut            Punto de corte: entre cut - 1 y cut.
    """

    strand: str
    protospacer: str
    pam: str
    start: int
    end: int
    pam_start: int
    pam_end: int
    cut: int

    @property
    def guide_rna(self) -> str:
        """El ARN guía tiene la misma secuencia que el protoespaciador, con U en vez de T."""
        return self.protospacer.replace("T", "U")


def find_targets(dna: str, guide_length: int = GUIDE_LENGTH) -> list[Target]:
    """Todas las dianas NGG de ambas hebras, ordenadas por punto de corte."""
    if guide_length <= CUT_OFFSET:
        raise ValueError(f"La guía debe tener más de {CUT_OFFSET} nucleótidos.")
    dna = validate_dna(dna)
    n = len(dna)
    targets = []

    # Hebra superior: protoespaciador en [i - L, i), PAM en [i, i + 3).
    for i in range(guide_length, n - PAM_LENGTH + 1):
        if _is_pam(dna[i:i + PAM_LENGTH]):
            targets.append(Target("+", dna[i - guide_length:i], dna[i:i + PAM_LENGTH],
                                  i - guide_length, i, i, i + PAM_LENGTH, i - CUT_OFFSET))

    # Hebra inferior: se busca igual en la complementaria inversa y se traducen
    # las posiciones. Una posición p de la inversa es la n - 1 - p de la superior,
    # y un límite entre bases c de la inversa es el límite n - c de la superior.
    rc = reverse_complement(dna)
    for i in range(guide_length, n - PAM_LENGTH + 1):
        if _is_pam(rc[i:i + PAM_LENGTH]):
            targets.append(Target("-", rc[i - guide_length:i], rc[i:i + PAM_LENGTH],
                                  n - i, n - i + guide_length,
                                  n - i - PAM_LENGTH, n - i, n - i + CUT_OFFSET))

    return sorted(targets, key=lambda t: (t.cut, t.strand))


@dataclass(frozen=True)
class Edit:
    """Un cambio en la secuencia: en 'position' se quitan 'deleted' y se ponen 'inserted'."""

    original: str
    position: int
    deleted: str
    inserted: str

    @property
    def edited(self) -> str:
        end = self.position + len(self.deleted)
        return self.original[:self.position] + self.inserted + self.original[end:]

    @property
    def length_change(self) -> int:
        return len(self.inserted) - len(self.deleted)

    @property
    def frameshift(self) -> bool:
        """True si desplaza el marco de lectura (el cambio no es múltiplo de 3)."""
        return self.length_change % 3 != 0

    @property
    def kind(self) -> str:
        """'insertion', 'deletion' o 'substitution'."""
        if not self.deleted:
            return "insertion"
        if not self.inserted:
            return "deletion"
        return "substitution"


def repair_nhej(dna: str, cut: int, rng: random.Random | None = None) -> Edit:
    """Reparación imprecisa: inserción de 1 base o deleción de 1 a 10 bases en el corte.

    Es una simplificación: en la célula los resultados más frecuentes son
    precisamente inserciones de 1 base y deleciones pequeñas.
    """
    dna = validate_dna(dna)
    _check_cut(dna, cut)
    rng = rng or random.Random()

    if rng.random() < 0.5:
        return Edit(dna, cut, "", rng.choice("ACGT"))

    size = rng.randint(1, min(MAX_NHEJ_DELETION, len(dna) - 1))
    left = rng.randint(0, size)                 # cuántas bases se pierden a la izquierda
    start = min(max(0, cut - left), len(dna) - size)
    return Edit(dna, start, dna[start:start + size], "")


def repair_hdr(dna: str, cut: int, insert: str) -> Edit:
    """Reparación con molde: inserta 'insert' exactamente en el punto de corte."""
    dna, insert = validate_dna(dna), validate_dna(insert)
    _check_cut(dna, cut)
    return Edit(dna, cut, "", insert)


def donor_template(dna: str, cut: int, insert: str, arm: int = HOMOLOGY_ARM) -> tuple[str, str, str]:
    """Molde para HDR: (brazo izquierdo, inserto, brazo derecho).
    Los brazos copian el ADN que rodea al corte para que la célula lo reconozca."""
    dna, insert = validate_dna(dna), validate_dna(insert)
    _check_cut(dna, cut)
    return dna[max(0, cut - arm):cut], insert, dna[cut:cut + arm]


def _check_cut(dna: str, cut: int) -> None:
    if not 0 < cut < len(dna):
        raise InvalidSequenceError(f"El punto de corte {cut} está fuera de la secuencia.")
