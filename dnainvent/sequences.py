"""Operaciones básicas sobre secuencias (opciones 1, 2, 3 y 5 del menú principal)."""

from collections import Counter
from dataclasses import dataclass

from .constants import (AMINO_ACID_NAMES, BASE_TO_BINARY, BINARY_TO_BASE, BITS_PER_BASE,
                        CODON_TABLE, COMPLEMENT, START_CODON, STOP_CODONS)
from .validation import validate_bits, validate_dna, validate_rna


def complement(dna: str) -> str:
    """Cadena complementaria base a base, escrita 3'->5'.

    Sirve para dibujar la doble hélice alineada: cada base queda debajo de su pareja.
        5'-ATGC-3'
        3'-TACG-5'   <- complement("ATGC")
    """
    dna = validate_dna(dna)
    return "".join(COMPLEMENT[base] for base in dna)


def reverse_complement(dna: str) -> str:
    """Cadena complementaria escrita 5'->3', que es la convención estándar.

    Las dos hebras del ADN son antiparalelas: la complementaria se lee en sentido
    contrario, así que hay que invertirla para escribirla también de 5' a 3'.
        reverse_complement("ATGC") == "GCAT"
    """
    return complement(dna)[::-1]


def transcribe(dna: str) -> str:
    """ADN (cadena codificante) -> ARN."""
    return validate_dna(dna).replace("T", "U")


@dataclass(frozen=True)
class Translation:
    """Resultado de traducir un ARN mensajero.

    protein     Aminoácidos en código de una letra, empezando por la metionina (M).
    start       Posición (desde 0) del codón de inicio AUG, o None si no hay.
    end         Posición justo después del último codón leído (útil para resaltar).
    stop_codon  Codón de parada encontrado (UAA, UAG o UGA), o None si no hay.
    """

    protein: str
    start: int | None
    end: int | None
    stop_codon: str | None

    @property
    def has_start(self) -> bool:
        return self.start is not None

    @property
    def is_complete(self) -> bool:
        """True si la traducción terminó en un codón de parada."""
        return self.stop_codon is not None


def translate(rna: str) -> Translation:
    """Traduce ARN a proteína como lo hace el ribosoma:
    empieza en el primer AUG, lee de 3 en 3 y se detiene en el primer codón de parada.
    Si no hay codón de parada, traduce hasta el último codón completo."""
    rna = validate_rna(rna)
    start = rna.find(START_CODON)
    if start == -1:
        return Translation(protein="", start=None, end=None, stop_codon=None)

    protein = []
    position = start
    while position + 3 <= len(rna):
        codon = rna[position:position + 3]
        position += 3
        if codon in STOP_CODONS:
            return Translation("".join(protein), start, position, codon)
        protein.append(CODON_TABLE[codon])
    return Translation("".join(protein), start, position, None)


def count_amino_acids(protein: str) -> Counter:
    return Counter(protein)


def amino_acid_names(protein: str) -> list[str]:
    return [AMINO_ACID_NAMES[aa] for aa in protein]


def to_binary(dna: str) -> str:
    """ADN -> binario, 2 bits por base: to_binary("ACGT") == "00011011"."""
    dna = validate_dna(dna)
    return "".join(BASE_TO_BINARY[base] for base in dna)


def from_binary(bits: str) -> str:
    """Binario -> ADN, leyendo de 2 en 2: from_binary("00011011") == "ACGT"."""
    bits = validate_bits(bits)
    return "".join(BINARY_TO_BASE[bits[i:i + BITS_PER_BASE]]
                   for i in range(0, len(bits), BITS_PER_BASE))
