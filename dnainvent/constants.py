"""Tablas de datos biológicos. Solo datos, sin lógica."""

DNA_BASES = frozenset("ACGT")
RNA_BASES = frozenset("ACGU")

# Pares de bases de Watson-Crick.
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

START_CODON = "AUG"
STOP_CODONS = frozenset({"UAA", "UAG", "UGA"})

# Código genético estándar: codón de ARN -> aminoácido (una letra). '_' = parada.
CODON_TABLE = {
    "AUA": "I", "AUC": "I", "AUU": "I", "AUG": "M",
    "ACA": "T", "ACC": "T", "ACG": "T", "ACU": "T",
    "AAC": "N", "AAU": "N", "AAA": "K", "AAG": "K",
    "AGC": "S", "AGU": "S", "AGA": "R", "AGG": "R",
    "CUA": "L", "CUC": "L", "CUG": "L", "CUU": "L",
    "CCA": "P", "CCC": "P", "CCG": "P", "CCU": "P",
    "CAC": "H", "CAU": "H", "CAA": "Q", "CAG": "Q",
    "CGA": "R", "CGC": "R", "CGG": "R", "CGU": "R",
    "GUA": "V", "GUC": "V", "GUG": "V", "GUU": "V",
    "GCA": "A", "GCC": "A", "GCG": "A", "GCU": "A",
    "GAC": "D", "GAU": "D", "GAA": "E", "GAG": "E",
    "GGA": "G", "GGC": "G", "GGG": "G", "GGU": "G",
    "UCA": "S", "UCC": "S", "UCG": "S", "UCU": "S",
    "UUC": "F", "UUU": "F", "UUA": "L", "UUG": "L",
    "UAC": "Y", "UAU": "Y", "UAA": "_", "UAG": "_",
    "UGC": "C", "UGU": "C", "UGA": "_", "UGG": "W",
}

AMINO_ACID_NAMES = {
    "A": "Alanine (Ala)", "R": "Arginine (Arg)", "N": "Asparagine (Asn)",
    "D": "Aspartate (Asp)", "C": "Cysteine (Cys)", "Q": "Glutamine (Gln)",
    "E": "Glutamate (Glu)", "G": "Glycine (Gly)", "H": "Histidine (His)",
    "I": "Isoleucine (Ile)", "L": "Leucine (Leu)", "K": "Lysine (Lys)",
    "M": "Methionine (Met)", "F": "Phenylalanine (Phe)", "P": "Proline (Pro)",
    "S": "Serine (Ser)", "T": "Threonine (Thr)", "W": "Tryptophan (Trp)",
    "Y": "Tyrosine (Tyr)", "V": "Valine (Val)",
}

# Codificación de 2 bits por base (4 bases = 2^2 combinaciones).
# El orden no es casual: cada base y su complementaria son opuestas bit a bit
# (A=00 <-> T=11, C=01 <-> G=10), así que invertir los bits da la complementaria.
BASE_TO_BINARY = {"A": "00", "C": "01", "G": "10", "T": "11"}
BINARY_TO_BASE = {bits: base for base, bits in BASE_TO_BINARY.items()}
BITS_PER_BASE = 2

BASE_NAMES = {"A": "Adenina", "C": "Citosina", "G": "Guanina", "T": "Timina"}
