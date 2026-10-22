"""Limpieza y validación de secuencias. Todas las funciones de lógica pasan por aquí."""

from .constants import BITS_PER_BASE, DNA_BASES, RNA_BASES


class InvalidSequenceError(ValueError):
    """La secuencia está vacía o contiene caracteres no permitidos."""


def clean(sequence: str) -> str:
    """Quita espacios y saltos de línea y pasa a mayúsculas."""
    return "".join(sequence.split()).upper()


def _validate(sequence: str, alphabet: frozenset, kind: str) -> str:
    sequence = clean(sequence)
    if not sequence:
        raise InvalidSequenceError(f"La secuencia de {kind} está vacía.")
    invalid = sorted(set(sequence) - alphabet)
    if invalid:
        raise InvalidSequenceError(
            f"Caracteres no válidos para {kind}: {', '.join(invalid)} "
            f"(permitidos: {', '.join(sorted(alphabet))})."
        )
    return sequence


def validate_dna(sequence: str) -> str:
    return _validate(sequence, DNA_BASES, "ADN")


def validate_rna(sequence: str) -> str:
    return _validate(sequence, RNA_BASES, "ARN")


def validate_bits(bits: str) -> str:
    """Cadena de bits que se pueda leer de 2 en 2 como bases."""
    bits = _validate(bits, frozenset("01"), "binario")
    if len(bits) % BITS_PER_BASE:
        raise InvalidSequenceError(
            f"La cadena binaria tiene {len(bits)} bits; debe ser múltiplo de {BITS_PER_BASE}."
        )
    return bits
