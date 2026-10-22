"""'Generar un ser a partir de bits' (opción 5 del submenú '*').

Con 2 bits por base, cualquier cadena de bits es una secuencia de ADN,
y las operaciones sobre bits se corresponden con operaciones biológicas:
    flip(bits)                        -> cadena complementaria
    reverse_blocks(bits)              -> secuencia leída al revés
    flip(reverse_blocks(bits))        -> complementaria inversa (5'->3')
"""

import random
import time

from .constants import BITS_PER_BASE
from .sequences import from_binary


def generate_bits(n: int, rng: random.Random | None = None) -> tuple[str, float]:
    """Genera n bits aleatorios. Devuelve la cadena y los segundos que ha tardado."""
    if n <= 0:
        raise ValueError("El número de bits debe ser mayor que 0.")
    rng = rng or random.Random()
    start = time.perf_counter()
    bits = "".join(rng.choice("01") for _ in range(n))
    return bits, time.perf_counter() - start


def generate_being(n_bases: int, rng: random.Random | None = None) -> tuple[str, str, float]:
    """Genera un ADN aleatorio de n_bases. Devuelve (bits, adn, segundos)."""
    bits, elapsed = generate_bits(n_bases * BITS_PER_BASE, rng)
    return bits, from_binary(bits), elapsed


def flip(bits: str) -> str:
    """Cambia ceros por unos y viceversa."""
    return bits.translate(str.maketrans("01", "10"))


def reverse_blocks(bits: str, size: int = BITS_PER_BASE) -> str:
    """Invierte el orden de los bloques sin desordenar los bits de dentro.

    reverse_blocks("0001") == "0100"   (A|C -> C|A)
    Ojo: invertir bit a bit ("0001"[::-1] == "1000") daría G|A, porque
    también da la vuelta a cada base y convierte C=01 en G=10.
    """
    return "".join(reversed(split_blocks(bits, size)))


def split_blocks(bits: str, size: int = BITS_PER_BASE) -> list[str]:
    return [bits[i:i + size] for i in range(0, len(bits), size)]
