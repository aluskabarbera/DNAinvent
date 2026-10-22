"""Herramientas del submenú '*': buscar, similitud y permutaciones."""

import math
from collections import Counter
from collections.abc import Iterator

from .alignment import align
from .validation import validate_dna


def find_positions(dna: str, target: str) -> list[int]:
    """Posiciones (desde 0) donde aparece target dentro de dna, sin solapamiento."""
    dna, target = validate_dna(dna), validate_dna(target)
    positions, start = [], dna.find(target)
    while start != -1:
        positions.append(start)
        start = dna.find(target, start + len(target))
    return positions


def similarity(first: str, second: str) -> float:
    """Porcentaje de identidad tras alinear las dos secuencias."""
    return align(first, second).identity


def count_positional_permutations(dna: str, r: int) -> int:
    """Ordenaciones de r posiciones, tratando cada base como distinta aunque se repita.
    Es lo que cuenta itertools.permutations: n! / (n - r)!"""
    _check_r(r)
    return math.perm(len(validate_dna(dna)), r)


def count_permutations(dna: str, r: int) -> int:
    """Cuántas secuencias DISTINTAS de r bases se pueden formar, sin generarlas.

    Se añade una base cada vez. Si ya hay 'ways[j]' secuencias de longitud j,
    añadir k copias de la base nueva da secuencias de longitud j + k, y hay
    C(j + k, k) formas de elegir en qué posiciones van esas k copias.
    Ejemplo: "AAC", r = 3 -> 3 (AAC, ACA, CAA), mientras que 3! = 6.
    """
    _check_r(r)
    counts = Counter(validate_dna(dna))
    ways = [1] + [0] * r  # ways[j]: secuencias distintas de longitud j
    for copies in counts.values():
        new_ways = [0] * (r + 1)
        for j, current in enumerate(ways):
            if current:
                for k in range(min(copies, r - j) + 1):
                    new_ways[j + k] += current * math.comb(j + k, k)
        ways = new_ways
    return ways[r]


def iter_permutations(dna: str, r: int) -> Iterator[str]:
    """Genera las permutaciones distintas de r bases, en orden alfabético.

    En lugar de ordenar posiciones (que repite resultados cuando hay bases
    iguales), en cada paso elige qué base poner entre las que aún quedan.
    Como cada base se prueba una sola vez por paso, no salen duplicados.
    """
    _check_r(r)
    remaining = Counter(validate_dna(dna))
    bases = sorted(remaining)
    current: list[str] = []

    def build() -> Iterator[str]:
        if len(current) == r:
            yield "".join(current)
            return
        for base in bases:
            if remaining[base]:
                remaining[base] -= 1
                current.append(base)
                yield from build()
                current.pop()
                remaining[base] += 1

    if r <= sum(remaining.values()):
        yield from build()


def _check_r(r: int) -> None:
    if r < 0:
        raise ValueError("El tamaño de la permutación no puede ser negativo.")
