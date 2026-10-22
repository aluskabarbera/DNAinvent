"""Alineamiento global de dos secuencias (algoritmo de Needleman-Wunsch, 1970).

Idea: buscar la forma de colocar una secuencia debajo de la otra, insertando
huecos ('-') donde haga falta, para que se parezcan lo máximo posible.

    Sin alinear:   ACGTACGTAC        Alineadas:   -ACGTACGTAC
                   TACGTACGTAC                    TACGTACGTAC

Cada posición del alineamiento suma o resta puntos:
    coincidencia (match)       misma base arriba y abajo   -> +1
    sustitución (mismatch)     bases distintas              -> -1
    hueco (gap)                inserción o deleción         -> -2

Probar todas las colocaciones posibles sería inabarcable, así que se usa
programación dinámica: se rellena una tabla donde cada casilla guarda la
mejor puntuación posible para alinear los primeros i caracteres de una
secuencia con los primeros j de la otra, reutilizando las casillas ya
calculadas. Coste: proporcional a len(a) * len(b).
"""

from dataclasses import dataclass

from .validation import validate_dna

GAP = "-"


@dataclass(frozen=True)
class Scoring:
    match: int = 1
    mismatch: int = -1
    gap: int = -2


@dataclass(frozen=True)
class Alignment:
    """Resultado de un alineamiento. Las dos cadenas tienen la misma longitud."""

    first: str
    second: str
    score: int

    @property
    def length(self) -> int:
        return len(self.first)

    @property
    def matches(self) -> int:
        return sum(1 for a, b in zip(self.first, self.second) if a == b)

    @property
    def gaps(self) -> int:
        return sum(1 for a, b in zip(self.first, self.second) if GAP in (a, b))

    @property
    def mismatches(self) -> int:
        return self.length - self.matches - self.gaps

    @property
    def identity(self) -> float:
        """Porcentaje de posiciones idénticas sobre la longitud del alineamiento."""
        return self.matches / self.length * 100

    @property
    def markup(self) -> str:
        """Línea intermedia: '|' coincidencia, '.' sustitución, ' ' hueco."""
        return "".join(
            " " if GAP in (a, b) else "|" if a == b else "."
            for a, b in zip(self.first, self.second)
        )


def align(first: str, second: str, scoring: Scoring = Scoring()) -> Alignment:
    """Alineamiento global óptimo de dos secuencias de ADN."""
    a, b = validate_dna(first), validate_dna(second)
    n, m = len(a), len(b)

    def pair_score(x: str, y: str) -> int:
        return scoring.match if x == y else scoring.mismatch

    # 1) Tabla de puntuaciones. La fila 0 y la columna 0 representan alinear
    #    con una secuencia vacía, es decir, poner solo huecos.
    table = [[0] * (m + 1) for _ in range(n + 1)]
    for i in range(1, n + 1):
        table[i][0] = i * scoring.gap
    for j in range(1, m + 1):
        table[0][j] = j * scoring.gap

    # 2) Cada casilla elige la mejor de tres opciones:
    #    diagonal  -> emparejar a[i-1] con b[j-1]
    #    arriba    -> a[i-1] frente a un hueco
    #    izquierda -> un hueco frente a b[j-1]
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            table[i][j] = max(
                table[i - 1][j - 1] + pair_score(a[i - 1], b[j - 1]),
                table[i - 1][j] + scoring.gap,
                table[i][j - 1] + scoring.gap,
            )

    # 3) Reconstrucción (traceback): desde la esquina final hacia atrás,
    #    averiguando de qué casilla vino cada puntuación.
    aligned_a, aligned_b = [], []
    i, j = n, m
    while i > 0 or j > 0:
        if i > 0 and j > 0 and table[i][j] == table[i - 1][j - 1] + pair_score(a[i - 1], b[j - 1]):
            aligned_a.append(a[i - 1])
            aligned_b.append(b[j - 1])
            i, j = i - 1, j - 1
        elif i > 0 and table[i][j] == table[i - 1][j] + scoring.gap:
            aligned_a.append(a[i - 1])
            aligned_b.append(GAP)
            i -= 1
        else:
            aligned_a.append(GAP)
            aligned_b.append(b[j - 1])
            j -= 1

    return Alignment(
        first="".join(reversed(aligned_a)),
        second="".join(reversed(aligned_b)),
        score=table[n][m],
    )
