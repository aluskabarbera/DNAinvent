"""DNAinvent: herramientas para trabajar con secuencias de ADN y ARN.

Este paquete solo contiene lógica: no pregunta nada ni imprime nada,
así cualquier interfaz (consola, ventana, web) puede usarlo.
"""

from .validation import InvalidSequenceError

__version__ = "2.0.0-dev"
__all__ = ["InvalidSequenceError"]
