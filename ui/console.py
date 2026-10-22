"""Interfaz de consola. Aquí (y solo aquí) hay input(), print() y colores."""

import sys

from colorama import Fore, Style, init

from dnainvent import InvalidSequenceError, alignment, bits, crispr, sequences, tools
from dnainvent.constants import AMINO_ACID_NAMES
from dnainvent.validation import validate_bits, validate_dna, validate_rna

BASE_COLORS = {"A": Fore.YELLOW, "C": Fore.BLUE, "G": Fore.GREEN, "T": Fore.RED}
PERMUTATION_WARNING = 10_000
ALIGNMENT_WARNING = 4_000_000   # unas 2.000 x 2.000 bases
LINE_WIDTH = 60


# ---------- Utilidades de entrada y salida ----------

def ask_sequence(prompt: str, validator) -> str:
    """Pide una secuencia hasta que sea válida."""
    while True:
        try:
            seq = validator(input(prompt))
            print(f"Esta secuencia tiene {len(seq)} caracteres.")
            return seq
        except InvalidSequenceError as error:
            print(Fore.RED + str(error) + Style.RESET_ALL)


def ask_dna(prompt: str = "Introduce la cadena de ADN (A, C, G, T): ") -> str:
    return ask_sequence(prompt, validate_dna)


def ask_rna(prompt: str = "Introduce la cadena de ARN (A, C, G, U): ") -> str:
    def validator(text: str) -> str:
        try:
            return validate_rna(text)
        except InvalidSequenceError as error:
            if "T" in text.upper():
                raise InvalidSequenceError(
                    f"{error} ¿Es ADN? Transcríbelo antes con la opción 2.") from error
            raise
    return ask_sequence(prompt, validator)


def ask_int(prompt: str, minimum: int = 0, maximum: int | None = None,
            default: int | None = None) -> int:
    limits = f"entre {minimum} y {maximum}" if maximum is not None else f"mayor o igual que {minimum}"
    while True:
        text = input(prompt).strip()
        if not text and default is not None:
            return default
        try:
            value = int(text)
            if value >= minimum and (maximum is None or value <= maximum):
                return value
        except ValueError:
            pass
        print(Fore.RED + f"Escribe un número entero {limits}." + Style.RESET_ALL)


def highlight(seq: str, positions: list[int], length: int) -> str:
    """Resalta en cian los fragmentos que empiezan en positions."""
    out, last = [], 0
    for pos in positions:
        out.append(seq[last:pos])
        out.append(Fore.CYAN + seq[pos:pos + length] + Style.RESET_ALL)
        last = pos + length
    out.append(seq[last:])
    return "".join(out)


def color_dna(dna: str) -> str:
    return "".join(BASE_COLORS[b] + b for b in dna) + Style.RESET_ALL


def color_bits(bit_string: str) -> str:
    """Bits agrupados de 2 en 2 con el color de la base que representan."""
    return " ".join(BASE_COLORS[base] + block
                    for block, base in zip(bits.split_blocks(bit_string),
                                           sequences.from_binary(bit_string))) + Style.RESET_ALL


# ---------- Menú principal ----------

def option_complement():
    dna = ask_dna()
    print("\nDoble hélice (hebras antiparalelas):")
    print(f"  5'-{dna}-3'")
    print(f"     {'|' * len(dna)}")
    print(f"  3'-{sequences.complement(dna)}-5'")
    print("\nComplementaria inversa (5'->3'):")
    print(Fore.CYAN + f"  5'-{sequences.reverse_complement(dna)}-3'" + Style.RESET_ALL)


def option_transcribe():
    print(sequences.transcribe(ask_dna()))


def option_translate():
    rna = ask_rna()
    result = sequences.translate(rna)

    if not result.has_start:
        print(Fore.YELLOW + "No hay codón de inicio (AUG): el ribosoma no empezaría a traducir." + Style.RESET_ALL)
        return

    # Marco de lectura en cian y codón de parada en rojo.
    stop_len = 3 if result.is_complete else 0
    orf_end = result.end - stop_len
    print(rna[:result.start]
          + Fore.CYAN + rna[result.start:orf_end]
          + Fore.RED + rna[orf_end:result.end] + Style.RESET_ALL
          + rna[result.end:])

    print(f"\nInicio en la posición {result.start + 1}.")
    if result.is_complete:
        print(f"Parada en el codón {result.stop_codon} (posición {orf_end + 1}).")
    else:
        print(Fore.YELLOW + "Sin codón de parada: la secuencia podría estar incompleta." + Style.RESET_ALL)

    print(f"\nProteína ({len(result.protein)} aminoácidos): {result.protein}")
    for aa, count in sequences.count_amino_acids(result.protein).most_common():
        print(f"  {AMINO_ACID_NAMES[aa]:<22} x{count}")


def paint_target(dna: str, target: "crispr.Target") -> str:
    """Protoespaciador en cian, PAM en magenta y el corte marcado con una barra roja."""
    out = []
    for index, base in enumerate(dna):
        if index == target.cut:
            out.append(Fore.RED + Style.BRIGHT + "|" + Style.RESET_ALL)
        if target.start <= index < target.end:
            out.append(Fore.CYAN + base)
        elif target.pam_start <= index < target.pam_end:
            out.append(Fore.MAGENTA + base)
        else:
            out.append(Style.RESET_ALL + base)
    return "".join(out) + Style.RESET_ALL


def show_protein_effect(before: str, after: str):
    """Traduce la hebra antes y después de editar para ver el efecto en la proteína."""
    old = sequences.translate(sequences.transcribe(before))
    new = sequences.translate(sequences.transcribe(after))
    if not old.has_start and not new.has_start:
        print("  (No hay codón de inicio ATG: esta secuencia no codifica proteína.)")
        return
    print(f"  Antes:   {old.protein or '-'}")
    print(f"  Después: {new.protein or '-'}")
    if old.protein == new.protein:
        print("  La proteína no cambia.")


def option_crispr():
    dna = ask_dna()
    guide_length = ask_int(f"Longitud de la guía (Enter = {crispr.GUIDE_LENGTH}): ",
                           minimum=crispr.CUT_OFFSET + 1, default=crispr.GUIDE_LENGTH)

    # 1) Buscar dianas
    targets = crispr.find_targets(dna, guide_length)
    if not targets:
        print(Fore.YELLOW + f"No hay ninguna PAM (NGG) con {guide_length} bases delante en "
              "ninguna de las dos hebras. Prueba una secuencia más larga o una guía más corta."
              + Style.RESET_ALL)
        return

    print(f"\nDianas encontradas: {len(targets)}")
    for number, t in enumerate(targets, start=1):
        print(f"  {number:>2}] hebra {t.strand}  corte tras la base {t.cut:>4}  "
              f"guía 5'-{t.guide_rna}-3'  PAM {t.pam}")
    choice = ask_int("Elige una diana: ", minimum=1, maximum=len(targets))
    target = targets[choice - 1]

    # 2) Cortar
    print("\nCas9 corta aquí (cian: diana, magenta: PAM, |: corte):")
    print("  " + paint_target(dna, target))

    # 3) Reparar
    print("\n¿Cómo repara la célula el corte?")
    print("  1] NHEJ: pega los extremos y suele equivocarse (inactiva el gen)")
    print("  2] HDR: copia un ADN molde que aportamos (inserta una secuencia)")
    repair = ask_int("Elige una opción: ", minimum=1, maximum=2)

    if repair == 1:
        edit = crispr.repair_nhej(dna, target.cut)
    else:
        insert = ask_dna("Secuencia que quieres insertar: ")
        edit = crispr.repair_hdr(dna, target.cut, insert)
        left, _, right = crispr.donor_template(dna, target.cut, insert)
        print("\nMolde de ADN donante (brazos de homología en gris, inserto en verde):")
        print("  " + Style.DIM + left + Style.RESET_ALL + Fore.GREEN + insert
              + Style.RESET_ALL + Style.DIM + right + Style.RESET_ALL)

    # Resultado
    # Posiciones contadas desde 1, como en biología.
    first, last = edit.position + 1, edit.position + len(edit.deleted)
    change = {
        "insertion": f"inserción de {len(edit.inserted)} base(s) ({edit.inserted}) "
                     f"entre las posiciones {edit.position} y {edit.position + 1}",
        "deletion": f"deleción de {len(edit.deleted)} base(s) ({edit.deleted}) "
                    f"en las posiciones {first}-{last}",
        "substitution": f"sustitución de {edit.deleted} por {edit.inserted} "
                        f"en las posiciones {first}-{last}",
    }[edit.kind]
    print(f"\nResultado: {change}.")
    show_alignment(alignment.align(edit.original, edit.edited))
    if edit.frameshift:
        print(Fore.YELLOW + "El cambio no es múltiplo de 3: desplaza el marco de lectura." + Style.RESET_ALL)
    else:
        print("El cambio es múltiplo de 3: el marco de lectura se conserva.")
    print("\nEfecto en la proteína:")
    show_protein_effect(edit.original, edit.edited)


def option_binary():
    dna = ask_dna()
    bit_string = sequences.to_binary(dna)
    print(color_bits(bit_string))
    print(f"\n{len(dna)} bases -> {len(bit_string)} bits ({len(bit_string) / 8:g} bytes).")
    print(f"Como texto normal (8 bits por letra) ocuparía {len(dna) * 8} bits: 4 veces más.")


def option_from_binary():
    bit_string = ask_sequence("Introduce la cadena binaria (0 y 1, 2 bits por base): ", validate_bits)
    print(color_dna(sequences.from_binary(bit_string)))


# ---------- Submenú de herramientas ----------

def tool_find():
    dna = ask_dna()
    target = ask_dna("Secuencia que quieres buscar: ")
    positions = tools.find_positions(dna, target)
    print(highlight(dna, positions, len(target)))
    print(f"Encontrada {len(positions)} veces en las posiciones: {[p + 1 for p in positions]}")


def tool_similarity():
    first = ask_dna("Primera cadena de ADN: ")
    second = ask_dna("Segunda cadena de ADN: ")
    cells = len(first) * len(second)
    if cells > ALIGNMENT_WARNING and input(
            f"El alineamiento necesita {cells:,} cálculos y puede tardar. ¿Continuar? (s/n): ").lower() != "s":
        return

    result = alignment.align(first, second)
    show_alignment(result)
    print(f"Identidad:     {result.identity:.2f} % ({result.matches}/{result.length})")
    print(f"Sustituciones: {result.mismatches}")
    print(f"Huecos:        {result.gaps}")
    print(f"Puntuación:    {result.score}")


def show_alignment(result: "alignment.Alignment"):
    print()
    for start in range(0, result.length, LINE_WIDTH):
        end = start + LINE_WIDTH
        top, mid, bottom = result.first[start:end], result.markup[start:end], result.second[start:end]
        print("  " + "".join(color_aligned(a, b) for a, b in zip(top, bottom)) + Style.RESET_ALL)
        print("  " + mid)
        print("  " + "".join(color_aligned(b, a) for a, b in zip(top, bottom)) + Style.RESET_ALL)
        print()


def color_aligned(base: str, other: str) -> str:
    """Verde si coincide, amarillo si es una sustitución, rojo si hay hueco."""
    if alignment.GAP in (base, other):
        return Fore.RED + base
    return (Fore.GREEN if base == other else Fore.YELLOW) + base


def tool_permutations():
    dna = ask_dna()
    r = ask_int("¿De cuántas bases quieres cada permutación?: ", minimum=1)
    if r > len(dna):
        print(f"No se pueden formar grupos de {r} bases con solo {len(dna)}.")
        return
    total = tools.count_permutations(dna, r)
    positional = tools.count_positional_permutations(dna, r)
    print(f"Ordenaciones de posiciones: {positional:,}")
    print(f"Secuencias distintas:       {total:,}")
    if positional > total:
        print(f"Se evitan {positional - total:,} repeticiones por las bases iguales.")
    if total > PERMUTATION_WARNING and input("Son muchas. ¿Mostrarlas igualmente? (s/n): ").lower() != "s":
        return
    for p in tools.iter_permutations(dna, r):
        print(p)


def show_being(title: str, bit_string: str):
    print(f"\n{title}")
    print("  Bits: " + color_bits(bit_string))
    print("  ADN:  " + color_dna(sequences.from_binary(bit_string)))


def tool_bits():
    n = ask_int("¿Cuántas bases quieres que tenga el ser?: ", minimum=1)
    original, _, elapsed = bits.generate_being(n)
    current = original  # estado actual: cada operación parte de aquí y lo actualiza
    show_being("Ser generado:", current)

    def apply(title: str, transform) -> str:
        """Aplica transform al estado actual, lo guarda y devuelve el estado anterior."""
        nonlocal current
        previous, current = current, transform(current)
        show_being(title, current)
        return previous

    def reverse_complement_step():
        previous = apply("Invertir ceros y unos + orden de las bases (complementaria inversa):",
                         lambda b: bits.flip(bits.reverse_blocks(b)))
        check = sequences.from_binary(current) == sequences.reverse_complement(
            sequences.from_binary(previous))
        print(f"  ¿Coincide con reverse_complement() del paso anterior? {'Sí' if check else 'No'}")

    def restore():
        nonlocal current
        current = original
        show_being("Ser original restaurado:", current)

    actions = {
        "1": lambda: print(f"\nGenerados {len(original)} bits en {elapsed * 1000:.3f} ms "
                           f"({elapsed / len(original) * 1e6:.3f} µs de media por bit)."),
        "2": lambda: apply("Invertir ceros y unos (complementaria):", bits.flip),
        "3": lambda: apply("Invertir el orden de las bases:", bits.reverse_blocks),
        "4": reverse_complement_step,
        "5": restore,
    }
    menu(
        "Ser a partir de bits",
        {"1": "Tiempo de generación", "2": "Invertir ceros y unos",
         "3": "Invertir el orden de las bases", "4": "Las dos cosas a la vez",
         "5": "Volver al ser original"},
        actions,
        back_label="Volver",
    )


def tools_menu():
    menu(
        "DNAtools",
        {"1": "Buscar una secuencia", "2": "Alinear y comparar dos secuencias",
         "3": "Permutaciones", "4": "Generar un ser a partir de bits"},
        {"1": tool_find, "2": tool_similarity, "3": tool_permutations, "4": tool_bits},
        back_label="Volver",
    )


# ---------- Motor de menús ----------

def menu(title: str, labels: dict[str, str], actions: dict, back_label: str) -> None:
    """Muestra un menú hasta que el usuario elige 0."""
    while True:
        print(f"\n        {title}")
        for key, label in labels.items():
            print(f"    {key}] {label}")
        print(f"    0] {back_label}")
        choice = input("Elige una opción: ").strip()
        if choice == "0":
            return
        action = actions.get(choice)
        if action is None:
            print("Opción no válida, prueba otra vez.")
            continue
        action()


def run() -> None:
    init()
    try:
        menu(
            "DNAinvent",
            {"1": "Doble hélice y complementaria inversa", "2": "Transcripción ADN -> ARN",
             "3": "Traducción ARN -> aminoácidos", "4": "CRISPR/Cas9",
             "5": "ADN a binario", "6": "Binario a ADN",
             "*": "Herramientas"},
            {"1": option_complement, "2": option_transcribe, "3": option_translate,
             "4": option_crispr, "5": option_binary, "6": option_from_binary,
             "*": tools_menu},
            back_label="Salir",
        )
    except (KeyboardInterrupt, EOFError):
        print()
    print("¡Hasta la próxima!")
    sys.exit(0)
