import itertools
import math
import random

import pytest

from dnainvent import InvalidSequenceError, alignment, bits, crispr, sequences, tools


def test_complement():
    assert sequences.complement("atgc") == "TACG"


def test_reverse_complement():
    assert sequences.reverse_complement("ATGC") == "GCAT"
    assert sequences.reverse_complement("AACG") == "CGTT"


def test_reverse_complement_twice_returns_original():
    dna = "ATGCGTTAGC"
    assert sequences.reverse_complement(sequences.reverse_complement(dna)) == dna


def test_palindromic_site():
    # EcoRI (GAATTC) es su propia complementaria inversa.
    assert sequences.reverse_complement("GAATTC") == "GAATTC"


def test_transcribe():
    assert sequences.transcribe("ATGC") == "AUGC"


def test_translate_starts_at_aug_and_stops():
    # Basura antes del AUG, y bases después de la parada que no se deben leer.
    result = sequences.translate("GGAUGUUUGGCUAAAUGCCC")
    assert result.protein == "MFG"
    assert result.start == 2
    assert result.stop_codon == "UAA"
    assert result.end == 14
    assert result.is_complete


def test_translate_without_start():
    result = sequences.translate("UUUCCCGGG")
    assert result.protein == ""
    assert not result.has_start


def test_translate_without_stop_is_incomplete():
    result = sequences.translate("AUGUUUGG")  # 'GG' sobrante se ignora
    assert result.protein == "MF"
    assert not result.is_complete


def test_translate_immediate_stop():
    assert sequences.translate("AUGUGA").protein == "M"


def test_translate_rejects_dna():
    with pytest.raises(InvalidSequenceError):
        sequences.translate("ATGTTT")


def test_amino_acid_names():
    assert sequences.amino_acid_names("MQ") == ["Methionine (Met)", "Glutamine (Gln)"]


def test_invalid_characters_raise():
    with pytest.raises(InvalidSequenceError):
        sequences.complement("ATGX")


def test_empty_sequence_raises():
    with pytest.raises(InvalidSequenceError):
        tools.find_positions("", "A")


def test_find_positions():
    assert tools.find_positions("ACGACG", "ACG") == [0, 3]




def test_count_permutations():
    assert tools.count_permutations("ACGT", 2) == 12  # sin repetidas: igual que math.perm


def test_permutations_without_duplicates():
    assert list(tools.iter_permutations("AAC", 3)) == ["AAC", "ACA", "CAA"]
    assert tools.count_permutations("AAC", 3) == 3
    assert tools.count_positional_permutations("AAC", 3) == 6


def test_partial_permutations_without_duplicates():
    assert list(tools.iter_permutations("AAC", 2)) == ["AA", "AC", "CA"]


def test_all_same_base():
    assert list(tools.iter_permutations("AAAA", 4)) == ["AAAA"]


def test_r_larger_than_sequence():
    assert tools.count_permutations("ACG", 5) == 0
    assert list(tools.iter_permutations("ACG", 5)) == []


def test_negative_r_raises():
    with pytest.raises(ValueError):
        tools.count_permutations("ACG", -1)


def test_count_is_fast_for_huge_numbers():
    # 50 A y 50 C: C(100, 50) secuencias distintas, imposible de generar.
    assert tools.count_permutations("A" * 50 + "C" * 50, 100) == math.comb(100, 50)


@pytest.mark.parametrize("seed", range(15))
def test_permutations_match_brute_force(seed):
    rng = random.Random(seed)
    dna = "".join(rng.choice("ACG") for _ in range(rng.randint(1, 7)))
    for r in range(len(dna) + 1):
        expected = sorted({"".join(p) for p in itertools.permutations(dna, r)})
        generated = list(tools.iter_permutations(dna, r))
        assert generated == expected           # mismas, sin duplicados y en orden
        assert tools.count_permutations(dna, r) == len(expected)


def test_to_binary():
    assert sequences.to_binary("acgt") == "00011011"


def test_from_binary():
    assert sequences.from_binary("00 01 10 11") == "ACGT"


def test_binary_round_trip():
    dna = "GATTACA"
    assert sequences.from_binary(sequences.to_binary(dna)) == dna


def test_from_binary_rejects_odd_length():
    with pytest.raises(InvalidSequenceError):
        sequences.from_binary("001")


def test_from_binary_rejects_other_digits():
    with pytest.raises(InvalidSequenceError):
        sequences.from_binary("0012")


def test_reverse_blocks_keeps_bases_intact():
    assert bits.reverse_blocks("0001") == "0100"  # AC -> CA, no GA


@pytest.mark.parametrize("seed", range(20))
def test_bit_operations_match_biology(seed):
    rng = random.Random(seed)
    bit_string, dna, _ = bits.generate_being(30, rng)
    assert sequences.from_binary(bits.flip(bit_string)) == sequences.complement(dna)
    assert sequences.from_binary(bits.reverse_blocks(bit_string)) == dna[::-1]
    assert (sequences.from_binary(bits.flip(bits.reverse_blocks(bit_string)))
            == sequences.reverse_complement(dna))


# ---------- Alineamiento ----------

def test_align_identical():
    result = alignment.align("GATTACA", "GATTACA")
    assert result.identity == 100
    assert result.gaps == 0


def test_align_handles_insertion():
    # El caso que antes daba 0 %: la misma secuencia con una T delante.
    result = alignment.align("ACGTACGTAC", "TACGTACGTAC")
    assert result.first == "-ACGTACGTAC"
    assert result.second == "TACGTACGTAC"
    assert result.gaps == 1
    assert result.matches == 10
    assert round(result.identity, 2) == 90.91


def test_align_substitution():
    result = alignment.align("ACGT", "AGGT")
    assert result.markup == "|.||"
    assert result.mismatches == 1
    assert result.score == 3 - 1


def test_aligned_strings_have_same_length_and_original_bases():
    first, second = "GATTACAGATTACA", "GCATGCTTACA"
    result = alignment.align(first, second)
    assert len(result.first) == len(result.second)
    assert result.first.replace("-", "") == first
    assert result.second.replace("-", "") == second


def test_align_score_is_symmetric():
    a, b = "ACGGTAC", "ACTTAC"
    assert alignment.align(a, b).score == alignment.align(b, a).score


def test_counts_add_up():
    r = alignment.align("AAGCTTGA", "AGCTAGGA")
    assert r.matches + r.mismatches + r.gaps == r.length


def test_custom_scoring_changes_result():
    # Si los huecos salen muy caros, prefiere sustituciones.
    expensive = alignment.align("ACGT", "AGT", alignment.Scoring(gap=-10))
    assert expensive.gaps == 1  # la diferencia de longitud obliga a 1 hueco
    cheap = alignment.align("AAAA", "TTTT", alignment.Scoring(mismatch=-5, gap=-1))
    assert cheap.mismatches == 0  # con huecos baratos, evita sustituciones


def test_similarity_uses_alignment():
    assert round(tools.similarity("ACGTACGTAC", "TACGTACGTAC"), 2) == 90.91


# ---------- CRISPR/Cas9 ----------

GUIDE = "GACGTTACCGATGCATGCAT"          # 20 nt
PLUS_DNA = GUIDE + "TGG" + "AAAA"       # diana en la hebra superior


def test_find_target_on_plus_strand():
    [target] = crispr.find_targets(PLUS_DNA)
    assert target.strand == "+"
    assert target.protospacer == GUIDE
    assert target.pam == "TGG"
    assert target.cut == 17             # 3 pb antes de la PAM (que empieza en 20)
    assert target.guide_rna == GUIDE.replace("T", "U")


def test_find_target_on_minus_strand():
    # La misma diana, pero escrita en la otra hebra: la PAM aparece como CCA.
    dna = sequences.reverse_complement(PLUS_DNA)
    [target] = crispr.find_targets(dna)
    assert target.strand == "-"
    assert target.protospacer == GUIDE
    assert target.pam == "TGG"
    assert target.cut == len(dna) - 17
    assert dna[target.pam_start:target.pam_end] == "CCA"


def test_no_targets_in_short_sequence():
    assert crispr.find_targets("ACGTTGG") == []


@pytest.mark.parametrize("seed", range(10))
def test_targets_are_consistent(seed):
    rng = random.Random(seed)
    dna = "".join(rng.choice("ACGT") for _ in range(120))
    rc = sequences.reverse_complement(dna)
    for t in crispr.find_targets(dna):
        assert t.pam[1:] == "GG"
        assert len(t.protospacer) == 20
        region = dna[t.start:t.end]
        assert region == (t.protospacer if t.strand == "+" else sequences.reverse_complement(t.protospacer))
        # El corte está dentro del protoespaciador, a 3 pb de la PAM.
        assert t.start < t.cut < t.end
        assert (t.pam_start - t.cut if t.strand == "+" else t.cut - t.pam_end) == 3


@pytest.mark.parametrize("seed", range(10))
def test_targets_mirror_between_strands(seed):
    # Buscar en una secuencia o en su complementaria inversa da las mismas
    # dianas, con la hebra cambiada y las posiciones reflejadas.
    rng = random.Random(seed)
    dna = "".join(rng.choice("ACGT") for _ in range(100))
    n = len(dna)
    direct = {(t.strand, t.protospacer, t.cut) for t in crispr.find_targets(dna)}
    mirrored = {("+" if t.strand == "-" else "-", t.protospacer, n - t.cut)
                for t in crispr.find_targets(sequences.reverse_complement(dna))}
    assert direct == mirrored


@pytest.mark.parametrize("seed", range(30))
def test_nhej_changes_sequence_at_cut(seed):
    edit = crispr.repair_nhej(PLUS_DNA, 17, random.Random(seed))
    assert edit.edited != PLUS_DNA
    assert len(edit.edited) == len(PLUS_DNA) + edit.length_change
    assert edit.position <= 17 <= edit.position + len(edit.deleted)


def test_hdr_inserts_at_cut():
    edit = crispr.repair_hdr(PLUS_DNA, 17, "ccc")
    assert edit.edited == PLUS_DNA[:17] + "CCC" + PLUS_DNA[17:]
    assert edit.kind == "insertion"
    assert not edit.frameshift


def test_frameshift_detection():
    assert crispr.Edit("ACGTACGT", 4, "", "A").frameshift
    assert crispr.Edit("ACGTACGT", 2, "GT", "").frameshift
    assert not crispr.Edit("ACGTACGT", 2, "GTA", "").frameshift


def test_donor_template_has_homology_arms():
    left, insert, right = crispr.donor_template(PLUS_DNA, 17, "CCC", arm=5)
    assert (left, insert, right) == (PLUS_DNA[12:17], "CCC", PLUS_DNA[17:22])


def test_cut_outside_sequence_raises():
    with pytest.raises(InvalidSequenceError):
        crispr.repair_hdr("ACGT", 10, "A")
