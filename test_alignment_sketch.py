"""
Tests for alignment_sketch.py — built around real Arabidopsis IGS variant sequences.

  data/variants.fa            — 4 unaligned IGS variants (different lengths)
  data/variants_clustalo.fa   — Clustal Omega alignment of the same 4 sequences
"""

import sys
import pytest
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
from alignment_sketch import (
    read_alignment,
    compute_borders,
    compute_entropy,
    write_entropy_txt,
    draw_alignment_svg,
    draw_entropy_svg,
    draw_consensus_svg,
    draw_logo_svg,
)

DATA_DIR  = Path(__file__).parent / "data"
ALIGNED   = DATA_DIR / "variants_clustalo.fa"   # aligned — all same length
UNALIGNED = DATA_DIR / "variants.fa"            # raw sequences — different lengths


# ── shared fixtures (computed once per session) ───────────────────────────────

@pytest.fixture(scope="module")
def records():
    with open(ALIGNED) as f:
        return read_alignment(f)


@pytest.fixture(scope="module")
def sequences(records):
    return [str(r.seq) for r in records]


@pytest.fixture(scope="module")
def borders(sequences):
    return compute_borders(sequences)


@pytest.fixture(scope="module")
def entropy(sequences):
    return compute_entropy(sequences)


# ── read_alignment ────────────────────────────────────────────────────────────

def test_sequence_count(records):
    """Alignment must contain 4 IGS variants."""
    assert len(records) == 4


def test_all_sequences_same_length(records):
    """Clustal Omega output must be a proper alignment — all rows equal length."""
    lengths = {len(r.seq) for r in records}
    assert lengths == {490}


def test_rejects_unaligned_input():
    """Raw sequences with different lengths must cause SystemExit with a clear error."""
    with open(UNALIGNED) as f:
        with pytest.raises(SystemExit):
            read_alignment(f)


# ── compute_borders ───────────────────────────────────────────────────────────

def test_borders_count(borders):
    """Known number of conserved/variable transitions in this alignment."""
    assert len(borders) == 88


def test_borders_last_equals_alignment_length(borders):
    """Last border must be the alignment length (sentinel value)."""
    assert borders[-1] == 490


def test_borders_are_sorted(borders):
    """Borders must be strictly increasing."""
    assert borders == sorted(borders)


def test_borders_within_range(borders, sequences):
    """Every border index must be within [0, alignment_length] (0-indexed columns)."""
    aln_len = len(sequences[0])
    assert all(0 <= b <= aln_len for b in borders)


def test_borders_first(borders):
    """First transition is at column 17 in this alignment."""
    assert borders[0] == 17


# ── compute_entropy ───────────────────────────────────────────────────────────

def test_entropy_length(entropy, sequences):
    """One entropy value per alignment column."""
    assert len(entropy) == len(sequences[0])


def test_entropy_minimum_is_zero(entropy):
    """Fully conserved columns must have entropy exactly 0."""
    assert min(entropy) == 0.0


def test_entropy_maximum_within_dna_range(entropy):
    """DNA entropy can never exceed 2.0 bits (4 equally frequent bases)."""
    assert max(entropy) <= 2.0


def test_entropy_maximum_value(entropy):
    """Maximum entropy in this alignment is 1.5 bits (3 equally frequent bases)."""
    assert abs(max(entropy) - 1.5) < 0.001


def test_entropy_fully_conserved_column_count(entropy):
    """381 columns are fully conserved among non-gap sequences (H = 0.0)."""
    assert sum(1 for h in entropy if h == 0.0) == 381


def test_entropy_excludes_gaps(entropy):
    """Column 18 has 1 gap and 3 × 'C' — gaps excluded so H must be 0.0."""
    assert entropy[17] == 0.0


def test_entropy_two_base_column(sequences, entropy):
    """Column 20 has 3 × A and 1 × C with no gaps → H ≈ 0.8113 bits."""
    col = [seq[19].upper() for seq in sequences]
    assert col.count("-") == 0, "test assumption: column 20 has no gaps"
    assert abs(entropy[19] - 0.8113) < 0.001


# ── write_entropy_txt ─────────────────────────────────────────────────────────

def test_entropy_txt_row_count(sequences, entropy, tmp_path):
    """Output file must have a header line plus one row per alignment column."""
    out = tmp_path / "entropy.txt"
    write_entropy_txt(sequences, entropy, str(out))
    lines = out.read_text().splitlines()
    assert lines[0].startswith("pos\tentropy")   # header present
    assert len(lines) == 490 + 1                # header + 490 data rows


def test_entropy_txt_column_count(sequences, entropy, tmp_path):
    """Every data row must have 8 tab-separated columns."""
    out = tmp_path / "entropy.txt"
    write_entropy_txt(sequences, entropy, str(out))
    lines = out.read_text().splitlines()
    for line in lines[1:]:
        assert len(line.split("\t")) == 8


def test_entropy_txt_conserved_column_values(sequences, entropy, tmp_path):
    """First column is conserved: position=1, entropy=0.0000, IC=2.0000."""
    out = tmp_path / "entropy.txt"
    write_entropy_txt(sequences, entropy, str(out))
    fields = out.read_text().splitlines()[1].split("\t")
    assert fields[0] == "1"        # 1-based position
    assert fields[1] == "0.0000"   # entropy
    assert fields[2] == "2.0000"   # information content


def test_entropy_txt_positions_are_sequential(sequences, entropy, tmp_path):
    """Position column must run 1 … 490 without gaps."""
    out = tmp_path / "entropy.txt"
    write_entropy_txt(sequences, entropy, str(out))
    lines = out.read_text().splitlines()[1:]
    positions = [int(line.split("\t")[0]) for line in lines]
    assert positions == list(range(1, 491))


# ── SVG output files ──────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def svg_dir(tmp_path_factory, records, sequences, borders, entropy):
    """Run the full drawing pipeline once; return the output directory."""
    d = tmp_path_factory.mktemp("svg")
    draw_alignment_svg(records, sequences, borders, output=str(d / "alignment.svg"))
    draw_entropy_svg(records, sequences, entropy,  output=str(d / "entropy.svg"))
    draw_consensus_svg(sequences, borders,          output=str(d / "consensus.svg"))
    draw_logo_svg(records, sequences, entropy,      output=str(d / "logo.svg"))
    return d


@pytest.mark.parametrize("filename", [
    "alignment.svg", "entropy.svg", "consensus.svg", "logo.svg"
])
def test_svg_file_created_and_non_empty(svg_dir, filename):
    """Each drawing function must produce a non-empty file."""
    assert (svg_dir / filename).stat().st_size > 0


@pytest.mark.parametrize("filename", [
    "alignment.svg", "entropy.svg", "consensus.svg", "logo.svg"
])
def test_svg_files_are_valid_xml(svg_dir, filename):
    """Each SVG must be well-formed XML — catches truncated or broken output."""
    import xml.etree.ElementTree as ET
    ET.parse(svg_dir / filename)   # raises ParseError if malformed
