# AlignmentSketch

[![Tests](https://github.com/katkah/AlignmentDrawer/actions/workflows/test.yml/badge.svg)](https://github.com/katkah/AlignmentDrawer/actions/workflows/test.yml)

AlignmentSketch is a command-line tool for visualizing multiple sequence alignments (MSA). It takes a pre-aligned FASTA file and produces four complementary SVG figures: a full alignment view, a per-sequence conservation heatmap, a single-row consensus, and a sequence logo. All outputs are scalable vector graphics suitable for publication.

---

## Background

AlignmentSketch was developed to visualise variation in the ribosomal DNA intergenic spacer (IGS) of *Arabidopsis thaliana*. The original figures produced by this tool were published in:

> Havlová, K., Dvořáčková, M., Peiro, R., Abia, D., Mozgová, I., Vansáčová, L., Gutierrez, C., & Fajkus, J. (2016).
> **Variation of 45S rDNA intergenic spacers in *Arabidopsis thaliana*.**
> *Plant Molecular Biology.* https://doi.org/10.1007/s11103-016-0543-y

<img width="919" height="511" alt="Original alignment figure from Havlová et al. 2016" src="https://github.com/user-attachments/assets/4a6b5052-bf40-4442-a3c5-b2bf3999c21c" />

---

## Requirements

- Python 3.8 or newer
- [Biopython](https://biopython.org/) — for reading FASTA files
- [svgwrite](https://svgwrite.readthedocs.io/) — for generating SVG output

```bash
pip install biopython svgwrite
```

---

## Input format

AlignmentSketch expects a **pre-aligned FASTA file** — all sequences must be the same length. Gap characters are represented by `-`.

```
>Mouse
AAACCCGTTTTATGCAAACCCGTTTACGT--TTTGG
>Human
AAACCCGTTTCAGCGAAACCCGTTTGCGTA-TTTGG
>Chimp
AAACCCGTTTCAGCGAAACCCGTTTGCGT--TTTGG
>Rat
AAACCCGTTTATGCGAAACCCGTTTACGTA-TTTGG
```

If sequences are not all the same length, the tool exits immediately with a clear error. Input can be a file or piped from stdin, which makes it easy to chain with an alignment tool:

```bash
mafft sequences.fasta | python alignment_sketch.py
```

---

## Usage

```bash
# minimal — reads alignment.fasta, writes four SVGs with default names
python alignment_sketch.py -i alignment.fasta

# custom output names and wider line width
python alignment_sketch.py \
  -i alignment.fasta \
  -o alignment.svg \
  -e entropy.svg \
  -c consensus.svg \
  -l logo.svg \
  -w 1200

# also write a tab-separated entropy table
python alignment_sketch.py -i alignment.fasta -t entropy.txt
```

---

## Options

| Flag | Short | Default | Description |
|------|-------|---------|-------------|
| `--input` | `-i` | stdin | Pre-aligned FASTA input file |
| `--output` | `-o` | `alignment.svg` | Alignment view output |
| `--entropy-output` | `-e` | `entropy.svg` | Entropy heatmap output |
| `--consensus-output` | `-c` | `consensus.svg` | Consensus view output |
| `--logo-output` | `-l` | `logo.svg` | Sequence logo output |
| `--entropy-txt` | `-t` | *(none)* | Per-position entropy table (TSV) |
| `--width` | `-w` | `800` | Max sequence line width in pixels before wrapping |

Run `python alignment_sketch.py -h` for the full help text.

---

## Output files

### alignment.svg — Alignment view

The main view showing all sequences side by side.

- **Coloured letters** at positions where sequences differ
- **Grey numbered boxes** at positions where all sequences are identical — the number is the length of that conserved block
- Coordinate ticks below the last sequence at the start of each variable region
- Long alignments wrap automatically; line width is controlled with `--width`

Base colours: A = green, C = dark blue, G = black, T = dark red.

### entropy.svg — Entropy heatmap

A compact fixed-width overview strip showing conservation across the full alignment at a glance.

- **Dark grey** = conserved (low Shannon entropy), **light grey** = variable (high entropy)
- **Thin black bar** = gap (`-`) in that sequence
- The whole alignment always fits within `--width` pixels — no line wrapping

#### How entropy is calculated

For each alignment column, Shannon entropy is computed as:

```
H = −Σ (freq × log₂ freq)
```

where `freq` is the frequency of each base (A, C, G, T) among the **non-gap** sequences in that column. Gaps are excluded so that a column where all present sequences agree is correctly scored as conserved (H = 0) even when some sequences have a gap there.

- `H = 0.0` — fully conserved, drawn as **dark grey**
- `H = 2.0` — all four bases equally frequent, drawn as **light grey**

### consensus.svg — Consensus view

The complement of the alignment view — a single row showing only what is shared across all sequences.

- **Coloured letters** at positions conserved across all sequences
- **Grey numbered boxes** at variable positions, labelled with the region length
- Coordinate ticks at the start of each conserved block

### logo.svg — Sequence logo

A classic sequence logo showing conservation and base composition per column.

- **Bar height** = information content (IC = 2 − entropy, in bits). A fully conserved column reaches maximum height (2 bits); a completely random column produces no bar.
- **Coloured segments** show base frequency within each bar — most frequent base at the bottom
- Colour legend and position ticks every 20 bp are included

### entropy table (`--entropy-txt`)

A tab-separated file with one row per alignment column:

```
pos   entropy   IC      A   C   G   T   gaps
1     0.0000    2.0000  4   0   0   0   0
11    1.5000    0.5000  1   2   0   1   0
```

| Column | Description |
|--------|-------------|
| `pos` | 1-based alignment position |
| `entropy` | Shannon entropy (0 = conserved, up to 2.0 for DNA) |
| `IC` | Information content in bits (2 − entropy) |
| `A C G T` | Raw base counts at that column |
| `gaps` | Number of gap characters (`-`) |

---

## Examples

The figures below are produced from a real Arabidopsis IGS variant alignment (`data/variants_clustalo.fa`, 4 sequences × 490 bp).

**Alignment view** — coloured letters at variable positions, conserved blocks as numbered grey boxes:
![Alignment view](examples/alignment.svg)

**Entropy heatmap** — compact overview strip, dark grey = conserved, light grey = variable:
![Entropy heatmap](examples/entropy.svg)

**Consensus view** — single row showing only conserved positions:
![Consensus view](examples/consensus.svg)

**Sequence logo** — bar height = information content, segments = base frequencies:
![Sequence logo](examples/logo.svg)

---

## Testing

The test suite is built around a real Arabidopsis IGS variant alignment and covers alignment parsing, border detection, entropy calculation, table output, and SVG generation.

```bash
pip install pytest
python -m pytest test_alignment_sketch.py -v
```
