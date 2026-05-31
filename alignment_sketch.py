import sys
import math
import argparse
import svgwrite
from Bio import SeqIO

FONT_SIZE = 14
CHAR_W = FONT_SIZE * 0.601   # monospace character width approximation
ROW_H = int(FONT_SIZE * 1.9) # vertical spacing between sequence rows
MAX_LINE_W = 800              # sequence area width (px) before line wrapping
HEAD_PAD = 12                 # gap between header text and first sequence column

BASE_COLORS = {
    'A': '#008000',  # green
    'C': '#3A5FCD',  # dark blue
    'G': '#000000',  # black
    'T': '#CC0000',  # dark red
}


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "AlignmentSketch — compact multiple sequence alignment visualizer.\n"
            "Reads a pre-aligned FASTA file (all sequences must be the same length,\n"
            "gaps represented as '-') and writes up to four SVG figures plus an\n"
            "optional tab-separated entropy table."
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "outputs\n"
            "-------\n"
            "  alignment.svg   Coloured letters at variable positions; conserved blocks\n"
            "                  compressed into grey numbered boxes showing block length.\n"
            "  entropy.svg     Per-sequence heatmap: dark grey = conserved, light grey =\n"
            "                  variable, thin black bar = gap in that sequence.\n"
            "  consensus.svg   Single-row view: coloured letters where all sequences agree,\n"
            "                  grey numbered boxes where they differ.\n"
            "  logo.svg        Sequence logo: bar height = information content (0-2 bits),\n"
            "                  coloured segments show base frequencies (A=green, C=blue,\n"
            "                  G=black, T=red). Tick marks every 20 bp.\n"
            "\n"
            "examples\n"
            "--------\n"
            "  python alignment_sketch.py -i alignment.fasta\n"
            "  python alignment_sketch.py -i alignment.fasta -w 1200 -t entropy.txt\n"
            "  cat alignment.fasta | python alignment_sketch.py -o out.svg\n"
        )
    )
    parser.add_argument(
        "--input", "-i", default=None, metavar="FILE",
        help="Pre-aligned FASTA input file. All sequences must be the same length. "
             "If omitted, reads from stdin.")
    parser.add_argument(
        "--output", "-o", default="alignment.svg", metavar="FILE",
        help="SVG output for the alignment view — coloured letters at variable positions, "
             "grey boxes at conserved blocks. (default: alignment.svg)")
    parser.add_argument(
        "--width", "-w", type=int, default=800, metavar="PX",
        help="Maximum width of the sequence area in pixels before the alignment wraps "
             "to a new line. Increase for wide monitors or long alignments. (default: 800)")
    parser.add_argument(
        "--entropy-output", "-e", default="entropy.svg", metavar="FILE",
        help="SVG output for the per-sequence entropy heatmap. Dark grey = conserved "
             "(low Shannon entropy), light grey = variable (high entropy), thin black "
             "bar = gap '-'. (default: entropy.svg)")
    parser.add_argument(
        "--consensus-output", "-c", default="consensus.svg", metavar="FILE",
        help="SVG output for the single-row consensus view. Positions conserved across "
             "all sequences are shown as coloured letters; variable regions are compressed "
             "into grey numbered boxes. (default: consensus.svg)")
    parser.add_argument(
        "--logo-output", "-l", default="logo.svg", metavar="FILE",
        help="SVG output for the sequence logo. Bar height represents information content "
             "(2 - Shannon entropy, in bits). Coloured segments within each bar show base "
             "frequencies. (default: logo.svg)")
    parser.add_argument(
        "--entropy-txt", "-t", default=None, metavar="FILE",
        help="Write a tab-separated entropy table to this file. Columns: position (1-based), "
             "Shannon entropy, information content, base counts (A C G T), gap count. "
             "Not written if this option is omitted.")
    return parser.parse_args()


def read_alignment(source):
    records = list(SeqIO.parse(source, "fasta"))
    if not records:
        print("Error: no sequences found in input.", file=sys.stderr)
        sys.exit(1)
    lengths = set(len(r.seq) for r in records)
    if len(lengths) > 1:
        print("Error: sequences are not all the same length.", file=sys.stderr)
        for r in records:
            print(f"  {r.id}: {len(r.seq)}", file=sys.stderr)
        sys.exit(1)
    return records


def compute_borders(sequences):
    borders = []
    prev_same = True
    for i in range(len(sequences[0])):
        this_same = all(seq[i] == sequences[0][i] for seq in sequences)
        if this_same != prev_same:
            borders.append(i)
            prev_same = this_same
    borders.append(len(sequences[0]))
    return borders


def compute_entropy(sequences):
    """Shannon entropy per alignment column.
    Gaps are excluded from the frequency calculation.
    Returns a list of floats (0.0 = fully conserved, up to 2.0 = max variable for DNA).
    """
    result = []
    base_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    for i in range(len(sequences[0])):
        counts = [0, 0, 0, 0]
        gaps = 0
        for seq in sequences:
            base = seq[i].upper()
            if base == '-':
                gaps += 1
            elif base in base_index:
                counts[base_index[base]] += 1
        total = len(sequences) - gaps
        h = 0.0
        if total > 0:
            for count in counts:
                if count > 0:
                    freq = count / total
                    h -= freq * math.log2(freq)
        result.append(h)
    return result


def write_entropy_txt(sequences, entropy, output):
    """Write a tab-separated table: position, entropy, information content, bases counts."""
    base_index = {'A': 0, 'C': 1, 'G': 2, 'T': 3}
    with open(output, 'w') as f:
        f.write("pos\tentropy\tIC\tA\tC\tG\tT\tgaps\n")
        for i, h in enumerate(entropy):
            counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
            gaps = 0
            for seq in sequences:
                b = seq[i].upper()
                if b == '-':
                    gaps += 1
                elif b in counts:
                    counts[b] += 1
            ic = max(0.0, 2.0 - h)
            f.write(f"{i+1}\t{h:.4f}\t{ic:.4f}\t"
                    f"{counts['A']}\t{counts['C']}\t{counts['G']}\t{counts['T']}\t{gaps}\n")
    print(f"Saved: {output}")


def draw_alignment_svg(records, sequences, borders, output="alignment.svg", max_line_w=MAX_LINE_W):
    n = len(sequences)
    head_w = max(len(r.id) for r in records) * CHAR_W + HEAD_PAD

    # Height of one wrapped line block: n sequence rows + space below for coordinates
    block_h = (n + 1) * ROW_H + ROW_H

    def y_seq(s, line_num):
        """Baseline Y for sequence row s on wrapped line line_num."""
        return line_num * block_h + (s + 1) * ROW_H

    # Build region list: (start_col, end_col, is_conserved)
    starts_conserved = borders[0] > 0
    regions = []
    is_con = starts_conserved
    prev = 0
    for b in borders:
        regions.append((prev, b, is_con))
        is_con = not is_con
        prev = b

    # --- First pass: collect drawing commands ---
    cmds = []
    x = head_w
    cur_line = 0
    lines_with_headers = set()

    def add_headers(ln):
        for s, r in enumerate(records):
            cmds.append(('text', 0, y_seq(s, ln), r.id, '#606060', FONT_SIZE))
        lines_with_headers.add(ln)

    def check_wrap(needed_w):
        nonlocal x, cur_line
        if x + needed_w > head_w + max_line_w and x > head_w:
            cur_line += 1
            x = head_w

    add_headers(0)

    for start, end, is_con in regions:
        if is_con:
            num = str(end - start)
            box_w = len(num) * CHAR_W + CHAR_W * 1.0  # text width + padding
            gap = CHAR_W * 0.1

            check_wrap(box_w + gap * 2)
            if cur_line not in lines_with_headers:
                add_headers(cur_line)

            x += gap
            text_x = x + (box_w - len(num) * CHAR_W) / 2
            for s in range(n):
                y_top = y_seq(s, cur_line) - FONT_SIZE - 2
                box_h = FONT_SIZE + 6
                cmds.append(('rect', x, y_top, box_w, box_h, '#EFEFEF', '#AAAAAA', 1))
                cmds.append(('text', text_x, y_seq(s, cur_line) + 2, num, '#555555', FONT_SIZE))

            x += box_w + gap

        else:
            for col_offset, col in enumerate(range(start, end)):
                check_wrap(CHAR_W)
                if cur_line not in lines_with_headers:
                    add_headers(cur_line)

                # coordinate tick at the first column of each variable block
                if col_offset == 0:
                    y_coord = y_seq(n - 1, cur_line) + ROW_H * 0.7
                    tick_x = x + CHAR_W / 2
                    cmds.append(('line', tick_x, y_coord - 5, tick_x, y_coord - 1, '#909090', 1))
                    cmds.append(('small_text', x, y_coord + 9, str(start + 1), '#909090', 9))

                for s, seq in enumerate(sequences):
                    base = seq[col].upper()
                    color = BASE_COLORS.get(base, '#000000')
                    cmds.append(('text', x, y_seq(s, cur_line), base, color, FONT_SIZE))

                x += CHAR_W

    # --- Canvas size ---
    svg_w = head_w + max_line_w + HEAD_PAD
    svg_h = (cur_line + 1) * block_h + ROW_H

    # --- Second pass: build SVG ---
    dwg = svgwrite.Drawing(output, size=(f'{int(svg_w)}px', f'{int(svg_h)}px'))

    for cmd in cmds:
        kind = cmd[0]
        if kind in ('text', 'small_text'):
            _, cx, cy, txt, fill, fs = cmd
            dwg.add(dwg.text(txt, insert=(cx, cy), fill=fill,
                             font_size=f'{fs}px',
                             font_family='Courier New, monospace'))
        elif kind == 'rect':
            _, cx, cy, w, h, fill, stroke, sw = cmd
            dwg.add(dwg.rect(insert=(cx, cy), size=(w, h),
                             fill=fill, stroke=stroke, stroke_width=sw))
        elif kind == 'line':
            _, x1, y1, x2, y2, stroke, sw = cmd
            dwg.add(dwg.line(start=(x1, y1), end=(x2, y2),
                             stroke=stroke, stroke_width=sw))

    dwg.save()
    print(f"Saved: {output}")


def count_equal_color_run(seq, entropy_vals, position):
    """Count consecutive non-gap positions from 'position' that share the same entropy color.
    Direct translation of the original C++ number_of_equal_elements."""
    count = 1
    i = position
    while i + 1 < len(seq) and seq[i + 1] != '-':
        c1 = min(255, int(256 - entropy_vals[i] * 100))
        c2 = min(255, int(256 - entropy_vals[i + 1] * 100))
        if c1 == c2:
            count += 1
            i += 1
        else:
            break
    return count


def find_deletion_end(seq, position):
    """Return the index just past the end of a gap run starting at position."""
    i = position
    while i < len(seq) and seq[i] == '-':
        i += 1
    return i


def draw_entropy_svg(records, sequences, entropy, output="entropy.svg", max_line_w=MAX_LINE_W):
    """Compact overview strip: the whole alignment fits in max_line_w — no line wrapping.
    COL_W is computed from the sequence length so every column gets a proportional slice."""
    n = len(sequences)
    seq_len = len(sequences[0])

    COL_W = max_line_w / seq_len  # float — may be well below 1 px for long alignments
    RECT_H = ROW_H - 2
    DEL_H = max(2, ROW_H // 8)
    head_w = int(max(len(r.id) for r in records) * CHAR_W + HEAD_PAD)

    def y_row(s):
        return (s + 1) * ROW_H

    def y_rect_top(s):
        return y_row(s) - FONT_SIZE // 2 - RECT_H // 2

    def y_del_top(s):
        return y_row(s) - FONT_SIZE // 2 - DEL_H // 2

    cmds = []

    # headers — drawn once on the left
    for s, r in enumerate(records):
        cmds.append(('text', 0, y_row(s), r.id, '#606060', FONT_SIZE))

    # one horizontal strip per sequence, no wrapping
    for s in range(n):
        x = head_w
        position = 0

        while position < seq_len:
            if sequences[s][position] != '-':
                run = count_equal_color_run(sequences[s], entropy, position)
                gray = min(220, int(150 + entropy[position] * 100))
                color = f'rgb({gray},{gray},{gray})'
                cmds.append(('rect', x, y_rect_top(s), COL_W * run, RECT_H, color, color, 0))
                x += COL_W * run
                position += run
            else:
                del_end = find_deletion_end(sequences[s], position)
                del_run = del_end - position
                cmds.append(('rect', x, y_del_top(s), COL_W * del_run, DEL_H,
                             '#000000', '#000000', 0))
                x += COL_W * del_run
                position = del_end

    svg_w = head_w + max_line_w + HEAD_PAD
    svg_h = (n + 1) * ROW_H

    dwg = svgwrite.Drawing(output, size=(f'{int(svg_w)}px', f'{int(svg_h)}px'))

    for cmd in cmds:
        kind = cmd[0]
        if kind == 'text':
            _, cx, cy, txt, fill, fs = cmd
            dwg.add(dwg.text(txt, insert=(cx, cy), fill=fill,
                             font_size=f'{fs}px',
                             font_family='Courier New, monospace'))
        elif kind == 'rect':
            _, cx, cy, w, h, fill, stroke, sw = cmd
            if w > 0 and h > 0:
                dwg.add(dwg.rect(insert=(cx, cy), size=(w, h),
                                 fill=fill, stroke=stroke, stroke_width=sw))

    dwg.save()
    print(f"Saved: {output}")


def draw_logo_svg(records, sequences, entropy, output="logo.svg", max_line_w=MAX_LINE_W):
    """Sequence logo: stacked coloured bars per column.
    Bar height = information content (2 - entropy in bits).
    Each segment height = base frequency × total bar height.
    Segments ordered most → least frequent bottom to top."""

    n = len(sequences)
    seq_len = len(sequences[0])

    COL_W = 6            # column width in pixels
    LOGO_MAX_H = 80      # height of a fully conserved column (2 bits)
    SCALE = LOGO_MAX_H / 2.0  # pixels per bit

    head_w = HEAD_PAD + 10  # logo only needs room for the "0"/"2" scale labels
    block_h = LOGO_MAX_H + 4 * ROW_H  # title row + legend row + bars + tick labels

    def y_baseline(ln):
        """Bottom edge of the logo bars for wrapped line ln."""
        return ln * block_h + 2 * ROW_H + LOGO_MAX_H

    cmds = []
    lines_used = {0}
    x = head_w
    cur_line = 0

    def check_wrap():
        nonlocal x, cur_line
        if x + COL_W > head_w + max_line_w and x > head_w:
            cur_line += 1
            x = head_w
            lines_used.add(cur_line)

    for col in range(seq_len):
        check_wrap()
        yb_cur = y_baseline(cur_line)

        # count bases, exclude gaps
        counts = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        gaps = sum(1 for seq in sequences if seq[col].upper() == '-')
        for seq in sequences:
            b = seq[col].upper()
            if b in counts:
                counts[b] += 1

        total = n - gaps

        if total > 0:
            ic = max(0.0, 2.0 - entropy[col])
            total_bar_h = ic * SCALE
            if total_bar_h >= 0.5:
                # draw segments bottom → top, most frequent first
                sorted_bases = sorted(counts, key=lambda b: counts[b], reverse=True)
                y = yb_cur
                for base in sorted_bases:
                    if counts[base] == 0:
                        continue
                    freq = counts[base] / total
                    seg_h = freq * total_bar_h
                    if seg_h < 0.5:
                        continue
                    color = BASE_COLORS.get(base, '#000000')
                    cmds.append(('rect', x, y - seg_h, COL_W, seg_h, color, color, 0))
                    y -= seg_h

        # tick mark every 20 bp (positions 20, 40, 60 … in 1-based coords)
        if (col + 1) % 20 == 0:
            tick_x = x + COL_W / 2
            cmds.append(('line', tick_x, yb_cur, tick_x, yb_cur + 5, '#AAAAAA', 1))
            cmds.append(('center_text', tick_x, yb_cur + 15, str(col + 1), '#AAAAAA', 9))

        x += COL_W

    # per-block decorations: title, baseline, scale ticks, legend
    sq = 8   # legend colour square size (px)
    for ln in sorted(lines_used):
        yb = y_baseline(ln)
        # gray baseline
        cmds.append(('line', head_w, yb, head_w + max_line_w, yb, '#CCCCCC', 1))
        # centered title at the very top of the block
        cmds.append(('center_text', head_w + max_line_w / 2,
                     ln * block_h + int(ROW_H * 0.85),
                     'Sequence Logo', '#606060', FONT_SIZE + 2))
        # scale ticks: "2" at top of bars, "0" at baseline
        cmds.append(('small_text', 0, ln * block_h + 2 * ROW_H, '2', '#AAAAAA', 9))
        cmds.append(('small_text', 0, yb + 10, '0', '#AAAAAA', 9))
        # legend: coloured square + letter for each base, between title and bars
        leg_x = head_w
        leg_y = ln * block_h + int(ROW_H * 1.75)
        for base in ['A', 'C', 'G', 'T']:
            color = BASE_COLORS[base]
            cmds.append(('rect', leg_x, leg_y - sq, sq, sq, color, color, 0))
            cmds.append(('small_text', leg_x + sq + 2, leg_y, base, color, 9))
            leg_x += sq + 2 + CHAR_W + 8   # square + gap + letter + spacing

    svg_w = head_w + max_line_w + HEAD_PAD
    svg_h = (max(lines_used) + 1) * block_h + ROW_H

    dwg = svgwrite.Drawing(output, size=(f'{int(svg_w)}px', f'{int(svg_h)}px'))

    for cmd in cmds:
        kind = cmd[0]
        if kind in ('text', 'small_text'):
            _, cx, cy, txt, fill, fs = cmd
            dwg.add(dwg.text(txt, insert=(cx, cy), fill=fill,
                             font_size=f'{fs}px',
                             font_family='Courier New, monospace'))
        elif kind == 'center_text':
            _, cx, cy, txt, fill, fs = cmd
            dwg.add(dwg.text(txt, insert=(cx, cy), fill=fill,
                             font_size=f'{fs}px',
                             font_family='Courier New, monospace',
                             **{'text-anchor': 'middle'}))
        elif kind == 'rect':
            _, cx, cy, w, h, fill, stroke, sw = cmd
            if h > 0:
                dwg.add(dwg.rect(insert=(cx, cy), size=(w, h),
                                 fill=fill, stroke=stroke, stroke_width=sw))
        elif kind == 'line':
            _, x1, y1, x2, y2, stroke, sw = cmd
            dwg.add(dwg.line(start=(x1, y1), end=(x2, y2),
                             stroke=stroke, stroke_width=sw))

    dwg.save()
    print(f"Saved: {output}")


def draw_consensus_svg(sequences, borders, output="consensus.svg", max_line_w=MAX_LINE_W):
    """Single-row consensus view: conserved positions as coloured letters,
    variable regions compressed into numbered grey boxes."""

    head_w = int(len("Consensus") * CHAR_W + HEAD_PAD)
    block_h = 3 * ROW_H  # single row + coordinate space + margin

    def y_row(ln):
        return ln * block_h + ROW_H

    # Build region list — same as draw_alignment_svg but roles are swapped:
    #   is_con=True  → conserved → draw letters
    #   is_con=False → variable  → draw numbered box
    starts_conserved = borders[0] > 0
    regions = []
    is_con = starts_conserved
    prev = 0
    for b in borders:
        regions.append((prev, b, is_con))
        is_con = not is_con
        prev = b

    cmds = []
    x = head_w
    cur_line = 0
    lines_with_headers = set()

    def add_header(ln):
        cmds.append(('text', 0, y_row(ln), 'Consensus', '#606060', FONT_SIZE))
        lines_with_headers.add(ln)

    def check_wrap(needed):
        nonlocal x, cur_line
        if x + needed > head_w + max_line_w and x > head_w:
            cur_line += 1
            x = head_w

    add_header(0)

    for start, end, is_con in regions:
        if is_con:
            # conserved region → individual coloured letters
            for col_offset, col in enumerate(range(start, end)):
                check_wrap(CHAR_W)
                if cur_line not in lines_with_headers:
                    add_header(cur_line)

                # coordinate tick at the start of each conserved block
                if col_offset == 0:
                    y_coord = y_row(cur_line) + ROW_H * 0.7
                    tick_x = x + CHAR_W / 2
                    cmds.append(('line', tick_x, y_coord - 5, tick_x, y_coord - 1, '#909090', 1))
                    cmds.append(('small_text', x, y_coord + 9, str(start + 1), '#909090', 9))

                base = sequences[0][col].upper()
                color = BASE_COLORS.get(base, '#000000')
                cmds.append(('text', x, y_row(cur_line), base, color, FONT_SIZE))
                x += CHAR_W

        else:
            # variable region → grey numbered box
            num = str(end - start)
            box_w = len(num) * CHAR_W + CHAR_W * 1.0
            gap = CHAR_W * 0.1

            check_wrap(box_w + gap * 2)
            if cur_line not in lines_with_headers:
                add_header(cur_line)

            x += gap
            text_x = x + (box_w - len(num) * CHAR_W) / 2
            y_top = y_row(cur_line) - FONT_SIZE - 2
            cmds.append(('rect', x, y_top, box_w, FONT_SIZE + 6, '#EFEFEF', '#AAAAAA', 1))
            cmds.append(('text', text_x, y_row(cur_line) + 2, num, '#555555', FONT_SIZE))
            x += box_w + gap

    svg_w = head_w + max_line_w + HEAD_PAD
    svg_h = (cur_line + 1) * block_h + ROW_H

    dwg = svgwrite.Drawing(output, size=(f'{int(svg_w)}px', f'{int(svg_h)}px'))

    for cmd in cmds:
        kind = cmd[0]
        if kind in ('text', 'small_text'):
            _, cx, cy, txt, fill, fs = cmd
            dwg.add(dwg.text(txt, insert=(cx, cy), fill=fill,
                             font_size=f'{fs}px',
                             font_family='Courier New, monospace'))
        elif kind == 'rect':
            _, cx, cy, w, h, fill, stroke, sw = cmd
            dwg.add(dwg.rect(insert=(cx, cy), size=(w, h),
                             fill=fill, stroke=stroke, stroke_width=sw))
        elif kind == 'line':
            _, x1, y1, x2, y2, stroke, sw = cmd
            dwg.add(dwg.line(start=(x1, y1), end=(x2, y2),
                             stroke=stroke, stroke_width=sw))

    dwg.save()
    print(f"Saved: {output}")


def main():
    args = parse_args()

    if args.input:
        with open(args.input) as f:
            records = read_alignment(f)
    else:
        records = read_alignment(sys.stdin)

    sequences = [str(r.seq) for r in records]
    borders = compute_borders(sequences)
    entropy = compute_entropy(sequences)
    draw_alignment_svg(records, sequences, borders, output=args.output, max_line_w=args.width)
    draw_entropy_svg(records, sequences, entropy, output=args.entropy_output, max_line_w=args.width)
    draw_consensus_svg(sequences, borders, output=args.consensus_output, max_line_w=args.width)
    draw_logo_svg(records, sequences, entropy, output=args.logo_output, max_line_w=args.width)
    if args.entropy_txt:
        write_entropy_txt(sequences, entropy, args.entropy_txt)


if __name__ == "__main__":
    main()
