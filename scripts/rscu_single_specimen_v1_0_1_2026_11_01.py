#!/usr/bin/env python3
# =============================================================================
# rscu_single_specimen_v1_0_1_2026_11_01.py
# -----------------------------------------------------------------------------
# Program   : rscu_single_specimen
# Version   : v1.0.1
# Created   : 2026-11-01   (this version's creation date; stamped into the
#                           header, the output folder name, and every output
#                           filename)
# Author    : Generated for L. Brian Patrick
# -----------------------------------------------------------------------------
# CHANGES IN v1.0.1 (2026-11-01)
#   * NEW internal-stop check: a TAA/TAG anywhere other than the gene's terminal
#     stop position is counted per gene (Internal_Stops column in the diagnostic
#     CSV; Internal_Stops total in the summary CSV) and triggers a console
#     warning naming the affected gene(s). RSCU values are unchanged.
#   * Sharp, Tuohy & Mosurski (1986) added as the original RSCU reference.
#   * Removed the unverified claim of byte-for-byte identity with
#     rscu_integrated v2.0.2.
# -----------------------------------------------------------------------------
# WHAT THIS PROGRAM DOES
#   Computes Relative Synonymous Codon Usage (RSCU; Sharp et al. 1986; Sharp &
#   Li 1987) for a SINGLE specimen (one FASTA file whose records are the
#   protein-coding genes of that one specimen) under NCBI Genetic Code 5
#   (invertebrate mitochondrial; Osawa et al. 1992), and draws the stacked-bar
#   RSCU graph for that specimen.
#
#   It is a single-specimen distillation of the batch/comparative engine
#   rscu_integrated (v2.0.2): the per-specimen RSCU math and the stacked-bar
#   figure are derived from that engine's Phase 3, but everything that
#   requires MORE THAN ONE specimen (master matrices, phylogenetic PCA,
#   Ward clustering, the family codon-usage tests, metadata grouping, and the
#   Newick tree) has been removed, because none of it is defined for n = 1.
#
# OUTPUTS  (written to  <input folder>/RSCU_SingleSpecimen_Output_v1.0.1/ )
#   <specimen>_rscu_data_v1.0.1_2026-11-01.csv     Family, Codon, Count,
#                                                  Family_Total, RSCU
#   <specimen>_diagnostic_v1.0.1_2026-11-01.csv    one row per gene: Length_bp,
#                                                  Sense_Codons, Stop_Codons,
#                                                  Internal_Stops, Ambiguous,
#                                                  Truncated_Stop, Partial
#   <specimen>_summary_v1.0.1_2026-11-01.csv       specimen-level totals
#   <specimen>_rscu_chart_v1.0.1_2026-11-01.png    the RSCU stacked-bar graph
#
# RUN ENVIRONMENTS
#   PRIMARY (default): Windows 11 + IDLE shell, Python 3.13. Save this script in
#   the same folder as your specimen FASTA and press F5. With exactly one FASTA
#   in the folder it runs with no typing; with several it lists them and asks
#   which one. No command line is needed.
#
#   SECONDARY (optional): Ubuntu / Linux, including WSL on Windows. Run from a
#   shell and pass paths on the command line, e.g.
#       python3 rscu_single_specimen_v1_0_1_2026_11_01.py \
#           --input-file /path/to/Lycosa_shansia_OK032619.fasta
#   Any option below can also be set on the command line (see --help). Command
#   line values override the in-file USER OPTIONS.
#
# DEPENDENCIES:  numpy, pandas, matplotlib
#   pip install numpy pandas matplotlib          (Windows: from a Command Prompt)
#   pip3 install numpy pandas matplotlib          (Ubuntu / WSL)
#   (Biopython, SciPy and scikit-learn are NOT needed here — no tree, no stats
#    tests, no ordination.)
# =============================================================================

import os
import sys
import argparse
import traceback
from collections import OrderedDict

import numpy as np
import pandas as pd

import matplotlib
matplotlib.use("Agg")            # headless-safe on both IDLE and WSL: the figure
                                 # is SAVED to a PNG, never popped up in a window,
                                 # so no GUI backend is required in either place.
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.patches import Rectangle


# ═════════════════════════════════════════════════════════════════════════════
# ★★★ USER OPTIONS — EDIT BELOW, OR OVERRIDE ANY OF THEM ON THE COMMAND LINE ★★★
# ═════════════════════════════════════════════════════════════════════════════

INPUT_FILE = None
# Path to the ONE specimen FASTA file. Leave as None (DEFAULT) to auto-discover
# a FASTA/TXT file in INPUT_DIR (see below). Set it to point straight at your
# file, e.g.  INPUT_FILE = r"C:\Users\you\Documents\Lycosa_shansia_OK032619.fasta"
# CLI: --input-file PATH

INPUT_DIR = "."
# Folder to scan for the specimen file when INPUT_FILE is None. "." means "the
# folder this script is run from". WINDOWS/IDLE TIP: save this script INTO the
# folder that holds your specimen FASTA and just press F5 — then leave this as
# ".". If exactly one FASTA/TXT is found it is used automatically; if several are
# found you are asked which one. CLI: --input-dir PATH

SPECIMEN_LABEL = None
# Optional name for the specimen, used for the output filenames and the chart
# title. Leave None (DEFAULT) to derive it automatically: if every record header
# starts with the same first token that token is used (e.g. the headers
# ">Lycosa_shansia_OK032619 ATP6", ">Lycosa_shansia_OK032619 ATP8", ... give
# "Lycosa_shansia_OK032619"); otherwise the input file-name stem is used.
# CLI: --specimen-label "Genus_species"

SHOW_TITLES = True
# Set to False to suppress the figure title (publication mode). CLI: --titles off

OUTPUT_DIR = None
# When None (DEFAULT) the output folder RSCU_SingleSpecimen_Output_v1.0.1/ is
# created next to the input file, so inputs and results stay together. Set an
# absolute path to send outputs elsewhere. CLI: --output-dir PATH

OPEN_CHART_WHEN_DONE = False
# If True, try to open the finished PNG in the operating system's default image
# viewer after saving it (os.startfile on Windows, xdg-open on Linux/WSL).
# Default False so IDLE runs stay quiet and WSL runs never error on a missing
# viewer. CLI: --open-chart on

# ═════════════════════════════════════════════════════════════════════════════


# ── INTERNAL CONFIGURATION (do not edit unless you know why) ──────────────────
VERSION = "v1.0.1"
VERSION_DATE = "2026-11-01"                         # this version's creation date
OUTPUT_DIRNAME = f"RSCU_SingleSpecimen_Output_{VERSION}"
FASTA_EXTS = {".fasta", ".fa", ".fna", ".fas", ".txt"}

RANK_COLORS = ["#648FFF", "#785EF0", "#DC267F", "#FE6100"]   # IBM colorblind-safe

# NCBI Genetic Code Table 5 — invertebrate mitochondrial (Osawa et al. 1992):
# ATA = Met, AGA/AGG = Ser, TGA = Trp; only TAA/TAG are stops. Feeding in any
# other genome (nuclear, vertebrate mtDNA, ...) will MIS-CLASSIFY codons.
CODON_FAMILIES = OrderedDict([
    ('Ala', ['GCT', 'GCC', 'GCA', 'GCG']),
    ('Arg', ['CGT', 'CGC', 'CGA', 'CGG']),
    ('Asn', ['AAT', 'AAC']),
    ('Asp', ['GAT', 'GAC']),
    ('Cys', ['TGT', 'TGC']),
    ('Gln', ['CAA', 'CAG']),
    ('Glu', ['GAA', 'GAG']),
    ('Gly', ['GGT', 'GGC', 'GGA', 'GGG']),
    ('His', ['CAT', 'CAC']),
    ('Ile', ['ATT', 'ATC']),
    ('Leu1', ['CTT', 'CTC', 'CTA', 'CTG']),
    ('Leu2', ['TTA', 'TTG']),
    ('Lys', ['AAA', 'AAG']),
    ('Met', ['ATA', 'ATG']),
    ('Phe', ['TTT', 'TTC']),
    ('Pro', ['CCT', 'CCC', 'CCA', 'CCG']),
    ('Ser1', ['TCT', 'TCC', 'TCA', 'TCG']),
    ('Ser2', ['AGT', 'AGC', 'AGA', 'AGG']),
    ('Thr', ['ACT', 'ACC', 'ACA', 'ACG']),
    ('Trp', ['TGA', 'TGG']),
    ('Tyr', ['TAT', 'TAC']),
    ('Val', ['GTT', 'GTC', 'GTA', 'GTG']),
])
STOP_CODONS = {'TAA', 'TAG'}
CODON_TO_FAMILY = {c: fam for fam, codons in CODON_FAMILIES.items() for c in codons}


# ═════════════════════════════════════════════════════════════════════════════
# COMMAND-LINE INTERFACE  (used mainly from Ubuntu / WSL; IDLE just presses F5)
# ═════════════════════════════════════════════════════════════════════════════

def build_arg_parser():
    p = argparse.ArgumentParser(
        prog="rscu_single_specimen",
        description="Single-specimen RSCU (Sharp et al. 1986; Sharp & Li 1987) and RSCU stacked-bar "
                    "graph for one invertebrate mitogenome FASTA under NCBI "
                    "Genetic Code 5. Any option omitted keeps the in-file "
                    "default. Default run environment is the Windows 11 IDLE "
                    "shell (press F5); this CLI is for Ubuntu / WSL.")
    p.add_argument("--input-file", default=None,
                   help="Path to the single specimen FASTA/TXT file. If given, "
                        "no folder scan or prompt happens.")
    p.add_argument("--input-dir", default=None,
                   help="Folder to scan for one FASTA/TXT file when --input-file "
                        "is not given (default: the INPUT_DIR set in the script).")
    p.add_argument("--specimen-label", default=None,
                   help="Name for the specimen (output filenames + chart title). "
                        "Default: derived from the headers, else the file stem.")
    p.add_argument("--output-dir", default=None,
                   help="Absolute output folder. Default: an output folder is "
                        "created next to the input file.")
    p.add_argument("--titles", choices=["on", "off"], default=None,
                   help="Show or hide the figure title.")
    p.add_argument("--open-chart", choices=["on", "off"], default=None,
                   help="Open the finished PNG in the OS image viewer.")
    return p


def apply_cli_args(args):
    """Override the module-level USER OPTIONS with any CLI values supplied."""
    g = globals()
    if args.input_file is not None:
        g["INPUT_FILE"] = args.input_file
    if args.input_dir is not None:
        g["INPUT_DIR"] = args.input_dir
    if args.specimen_label is not None:
        g["SPECIMEN_LABEL"] = args.specimen_label
    if args.output_dir is not None:
        g["OUTPUT_DIR"] = args.output_dir
    if args.titles is not None:
        g["SHOW_TITLES"] = (args.titles == "on")
    if args.open_chart is not None:
        g["OPEN_CHART_WHEN_DONE"] = (args.open_chart == "on")


# ═════════════════════════════════════════════════════════════════════════════
# INPUT DISCOVERY AND PARSING
# ═════════════════════════════════════════════════════════════════════════════

def parse_fasta(filepath):
    """
    Yield (header, sequence) pairs. utf-8-sig transparently strips a UTF-8
    byte-order mark if present; U is converted to T and the sequence is
    upper-cased. Windows CRLF line endings are handled by .strip().
    """
    with open(filepath, 'r', encoding='utf-8-sig') as f:
        header, seq = None, []
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if header:
                    yield header, ''.join(seq).upper().replace('U', 'T')
                header = line[1:]
                seq = []
            else:
                seq.append(line)
        if header:
            yield header, ''.join(seq).upper().replace('U', 'T')


def detect_fasta_content(filepath):
    """A file is treated as sequence input only if its first non-empty line
    begins with '>' (so a stray .txt that is not FASTA is skipped)."""
    try:
        with open(filepath, 'r', encoding='utf-8-sig') as f:
            for line in f:
                line = line.strip()
                if line:
                    return line.startswith('>')
    except Exception:
        return False
    return False


def discover_input_file(input_file, input_dir):
    """
    Resolve which single FASTA/TXT file to analyze.
      * If input_file is given, validate and use it.
      * Else scan input_dir. 0 files -> clear error; 1 file -> use it;
        >1 files -> interactive numbered prompt (IDLE-friendly). In a headless
        run (no console) the prompt raises EOFError and we instruct the user to
        pass --input-file.
    Returns an absolute path.
    """
    if input_file:
        if not os.path.isfile(input_file):
            print(f"  ✗ INPUT_FILE not found: {input_file}")
            sys.exit(1)
        if not detect_fasta_content(input_file):
            print(f"  ✗ INPUT_FILE does not look like FASTA (no leading '>'): "
                  f"{input_file}")
            sys.exit(1)
        return os.path.abspath(input_file)

    if not os.path.isdir(input_dir):
        print(f"  ✗ INPUT_DIR is not a folder: {os.path.abspath(input_dir)}")
        sys.exit(1)

    candidates = []
    for name in sorted(os.listdir(input_dir)):
        ext = os.path.splitext(name)[1].lower()
        if ext in FASTA_EXTS:
            path = os.path.join(input_dir, name)
            if os.path.isfile(path) and detect_fasta_content(path):
                candidates.append(path)

    if not candidates:
        print("  ✗ No FASTA/TXT specimen file found in:")
        print(f"      {os.path.abspath(input_dir)}")
        print("    Fix: put your specimen FASTA there, set INPUT_FILE / INPUT_DIR "
              "near the top of the script, or pass --input-file PATH.")
        sys.exit(1)

    if len(candidates) == 1:
        return os.path.abspath(candidates[0])

    # More than one candidate: ask which specimen to analyze.
    print("  Multiple FASTA/TXT files found — this program analyzes ONE specimen:")
    for i, path in enumerate(candidates, start=1):
        print(f"    [{i}] {os.path.basename(path)}")
    try:
        sel = input(f"  Select the specimen file [1-{len(candidates)}]: ").strip()
        idx = int(sel) - 1
        if not (0 <= idx < len(candidates)):
            print("  ✗ Invalid selection.")
            sys.exit(1)
        return os.path.abspath(candidates[idx])
    except (ValueError, EOFError):
        print("  ✗ No valid selection. In a non-interactive run, pass the file "
              "explicitly with --input-file PATH.")
        sys.exit(1)


def read_specimen(path):
    """
    Read one specimen file into (specimen_id, [(gene_name, sequence), ...]).
      * gene_name = the record header text after the first whitespace-delimited
        token ("ID gene" -> gene); a single-token header -> that token.
      * specimen_id = SPECIMEN_LABEL if set; else the shared first header token
        when all records agree on it; else the input file-name stem.
    """
    records, first_tokens = [], []
    for header, seq in parse_fasta(path):
        parts = header.split(None, 1)
        if len(parts) > 1:
            first_tokens.append(parts[0])
            gene = parts[1].strip()
        elif parts:
            first_tokens.append(None)
            gene = parts[0]
        else:
            first_tokens.append(None)
            gene = "unknown"
        records.append((gene, seq))

    if not records:
        print(f"  ✗ No sequence records found in {os.path.basename(path)}.")
        sys.exit(1)

    stem = os.path.splitext(os.path.basename(path))[0]
    if SPECIMEN_LABEL:
        specimen_id = SPECIMEN_LABEL
    elif first_tokens and all(t is not None and t == first_tokens[0]
                              for t in first_tokens):
        specimen_id = first_tokens[0]
    else:
        specimen_id = stem
    return specimen_id, records


# ═════════════════════════════════════════════════════════════════════════════
# RSCU COMPUTATION  (derived from rscu_integrated v2.0.2, Phase 3)
# ═════════════════════════════════════════════════════════════════════════════

def compute_specimen_rscu(specimen_id, gene_records):
    """
    Pool all of the specimen's gene records and compute per-codon RSCU under
    Genetic Code 5. Returns (rscu_df, diag_df, summary).

    RSCU (Sharp et al. 1986; Sharp & Li 1987): for an amino-acid family with n
    synonymous codons and family total count F,
    RSCU(codon) = count(codon) * n / F   (0 if F=0).

    Internal stop codons (v1.0.1): a TAA/TAG is "internal" if it is NOT at the
    gene's terminal stop position. The terminal position is the last complete
    codon, UNLESS the gene ends in a truncated stop ("T"/"TA"), in which case
    the truncated stop is the terminal stop and EVERY complete TAA/TAG in the
    gene is internal. An internal stop means a frameshift, mis-annotation, or
    sequencing error; it is counted in Internal_Stops and warned about in the
    console. The codon is still tallied in Stop_Codons (not in RSCU), exactly
    as before, so RSCU values are unchanged.

    Truncated terminal stop codons (Ojala, Montoya & Attardi 1981): in
    invertebrate mtDNA many protein-coding genes end in an INCOMPLETE stop
    codon — a trailing "T" or "TA" — completed to a full TAA post-
    transcriptionally by polyadenylation of the mRNA 3' end. A trailing "T"/"TA"
    is counted as one complete terminal stop, recorded in Truncated_Stop, and
    does NOT flag the gene "Partial". Stop codons are not part of any sense
    family, so this never changes an RSCU value; it only keeps the diagnostics
    biologically correct.
    """
    codon_counts = {c: 0 for c in CODON_TO_FAMILY}
    diagnostic_rows = []
    n_genes = total_bp = total_ambig = partial_genes = total_trunc_stop = 0
    total_internal_stop = 0

    for gene, seq in gene_records:
        n_genes += 1
        total_bp += len(seq)

        full_len = (len(seq) // 3) * 3
        remainder = seq[full_len:]
        truncated_stop = remainder in ("T", "TA")
        # Start index of the codon allowed to be the terminal stop; with a
        # truncated stop there is none among the complete codons (-1).
        terminal_idx = -1 if truncated_stop else full_len - 3

        gene_sense = gene_ambig = gene_stop = gene_internal_stop = 0
        for i in range(0, full_len, 3):
            codon = seq[i:i + 3]
            if codon in CODON_TO_FAMILY:
                codon_counts[codon] += 1
                gene_sense += 1
            elif codon in STOP_CODONS:
                gene_stop += 1
                if i != terminal_idx:
                    gene_internal_stop += 1
            else:
                gene_ambig += 1
                total_ambig += 1
        total_internal_stop += gene_internal_stop

        if truncated_stop:
            gene_stop += 1
            total_trunc_stop += 1
        is_partial = bool(remainder) and not truncated_stop
        if is_partial:
            partial_genes += 1

        diagnostic_rows.append({
            'Gene': gene,
            'Length_bp': len(seq),
            'Sense_Codons': gene_sense,
            'Stop_Codons': gene_stop,
            'Internal_Stops': gene_internal_stop,
            'Ambiguous': gene_ambig,
            'Truncated_Stop': truncated_stop,
            'Partial': is_partial,
        })

    rscu_rows = []
    for fam, codons in CODON_FAMILIES.items():
        n = len(codons)
        fam_total = sum(codon_counts[c] for c in codons)
        for codon in codons:
            count = codon_counts[codon]
            rscu = (count * n) / fam_total if fam_total > 0 else 0.0
            rscu_rows.append({
                'Family': fam, 'Codon': codon, 'Count': count,
                'Family_Total': fam_total, 'RSCU': round(rscu, 4),
            })

    rscu_df = pd.DataFrame(rscu_rows)
    diag_df = pd.DataFrame(diagnostic_rows)
    summary = {
        'Specimen': specimen_id, 'N_Genes': n_genes, 'Total_bp': total_bp,
        'Total_Sense': sum(codon_counts.values()),
        'Total_Ambiguous': total_ambig,
        'Truncated_Stops': total_trunc_stop,
        'Internal_Stops': total_internal_stop,
        'Partial_Genes': partial_genes,
    }
    return rscu_df, diag_df, summary


# ═════════════════════════════════════════════════════════════════════════════
# RSCU GRAPH  (3-panel stacked bar: chart -> family labels -> codon key)
# ═════════════════════════════════════════════════════════════════════════════

def get_rank_colors(n):
    """Rank colors that never IndexError, even for families with >4 codons."""
    colors = list(RANK_COLORS)
    if n > len(colors):
        cmap = plt.get_cmap('viridis')
        extra = [cmap(x) for x in np.linspace(0, 1, n - len(colors))]
        colors = colors + extra
    return colors[:n]


def plot_stacked_rscu(rscu_df, specimen_id, output_path):
    """
    Stacked-bar RSCU chart: one bar per amino-acid family, segments stacked by
    codon rank (most- to least-used synonymous codon), a dashed reference line
    at RSCU = 1.0, and a codon/rank key. Panels top->bottom:
    chart (5) -> family labels (0.5) -> codon key (1.6).
    """
    fam_totals = rscu_df.groupby('Family')['Count'].sum().sort_values(ascending=False)
    sorted_fams = fam_totals.index.tolist()
    n_fams = len(sorted_fams)

    max_codons = int(rscu_df.groupby('Family').size().max())
    max_codons = max(max_codons, 4)
    rank_colors = get_rank_colors(max_codons)

    # Codons are ranked within each family by Count (ties broken arbitrarily).
    # Each rank has ONE fixed color, the bars are stacked rank 1 (bottom) to
    # rank 4 (top), and the codon key uses the same rank -> color -> row
    # mapping, so the stacking always matches the legend. Codons with Count 0
    # still receive a rank and appear in the codon key (zero-height segment).
    rank_lookup = {}
    for fam in sorted_fams:
        fam_data = (rscu_df[rscu_df['Family'] == fam]
                    .sort_values('Count', ascending=False)
                    .reset_index(drop=True))
        rank_lookup[fam] = {}
        for rank in range(max_codons):
            if rank < len(fam_data):
                rank_lookup[fam][rank] = (fam_data.iloc[rank]['Codon'],
                                          fam_data.iloc[rank]['RSCU'])
            else:
                rank_lookup[fam][rank] = ("", 0.0)

    fig = plt.figure(figsize=(16, 10))
    gs = GridSpec(3, 1, figure=fig, height_ratios=[5, 0.5, 1.6], hspace=0.08)

    # --- Panel 1: the stacked-bar chart -------------------------------------
    ax_chart = fig.add_subplot(gs[0])
    bottom = np.zeros(n_fams)
    for rank in range(max_codons):
        heights = [rank_lookup[fam][rank][1] for fam in sorted_fams]
        ax_chart.bar(range(n_fams), heights, bottom=bottom,
                     color=rank_colors[rank], edgecolor='white', linewidth=0.5,
                     label=f"Rank {rank + 1}")
        bottom += np.array(heights)

    ax_chart.axhline(y=1.0, color='black', linestyle='--', alpha=0.5, linewidth=0.8)
    max_stack = float(bottom.max())
    ax_chart.set_ylim(0, max(4.2, max_stack * 1.05))
    ax_chart.spines['top'].set_visible(False)
    ax_chart.spines['right'].set_visible(False)

    if SHOW_TITLES:
        ax_chart.set_title(f"RSCU Analysis: {specimen_id.replace('_', ' ')}",
                           fontweight='bold', fontsize=14, fontstyle='italic')
    ax_chart.set_ylabel("RSCU (stacked by rank)", fontsize=12, fontweight='bold')
    ax_chart.set_xticks(range(n_fams))
    ax_chart.set_xticklabels([''] * n_fams)
    ax_chart.set_xlim(-0.5, n_fams - 0.5)

    handles, labels = ax_chart.get_legend_handles_labels()
    ax_chart.legend(handles[::-1], labels[::-1], loc='upper right', fontsize=9,
                    framealpha=0.9, title='Codon Rank', title_fontsize=9)

    # --- Panel 3: the codon key (bottom) ------------------------------------
    ax_legend = fig.add_subplot(gs[2])
    ax_legend.set_xlim(-0.5, n_fams - 0.5)
    ax_legend.set_ylim(-0.5, max_codons - 0.5)
    box_w, box_h = 0.85, 0.85
    for col, fam in enumerate(sorted_fams):
        for rank in range(max_codons):
            codon_str, _ = rank_lookup[fam][rank]
            if codon_str:
                rect = Rectangle((col - box_w / 2, rank - box_h / 2), box_w, box_h,
                                 facecolor=rank_colors[rank],
                                 edgecolor='white', linewidth=0.8)
                ax_legend.add_patch(rect)
                ax_legend.text(col, rank, codon_str, ha='center', va='center',
                               fontsize=8.5, fontweight='bold',
                               fontfamily='monospace', color='white')
    ax_legend.set_xticks([])
    ax_legend.set_yticks([])
    for spine in ax_legend.spines.values():
        spine.set_visible(False)

    # --- Panel 2: family labels (middle, under the chart) -------------------
    ax_xlabel = fig.add_subplot(gs[1])
    ax_xlabel.set_xlim(-0.5, n_fams - 0.5)
    ax_xlabel.set_ylim(0, 1)
    for col, fam in enumerate(sorted_fams):
        ax_xlabel.text(col, 0.5, fam, ha='center', va='center',
                       fontsize=10, fontweight='bold')
    ax_xlabel.text(-1.0, 0.5, "Family:", ha='right', va='center',
                   fontsize=10, fontweight='bold', color='dimgray')
    ax_xlabel.set_xticks([])
    ax_xlabel.set_yticks([])
    for spine in ax_xlabel.spines.values():
        spine.set_visible(False)

    plt.savefig(output_path, dpi=300, facecolor='white', bbox_inches='tight')
    plt.close()


def maybe_open_chart(png_path):
    """Optionally open the PNG in the OS default viewer (guarded, never fatal)."""
    if not OPEN_CHART_WHEN_DONE:
        return
    try:
        if sys.platform.startswith("win"):
            os.startfile(png_path)                        # noqa: S606 (Windows)
        elif sys.platform == "darwin":
            os.system(f'open "{png_path}"')
        else:
            os.system(f'xdg-open "{png_path}" >/dev/null 2>&1 &')
    except Exception as e:
        print(f"  (Could not auto-open the chart: {e})")


# ═════════════════════════════════════════════════════════════════════════════
# MAIN
# ═════════════════════════════════════════════════════════════════════════════

def main():
    apply_cli_args(build_arg_parser().parse_args())

    print(f"\n{'=' * 72}")
    print(f"  rscu_single_specimen {VERSION}   (version date: {VERSION_DATE})")
    print(f"  Single-specimen RSCU + RSCU graph  |  NCBI Genetic Code 5")
    print(f"  Titles: {'ON' if SHOW_TITLES else 'OFF'}")
    print(f"{'=' * 72}\n")

    input_dir = INPUT_DIR if INPUT_DIR else "."
    input_path = discover_input_file(INPUT_FILE, input_dir)
    print(f"  Input file: {input_path}")

    specimen_id, gene_records = read_specimen(input_path)
    print(f"  Specimen:   {specimen_id}")
    print(f"  Genes read: {len(gene_records)} "
          f"({', '.join(g for g, _ in gene_records)})\n")

    # Output folder: next to the input file unless OUTPUT_DIR overrides it.
    if OUTPUT_DIR:
        out_dir = OUTPUT_DIR
    else:
        out_dir = os.path.join(os.path.dirname(input_path), OUTPUT_DIRNAME)
    os.makedirs(out_dir, exist_ok=True)

    try:
        rscu_df, diag_df, summary = compute_specimen_rscu(specimen_id, gene_records)
    except Exception as e:
        print(f"  ✗ ERROR computing RSCU for {specimen_id}: {e}")
        traceback.print_exc(limit=2)
        sys.exit(1)

    tag = f"{VERSION}_{VERSION_DATE}"            # e.g. v1.0.1_2026-11-01
    data_path = os.path.join(out_dir, f"{specimen_id}_rscu_data_{tag}.csv")
    diag_path = os.path.join(out_dir, f"{specimen_id}_diagnostic_{tag}.csv")
    summ_path = os.path.join(out_dir, f"{specimen_id}_summary_{tag}.csv")
    chart_path = os.path.join(out_dir, f"{specimen_id}_rscu_chart_{tag}.png")

    rscu_df.to_csv(data_path, index=False)
    diag_df.to_csv(diag_path, index=False)
    pd.DataFrame([summary]).to_csv(summ_path, index=False)
    plot_stacked_rscu(rscu_df, specimen_id, chart_path)

    # --- Console report ------------------------------------------------------
    print("─" * 72)
    print("  RSCU COMPUTED (Sharp et al. 1986; Sharp & Li 1987; NCBI Genetic Code 5)")
    print("─" * 72)
    print(f"  Genes pooled     : {summary['N_Genes']}")
    print(f"  Total bp         : {summary['Total_bp']}")
    print(f"  Sense codons     : {summary['Total_Sense']}")
    print(f"  Ambiguous codons : {summary['Total_Ambiguous']}")
    print(f"  Truncated stops  : {summary['Truncated_Stops']}  "
          f"(poly-A-completed; Ojala et al. 1981)")
    print(f"  Internal stops   : {summary['Internal_Stops']}")
    print(f"  Partial genes    : {summary['Partial_Genes']}")

    if summary['Total_Sense'] < 300:
        print("\n  ⚠ NOTE: fewer than 300 pooled sense codons — RSCU for rare "
              "amino-acid families (often Cys, Gln, Trp) rests on very few "
              "codons and should be read with caution.")
    if summary['Internal_Stops'] > 0:
        bad = diag_df[diag_df['Internal_Stops'] > 0]
        detail = ", ".join(f"{g} ({n})" for g, n in
                           zip(bad['Gene'], bad['Internal_Stops']))
        print(f"\n  ⚠ WARNING: {summary['Internal_Stops']} internal stop "
              f"codon(s) (TAA/TAG before the terminal stop) in "
              f"{len(bad)} gene(s): {detail}. This indicates a frameshift, "
              "mis-annotation, or sequencing error; the RSCU values for this "
              "specimen are affected. Check the Internal_Stops column in the "
              "diagnostic CSV.")
    if summary['Partial_Genes'] > 0:
        print(f"\n  ⚠ NOTE: {summary['Partial_Genes']} gene(s) had leftover bases "
              "that were NOT a recognized truncated stop (a frame/annotation "
              "issue). Check the Partial column in the diagnostic CSV.")

    print(f"\n  Output folder: {os.path.abspath(out_dir)}")
    print(f"    • {os.path.basename(data_path)}")
    print(f"    • {os.path.basename(diag_path)}")
    print(f"    • {os.path.basename(summ_path)}")
    print(f"    • {os.path.basename(chart_path)}")

    maybe_open_chart(chart_path)

    print(f"\n{'=' * 72}")
    print("  DONE")
    print(f"{'=' * 72}\n")


if __name__ == "__main__":
    main()
