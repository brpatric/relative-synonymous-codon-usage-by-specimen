===============================================================================
  SINGLE-SPECIMEN RSCU
  Relative Synonymous Codon Usage values + RSCU graph for ONE specimen
  README
===============================================================================

  README file  : README_rscu_single_specimen_v1_0_1_2026_11_01.txt
                 (README v1.0.1, created 2026-11-01; tracks the program version
                 it documents)
  Program      : rscu_single_specimen_v1_0_1_2026_11_01.py   (v1.0.1)
  Example data : example/Lycosa_shansia_OK032619.fasta   (the 13 protein-coding
                 genes of the published Lycosa shansia mitogenome, GenBank
                 OK032619; Ye et al. 2022; used throughout Section 8)
  Author       : L. Brian Patrick, Department of Biological Sciences, Dakota
                 Wesleyan University, Mitchell, SD, USA (ORCID 0000-0003-1969-9600)
  Repository   : https://github.com/brpatric/relative-synonymous-codon-usage-by-specimen
  Archive/DOI  : Zenodo, doi:10.5281/zenodo.24182022 (release v1.0.1)
  Citation     : Patrick LB. 2026. rscu_single_specimen, version 1.0.1.
                 Zenodo. doi:10.5281/zenodo.24182022.
  License      : none specified in the repository

  FILE-NAMING CONVENTION. Every delivered file carries this version and the date
  that version was created. The program and this README use the fully-underscored
  form  <name>_v<major>_<minor>_<patch>_<YYYY>_<MM>_<DD>  (e.g.
  rscu_single_specimen_v1_0_1_2026_11_01.py). The data/graph files the program
  writes use the dotted-version, hyphen-date form  <name>_v1.0.1_2026-11-01.<ext>
  so the version identity travels with every output. The date component is this
  version's CREATION date (2026-11-01), not the wall-clock run date, so re-running
  the same version reproduces the same filenames.

  Genetic code : NCBI Translation Table 5 (invertebrate mitochondrial), fixed
                 (Osawa et al. 1992).
  Language     : Python 3 (developed/tested on 3.13; runs from the IDLE shell).
  Environment  : Windows 11 + IDLE (PRIMARY, recommended, the default). Ubuntu /
                 Linux, including WSL on Windows, is a documented SECONDARY
                 option. See "RUN ENVIRONMENTS" below.
  Dependencies : numpy, pandas, matplotlib   (no Biopython, SciPy, or
                 scikit-learn — there is no tree, no ordination, and no
                 statistical test in a single-specimen analysis).


-------------------------------------------------------------------------------
  WHAT THIS VERSION IS  (v1.0.1, 2026-11-01)
-------------------------------------------------------------------------------
  CHANGES IN v1.0.1 (from v1.0.0, 2026-08-20):
     * NEW internal-stop check. A TAA/TAG codon anywhere other than the gene's
       terminal stop position is now counted per gene (new Internal_Stops
       column in the diagnostic CSV, new Internal_Stops total in the summary
       CSV), and the console prints a WARNING naming each affected gene. See
       Sections 2.2 and 6. RSCU values are unchanged.
     * Sharp, Tuohy & Mosurski (1986), the paper that first defined RSCU, is
       now cited alongside Sharp & Li (1987).
     * The claim that this tool reproduces rscu_integrated's numbers exactly
       has been removed (it was not independently verified).
     * Section 8 wording corrected (see the RSCU sanity check there), and the
       worked example now uses the published Lycosa shansia mitogenome
       (GenBank OK032619; Ye et al. 2022), included in example/ with its
       output files.

  This is a SINGLE-SPECIMEN tool: it takes ONE FASTA file whose records are the
  protein-coding genes of ONE specimen, computes the Relative Synonymous Codon
  Usage (RSCU) values for that specimen (Sharp et al. 1986; Sharp & Li 1987)
  under NCBI Genetic Code 5 (Osawa et al. 1992), and draws the RSCU stacked-bar
  graph.

  It is derived from the per-specimen "Phase 3" of the batch / comparative
  engine rscu_integrated (v2.0.2, 2026-07-24): the codon-counting, the RSCU
  formula, the truncated-stop handling, and the stacked-bar figure were taken
  from that engine. Everything in that engine that requires MORE THAN ONE
  specimen has been removed on purpose, because none of it is defined for a
  sample size of one:

     * master species x codon matrices                (needs >= 2 specimens)
     * phylogenetic PCA ordination                    (needs a tree + >= 3 taxa)
     * phylogenetically-corrected Ward clustering     (needs a tree + >= 3 taxa)
     * the RSCU heatmap across specimens              (needs >= 2 specimens)
     * the per-family codon-usage tests (Friedman /   (a blocked design across
       Wilcoxon with Benjamini-Hochberg FDR)           specimens; n=1 has no test)
     * metadata grouping / color coding               (compares groups of taxa)

  Removing those also removes their dependencies: this tool needs only numpy,
  pandas, and matplotlib.


===============================================================================
  RUN ENVIRONMENTS  (Windows 11 / IDLE is the primary, recommended, default;
                     Ubuntu / Linux, incl. WSL, is a documented secondary option)
===============================================================================

  This program is developed and run primarily on WINDOWS 11 using the IDLE shell
  (Python 3.13). That is the default path and everything below is written for it
  first. Ubuntu / Linux (including WSL on a Windows machine) is fully supported
  as a SECONDARY option and is documented in ENVIRONMENT 2. The science is
  identical in both; only how you launch the program differs.

  -----------------------------------------------------------------------------
  ENVIRONMENT 1 (PRIMARY, RECOMMENDED, DEFAULT): WINDOWS 11 + IDLE
  -----------------------------------------------------------------------------
  1. Install Python 3.13. During install, tick "Add Python to PATH" (convenient;
     not strictly required for this program).

  2. Install the three packages once, from a Command Prompt:
         pip install numpy pandas matplotlib

  3. Put this script IN THE SAME FOLDER as your specimen FASTA, for example:
         rscu_single_specimen_v1_0_1_2026_11_01.py
         Lycosa_shansia_OK032619.fasta

  4. Open rscu_single_specimen_v1_0_1_2026_11_01.py in IDLE and press F5.
       * If the folder holds EXACTLY ONE FASTA/TXT file, it is analyzed with no
         typing at all.
       * If the folder holds SEVERAL FASTA/TXT files, IDLE lists them in the
         shell and asks which single specimen to analyze (type the number).
     The RSCU values, diagnostics, summary, and the RSCU graph are written to a
     new sub-folder RSCU_SingleSpecimen_Output_v1.0.1/ beside your FASTA, and the
     four output paths are printed in the shell.

  5. (Optional) Instead of relying on the folder scan, you can point the script
     straight at a file by editing the USER OPTIONS block near the top:
         INPUT_FILE = r"C:\Users\you\Documents\Lycosa_shansia_OK032619.fasta"
     Use a raw string (the r"...") so Windows backslashes are taken literally.

  Notes for IDLE:
     * The figure is SAVED as a 300-dpi PNG; it is not popped up in a window
       (the program uses matplotlib's file-only "Agg" backend, which is the
       robust choice in both IDLE and WSL). To have Windows open the PNG in your
       default image viewer when the run finishes, set OPEN_CHART_WHEN_DONE=True
       in the USER OPTIONS block.
     * The interactive file-selection prompt (step 4) appears in the IDLE shell
       and works normally there.

  -----------------------------------------------------------------------------
  ENVIRONMENT 2 (SECONDARY): UBUNTU / LINUX COMMAND LINE  (incl. WSL on Windows)
  -----------------------------------------------------------------------------
  Use this when you want to run on a Linux box or inside WSL2 on a Windows
  machine. On Windows, install WSL once with:  wsl --install -d Ubuntu  (then
  open "Ubuntu").

  1. Python + packages (Ubuntu/WSL shell):
         sudo apt update
         sudo apt install -y python3 python3-pip
         pip3 install numpy pandas matplotlib
     (If pip refuses on a system-managed Python, either use a virtual
      environment, or add --break-system-packages.)

  2. Run it, passing the specimen file explicitly (recommended in a shell so no
     prompt is needed):
         python3 rscu_single_specimen_v1_0_1_2026_11_01.py \
             --input-file /path/to/Lycosa_shansia_OK032619.fasta

     Other options (each mirrors an in-file USER OPTION; anything omitted keeps
     the in-file default):
         --input-dir DIR         scan DIR for one FASTA/TXT instead of naming a
                                 file (0 found -> error; 1 -> used; >1 -> prompt)
         --specimen-label NAME   name used in output filenames and chart title
         --output-dir DIR        write outputs to DIR instead of beside the input
         --titles on|off         show or hide the figure title (off = publication)
         --open-chart on|off     open the finished PNG in the OS image viewer
         -h, --help              full option help

  3. WSL path tip: you can read a file that lives on the Windows side directly,
     e.g. --input-file /mnt/c/Users/you/Documents/specimen.fasta . Outputs land
     beside the input (so also on the Windows drive) unless --output-dir says
     otherwise.

  NB: the IDLE (native-Windows) and WSL Pythons are separate environments; run in
  whichever one you installed the three packages into. The numbers and the figure
  are identical either way.


===============================================================================
  TABLE OF CONTENTS
===============================================================================
  1. OVERVIEW — what the program does and what RSCU means
  2. WARNINGS AND CAVEATS (read before using results)
  3. INPUTS
  4. USER OPTIONS (configuration)
  5. HOW TO RUN (quick reference)
  6. METHOD DETAIL
  7. OUTPUT FILES
  8. WORKED EXAMPLE (Lycosa shansia, GenBank OK032619)
  9. IN-TEXT CITATION MAP
  10. REFERENCES (with DOIs)


-------------------------------------------------------------------------------
1. OVERVIEW
-------------------------------------------------------------------------------
  The program reads one specimen's protein-coding genes from a single FASTA
  file, pools them, counts codons under NCBI Genetic Code 5 (invertebrate
  mitochondrial; Osawa et al. 1992), computes per-codon RSCU (Sharp et al. 1986;
  Sharp & Li 1987), and draws a stacked-bar RSCU graph.

  RSCU measures whether a given synonymous codon is used more or less often than
  expected if all synonymous codons for an amino acid were used equally (Sharp
  et al. 1986; Sharp & Li 1987). An RSCU of 1.0 means "used exactly as
  expected"; > 1.0 means the codon is over-represented (preferred) within its
  amino-acid family, and < 1.0 means it is under-represented. Because RSCU
  normalizes raw counts by both amino-acid usage and family size, it isolates
  synonymous-codon preference from amino-acid composition and is comparable
  across genes and families within the specimen.

  A single specimen is a complete, self-contained unit of analysis for RSCU: the
  measure is defined per specimen from that specimen's own pooled codon counts.
  It does NOT need other specimens, a phylogeny, or any grouping information —
  which is exactly why this stand-alone tool exists.


-------------------------------------------------------------------------------
2. WARNINGS AND CAVEATS  (READ BEFORE USING RESULTS)
-------------------------------------------------------------------------------
  2.1  GENETIC CODE 5 ONLY. Codons are classified using the invertebrate
       mitochondrial code (Table 5; Osawa et al. 1992), where ATA = Met,
       AGA/AGG = Ser, TGA = Trp, and only TAA/TAG are stops. Feeding in nuclear
       genes, vertebrate mtDNA, or any other genome will MIS-CLASSIFY codons and
       produce wrong RSCU values. The program is not configurable to other codes
       without editing the CODON_FAMILIES / STOP_CODONS tables in the source.

  2.2  INPUT MUST BE CODING SEQUENCE (CDS), IN FRAME. Reading always starts at
       the first base of each record. A record that is not a protein-coding gene,
       or that does not begin in reading frame, will inflate the "ambiguous" and
       "stop" counts and distort RSCU. A record whose length is not a multiple of
       3 is flagged "Partial" UNLESS the leftover bases are a recognized
       truncated stop (see 2.3). A TAA/TAG codon anywhere other than the gene's
       terminal stop position is an INTERNAL STOP (new in v1.0.1): it is counted
       in the Internal_Stops diagnostic column and the console prints a WARNING
       naming the gene. An internal stop means a frameshift, mis-annotation, or
       sequencing error in that gene, and the specimen's RSCU values are affected
       until it is fixed.

  2.3  TRUNCATED / ABBREVIATED TERMINAL STOP CODONS ARE RECOGNIZED. In
       invertebrate mtDNA, many protein-coding genes terminate in an incomplete
       stop codon — a trailing "T" or "TA" — completed to a full TAA after
       transcription by polyadenylation of the mRNA 3' end (Ojala, Montoya &
       Attardi 1981). The program treats a trailing "T" or "TA" (the 1-2 leftover
       bases after the last complete codon) as ONE complete terminal stop codon:
       it is added to the gene's stop count, recorded in the Truncated_Stop
       diagnostic column, and the gene is NOT flagged "Partial". Because stop
       codons are not part of any sense family, this does not change any RSCU
       value; it only makes the stop/partial diagnostics biologically correct.
       Leftover bases that are anything OTHER than "T"/"TA" still flag the gene
       "Partial" (a genuine frame or annotation issue).

  2.4  SINGLE SPECIMEN ONLY — NO COMPARATIVE STATISTICS. This tool describes ONE
       specimen. There is no significance test, ordination, clustering, or
       heatmap here, because all of those compare specimens and are undefined for
       n = 1. RSCU here is a descriptive profile of that specimen's codon usage.
       Comparing several specimens (with a phylogeny, family tests, pPCA and
       Ward clustering) is outside the scope of this tool.

  2.5  LOW-COUNT AMINO-ACID FAMILIES. RSCU from a family with only a few codons
       is dominated by sampling noise. Rare amino acids (often Cys, Gln, Trp) may
       have very few codons in a single specimen, so their RSCU can swing widely
       on one or two counts. Read the Count and Family_Total columns of the RSCU
       CSV alongside the RSCU value. The program prints a caution if the pooled
       sense-codon total is under 300; adding more genes raises the counts and
       stabilizes every RSCU value.

  2.6  NOT MEDICAL/LEGAL/CLINICAL SOFTWARE. This is a research analysis tool.
       Validate results independently before drawing biological conclusions.


-------------------------------------------------------------------------------
3. INPUTS
-------------------------------------------------------------------------------
  ONE FASTA file = ONE specimen. Each ">" record is one protein-coding gene of
  that specimen. Up to 13 mitochondrial protein-coding genes is typical, but the
  program is agnostic to how many are supplied (2 genes or 13 — it pools whatever
  is there).

  Accepted extensions : .fasta, .fa, .fna, .fas, .txt (case-insensitive).
  Format detection    : a file is treated as sequence input only if its first
                        non-empty line begins with ">".
  Encoding            : UTF-8; a UTF-8 byte-order mark (BOM) is tolerated, and
                        Windows CRLF line endings are handled.
  Bases               : upper/lower case both fine; U is converted to T.

  Record headers. Each header is read as "ID gene":
     >Lycosa_shansia_OK032619 ATP6
      ^^^^^^^^^^^^^^^^^^^^^^^^  ^^^^
      first token (specimen)   rest of header (gene name)
  The GENE NAME is the header text after the first whitespace-delimited token; a
  single-token header is used as the gene name itself. The SPECIMEN NAME is
  chosen as follows (used for the output filenames and the chart title):
     (a) SPECIMEN_LABEL / --specimen-label, if you set one; otherwise
     (b) the shared first token, if every record header starts with the same one
         (the example headers all start with "Lycosa_shansia_OK032619", so that
         is used); otherwise
     (c) the input file-name stem.


-------------------------------------------------------------------------------
4. USER OPTIONS  (top of the .py; each has a matching command-line flag)
-------------------------------------------------------------------------------
  INPUT_FILE            (--input-file)   Path to the one specimen FASTA; None =
                                         auto-discover in INPUT_DIR.
  INPUT_DIR             (--input-dir)    Folder to scan when INPUT_FILE is None
                                         ("." = the folder you run from).
  SPECIMEN_LABEL        (--specimen-label) Name for the specimen; None = derive
                                         it (Section 3).
  SHOW_TITLES           (--titles on|off)  Figure title on/off.
  OUTPUT_DIR            (--output-dir)   Where outputs go; None = beside the
                                         input file.
  OPEN_CHART_WHEN_DONE  (--open-chart on|off)  Open the PNG in the OS viewer when
                                         finished.

  Command-line values override the in-file defaults. Under IDLE you normally set
  the in-file options (or just rely on the folder scan) and press F5; under
  Ubuntu/WSL you normally pass flags.


-------------------------------------------------------------------------------
5. HOW TO RUN  (quick reference)
-------------------------------------------------------------------------------
  Windows 11 / IDLE (default):
     Put the script beside your FASTA, open it in IDLE, press F5. Answer the
     one-line file prompt only if the folder has more than one FASTA.

  Ubuntu / WSL:
     python3 rscu_single_specimen_v1_0_1_2026_11_01.py \
         --input-file /path/to/specimen.fasta
     Add --titles off for a title-free publication figure; --help for all flags.


-------------------------------------------------------------------------------
6. METHOD DETAIL
-------------------------------------------------------------------------------
  STEP 1 — Parse the specimen and count codons.
     All of the specimen's gene records are pooled. Each record is read as a
     gene; U is converted to T and the sequence is upper-cased. The sequence is
     walked in non-overlapping triplets from position 0 up to the last complete
     codon, and each codon is classified as:
        SENSE     : belongs to one of the 22 Genetic-Code-5 families
                    (Osawa et al. 1992)
        STOP      : TAA or TAG
        AMBIGUOUS : anything else (e.g. codons containing N, or gaps)
     The 1-2 leftover bases after the last complete codon are then examined: a
     trailing "T" or "TA" is counted as one complete TRUNCATED terminal stop
     codon (completed to TAA by polyadenylation; Ojala, Montoya & Attardi 1981)
     and recorded in the Truncated_Stop column; any other leftover flags the gene
     "Partial". Every TAA/TAG is then checked for its position: the terminal
     stop position is the last complete codon, unless the gene ends in a
     truncated stop, in which case the truncated stop is the terminal stop and
     every complete TAA/TAG in the gene is internal. A TAA/TAG anywhere else is
     counted as an INTERNAL stop (Internal_Stops column) and triggers a console
     warning (new in v1.0.1). Internal stops are still tallied in Stop_Codons
     and never enter RSCU. This sense/stop/internal/ambiguous/truncated/partial
     bookkeeping is a data-quality diagnostic that lets you spot frame or
     contamination problems while correctly handling the abbreviated stop
     codons typical of mtDNA.

  STEP 2 — Compute RSCU per codon (Sharp et al. 1986; Sharp & Li 1987).
     For each amino-acid family with n synonymous codons and family total count
     F:
            RSCU(codon) = count(codon) * n / F         (0 if F = 0)
     RSCU normalizes raw counts by amino-acid usage and family size, isolating
     synonymous-codon preference from amino-acid composition (Sharp et al. 1986;
     Sharp & Li 1987). Stop codons — including truncated ones — are not part
     of any family and do not enter RSCU. Within each family the RSCU values sum
     to n (the number of synonymous codons), which is a handy internal check
     on the table.

  STEP 3 — Draw the RSCU graph.
     A stacked-bar chart with one bar per amino-acid family, families ordered by
     total codon count (most-used family first). Within each bar the synonymous
     codons are stacked by rank (most- to least-used codon of that family), each
     rank in a fixed colorblind-safe color. Rank 1 is always at the bottom of
     the bar and the stacking order always matches the figure legend. When two
     codons of a family have the same count, which one takes the higher rank
     (and color) is arbitrary. Codons that were never used (count 0) still
     appear in the codon key, as a zero-height segment in the bar. A dashed
     reference line marks RSCU = 1.0 (the "used-as-expected" level). Below the
     bars a family-label row names each family, and a codon key maps every
     colored segment to its codon.
     The figure is rendered with matplotlib (Hunter 2007) at 300 dpi. (The rank
     color palette auto-extends if a family ever has more than four codons.)

  Underlying numerical work uses NumPy (Harris et al. 2020) and pandas (McKinney
  2010).


-------------------------------------------------------------------------------
7. OUTPUT FILES  (in  <input folder>/RSCU_SingleSpecimen_Output_v1.0.1/ )
-------------------------------------------------------------------------------
  All names are version+date stamped ( {specimen}_..._v1.0.1_2026-11-01.<ext> ).

     {specimen}_rscu_data_v1.0.1_2026-11-01.csv
         The RSCU VALUES. Columns: Family, Codon, Count, Family_Total, RSCU.

     {specimen}_diagnostic_v1.0.1_2026-11-01.csv    (one row per gene)
         Columns: Gene, Length_bp, Sense_Codons, Stop_Codons, Internal_Stops,
         Ambiguous, Truncated_Stop, Partial. Your data-quality record.
         Internal_Stops should be 0 for every gene (Section 2.2).

     {specimen}_summary_v1.0.1_2026-11-01.csv       (one row)
         Columns: Specimen, N_Genes, Total_bp, Total_Sense, Total_Ambiguous,
         Truncated_Stops, Internal_Stops, Partial_Genes.

     {specimen}_rscu_chart_v1.0.1_2026-11-01.png
         The RSCU GRAPH (300-dpi stacked-bar chart described in Step 3).

  HOW TO READ THEM. The rscu_data CSV is the primary result; sort or filter it to
  see which synonymous codons a family prefers (RSCU > 1) or avoids (RSCU < 1),
  and always glance at Count / Family_Total to judge how many codons a value
  rests on (Section 2.5). The chart is the same information at a glance. The
  diagnostic and summary CSVs confirm the input was clean, in-frame CDS.


-------------------------------------------------------------------------------
8. WORKED EXAMPLE  (Lycosa shansia, GenBank OK032619)
-------------------------------------------------------------------------------
  The bundled file example/Lycosa_shansia_OK032619.fasta holds the 13
  mitochondrial protein-coding genes (ATP6, ATP8, COX1-3, CytB, ND1-6, ND4L) of
  the wolf spider Lycosa shansia, from the published complete mitogenome GenBank
  OK032619 (Ye et al. 2022). Every header has the form
  ">Lycosa_shansia_OK032619 GENE". Running the program on it (default settings)
  reports:

     Specimen  : Lycosa_shansia_OK032619   (derived from the shared header token)
     Genes     : 13
     Total bp  : 10,806
     Sense     : 3,590 codons
     Ambiguous : 0
     Trunc.stop: 2      (ND4 and ND4L end in an abbreviated stop; Ojala et al.
                         1981)
     Internal  : 0
     Partial   : 0

  The RSCU values land in Lycosa_shansia_OK032619_rscu_data_v1.0.1_2026-11-01.csv
  and the graph in Lycosa_shansia_OK032619_rscu_chart_v1.0.1_2026-11-01.png. As a
  sanity check, each amino-acid family's RSCU values sum to the number of
  synonymous codons in that family (e.g. the four Ala codons sum to 4.0, and the
  two Phe codons to 2.0), and the family totals sum to the 3,590 sense codons
  reported above. A clean run of a good-quality mitogenome looks like this: zero
  ambiguous codons, zero internal stops, zero partial genes, and a small number
  of truncated (poly-A-completed) terminal stops — here two. The four output
  files from this run are included in the repository under
  example/RSCU_SingleSpecimen_Output_v1.0.1/ so you can compare your own run.


-------------------------------------------------------------------------------
9. IN-TEXT CITATION MAP  (method -> reference)
-------------------------------------------------------------------------------
  * RSCU definition and computation ......... Sharp, Tuohy & Mosurski (1986);
                                               Sharp & Li (1987)
  * Invertebrate mitochondrial code (code 5)  Osawa et al. (1992)
  * Truncated (poly-A-completed) stop codons  Ojala, Montoya & Attardi (1981)
  Software:
  * Array / numerical computation ........... Harris et al. (2020)   [NumPy]
  * Data frames / CSV tables ................ McKinney (2010)        [pandas]
  * Figure rendering ........................ Hunter (2007)          [Matplotlib]
  Data standard:
  * Genetic Code Table 5 (reference table) .. NCBI (web resource; no DOI)
  Example data:
  * Lycosa shansia mitogenome (OK032619) .... Ye et al. (2022)


-------------------------------------------------------------------------------
10. REFERENCES  (with DOIs; every DOI below was checked and resolves)
-------------------------------------------------------------------------------

--- Codon usage (RSCU) ---

  Sharp PM, Tuohy TMF, Mosurski KR. 1986. Codon usage in yeast: cluster
     analysis clearly differentiates highly and lowly expressed genes. Nucleic
     Acids Research 14(13):5125-5143. doi:10.1093/nar/14.13.5125

  Sharp PM, Li W-H. 1987. The codon adaptation index — a measure of directional
     synonymous codon usage bias, and its potential applications. Nucleic Acids
     Research 15(3):1281-1295. doi:10.1093/nar/15.3.1281

--- Genetic code (invertebrate mitochondrial, NCBI Table 5) ---

  Osawa S, Jukes TH, Watanabe K, Muto A. 1992. Recent evidence for evolution of
     the genetic code. Microbiological Reviews 56(1):229-264.
     doi:10.1128/mr.56.1.229-264.1992

--- Molecular biology / mtDNA RNA processing ---

  Ojala D, Montoya J, Attardi G. 1981. tRNA punctuation model of RNA processing
     in human mitochondria. Nature 290(5806):470-474. doi:10.1038/290470a0

--- Example data ---

  Ye W, Wang J, Zhao X, Liu H, Zhu S. 2022. Mitochondrial genomes of two Lycosa
     spiders (Araneae, Lycosidae): genome description and phylogenetic
     implications. Diversity 14(7):538. doi:10.3390/d14070538
     (GenBank OK032619, Lycosa shansia)

--- Software ---

  Harris CR, Millman KJ, van der Walt SJ, Gommers R, Virtanen P, Cournapeau D,
     et al. 2020. Array programming with NumPy. Nature 585(7825):357-362.
     doi:10.1038/s41586-020-2649-2

  Hunter JD. 2007. Matplotlib: a 2D graphics environment. Computing in Science &
     Engineering 9(3):90-95. doi:10.1109/MCSE.2007.55

  McKinney W. 2010. Data structures for statistical computing in Python. In:
     Proceedings of the 9th Python in Science Conference, pp. 56-61.
     doi:10.25080/Majora-92bf1922-00a

--- Reference data standard (web resource; no DOI) ---

  NCBI. The Genetic Codes — Translation Table 5 (Invertebrate Mitochondrial
     Code). https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi


===============================================================================
  END OF README
===============================================================================
