# Changelog

## Documentation updates after 1.0.1 (November 2026)

No change to the program's calculations or outputs.

- Full reference README (`docs/`): replaced template placeholders with the author, repository, DOI and citation; removed a citation (Revell 2009) that had no reference entry; removed a pointer to a batch engine that is not part of this repository.
- Script header: author, repository and DOI added.
- Example Methods Word file: author and title metadata set.

## 1.0.1 (November 2026)

- Internal-stop check: a TAA/TAG before the gene's terminal stop is counted per gene (`Internal_Stops` column in the diagnostic CSV and a total in the summary CSV) and triggers a console warning naming the affected genes.
- Sharp, Tuohy & Mosurski (1986) cited as the original RSCU reference, alongside Sharp & Li (1987).
- Removed the unverified claim that results are byte-for-byte identical to the batch engine `rscu_integrated`.
- Documentation: corrected the RSCU sanity check in the worked example (each family's RSCU values sum to its number of synonymous codons); described chart stacking order, tie ranking, and unused codons in the codon key.
- Worked example: the 13 protein-coding genes of *Lycosa shansia* (GenBank OK032619; Ye et al. 2022) added in `example/`, with the output files from running the program on them.
- Example Methods paragraph for the worked example, on the front page and as `docs/RSCU_Methods_example_Lycosa_shansia.docx`.
- "How to cite this" section and `CITATION.cff` (GitHub's "Cite this repository" button).
- Archived on Zenodo with DOI [10.5281/zenodo.24182022](https://doi.org/10.5281/zenodo.24182022); DOI added to the README, citation and `CITATION.cff`.
- RSCU values are unchanged from 1.0.0.

## 1.0.0 (August 2026)

- Initial release: single-specimen RSCU under NCBI Genetic Code 5, truncated-stop handling, diagnostic and summary tables, and the stacked-bar RSCU graph.
