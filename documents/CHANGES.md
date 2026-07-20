# Manuscript text corrections — ECS25-0887.R1

Source: `Main Document.docx` (accepted version)
Output: `Main_Document_corrected.docx` (+ `Main_Document_corrected.pdf` render)

All edits are **surgical run-level replacements**. Superscript/subscript units,
paragraph styles, and the 86 Word CITATION reference fields are untouched
(verified: 86 fields before and after).

## Driver: figures were flattened to sequential numbering (Figures 1–10)

| Old figure | New figure(s) |
|---|---|
| Fig 1 study area | Figure 1 (unchanged) |
| Fig 2 schematic | Figure 2 (unchanged) |
| Fig 3 temporal traits | Figure 3 (unchanged) |
| Fig 4 seasonal histograms (Early/Mid/Late in one figure) | Figure 4 (Early), Figure 5 (Mid), Figure 6 (Late) |
| Fig 5 spatial trait maps (trait \| PFT \| difference) | Figure 7 (trait + PFT), Figure 8 (difference) |
| Fig 6 GPP/SIF diff maps (a,b) + time series (c,d) | Figure 9 (maps a,b), Figure 10 (time series a,b) |

## Body callout edits

- §3.1 LWC para: added cross-reference **(Figure 3d)** for parallelism with LAI (3a) and CHL/LMA (3b,3c).
- "Figure 4 shows seasonal histograms…" → "**Figures 4–6** show…", with each season tagged (Figure 4 / 5 / 6).
- "Figure 5 shows the spatial distribution…" → "**Figures 7 and 8** show…".
- "left panel" → "**left column of Figure 7**"; "middle panels" → "**right column of Figure 7**".
- "The right panels, quantifying…offer" → "**Figure 8**, quantifying…offers".
- "Figure 6 offers" → "**Figures 9 and 10** offer".
- "Figure 6a and Figure 6b" → "**Figure 9a and Figure 9b**"; "Similarly to Figure 5" → "**Similarly to Figure 8**".
- "time series graphs below the spatial maps in Figure 6c and Figure 6d" → "**time series in Figure 10a and Figure 10b**".
- "Figure 6d serves" → "**Figure 10b** serves"; "(Figure 6d)" → "**(Figure 10b)**".
- "overlapping areas of Figure 4" → "**overlapping areas of Figures 4–6**".
- "(Fig.4 a,d,g)" → "**(Figs. 4a, 5a, 6a)**"; "(Fig. 4b,e,h)" → "**(Figs. 4b, 5b, 6b)**"; "(Fig. 4c,f,i)" → "**(Figs. 4c, 5c, 6c)**".
- "trait maps (Figure 5)" → "**trait maps (Figure 7)**".

Unchanged figure references (correctly point to same figures): Figure 1 (study area, ×3 + topography callout), Figure 2 (schematic ×2), Figure 3a/3b/3c.

## Caption edits (split to match new numbering)

- **Figure 4** caption → split into Figure 4 (Early), Figure 5 (Mid), Figure 6 (Late), each with panels (a) CHL, (b) LMA, (c) LWC.
- **Figure 5** caption → Figure 7 (trait + PFT columns) and Figure 8 (trait − PFT differences).
- **Figure 6** caption → Figure 9 (GPP/SIF difference maps a,b) and Figure 10 (GPP/SIF time series a,b; old panels c,d relettered to a,b).

## Deliberately NOT changed

- **Citations.** The duplicated citations / stray parentheses seen in the journal's typeset proof
  (e.g. "Nature Conservancy 2022" twice, "Clark 2016" twice, a stray "(") **do not exist in this
  .docx** — they were field-rendering artifacts of the typeset proof. Verified: each citation
  appears once. No citation text was edited.
- **§3.1 LWC absolute values** (0.0057, 0.0012, 0.0053, 0.00092, 0.0100, 0.0026 g cm⁻²) left as
  published — they are numerically correct; only figure axes use the ×10⁻³ compact notation.
  Say the word if you want the prose converted to ×10⁻³ for notation consistency.

## Editorial-email fixes applied (second pass)

- **Item #5 — Table 1 relocated.** The Table 1 caption + table object were moved from inline
  in Methods to **after the References list** (immediately before the Figure captions block).
  The in-text "(Table 1)" callout stays in Methods.
- **Item #12 — Supporting Information references reformatted** to the Ecosphere form
  "**Appendix S1: Figure S#**" / "Appendix S1: Section S1" / "Appendix S1" (17 references
  across the body). No residual "Fig. S#", "Supporting Information", or "SI (…)" strings remain.

## Still outstanding for you (from the editorial email)

- **#1** Zenodo DOI for the GitHub code (statement currently cites a release tag), and add
  reference-list citations for the CaltechDATA dataset (`10.22002/7xgrn-qtc49`, confirm final/public)
  and the code. *(Hard blocker for file conveyance.)*
- **#2** Department name for affiliation 1 (Caltech); optional: add "CA" to affiliation 2.
- **#3** Number the four equations as Equation 1–4 (whole numbers) and add in-text references
  — best done in Word (they are math objects).
- **#7 / #9** Figure-file legibility: Fig 2 (Trait/Flux panels), Fig 7 scale numbers, Fig 9 insets
  — verify readable at 18 cm width.
- **#4** Final references ↔ in-text citation completeness pass (86 citation fields).
- **#11** Assemble Appendix S1 as a single PDF (`Appendix_S1`), fix running head, add
  journal/title/authors on page 1.

## Second review pass (editorial polish)

- **Removed 3 duplicated citations** that were hidden inside the Mendeley fields (only visible
  after flattening): "(The Nature Conservancy, 2022)", "(Kelley, 1999; Nocedal & Wright, 1999)",
  "(Clark, 2016)".
- **Stray symbols fixed:** "(Farquhar et al., 1980)(," → "),"; "G(θ) ((Ross, 1981)" → "(Ross, 1981)".
- **Double spaces** collapsed throughout (equation-alignment spacing preserved).
- **Affiliation 1:** added "Division of Geological and Planetary Sciences".
- **Figure captions** now each define their abbreviations (Fig 1 PFT; Figs 5/6/8 CHL/LMA/LWC/PFT;
  Fig 9 PFT; Fig 10 GPP/SIF740/TROPOMI). Superscripts/subscripts verified intact.

Still pending your decision:
- **Equations** are numbered (1.1)/(1.2)/(2.1)/(2.2); editor wants whole numbers (1)–(4), plus
  update the one in-text reference "Eq. (2.2)" → "Eq. (4)".
- **Affiliation 2** (JPL) is missing "CA".

## Tracked-changes version

This is a clean corrected file. For a tracked-changes copy, open Word →
**Review ▸ Compare** with `Main Document.docx` as original and
`Main_Document_corrected.docx` as revised — this yields accurate, reviewable markup.
