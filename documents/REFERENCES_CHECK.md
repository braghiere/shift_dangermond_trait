# Reference & citation check + corrections — ECS25-0887.R1 (editorial item #4)

Status: **DONE and independently verified to consensus.**
Every factual correction is grounded in **CrossRef** authoritative metadata (not model
memory). An independent verification agent re-checked all 104 entries against fresh CrossRef
pulls over two rounds and **concurs the list is correct**.

Mechanics: the bibliography and in-text citations were **Mendeley fields**. They have been
**flattened to static text** (formatting preserved; text was byte-identical before edits), so
the citations no longer auto-update — appropriate for final submission. Do NOT re-open/refresh
in Mendeley after this, or it may regenerate and undo the fixes.

## 1. In-text ↔ list matching — PASS
104 list entries; 104 distinct works cited; every in-text citation matches a list entry and
vice versa. No orphans, no uncited entries. (Author-name fixes below were applied to both the
list and the in-text citations, so matching is preserved.)

## 2. Corrections applied (all CrossRef-verified)

### Author-name errors (caught by the independent verifier / particle scan)
- **Cranko Page, J.** (2022 & 2024) — Mendeley mis-parsed the compound surname as "Page, J. C."
  Fixed in the list **and** in-text ("Page et al." → "Cranko Page et al."); also "De Kauwe, M. G."
  restored. The two entries were re-alphabetized into the C section (Cornwell < Cranko Page < Croft).
- **Farquhar, G. D., von Caemmerer, S., & Berry, J. A.** (1980) — "von" had been dropped.
- **Hakkenberg, C. R., & Goetz, S. J.** (2021) — removed a garbled duplicate 3rd author
  ("Christopher Hakkenberg, C. R."); CrossRef confirms only 2 authors.

### Wrong/missing bibliographic data
- **Enquist et al. (2019)** — journal was "Ary T. Oliveira-Filho"; corrected to
  *Science Advances*, 5(11), eaaz0414, https://doi.org/10.1126/sciadv.aaz0414.
- **Kampe et al. (2010)** — journal was a DOI URL; corrected to *Journal of Applied Remote
  Sensing*, 4(1), 043510.
- **Missing article numbers added:** Green 2022 → 20220071; Wang, H. 2023 → eadd5667;
  O'Sullivan 2022 → 4781; Wang, Y. 2025 → 4968; Wang, Y. 2022 → 258; Randerson 2025 → eadr5489.
- **Miller (1967)** — added pages **141–144** (CrossRef; my earlier estimate of 141–150 was wrong).
- **Alton (2011)** — citation number "1030" → **G01030**.

### Book/report publishers
- **Kelley (1999)** — de-duplicated title; added **Society for Industrial and Applied
  Mathematics, Philadelphia, PA**.
- **Nocedal & Wright (1999)** — **Springer-Verlag, New York** (1st ed., DOI 10.1007/b98874;
  the "Second Edition" label was removed as inconsistent with the 1999 date — confirm if you
  meant the 2006 2nd edition).
- **Ross (1981)** — **Dr. W. Junk Publishers, The Hague, Netherlands** (was "Junk").

### Journal-name cleanup (Mendeley export artifacts)
Blinn → *Forests*; Guanter → *Remote Sensing*; Wright → *Nature*; Cardoso → *npj Biodiversity*
(embedded "YYYY, Vol. X, Page Y" removed).

### DOI cleanup
Stripped Mendeley junk suffixes (Deck, Garnier, Kattge 2020, Randerson, Singh),
supplement-file links (Harrison 2020, Reichstein, Osnas), `/BIBTEX` (Sousa), and the
doubled prefix (Brodrick data DOI 10.3334/ORNLDAAC/2183).

### Corporate author
- **California Dept. of Forestry and Fire Protection (CAL FIRE), Fire and Resource Assessment
  Program (FRAP)** (2015 & 2018) — was garbled.

### Minor
- **Best et al. (2011)** — repository URL replaced with DOI 10.5194/gmd-4-677-2011.

## 3. Not errors — left as-is (verified defensible)
- **Online-first vs print year:** Cai 2023, Kattge 2020, Piao 2020, Cranko Page 2024 — the entry
  year matches the print/volume year and the in-text citation; CrossRef's earlier "issued"
  (online) date is not used.
- **Hancock (2019)** uses e-locator 2018EA000506 (valid; CrossRef also lists pp. 294–310).

## 4. Still needs your input (non-CrossRef sources — not blockers)
- **Tang & Armston (2019)** GEDI L2B ATBD — grey literature; add a stable URL
  (GEDI ATBDs are distributed via gedi.umd.edu / LP DAAC). Institution added tentatively.
- **CAL FIRE/FRAP (2015, 2018)** and **The Nature Conservancy (2025)** — agency/web sources;
  URLs are live but not independently cross-checkable.
- **Data/code citations for Open Research (item #1)** still to be added once the final
  CaltechDATA + Zenodo DOIs are confirmed (see CHANGES.md).

## 5. Verification record
- Round 1: 29 corrections applied → independent agent found 4 author errors I missed.
- Round 2: those + Farquhar fixed and re-alphabetized → independent agent **concurs, no new
  errors**; co-author spot-check on 14 particle-heavy entries all match CrossRef.
- Automated dropped-particle scan across all 104 entries: **0 flags**.
