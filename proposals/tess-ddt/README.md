# TESS DDT request: TIC 402710280, Sectors 113–114

This is a Director's Discretionary Targets (DDT) proposal for two-minute cadence on TIC 402710280. It follows the official templates at <https://tess.mit.edu/science/ddt/>.

| File | What it is |
|---|---|
| `proposal.tex` → `Vora_DDT_TIC402710280_S113-114.pdf` | Two-page proposal, the template's limit (text, figure and references). |
| `Vora_DDT_TIC402710280_S113-114_targets.csv` | Target list in the template's CSV format. |

## Facts behind the request

- **Not already scheduled.** TIC 402710280 is not on the GI 2-minute or 20-second target lists for S113 or S114 (`GI_S113.pri`, `GI_S114.pri`, `GI-20s_S11{3,4}.pri`, checked 2026-10-18).
- **On silicon.** tess-point places the star on camera 4 / CCD 4 (col 385, row 209) in S113 and camera 4 / CCD 1 (col 646, row 475) in S114.
- **No later TESS coverage.** The star is not observed again through S134, the last sector tess-point knows.
- **Transit timing.** Transit 57 is predicted for 2027-02-19 11:29 UTC (BJD_TDB 2461455.9785 ± 30 min), inside S114 (13 Feb–14 Mar 2027). Transit 56 (14 Jan) falls before S113 starts, and transit 58 (27 Mar) after S114 ends.
- **Deadlines.** The DDT rule is six weeks before the sector starts. For S113 (starts 17 Jan 2027) that is **6 Dec 2026**; for S114 (starts 13 Feb 2027) it is **2 Jan 2027**. Submit before 6 Dec so that both sectors are covered.

## Email

**To:** tess-ddt@mit.edu
**Subject:** DDT request: TIC 402710280, two-minute cadence in Sectors 113–114
**Attach:** `Vora_DDT_TIC402710280_S113-114.pdf`, `Vora_DDT_TIC402710280_S113-114_targets.csv`

> Dear TESS Science Office,
>
> I would like to request two-minute cadence for TIC 402710280 in Sectors 113 and 114. The proposal and target list are attached.
>
> The star, a V = 10 K dwarf at 67 pc, shows four transits of a 3.75 Earth-radius candidate with a period of 36.34 days, all in full-frame images. The period became clear only with the Sector 102 transit on 29 March 2026, after the Cycle 9 deadline, so the star could not be included in a GI proposal.
>
> The next transit is due on 19 February 2027 (BJD_TDB 2461455.9785 ± 30 min), in Sector 114. On 15 October 2026 I registered this prediction on Zenodo, together with the data to be used and the criteria for judging the result. The test was designed around SPOC two-minute light curves, which will exist only if the star is a two-minute target. tess-point shows no further TESS coverage through Sector 134, so this is the one transit TESS can record before late 2028.
>
> Preprint: https://doi.org/10.5281/zenodo.23390317
> Data and code: https://doi.org/10.5281/zenodo.23390274
> Pre-registered test: https://doi.org/10.5281/zenodo.23390325
>
> Thank you for considering the request. I would be glad to answer any questions.
>
> Best regards,
> Rishit Vora
> Independent researcher · rishitv@protonmail.com · ORCID 0009-0003-8778-9883

## After submission

- The review team assigns a DDT number and emails the proposer. Approved targets appear at <https://tess.mit.edu/public/DDTs/DDT_list.html>.
- If the request is approved, SPOC two-minute data for S114 should reach MAST about a month after the sector ends (mid-April 2027). The pre-registered test is then run on those data first, as the protocol requires.
