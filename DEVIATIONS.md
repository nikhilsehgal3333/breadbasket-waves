# Deviations from the plan

The method was fixed in the study plan (approved 17 Nov 2026) before any result was computed.
Every later change is listed here with its date and reason. Changes made before any result was
computed are plan refinements; changes made after are deviations.

| Date | Change | Reason | Before any result? |
|---|---|---|---|
| 2026-11-17 | Read NCEP/NCAR R1 over OPeNDAP with the `pydap` engine instead of `netCDF4` | The netCDF4 1.7.5 Windows wheel has DAP support, but its built-in curl cannot find an SSL certificate bundle ("Problem with the SSL CA cert"). pydap was the plan's first fallback; it returned the expected 92 x 9 x 144 JJA subset for the leap year 1980 | Yes |
| 2026-11-17 | Method refinements from an independent review: wave thresholds over 1960-2014; China = FAO 41 (mainland); EEU trend split at 1992; exact counts, Katz intervals and hypergeometric p-values per LMF cell; whole-year permutation test with four pre-set primary tests; manifest-checked numbers | Listed in full under "Plan refinements" in the study plan (`prototype/PLAN.md` in the application folder) | Yes |
| 2026-11-17 | SSA trend = first two reconstructed components (k = 2), not one; k = 1 and an 11-year moving average become sensitivities | A synthetic test showed that one component leaves up to 2.5 % of a plain linear trend in the residual, as large as the effects being measured, which would bias "below trend" at both ends of each series. Two components recover a linear trend exactly. The paper's SSA reference (Vautard, Yiou and Ghil 1992) gives no setting for yields | Yes |
| 2026-11-17 | Natural Earth outlines taken as GeoJSON from the project's GitHub release v5.1.2 instead of the shapefile zip | GeoJSON reads with the standard library; the shapefile would need an extra package. Same data | Yes |
