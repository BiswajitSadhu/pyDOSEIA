# pyDOSEIA

## 📘 Summary of Reference Data Sources Used in pyDOSEIA

| **Category**    | **Source Description**                                                                                     | **Reference**                     |
|----------------|-------------------------------------------------------------------------------------------------------------|----------------------------------|
| Inhalation      | Effective Dose Coefficients for Inhalation of Radionuclides for Members of the Public                      | ICRP Publication 119 (2012)      |
| Inhalation      | Dose Coefficients for Soluble or Reactive Gases and Vapours                                                | ICRP Publication 119 (2012)      |
| Inhalation      | Effective Dose Coefficients from Inhaled Air                                                                | DOE-STD-1196-2011                |
| Inhalation      | Inhalation Dose Coefficients to age 70 years (Sv/Bq)                                                       | JAERI-Data/Code 2002-013         |
| Inhalation      | Dose Coefficients for Soluble or Reactive Gases and Vapours (Class SR-1 and SR-2)                          | JAERI-Data/Code 2002-013         |
| Ground Shine    | Dose Rate Coefficients for Ground Surface (Sv/Bq/s/m²)                                                      | FGR-15 (USDOE)                   |
| Submersion      | Dose Rate Coefficients for Air Submersion (Sv/Bq/s/m³)                                                      | FGR-15 (USDOE)                   |
| Ingestion       | Effective Dose Coefficients for Ingestion to 70 years (Sv/Bq)                                               | ICRP Publication 119 (2012)      |
| Ingestion       | Ingestion Dose Coefficients to age 70 years (Sv/Bq)                                                         | JAERI-Data/Code 2002-013         |
| Half-Life       | Radionuclide Properties: ICRP-07 Collection                                                                 | ICRP Publication 107 (2008)      |
| Half-Life       | Dose Coefficient DB: Ingestion and Inhalation of Particulates                                              | JAERI-Data/Code 2002-013         |
| Half-Life       | Dose Coefficient DB: Inert Gases                                                                            | JAERI-Data/Code 2002-013         |
| Gamma Emission  | Recommended Gamma Ray Energies and Emission Probabilities by Radionuclide                                  | IAEA (2007) Decay Data           |

> **Note:** *pyDOSEIA optionally includes progeny contribution via the defined equation; otherwise, it returns unmodified dose coefficients.*


ℹ️ pyDOSEIA allows incorporating the contribution from progeny (see Eq. #ref:progeny) if the user chooses to do so; otherwise, it returns unmodified dose coefficients.

📚 References
<a name="ref1">[1]</a> ICRP Publication 119, 2012.
<a name="ref2">[2]</a> DOE-STD-1196-2011, U.S. Department of Energy.
<a name="ref3">[3]</a> JAERI-Data/Code 2002-013, Japan Atomic Energy Research Institute.
<a name="ref4">[4]</a> FGR-15, U.S. Department of Energy.
<a name="ref5">[5]</a> ICRP Publication 107, 2008.
<a name="ref6">[6]</a> IAEA, 2007. Decay Data for Radionuclides.


### Ingestion Data Tables Used in pyDOSEIA

| **Ingestion Data Used in pyDOSEIA** | **Reference** |
|-------------------------------------|---------------|
| Table VII: Conservative Values for Mass Interception and Environmental Removal Rates from Plant Surfaces | IAEA SRS 19 [^1] |
| Table VIII: Conservative Values for Crop and Soil Exposure Periods and Delay Times | IAEA SRS 19 [^1] |
| Table IX: Effective Surface Soil Density for Screening Purposes | IAEA SRS 19 [^1] |
| Table X: Loss Rate Constant Values for Screening Purposes | IAEA SRS 19 [^1] |
| Table XI: Element-Specific Transfer Factors for Terrestrial Foods for Screening Purposes | IAEA SRS 19 [^1] |

[^1]: IAEA SRS 19


