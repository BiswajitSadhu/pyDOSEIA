📘 Summary of Reference Data Sources Used in pyDOSEIA
Category	Description	Reference
Inhalation	Effective dose coefficients for inhalation of radionuclides	ICRP Publication 119 [1]
Inhalation	Dose coefficients for soluble/reactive gases and vapours	—
Inhalation	Effective dose from inhaled air	DOE-STD-1196-2011 [2]
Inhalation	Dose coefficients to age 70 y (Sv/Bq)	JAERI 2002 [3]
Inhalation	Dose coefficients for soluble/reactive gases (SR-1/SR-2)	—
Ground Shine	Dose rate for ground surface (Sv/Bq/s/m²)	FGR-15 [4]
Submersion	Dose rate for air submersion (Sv/Bq/s/m³)	FGR-15 [4]
Ingestion	Ingestion dose coefficients to age 70 y	ICRP 119 [1], JAERI 2002 [3]
Half-Life	Radionuclide properties (ICRP-07)	ICRP 107 [5]
Half-Life	Radionuclides in dose coefficient DB	JAERI 2002 [3]
Gamma	Gamma ray energies & emission probabilities	IAEA 2007 [6]

ℹ️ pyDOSEIA allows incorporating the contribution from progeny (see Eq. #ref:progeny) if the user chooses to do so; otherwise, it returns unmodified dose coefficients.

📚 References
<a name="ref1">[1]</a> ICRP Publication 119, 2012.
<a name="ref2">[2]</a> DOE-STD-1196-2011, U.S. Department of Energy.
<a name="ref3">[3]</a> JAERI-Data/Code 2002-013, Japan Atomic Energy Research Institute.
<a name="ref4">[4]</a> FGR-15, U.S. Department of Energy.
<a name="ref5">[5]</a> ICRP Publication 107, 2008.
<a name="ref6">[6]</a> IAEA, 2007. Decay Data for Radionuclides.


