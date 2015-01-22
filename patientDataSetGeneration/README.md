matRad uses open CT datasets provided by Craft et. al (2014). There are four patients  avaible: A prostate case, a liver case, a head and neck case, and a standart IMRT phantom. All information regarding to this data it is avaible on http://www.gigasciencejournal.com/content/3/1/37

---

# INSTRUCTIONS

The following instructions it is for generate patient data for matRad using the open dataset.

1. Create a empty folder to generate patient data.
2. Download LIVER.zip, PROSTATE.zip, TG119.zip and HN\_withoutDij.zip files from http://gigadb.org/dataset/100110 server.
3. Extract like a folder every downloaded data set.
4. Download matRad code for generate patient: [matRad_genDataSet.m](matRad_genDataSet.m), [matRad_genCT.m](matRad_genCT.m) and [matRad_genCst.m](matRad_genCst.m).
5. Download CST input files: [LIVER_CST_Input.m](LIVER_CST_Input.m), [PROSTATE_CST_Input.m](PROSTATE_CST_Input.m), [TG119_CST_Input.m](TG119_CST_Input.m) and [HN_withoutDij_CST_Input.m](HN_withoutDij_CST_Input.m).
6. Edit CST input files. Column 3: organ class: OAR, TARGET or IGNORED. Use only Upper case. Column 4: maximum dose [Grays]. Column 5: minimum dose [Grays]. Column 6: maximum penalty, and Column 7: minimum penalty.
7. Download Hounsfield units to water equivalent look up table [HU2waterEqT.mat](HU2waterEqT.mat).
8. Run matRad_genDataSet.m code with Matlab to generate data for matRad.
