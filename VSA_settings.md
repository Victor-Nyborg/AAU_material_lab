# Aqualab Vapor Sorption Analyzer - Settings to use for measurements
**Last updated: 09-04-2026**

This document is to standardize the procedure for performing isotherm sorption measurements with the Aqualab Vapor sorption analyzer, VSA.

It tells the reader which settings to use for performing either the DVS or DDI method. In the future there might be suggested setting depending on type of material. 

For both methods the initial moisture content should be set to *Dry Basis* and *0.0%*, the true moisture content will be calculated in the data treatment script, [VSA.py](VSA.py).
## DVS
The DVS method, measure the moisture content for a given ambient water activity when equilibrium has been achieved.


## DDI
The DDI method continuously measure the moisture content at a given water activity without seeking equilibrium with the ambient relative humidity. This method is specific for the VSA.

### Recommended general settings

| Setting                  | Value |
|--------------------------|-------|
| Start (aw)               | 0.05  |
| Final (aw)               | 0.95  |
| Temperature ($\text{\textdegree}$C) | 23    |
| Resolution (aw)          | 0.01  |
| Flow (ml/min)            | 100.0 |
| Timeout                  | Off   |
| Loop                     | On    |
