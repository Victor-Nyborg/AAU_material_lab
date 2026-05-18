# Aqualab Vapor Sorption Analyzer - Settings to use for measurements
**Last updated: 12-05-2026**

This document is to standardize the procedure for performing isotherm sorption measurements with the Aqualab Vapor sorption analyzer, VSA.

It tells the reader which settings to use for performing either the DVS or DDI method. In the future there might be suggested setting depending on type of material. 

For both methods the initial moisture content should be set to *Dry Basis* and *0.0%*, the true moisture content will be calculated in the data treatment script, [VSA.py](VSA.py).
The mass of the sample should preferably be above 1 gram if possible to reduce the affect from scale-noice. Light materials such as fibrous isolation materials benefits of being cutted into small pieces.

## DVS
The DVS method, measure the moisture content for a given ambient water activity when equilibrium has been achieved. It more or less follows the ISO 12571 standard.

Temperature for testing should be set to 23 C

## DDI
The DDI method continuously measure the moisture content at a given water activity without seeking equilibrium with the ambient relative humidity. This method is specific for the VSA.

### Recommended general settings

| Setting                  | Value |
|--------------------------|-------|
| Start (aw)               | 0.05  |
| Final (aw)               | 0.95  |
| Temperature (C) | 23    |
| Resolution (aw)          | 0.01  |
| Flow (ml/min)            | 100.0 |
| Timeout                  | Off   |
| Loop                     | On    |
