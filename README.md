[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![AERONET](https://img.shields.io/badge/AERONET-blue)](https://aeronet.gsfc.nasa.gov/cgi-bin/webtool_aod_v3)
[![Google Earth Enginer](https://img.shields.io/badge/Google_Earth_Engine-orange)](https://developers.google.com/earth-engine/datasets/catalog)

## Introduction

This python module downloads the AERONET AOD data and the reflectance data from GEE for a given time period and a given bounding box.

## How to use
```python
    aeronet_df = extract_aeronet_and_reflectance(gee_product_id = 'LANDSAT/LC08/C02/T1_TOA', 
                                             start_date = '2021-10-01', 
                                             end_date = '2021-12-31',
                                             spectral_bands = ['B1', 'B5'],
                                             scale = 10, 
                                             aeronet_level = 1.5, 
                                             temporal_average = 30,
                                             bbox = [10.0, 60.0, 35.0, 90.0], 
                                             dest_folder = os.getcwd())
```
