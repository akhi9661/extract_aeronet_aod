[![License](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![AERONET](https://img.shields.io/badge/AERONET-blue)](https://aeronet.gsfc.nasa.gov/cgi-bin/webtool_aod_v3)
[![Google Earth Enginer](https://img.shields.io/badge/Google_Earth_Engine-orange)](https://developers.google.com/earth-engine/datasets/catalog)

## Introduction

This python module downloads the AERONET AOD data and the reflectance data from GEE for a given time period and a given bounding box.

## How to use
```python
path = r'D:\Aerosol Modelling\Aerosol\Scripts\AERONET_module'
aeronet_df = extract_aeronet_and_reflectance(gee_product_id = 'COPERNICUS/S2_SR_HARMONIZED', 
                                             start_date = '2018-01-01', 
                                             end_date = '2018-12-31',
                                             spectral_bands = ['AOT'],
                                             scale = 10, 
                                             aeronet_level = 1.5, 
                                             temporal_average = 30,
                                             bbox = [-60, -180, 60, 180], 
                                             dest_folder = os.getcwd(), 
                                             chunk_size = 1000,
                                             verbose = False, 
                                             aeronet_site = pd.read_csv(os.path.join(path, 'aeronet_sites_2018.csv')), 
                                             reflectance_df = pd.read_csv(os.path.join (path, 'S2_SR_HARMONIZED_2018-01-01_2018-12-31.csv')))
```
