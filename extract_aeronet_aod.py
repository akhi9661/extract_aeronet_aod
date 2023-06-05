def process_bitmask(dataframe):
    
    import pandas as pd

    '''
    This function takes a dataframe with a column named 'QA_PIXEL' or 'QA60' and creates a new column named 'Cloud/Snow' 
    depending on the binary values of the QA_PIXEL or QA60 column. The function returns a dataframe with the new column.

    Parameters:
        dataframe (pandas dataframe): A dataframe with a column named 'QA_PIXEL' or 'QA60'.

    Returns:
        updated_dataframe (pandas dataframe): A dataframe with a new column named 'Cloud/Snow'.
    
    '''
    
    def convert_to_dictionary(qa_pixel):

        '''
        This function takes a QA_PIXEL value and converts it to a dictionary with the binary values of the QA_PIXEL as keys 
        and the binary values as values.

        Parameters:
            qa_pixel (int): A QA_PIXEL value.

        Returns:
            binary_dict (dictionary): A dictionary with the binary values of the QA_PIXEL as keys and the binary values as values.
        
        '''

        binary_str = bin(qa_pixel)[2:].zfill(16)  # Convert to binary string and pad with zeros to 16 bits
        binary_dict = {f'Bit {i}': bit for i, bit in enumerate(binary_str[::-1])}
        return binary_dict

    rows_with_bits = []
    headers = dataframe.columns.tolist()
    headers.append('Cloud/Snow')

    for index, row in dataframe.iterrows():
        if 'QA_PIXEL' in headers:
            mask = 'QA_PIXEL'
            qa_pixel = int(row['QA_PIXEL'])
            binary_dict = convert_to_dictionary(qa_pixel)
            cloud_snow = 'Cloud' if binary_dict['Bit 3'] == '1' else 'Snow' if binary_dict['Bit 5'] == '1' else ''
            
        elif 'QA60' in headers:
            mask = 'QA60'
            qa_60 = int(row['QA60'])
            binary_dict = convert_to_dictionary(qa_60)
            cloud_snow = 'Cloud' if binary_dict['Bit 10'] == '1' or binary_dict['Bit 11'] == '1' else ''
            
        else:
            cloud_snow = ''
        
        row['Cloud/Snow'] = cloud_snow
        rows_with_bits.append(row)
    
    updated_dataframe = pd.DataFrame(rows_with_bits, columns = headers)
    return updated_dataframe

def gee_point_extract(point_filename, product = 'LANDSAT/LC08/C02/T1_TOA', start_date = '2020-12-01', end_date = '2020-12-31', id_col = None, 
                      bands = ['B1', 'B2', 'B3', 'B4'], scale = 30, dest_folder = None):
    
    '''
    This function takes a point shapefile and extracts the pixel values for the specified bands from the specified sensor from 
    Google Earth Engine. The function returns a dataframe with the pixel values for each band and the latitude and longitude of
    the point. The function also returns a dataframe with the metadata [angle, QA_mask] for each point. 

    Parameters:
        point_filename (string): The filename of the point shapefile.
        product (string): The product ID of the sensor. Default is 'LANDSAT/LC08/C02/T1_TOA'. 
        start_date (string): The start date of the time period of interest. Format is 'YYYY-MM-DD'. Default is '2020-12-01'
        end_date (string): The end date of the time period of interest. Format is 'YYYY-MM-DD'. Default is '2020-12-31'
        id_col (string): The name of the column in the point shapefile that contains the unique ID for each point. Default is 'FID'.
        bands (list): A list of the bands to extract from the sensor. Default is ['B1', 'B2', 'B3', 'B4'].
        scale (int): The scale of the pixel values. Default is 30.
        dest_folder (string): The destination folder to save the output csv file. Default is None.

    Returns:
        final_df (pandas dataframe): A dataframe with the pixel, angle, and QA_mask values for each band
    
    '''
    
    import ee
    import os
    import re
    from gee_subset import gee_subset
    import pandas as pd
    from datetime import datetime
    import geopandas as gpd
    
    if not ee.data._credentials:
        ee.Authenticate()

    if not ee.data._initialized:
        ee.Initialize()
        
    opf = os.path.join(dest_folder, f'{product.split("/")[-1]}_{start_date}_{end_date}.csv')
    print(f'Processing: Fetching satellite data from "{product}" for time period: [{start_date, end_date}]\n')

    if product == 'LANDSAT/LC08/C02/T1_TOA':
        extra_bands = ['SAA', 'SZA', 'VAA', 'VZA', 'QA_PIXEL']
        for band in extra_bands:
            if band not in bands:
                bands.append(band)
        scale = 30

    if product == 'COPERNICUS/S2_HARMONIZED':
        extra_bands = ['QA60']
        for band in extra_bands:
            if band not in bands:
                bands.append('QA60')
        scale = 10

    if product is None:
        print('Enter a valid product ID')

    if isinstance(point_filename, pd.DataFrame):
        points = point_filename
        if dest_folder is None:
            opf = os.path.join(os.getcwd(), f'{product.split("/")[-1]}_{start_date}_{end_date}.csv')
    elif isinstance(point_filename, str) and point_filename.endswith('.shp'):
        points = gpd.read_file(point_filename)
        if dest_folder is None:
            opf = os.path.join(os.path.dirname(point_filename), f'{product.split("/")[-1]}_{start_date}_{end_date}.csv')
    elif isinstance(point_filename, str) and point_filename.endswith('.csv'):
        points = gpd.read_csv(point_filename)
        if dest_folder is None:
            opf = os.path.join(os.path.dirname(point_filename), f'{product.split("/")[-1]}_{start_date}_{end_date}.csv')
    else:
        print("Invalid input. Expected either a pandas dataframe or csv/shapefile path.")

    count = len(points)
    site = list(range(0, count, 1))    

    values = []

    for i in site:
        print(start_date, end_date)
        
        print(f"Extracting for {id_col}: {points.iloc[i, points.columns.get_loc(id_col)]}")
        df = gee_subset.gee_subset(product = product,
                                   bands = bands,
                                   start_date = start_date,
                                   end_date = end_date,
                                   latitude = points.iloc[i, points.columns.get_loc('lat')],
                                   longitude = points.iloc[i, points.columns.get_loc('lon')], 
                                   scale = scale)
    
        sid =  str(points.iloc[i, points.columns.get_loc(id_col)])
        df[id_col] = sid
        values.append(df)
        
    df1 = pd.concat(values, ignore_index = True)        
    if 'QA_PIXEL' in bands or 'QA60' in bands:
        final_df = process_bitmask(df1)
    else:
        final_df = df1.copy()
    
    if product == 'COPERNICUS/S2_HARMONIZED':
        
        print(f'Processing: Fetching azimuth and zenith (solar & sensor) from "{product}"\n')
        for band in bands:
            if band == 'QA60':
                continue  
            final_df[band] = final_df[band] * 0.0001
        
        collection = ee.ImageCollection('COPERNICUS/S2')
        cols = [id_col, 'latitude', 'longitude', 'SAA', 'SZA']
            
        df_meta = pd.DataFrame(columns = cols)
        for index, row in points.iterrows():
            latitude = row['lat']
            longitude = row['lon']
            id_var = row[id_col]
            point = ee.Geometry.Point(longitude, latitude)
            filteredCollection = collection.filterBounds(point).filterDate(start_date, end_date)
            image = ee.Image(filteredCollection.sort('system:time_start').first())

            SAA = image.get('MEAN_SOLAR_AZIMUTH_ANGLE')
            SZA = image.get('MEAN_SOLAR_ZENITH_ANGLE')

            df_meta = pd.concat([df_meta, pd.DataFrame({id_col: [id_var], 'latitude': [latitude], 'longitude': [longitude],
                                            'SAA': [SAA.getInfo()], 'SZA': [SZA.getInfo()]},
                                            index=[len(df_meta)])],
                    ignore_index=True)

            for band in bands:
                if band == 'QA60':
                    continue
                azimuth = image.get(f'MEAN_INCIDENCE_AZIMUTH_ANGLE_{band}')
                zenith = image.get(f'MEAN_INCIDENCE_ZENITH_ANGLE_{band}')

                df_meta.at[len(df_meta) - 1, f'VAA_{band}'] = azimuth.getInfo()
                df_meta.at[len(df_meta) - 1, f'VZA_{band}'] = zenith.getInfo()
                
        merged_df = final_df.merge(df_meta, on = id_col)
        merged_df.to_csv(opf, index = False)
        return merged_df
        
    elif product == 'LANDSAT/LC08/C02/T1_TOA':
        final_df['SAA'] = final_df['SAA'] * 0.01
        final_df['SZA'] = final_df['SZA'] * 0.01 
        final_df['VAA'] = final_df['VAA'] * 0.01
        final_df['VZA'] = final_df['VZA'] * 0.01
        
        final_df.to_csv(opf, index = False)
        return final_df
    
    else:
        final_df.to_csv(opf, index = False)
        return final_df
    
def download_aeronet_sites(year = 2021, level = 1.5, bbox = [2.0, 65.0, 40.0, 100.0], dest_folder = os.getcwd()):

    '''
    This function downloads the AERONET sites for a given year, level and bounding box. 
    The format of the bounding box is [min_lat, min_lon, max_lat, max_lon]. 

    Parameters:
        year (int): The year for which the AERONET sites are to be downloaded. Default is 2021.
        level (float): The level of data to be downloaded. Default is 1.5. Available options are 1.0, 1.5, 2.0.
        bbox (list): The bounding box for which the AERONET sites are to be downloaded. Default is [2.0, 65.0, 40.0, 100.0].
        dest_folder (str): The destination folder where the AERONET sites are to be saved. Default is the current working directory.

    Returns:
        None
    
    '''

    import os, pandas as pd

    print(f'\nProcessing: Fetching AERONET sites for year {year}')
    op_name = os.path.join(dest_folder, f'aeronet_sites_{year}.csv')
    
    if os.path.exists(op_name):
        site_list = pd.read_csv(op_name)
    else: 
        url = f'https://aeronet.gsfc.nasa.gov/Site_Lists_V3/aeronet_locations_v3_{year}_lev{int(level*10)}.txt'
        site_list = pd.read_csv(url, skiprows = 1, sep = ',')
        site_list.to_csv(op_name, index = False)
        
    site_list = site_list.iloc[:, :4].rename(columns = {'Longitude(decimal_degrees)': 'lon', 'Latitude(decimal_degrees)': 'lat'})
    print(f'Processing: Fetching AERONET sites for year {year} within {bbox}\n')
    
    min_lat, min_lon, max_lat, max_lon = bbox[0], bbox[1], bbox[2], bbox[3]
    site_subset = site_list[(site_list['lon'] > min_lon) & (site_list['lon'] < max_lon) & 
                            (site_list['lat'] > min_lat) & (site_list['lat'] < max_lat)]
    site_subset.to_csv(os.path.join(dest_folder, f'selected_aeronet_{year}.csv'), index = False)
    
    print(f'Number of AERONET sites available: {len(site_subset)}')
    return site_subset

def average_aeronet(df):

    '''
    This function averages the AERONET data for a given dataframe. Inherites the dataframe from the function `download_aeronet_data()`.
    Also excludes averaging of certain columns like `AERONET_Site`, `Date(dd:mm:yyyy)`, `Time(hh:mm:ss)`, etc.
    
    Parameters:
        df (pandas.DataFrame): The dataframe containing the AERONET data.

    Returns:
        pandas.DataFrame: The dataframe containing the averaged AERONET data.

    '''
    
    import numpy as np
    
    exclude_patterns = [
    "AERONET_Site",
    "Date(dd:mm:yyyy)",
    "Time(hh:mm:ss)",
    "Day_of_Year",
    "Day_of_Year(Fraction)",
    "Data_Quality_Level",
    "AERONET_Instrument_Number",
    "AERONET_Site_Name",
    "Site_Latitude(Degrees)",
    "Site_Longitude(Degrees)",
    "Site_Elevation(m)",
    "Solar_Zenith_Angle(Degrees)",
    "Sensor_Temperature(Degrees_C)",
    "Last_Date_Processed",
    "Number_of_Wavelengths",
    "Exact_Wavelengths_of_PW(um)_935nm",
    "Exact_Wavelengths_of_AOD(um)_681nm",
    "Exact_Wavelengths_of_AOD(um)_709nm",
    "Exact_Wavelengths_of_AOD(um)_Empty"
    ]
    
    column_list = df.columns.tolist()
    exclude_columns = [column for column in column_list if not any(pattern in column for pattern in exclude_patterns)]

    df.replace(-999.000000, np.nan, inplace=True)

    average_values = {}

    for column in df.columns:
        if column in exclude_columns:
            unique_values = df[column].unique()
            if len(unique_values) > 1:
                average_values[column] = unique_values[0]
            else:
                average_values[column] = unique_values[0] if len(unique_values) > 0 else np.nan
        elif pd.api.types.is_numeric_dtype(df[column]):
            average_values[column] = df[column].mean()
        else:
            average_values[column] = df[column].unique()[0]

    new_df = pd.DataFrame([average_values], columns=df.columns)
    return new_df

def download_aeronet_aod(site, temporal_scale=60, level=1.5, id_col='Site_Name', verbose = True):

    '''
    This function downloads the AERONET AOD data for a given site, temporal scale and level.

    Parameters:
        site (pandas.DataFrame): The dataframe containing the AERONET site information.
        temporal_scale (int): The temporal scale for which the AERONET data is to be downloaded (in minutes). Default is 60.
        level (float): The level of data to be downloaded. Default is 1.5. Available options are 1.0, 1.5, 2.0.
        id_col (str): The column name containing the site name. Default is `Site_Name`.
        verbose (bool): Whether to print the progress of the download. Default is True.

    Returns:
        pandas.DataFrame: The dataframe containing the AERONET AOD data.

    '''
    
    from datetime import datetime, timedelta
    
    print(f'\nProcessing: Fetching and averaging AERONET Level{level} values at \u00B1{temporal_scale}min\n')
    site['date'] = pd.to_datetime(site['date'])
    site['hour'] = pd.to_datetime(site['date'].dt.strftime('%H:%M:%S.%f'), format='%H:%M:%S.%f').dt.time

    for i in range(len(site)):
        site_name = site.loc[i, id_col]
        date_obj = site.loc[i, 'date']
        year1 = year2 = int(date_obj.year)
        month1 = month2 = int(date_obj.month)
        day1 = day2 = int(date_obj.day)
        time_obj = site.loc[i, 'hour']

        if time_obj.minute >= (60 - temporal_scale):
            hour1 = int(time_obj.hour)
            hour2 = int((datetime.combine(date_obj, time_obj).replace(microsecond=0) + timedelta(minutes=temporal_scale)).hour)
        else:
            hour1 = int(time_obj.hour)
            hour2 = int((datetime.combine(date_obj, time_obj).replace(microsecond=0) + timedelta(minutes=temporal_scale)).hour)
         
        url = f'https://aeronet.gsfc.nasa.gov/cgi-bin/print_web_data_v3?site={site_name}&year={year1}&month={month1}&day={day1}&hour={hour1}&year2={year2}&month2={month2}&day2={day2}&hour2={hour2}&AOD{int(level*10)}=1&AVG=10'
        
        try:
            df = pd.read_csv(url, skiprows=7, sep=',').dropna()
            if verbose:
                print("File downloaded successfully from", url)
            if len(df) > 1:
                new_df = average_aeronet(df)
                site.loc[i, new_df.columns] = new_df.iloc[0]  # Assign new_df to the corresponding row in site
            else:
                if verbose: 
                    print("No data in file from", url)
        except (IOError, pd.errors.ParserError, pd.errors.EmptyDataError) as e:
            print("Error occurred while downloading or parsing the file from", url)
            print("Skipping to the next site.")
            continue
        
    return site

def extract_aeronet_and_reflectance(gee_product_id = 'LANDSAT/LC08/C02/T1_TOA', 
                                    start_date = '2021-04-01', 
                                    end_date = '2021-04-30',
                                    spectral_bands = ['B1', 'B2', 'B3', 'B4'], 
                                    scale = 30, 
                                    aeronet_level = 1.5, 
                                    temporal_average = 30,
                                    bbox = [2.0, 65.0, 40.0, 100.0], 
                                    dest_folder = os.getcwd()):
    
    '''
    This is the caller function for the AERONET and reflectance extraction process. 
    It calls the functions `download_aeronet_aod` and `extract_reflectance` to download the AERONET AOD data and 
    extract the reflectance data from GEE, respectively.

    Parameters:
        gee_product_id (str): The GEE product ID. Default is `LANDSAT/LC08/C02/T1_TOA`.
        start_date (str): The start date in the format `YYYY-MM-DD`. Default is `2021-04-01`.
        end_date (str): The end date in the format `YYYY-MM-DD`. Default is `2021-04-30`.
        spectral_bands (list): The list of spectral bands to be extracted. Default is `['B1', 'B2', 'B3', 'B4']`.
        scale (int): The scale at which the reflectance data is to be extracted. Default is 30.
        aeronet_level (float): The level of AERONET data to be downloaded. Default is 1.5. Available options are 1.0, 1.5, 2.0.
        temporal_average (int): The temporal scale for which the AERONET data is to be downloaded (in minutes). Default is 60.
        bbox (list): The bounding box for which the reflectance data is to be extracted. Default is `[2.0, 65.0, 40.0, 100.0]`.
        dest_folder (str): The destination folder where the AERONET and reflectance data is to be saved. Default is the current working directory.

    Returns:
        pandas.DataFrame: The dataframe containing the AERONET AOD and reflectance data.
    
    '''
    
    opf = os.path.join(dest_folder, 'AERONET_module')
    os.makedirs(opf, exist_ok = True)

    start_year = int(start_date.split('-')[0])
    end_year = int(end_date.split('-')[0])
    start_month = int(start_date.split('-')[1])
    start_day = int(start_date.split('-')[-1])

    year_list = list(range(start_year, end_year + 1))
    for i, year in enumerate(year_list):
        if i == len(year_list) - 1:
            end_month = int(end_date.split('-')[1])
            end_day = int(end_date.split('-')[-1])
        else:
            end_month = 12
            end_day = 31

        if year == start_year:
            start_month = int(start_date.split('-')[1])
            start_day = int(start_date.split('-')[-1])
        else:
            start_month = 1
            start_day = 1

        start_date_str = f'{year}-{start_month:02d}-{start_day:02d}'
        end_date_str = f'{year}-{end_month:02d}-{end_day:02d}'

        site_subset = download_aeronet_sites(year = year, level=aeronet_level, bbox = bbox, dest_folder = opf)
        site = gee_point_extract(site_subset, product = gee_product_id, start_date = start_date_str, 
                                 end_date = end_date_str, id_col = 'Site_Name',
                                 bands = spectral_bands, scale = scale, dest_folder = opf)
        aeronet_df = download_aeronet_aod(site, temporal_scale = temporal_average, level = aeronet_level, 
                                          id_col = 'Site_Name', verbose = False)
        
        aeronet_df.to_csv(os.path.join(opf, f'aeronet_{aeronet_level}_{year}_{gee_product_id.split("/")[-1]}.csv'))
        
    return aeronet_df
