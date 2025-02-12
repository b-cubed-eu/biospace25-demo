<a href="https://colab.research.google.com/github/AgentschapPlantentuinMeise/DEMO_BioSpace25/blob/main/DEMO_BioSpace25.ipynb" target="_parent"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/></a>

# B-Cubed DEMO BioSpace25
## Building data cubes



### Install missing packages


```python
%pip install pygbif
%pip install dask
```
    

### Only execute the following block when using the TPU kernel


```python
%pip install geopandas
%pip install pydrive
%pip install ee
%pip install eerepr
%pip install geemap
```

### Loading packages


```python
from pygbif import occurrences as occ
import pandas as pd
import geopandas as gpd
from pyproj import Proj, Transformer
from shapely.geometry import mapping
from shapely.geometry import Polygon
import matplotlib.pyplot as plt

from pydrive.auth import GoogleAuth
from pydrive.drive import GoogleDrive
from google.colab import auth
from google.colab import drive
from oauth2client.client import GoogleCredentials
import io
from io import StringIO
import zipfile
import math
```

    WARNING:root:pydrive is deprecated and no longer maintained. We recommend that you migrate your projects to pydrive2, the maintained fork of pydrive


### Loading Earth Engine


```python
import ee
import eerepr
import geemap

ee.Authenticate(force=True)
ee.Initialize(project='nithecs-436810')

LANDSAT_ID = "LANDSAT/LC08/C02/T1_L2"
BOUNDARIES_ID = 'FAO/GAUL/2015/level1'
WDPA_ID = 'WCMC/WDPA/current/polygons'
SENTINEL_ID = "COPERNICUS/S2_SR_HARMONIZED"


dataset = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2').filterDate('2021-05-01', '2021-06-01')
sa = ee.FeatureCollection(BOUNDARIES_ID).filter(
    'ADM0_NAME == "South Africa"')

dataset_eo = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED').filterDate('2020-01-01', '2020-01-30')

protected_areas = ee.FeatureCollection(WDPA_ID)


sa_landsat = dataset.filterBounds(sa)
sa_sentinel = dataset_eo.filterBounds(sa)

```

### Example of the GBIF API through pygbif


```python
from pygbif import occurrences
data = occurrences.search(speciesKey=5229490, limit=10)

print(data['results'])
```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>


## GBIF data Cubes

### Generating the Cube

#### Exemplar JSON query for generating a data cube


```
# DON'T EXECUTE THIS CELL. FOR DEMO PUPOSE ONLY #
{
  "sendNotification": true,
  "notificationAddresses": [
    "maarten.trekels@plantentuinmeise.be"
  ],
  "format": "SQL_TSV_ZIP",
  "sql": "SELECT  PRINTF('%04d-%02d', \"year\", \"month\") AS yearMonth,
   GBIF_EEARGCode(10000, decimalLatitude,  decimalLongitude,  COALESCE(coordinateUncertaintyInMeters, 1000) ) AS eeaCellCode,
   speciesKey,
   species,
   establishmentMeans,
   degreeOfEstablishment,
   pathway,
   COUNT(*) AS occurrences,
   COUNT(DISTINCT recordedBy) AS distinctObservers
   FROM  occurrence
   WHERE occurrenceStatus = 'PRESENT'
   AND countryCode = 'BE'
   AND hasCoordinate = TRUE
   AND NOT ARRAY_CONTAINS(issue, 'ZERO_COORDINATE')
   AND NOT ARRAY_CONTAINS(issue, 'COORDINATE_OUT_OF_RANGE')
   AND NOT ARRAY_CONTAINS(issue, 'COORDINATE_INVALID')
   AND NOT ARRAY_CONTAINS(issue, 'COUNTRY_COORDINATE_MISMATCH')
   AND \"month\" IS NOT NULL
   GROUP BY yearMonth,
   eeaCellCode,
   speciesKey,
   species,
   establishmentMeans,
   degreeOfEstablishment,
   pathway
   ORDER BY  yearMonth DESC,
   eeaCellCode ASC,
   speciesKey ASC"
}


```



## Loading the Data cube in pandas



You can download a pre generated data cube from GitHub or any other online resource

#### Download from Google Drive


```python
drive.mount('/content/drive')

def convert_to_int(x):
       try:
           return int(x)
       except ValueError:
           return pd.NA  # or np.nan if you prefer NumPy NaNs

data = pd.read_csv('/content/drive/Shareddrives/BioSpace25/supporting_data/Cube_ZA_QDGC_l3.csv', sep='\t', converters={'familykey': convert_to_int, 'specieskey': convert_to_int})

data['familykey'] = pd.to_numeric(data['familykey'], errors='coerce').astype('Int64')
data['specieskey'] = pd.to_numeric(data['specieskey'], errors='coerce').astype('Int64')

```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>



    Mounted at /content/drive



```python
print(data)
```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>



             yearmonth    qdgccode  familykey        family  specieskey  \
    0          2024-09  E016S28ADD       2406  Crassulaceae     7716880   
    1          2024-09  E016S28BCC       2406  Crassulaceae     7716880   
    2          2024-09  E016S28BDD       4676   Geraniaceae     3826148   
    3          2024-09  E016S28BDD       2406  Crassulaceae     7334236   
    4          2024-09  E016S28CBB       6752     Aizoaceae     8003531   
    ...            ...         ...        ...           ...         ...   
    17876161   1772-04  E025S28CCC       7689   Orchidaceae     2783834   
    17876162   1694-10  E027S32DBB       2430    Onagraceae     3188875   
    17876163   1678-02  E028S25DCB       7359     Araneidae        <NA>   
    17876164   1645-12  E030S30BBC       7016  Notodontidae     1824935   
    17876165   1608-10  E023S33ABC       4334        Apidae     5040145   
    
                                 species  occurrences  familycount  
    0                  Crassula sladenii            1      44434.0  
    1                  Crassula sladenii            1      44434.0  
    2         Pelargonium klinghardtense            1      46821.0  
    3                   Crassula elegans            1      44434.0  
    4               Conophytum saxetanum            1      57598.0  
    ...                              ...          ...          ...  
    17876161               Disa uniflora            1      54636.0  
    17876162             Oenothera rosea            1       2688.0  
    17876163                         NaN            1       9699.0  
    17876164            Antheua tricolor            1        874.0  
    17876165          Anthophora praecox            1      31913.0  
    
    [17876166 rows x 8 columns]


## Getting a Geopackage file from the Grid that you use


```python
# Load QDGC code

input_file = "/content/drive/Shareddrives/BioSpace25/supporting_data/qdgc_south_africa.gpkg"

qdgc_ref = gpd.read_file(input_file, layer='tbl_qdgc_03')
```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>




```python
print(qdgc_ref)
```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>



                 qdgc  level_qdgc  cellsize_degrees  lon_center  lat_center  \
    0      E016S46CDD           3             0.125     16.4375    -46.9375   
    1      E016S46CDB           3             0.125     16.4375    -46.8125   
    2      E016S46CBD           3             0.125     16.4375    -46.6875   
    3      E016S46CBB           3             0.125     16.4375    -46.5625   
    4      E016S46ADD           3             0.125     16.4375    -46.4375   
    ...           ...         ...               ...         ...         ...   
    34422  E037S22DBD           3             0.125     37.9375    -22.6875   
    34423  E037S22DBB           3             0.125     37.9375    -22.5625   
    34424  E037S22BDD           3             0.125     37.9375    -22.4375   
    34425  E037S22BDB           3             0.125     37.9375    -22.3125   
    34426  E037S22BBD           3             0.125     37.9375    -22.1875   
    
             area_km2                                           geometry  
    0      132.265148  MULTIPOLYGON (((16.375 -47, 16.375 -46.875, 16...  
    1      132.569719  MULTIPOLYGON (((16.375 -46.875, 16.375 -46.75,...  
    2      132.873641  MULTIPOLYGON (((16.375 -46.75, 16.375 -46.625,...  
    3      133.176911  MULTIPOLYGON (((16.375 -46.625, 16.375 -46.5, ...  
    4      133.479529  MULTIPOLYGON (((16.375 -46.5, 16.375 -46.375, ...  
    ...           ...                                                ...  
    34422  177.801185  MULTIPOLYGON (((37.875 -22.75, 37.875 -22.625,...  
    34423  177.959230  MULTIPOLYGON (((37.875 -22.625, 37.875 -22.5, ...  
    34424  178.116437  MULTIPOLYGON (((37.875 -22.5, 37.875 -22.375, ...  
    34425  178.272807  MULTIPOLYGON (((37.875 -22.375, 37.875 -22.25,...  
    34426  178.428337  MULTIPOLYGON (((37.875 -22.25, 37.875 -22.125,...  
    
    [34427 rows x 7 columns]


## Merging the Data cube with the grid


```python
test_merge = pd.merge(data, qdgc_ref, left_on='qdgccode', right_on='qdgc')
```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>




```python
# Convert to GeoDataFrame

gdf = gpd.GeoDataFrame(test_merge, geometry='geometry')

gdf = gdf.set_crs(epsg=4326, inplace=False)
```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>



## Filtering data (e.g. on species)
Check for a single species (Loxodonta africana: https://www.gbif.org/species/2435350)


```python
filtered_gdf = gdf[gdf['specieskey'].eq(2435350)]

print(filtered_gdf)

```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>



             yearmonth    qdgccode  familykey        family  specieskey  \
    6401       2024-09  E025S33BAD       9427  Elephantidae     2435350   
    6402       2024-09  E025S33BBC       9427  Elephantidae     2435350   
    6403       2024-09  E025S33BCA       9427  Elephantidae     2435350   
    6405       2024-09  E025S33BCB       9427  Elephantidae     2435350   
    6406       2024-09  E025S33BCC       9427  Elephantidae     2435350   
    ...            ...         ...        ...           ...         ...   
    14708306   1990-01  E025S33BCB       9427  Elephantidae     2435350   
    14708309   1990-01  E025S33BDA       9427  Elephantidae     2435350   
    16222304   1987-01  E031S24BBC       9427  Elephantidae     2435350   
    16427602   1985-09  E031S23ABA       9427  Elephantidae     2435350   
    16427611   1985-09  E031S23ADD       9427  Elephantidae     2435350   
    
                         species  occurrences  familycount        qdgc  \
    6401      Loxodonta africana           14       5301.0  E025S33BAD   
    6402      Loxodonta africana            3       5301.0  E025S33BBC   
    6403      Loxodonta africana            1       5301.0  E025S33BCA   
    6405      Loxodonta africana           21       5301.0  E025S33BCB   
    6406      Loxodonta africana            2       5301.0  E025S33BCC   
    ...                      ...          ...          ...         ...   
    14708306  Loxodonta africana            2       5301.0  E025S33BCB   
    14708309  Loxodonta africana            2       5301.0  E025S33BDA   
    16222304  Loxodonta africana            1       5301.0  E031S24BBC   
    16427602  Loxodonta africana            1       5301.0  E031S23ABA   
    16427611  Loxodonta africana            1       5301.0  E031S23ADD   
    
              level_qdgc  cellsize_degrees  lon_center  lat_center    area_km2  \
    6401               3             0.125     25.6875    -33.1875  161.604840   
    6402               3             0.125     25.8125    -33.1875  161.604840   
    6403               3             0.125     25.5625    -33.3125  161.378183   
    6405               3             0.125     25.6875    -33.3125  161.378183   
    6406               3             0.125     25.5625    -33.4375  161.150755   
    ...              ...               ...         ...         ...         ...   
    14708306           3             0.125     25.6875    -33.3125  161.378183   
    14708309           3             0.125     25.8125    -33.3125  161.378183   
    16222304           3             0.125     31.8125    -24.1875  175.839539   
    16427602           3             0.125     31.3125    -23.0625  177.322031   
    16427611           3             0.125     31.4375    -23.4375  176.835361   
    
                                                       geometry  
    6401      MULTIPOLYGON (((25.625 -33.25, 25.625 -33.125,...  
    6402      MULTIPOLYGON (((25.75 -33.25, 25.75 -33.125, 2...  
    6403      MULTIPOLYGON (((25.5 -33.375, 25.5 -33.25, 25....  
    6405      MULTIPOLYGON (((25.625 -33.375, 25.625 -33.25,...  
    6406      MULTIPOLYGON (((25.5 -33.5, 25.5 -33.375, 25.6...  
    ...                                                     ...  
    14708306  MULTIPOLYGON (((25.625 -33.375, 25.625 -33.25,...  
    14708309  MULTIPOLYGON (((25.75 -33.375, 25.75 -33.25, 2...  
    16222304  MULTIPOLYGON (((31.75 -24.25, 31.75 -24.125, 3...  
    16427602  MULTIPOLYGON (((31.25 -23.125, 31.25 -23, 31.3...  
    16427611  MULTIPOLYGON (((31.375 -23.5, 31.375 -23.375, ...  
    
    [3851 rows x 15 columns]


## Apply the function to create a list of features in GEE


```python
filtered_gdf = filtered_gdf.set_crs(epsg=4326, inplace=False)

data_raw = geemap.geopandas_to_ee(filtered_gdf)
```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>



## Visualization of the data cubes on a map with different layers


```python
from google.colab import output
output.enable_custom_widget_manager()
```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>




```python
Map = geemap.Map(layout={"height": "400px", "width": "800px"})


# Add the original data layer in blue
Map.addLayer(data_raw, {"color": "blue"}, "Original data")

visualization = {
    'min': 0.0,
    'max': 0.3,
    'bands': ['B2', 'B3', 'B4'],
}

Map.addLayer(sa_sentinel, visualization, 'RGB')

Map.addLayer(protected_areas)


# Set the center of the map to the coordinates
Map.setCenter(-28.50, 29.41)
Map
```



<style>
    .geemap-dark {
        --jp-widgets-color: white;
        --jp-widgets-label-color: white;
        --jp-ui-font-color1: white;
        --jp-layout-color2: #454545;
        background-color: #383838;
    }

    .geemap-dark .jupyter-button {
        --jp-layout-color3: #383838;
    }

    .geemap-colab {
        background-color: var(--colab-primary-surface-color, white);
    }

    .geemap-colab .jupyter-button {
        --jp-layout-color3: var(--colab-primary-surface-color, white);
    }
</style>




    Map(center=[29.41, -28.5], controls=(WidgetControl(options=['position', 'transparent_bg'], widget=SearchDataGU…


## Export the data to GeoParquet


```python
#gdf.to_parquet('/content/drive/Shareddrives/BioSpace25/supporting_data/data_ZA.parquet')
```

# EBV data cubes in NetCDF format


```python
%pip install netCDF4
```

    Requirement already satisfied: netCDF4 in /usr/local/lib/python3.11/dist-packages (1.7.2)
    Requirement already satisfied: cftime in /usr/local/lib/python3.11/dist-packages (from netCDF4) (1.6.4.post1)
    Requirement already satisfied: certifi in /usr/local/lib/python3.11/dist-packages (from netCDF4) (2024.12.14)
    Requirement already satisfied: numpy in /usr/local/lib/python3.11/dist-packages (from netCDF4) (1.26.4)



```python
%pip install rioxarray
%pip install cartopy
%pip install basemap
```

    Requirement already satisfied: rioxarray in /usr/local/lib/python3.11/dist-packages (0.18.2)
    Requirement already satisfied: packaging in /usr/local/lib/python3.11/dist-packages (from rioxarray) (23.2)
    Requirement already satisfied: rasterio>=1.3.7 in /usr/local/lib/python3.11/dist-packages (from rioxarray) (1.4.3)
    Requirement already satisfied: xarray>=2024.7.0 in /usr/local/lib/python3.11/dist-packages (from rioxarray) (2025.1.1)
    Requirement already satisfied: pyproj>=3.3 in /usr/local/lib/python3.11/dist-packages (from rioxarray) (3.6.1)
    Requirement already satisfied: numpy>=1.23 in /usr/local/lib/python3.11/dist-packages (from rioxarray) (1.26.4)
    Requirement already satisfied: certifi in /usr/local/lib/python3.11/dist-packages (from pyproj>=3.3->rioxarray) (2024.12.14)
    Requirement already satisfied: affine in /usr/local/lib/python3.11/dist-packages (from rasterio>=1.3.7->rioxarray) (2.4.0)
    Requirement already satisfied: attrs in /usr/local/lib/python3.11/dist-packages (from rasterio>=1.3.7->rioxarray) (24.3.0)
    Requirement already satisfied: click>=4.0 in /usr/local/lib/python3.11/dist-packages (from rasterio>=1.3.7->rioxarray) (8.1.8)
    Requirement already satisfied: cligj>=0.5 in /usr/local/lib/python3.11/dist-packages (from rasterio>=1.3.7->rioxarray) (0.7.2)
    Requirement already satisfied: click-plugins in /usr/local/lib/python3.11/dist-packages (from rasterio>=1.3.7->rioxarray) (1.1.1)
    Requirement already satisfied: pyparsing in /usr/local/lib/python3.11/dist-packages (from rasterio>=1.3.7->rioxarray) (3.2.1)
    Requirement already satisfied: pandas>=2.1 in /usr/local/lib/python3.11/dist-packages (from xarray>=2024.7.0->rioxarray) (2.2.2)
    Requirement already satisfied: python-dateutil>=2.8.2 in /usr/local/lib/python3.11/dist-packages (from pandas>=2.1->xarray>=2024.7.0->rioxarray) (2.8.2)
    Requirement already satisfied: pytz>=2020.1 in /usr/local/lib/python3.11/dist-packages (from pandas>=2.1->xarray>=2024.7.0->rioxarray) (2024.2)
    Requirement already satisfied: tzdata>=2022.7 in /usr/local/lib/python3.11/dist-packages (from pandas>=2.1->xarray>=2024.7.0->rioxarray) (2025.1)
    Requirement already satisfied: six>=1.5 in /usr/local/lib/python3.11/dist-packages (from python-dateutil>=2.8.2->pandas>=2.1->xarray>=2024.7.0->rioxarray) (1.17.0)
    Requirement already satisfied: cartopy in /usr/local/lib/python3.11/dist-packages (0.24.1)
    Requirement already satisfied: numpy>=1.23 in /usr/local/lib/python3.11/dist-packages (from cartopy) (1.26.4)
    Requirement already satisfied: matplotlib>=3.6 in /usr/local/lib/python3.11/dist-packages (from cartopy) (3.8.4)
    Requirement already satisfied: shapely>=1.8 in /usr/local/lib/python3.11/dist-packages (from cartopy) (2.0.6)
    Requirement already satisfied: packaging>=21 in /usr/local/lib/python3.11/dist-packages (from cartopy) (23.2)
    Requirement already satisfied: pyshp>=2.3 in /usr/local/lib/python3.11/dist-packages (from cartopy) (2.3.1)
    Requirement already satisfied: pyproj>=3.3.1 in /usr/local/lib/python3.11/dist-packages (from cartopy) (3.6.1)
    Requirement already satisfied: contourpy>=1.0.1 in /usr/local/lib/python3.11/dist-packages (from matplotlib>=3.6->cartopy) (1.3.1)
    Requirement already satisfied: cycler>=0.10 in /usr/local/lib/python3.11/dist-packages (from matplotlib>=3.6->cartopy) (0.12.1)
    Requirement already satisfied: fonttools>=4.22.0 in /usr/local/lib/python3.11/dist-packages (from matplotlib>=3.6->cartopy) (4.55.6)
    Requirement already satisfied: kiwisolver>=1.3.1 in /usr/local/lib/python3.11/dist-packages (from matplotlib>=3.6->cartopy) (1.4.8)
    Requirement already satisfied: pillow>=8 in /usr/local/lib/python3.11/dist-packages (from matplotlib>=3.6->cartopy) (11.1.0)
    Requirement already satisfied: pyparsing>=2.3.1 in /usr/local/lib/python3.11/dist-packages (from matplotlib>=3.6->cartopy) (3.2.1)
    Requirement already satisfied: python-dateutil>=2.7 in /usr/local/lib/python3.11/dist-packages (from matplotlib>=3.6->cartopy) (2.8.2)
    Requirement already satisfied: certifi in /usr/local/lib/python3.11/dist-packages (from pyproj>=3.3.1->cartopy) (2024.12.14)
    Requirement already satisfied: six>=1.5 in /usr/local/lib/python3.11/dist-packages (from python-dateutil>=2.7->matplotlib>=3.6->cartopy) (1.17.0)
    Requirement already satisfied: basemap in /usr/local/lib/python3.11/dist-packages (1.4.1)
    Requirement already satisfied: basemap-data<1.4,>=1.3.2 in /usr/local/lib/python3.11/dist-packages (from basemap) (1.3.2)
    Requirement already satisfied: pyshp<2.4,>=1.2 in /usr/local/lib/python3.11/dist-packages (from basemap) (2.3.1)
    Requirement already satisfied: matplotlib<3.9,>=1.5 in /usr/local/lib/python3.11/dist-packages (from basemap) (3.8.4)
    Requirement already satisfied: pyproj<3.7.0,>=1.9.3 in /usr/local/lib/python3.11/dist-packages (from basemap) (3.6.1)
    Requirement already satisfied: packaging<24.0,>=16.0 in /usr/local/lib/python3.11/dist-packages (from basemap) (23.2)
    Requirement already satisfied: numpy<1.27,>=1.21 in /usr/local/lib/python3.11/dist-packages (from basemap) (1.26.4)
    Requirement already satisfied: contourpy>=1.0.1 in /usr/local/lib/python3.11/dist-packages (from matplotlib<3.9,>=1.5->basemap) (1.3.1)
    Requirement already satisfied: cycler>=0.10 in /usr/local/lib/python3.11/dist-packages (from matplotlib<3.9,>=1.5->basemap) (0.12.1)
    Requirement already satisfied: fonttools>=4.22.0 in /usr/local/lib/python3.11/dist-packages (from matplotlib<3.9,>=1.5->basemap) (4.55.6)
    Requirement already satisfied: kiwisolver>=1.3.1 in /usr/local/lib/python3.11/dist-packages (from matplotlib<3.9,>=1.5->basemap) (1.4.8)
    Requirement already satisfied: pillow>=8 in /usr/local/lib/python3.11/dist-packages (from matplotlib<3.9,>=1.5->basemap) (11.1.0)
    Requirement already satisfied: pyparsing>=2.3.1 in /usr/local/lib/python3.11/dist-packages (from matplotlib<3.9,>=1.5->basemap) (3.2.1)
    Requirement already satisfied: python-dateutil>=2.7 in /usr/local/lib/python3.11/dist-packages (from matplotlib<3.9,>=1.5->basemap) (2.8.2)
    Requirement already satisfied: certifi in /usr/local/lib/python3.11/dist-packages (from pyproj<3.7.0,>=1.9.3->basemap) (2024.12.14)
    Requirement already satisfied: six>=1.5 in /usr/local/lib/python3.11/dist-packages (from python-dateutil>=2.7->matplotlib<3.9,>=1.5->basemap) (1.17.0)



```python
import netCDF4 as nc
import xarray as xr


birds_file = xr.open_dataset('/content/drive/Shareddrives/BioSpace25/supporting_data/viti_spepop_id77_20240206_v1.nc')

print(birds_file)
```

    <xarray.Dataset> Size: 26kB
    Dimensions:  (lon: 559, lat: 437, time: 1, entity: 486)
    Coordinates:
      * lon      (lon) float64 4kB 9.45e+05 9.55e+05 ... 6.515e+06 6.525e+06
      * lat      (lat) float64 3kB 5.305e+06 5.295e+06 ... 9.55e+05 9.45e+05
      * time     (time) datetime64[ns] 8B 2018-01-01
      * entity   (entity) |S37 18kB b'Gavia stellata                       ' ... ...
    Data variables:
        crs      |S1 1B ...
    Attributes: (12/38)
        Conventions:                CF-1.8, ACDD-1.3, EBV-1.0
        naming_authority:           The German Centre for Integrative Biodiversit...
        history:                    EBV netCDF created using ebvcube, 2024-02-06
        ebv_vocabulary:             https://portal.geobon.org/api/v1/ebv
        ebv_cube_dimensions:        lon, lat, time, entity
        keywords:                   ebv_class: Species populations, ebv_name: Spe...
        ...                         ...
        geospatial_lat_units:       meter
        time_coverage_start:        2013-01-01
        time_coverage_end:          2018-12-31
        time_coverage_resolution:   P0000-00-00
        date_issued:                2024-02-12
        comment:                    List of species: https://cdr.eionet.europa.eu...



```python
def print_netcdf_structure(nc_file_path):
  """Prints the structure (groups, variables, and their paths) of a NetCDF file.

  Args:
    nc_file_path: Path to the NetCDF file.
  """
  def print_group_structure(group, path=""):
    """Recursively prints the structure of a group within the NetCDF file."""
    for var_name in group.variables:
      print(f"{path}/{var_name}")  # Print variable path
    for group_name in group.groups:
      subgroup = group.groups[group_name]
      print_group_structure(subgroup, f"{path}/{group_name}")  # Recursively explore subgroups

  with nc.Dataset(nc_file_path, 'r') as nc_file:
    print_group_structure(nc_file)  # Start with the root group

# Example usage:
nc_file_path = '/content/drive/Shareddrives/BioSpace25/supporting_data/viti_spepop_id77_20240206_v1.nc'
print_netcdf_structure(nc_file_path)
```

    /lon
    /lat
    /time
    /crs
    /entity
    /metric_1/ebv_cube



```python
print(birds_file.variables)
```

    Frozen({'lon': <xarray.IndexVariable 'lon' (lon: 559)> Size: 4kB
    array([ 945000.,  955000.,  965000., ..., 6505000., 6515000., 6525000.])
    Attributes:
        long_name:      lon
        standard_name:  projection_x_coordinate
        axis:           X
        units:          meter, 'lat': <xarray.IndexVariable 'lat' (lat: 437)> Size: 3kB
    array([5305000., 5295000., 5285000., ...,  965000.,  955000.,  945000.])
    Attributes:
        long_name:      lat
        standard_name:  projection_y_coordinate
        axis:           Y
        units:          meter, 'time': <xarray.IndexVariable 'time' (time: 1)> Size: 8B
    array(['2018-01-01T00:00:00.000000000'], dtype='datetime64[ns]')
    Attributes:
        long_name:  time
        axis:       T, 'crs': <xarray.Variable ()> Size: 1B
    [1 values with dtype=|S1]
    Attributes:
        spatial_ref:                     PROJCRS["ETRS89-extended / LAEA Europe",...
        GeoTransform:                    940000 10000 0.0 5310000 0.0 -10000
        grid_mapping_name:               lambert_azimuthal_equal_area
        latitude_of_projection_origin:   52.0
        longitude_of_projection_origin:  10.0
        false_easting:                   4321000.0
        false_northing:                  3210000.0
        semi_major_axis:                 6378137.0
        inverse_flattening:              298.257223563
        longitude_of_prime_meridian:     0.0
        long_name:                       CRS definition, 'entity': <xarray.IndexVariable 'entity' (entity: 486)> Size: 18kB
    array([b'Gavia stellata                       ',
           b'Gavia arctica                        ',
           b'Tachybaptus ruficollis               ', ...,
           b'Accipiter gentilis all others        ',
           b'Melanitta nigra s. str.              ',
           b'Sylvia subalpina                     '], dtype='|S37')
    Attributes:
        units:                           1
        ebv_entity_type:                 Species
        ebv_entity_scope:                Bird species listed under the Art. 12 of...
        ebv_entity_classification_name:  Species names as accepted by the Birds D...
        ebv_entity_classification_url:   https://cdr.eionet.europa.eu/help/birds_...
        long_name:                       entity})



```python
time = birds_file.variables['time']
print(time)

print(birds_file['entity'])
```

    <xarray.IndexVariable 'time' (time: 1)> Size: 8B
    array(['2018-01-01T00:00:00.000000000'], dtype='datetime64[ns]')
    Attributes:
        long_name:  time
        axis:       T
    <xarray.DataArray 'entity' (entity: 486)> Size: 18kB
    array([b'Gavia stellata                       ',
           b'Gavia arctica                        ',
           b'Tachybaptus ruficollis               ', ...,
           b'Accipiter gentilis all others        ',
           b'Melanitta nigra s. str.              ',
           b'Sylvia subalpina                     '], dtype='|S37')
    Coordinates:
      * entity   (entity) |S37 18kB b'Gavia stellata                       ' ... ...
    Attributes:
        units:                           1
        ebv_entity_type:                 Species
        ebv_entity_scope:                Bird species listed under the Art. 12 of...
        ebv_entity_classification_name:  Species names as accepted by the Birds D...
        ebv_entity_classification_url:   https://cdr.eionet.europa.eu/help/birds_...
        long_name:                       entity



```python
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import xarray as xr
import numpy as np
from pyproj import Transformer
from matplotlib.colors import ListedColormap, BoundaryNorm

# --- Load Dataset Efficiently ---
birds_file = xr.open_dataset(
    '/content/drive/Shareddrives/BioSpace25/supporting_data/viti_spepop_id77_20240206_v1.nc',
    group="metric_1",
    chunks={'entity': 1, 'time': 1}  # Load only one entity/time slice at a time
)

# --- Select Target Species and Time ---
species_index = 150  # Change to the species index you need
time_index = 0  # Change to the desired time index

# Extract only the required slice
species_data_subset = birds_file['ebv_cube'].sel(entity=species_index, time=time_index).compute()

# Convert to a 2D array
species_distribution_2d = np.squeeze(species_data_subset)

# --- Load Longitude and Latitude (Only Once) ---
with xr.open_dataset('/content/drive/Shareddrives/BioSpace25/supporting_data/viti_spepop_id77_20240206_v1.nc') as ds:
    lon = ds['lon'].values  # 1D array (size: 559)
    lat = ds['lat'].values  # 1D array (size: 437)

# --- Create Meshgrid Efficiently ---
lon_grid, lat_grid = np.meshgrid(lon, lat)

# --- Efficient Coordinate Transformation ---
transformer = Transformer.from_crs("epsg:3035", "epsg:4326", always_xy=True)

# Transform the entire 2D meshgrid
lon_deg, lat_deg = transformer.transform(lon_grid, lat_grid)

# --- Create the Plot ---
fig = plt.figure(figsize=(10, 6))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree())

# Add map features
ax.coastlines()
ax.add_feature(cfeature.LAND, edgecolor='black', facecolor='lightgray')

# --- 🔹 Fix: Use Discrete Colormap Without Color Bar ---
unique_values = np.unique(species_distribution_2d)

# If only one value, choose a single solid color
if len(unique_values) == 1:
    cmap = ListedColormap(["red"])  # Single-color for uniform data
    norm = None
else:
    cmap = ListedColormap(["white", "blue"])  # Adjust colors as needed
    norm = BoundaryNorm([0, 0.5, 1], cmap.N)

# --- 🔹 Fix: Use `shading="nearest"` to Ensure Correct Grid Alignment ---
cs = ax.pcolormesh(
    lon_deg, lat_deg, species_distribution_2d,
    transform=ccrs.PlateCarree(),
    cmap=cmap,
    norm=norm,
    shading='nearest'  # Prevents visual distortion
)

# --- 🔹 Completely Remove Color Bar ---
# No `fig.colorbar(cs)`, so no scale bar will be shown

# --- Get the Species Name Efficiently ---
with xr.open_dataset('/content/drive/Shareddrives/BioSpace25/supporting_data/viti_spepop_id77_20240206_v1.nc') as ds:
    species_name = ds['entity'].values[species_index].decode('utf-8').strip()

# --- Final Plot Customization ---
ax.set_title(f"Species {species_name} Distribution (10x10 km grid) at time {birds_file['time'].values[time_index].item()}")
ax.set_xlabel("Longitude")
ax.set_ylabel("Latitude")

# --- 🔹 Fix: Add Grid Lines to Show Exact 10x10 km Cells ---
gridlines = ax.gridlines(draw_labels=True, linestyle="--", linewidth=0.5, color="black", alpha=0.5)

# Show the plot
plt.show()

```

    /usr/local/lib/python3.11/dist-packages/cartopy/io/__init__.py:241: DownloadWarning: Downloading: https://naturalearth.s3.amazonaws.com/50m_physical/ne_50m_land.zip
      warnings.warn(f'Downloading: {url}', DownloadWarning)
    /usr/local/lib/python3.11/dist-packages/cartopy/io/__init__.py:241: DownloadWarning: Downloading: https://naturalearth.s3.amazonaws.com/50m_physical/ne_50m_coastline.zip
      warnings.warn(f'Downloading: {url}', DownloadWarning)



    
![png](DEMO_BioSpace25_files/DEMO_BioSpace25_43_1.png)
    



```python
from google.colab import output
output.disable_custom_widget_manager()
```
