# GeoLifeCLEF 2020
Automatically predicting the list of species that are the most likely to be observed at a given location is useful for 
many scenarios in biodiversity informatics. First of all, it could improve species identification processes and tools 
by reducing the list of candidate species that are observable at a given location (be they automated, semi-automated 
or based on classical field guides or flora). More generally, it could facilitate biodiversity inventories through the 
development of location-based recommendation services (typically on mobile phones) as well as the involvement of 
non-expert nature observers. Last but not least, it might serve educational purposes thanks to biodiversity discovery 
applications providing functionalities such as contextualized educational pathways.

The rest of this documents presents (1) the data, and (2) the python code.

## 1. Data
The data are composed of two parts: the environmental rasters and the actual dataset containing all the occurrences. All the data is downloadable on the [AIcrowd page](https://www.aicrowd.com/challenges/lifeclef-2020-geo).
This section will describe both. You can check the 
[Protocol note](https://docs.google.com/document/d/19PF68B30HNSXq6_Rp6-Rd9GzOtGTsnHF_js4SkxqW3g/edit) for more 
details.
### Dataset of Occurrences 

The dataset is composed in multiple files:
- occurrences_fr_train.csv
- occurrences_fr_test.csv
- occurrences_us_train.csv
- occurrences_us_test.csv
- species_metadata.csv

More details about the dataset are given in the protocol note. The datasets columns include :

| Name        | Description   |
| ------------- |:-------------|
|id| The GLC20 reference identifier for the occurrence.|
|lat| Decimal latitude in the WGS84 coordinate system.|
|lon| Decimal longitude in the WGS84 coordinate system.|
|species_id| The GLC20 reference identifier for the species.|

Notice that the most important fields are Latitude and Longitude in order to extract the environmental patch and 
glc19SpId which contains the species ID. 
### High resolution tensors

The data contains also a tensor of high spatial resolution variables for each occurrences. The variables are the satelite images (in 4 chanels: Red, Green, Blue, Near Infra-Red), the altitude and the land cover. All the details on the extraction of these tensors and the manipulation of their original data sources are given in the [Protocol note](https://docs.google.com/document/d/19PF68B30HNSXq6_Rp6-Rd9GzOtGTsnHF_js4SkxqW3g/edit). Tensors are stored given the occurrences ids. The tensor of an occurrence with the id XXXXXABCD is at the location /CD/AB/XXXXXABCD.npy and /CD/AB/XXXXXABCD_alti.npy
### Environmental Rasters
The rasters are available directly on 
[AICrowd](https://www.aicrowd.com). The following variables are
available:

| Name        | Description           | Nature  | Values |
| ------------- |:-------------| :-----:|-----:|
| BIO_1      | Annual Mean Temp. (mean of monthly) | quanti. |[-116, 269]|
| BIO_2      | Max-temp - min-temp | quanti. |[-53, 361]|
| BIO_3      | Isothermality (100*2/7) | quanti. |[19, 69]|
| BIO_4      | Temp. seasonality (std.dev*100) | quanti. |[1624, 13302]|
| BIO_5      | Max Temp of warmest month | quanti. |[-25, 457]|
| BIO_6      | Min Temp of coldest month | quanti. |[-276, 183]|
| BIO_7      | Temp. annual range | quanti. |[117, 515]|
| BIO_8      | Mean temp. of wettest quarter | quanti. |[-169, 332]|
| BIO_9      | Mean temp. of driest quarter | quanti. |[-181, 331]|
| BIO_10      | Mean temp. of warmest quarter | quanti. |[-53, 361]|
| BIO_11      | Mean temp. of coldest quarter | quanti. |[-186, 220]|
| BIO_12      | Annual precipitations | quanti. |[-35, 3385]|
| BIO_13      | Precipitations of wettest month | quanti. |[7, 570]| 
| BIO_14      | Precipitations of driest month | quanti. |[0, 184]|
| BIO_15      | Precipitations seasonality (coef. of var.) | quanti. |[5, 140]|
| BIO_16      | Precipitations of wettest quarter | quanti. |[19, 1546]|
| BIO_17      | Precipitations of driest quarter | quanti. |[0, 612]|
| BIO_18      | Precipitations of warmest quarter | quanti. |[1, 777]|
| BIO_19      | Precipitations of coldest quarter | quanti. |[5, 1485]|
| BDTICM      | Absolute depth to bedrock in cm | quanti. |[0, 112467]|
| BLDFIE      | Bulk density in kg/m3 at 15cm depth | quanti. |[93, 1829]|
| CECSOL      | Cation exchange capacity of soil in cmolc/kg 15cm depth | quanti. |[0, 385]|
| CLYPPT      | Clay (0-2 micro meter) mass fraction at 15cm depth | quanti. |[0, 81]|
| ORCDRC      | Soil organic carbon content (fine earth fraction) in g per kg 15cm depth | quanti. |[0, 524]|
| PHIHOX      | Ph x 10 in H20 15cm depth  | quanti. |[32, 98]|
| SLTPPT      | Silt mass fraction at 15cm depth | quanti. |[0, 86]|
| SNDPPT      | Sand mass fraction at 15cm depth | quanti. |[0, 99]|


More details about each raster are available within the archive.
## 2. Python3
The file ```environmental_raster_glc.py``` provides to the participant of the GLC19 challenge a mean to extract 
environmental patches or vectors given the provided rasters. Providing a set of input rasters, it enables the online (in memory) extraction of environmental patches at a given spatial position OR of the offline construction (on disk) of all the patches of a set of spatial positions. 

The following examples are for *python3* but the code should work with *python2*.

The ```environmental_raster_glc.py``` follows two goals: 
* in code use,
* command line use.

In code use enables to extract environmental tensors on the go, for instance within a **Pytorch** dataset, thus reducing IO 
and improving training performances.

The command line use enables to export the dataset on disk.

### In code

In addition to the standard libraries, this code requires the following ones:
```python
import rasterio
import pandas
import numpy
import matplotlib
```
#### Constructing the Extractor
The core object to manipulate is the ```PatchExtractor``` which will manage the multiple available rasters.
Constructing the extractor only requires to set up the ```root_path``` of the rasters data:

```python
# constructing the extractor
extractor = PatchExtractor(root_path='/home/test/rasters')
```
By default, the extractor will return nx64x64 (where n depends on the rasters) patches. For custom size (other then 64),
 the constructor also accept an additional ```size``` parameter:

```python
# constructing the extractor
extractor = PatchExtractor(root_path='/home/test/rasters', size=256)
```
Attention, if size is too big, some patches from the dataset will be smaller due to an overflow in the raster map.

If size equals 1, then the extractor will return an environmental vector instead of the environmental tensor. 

Once the extractor is available, the rasters can be added. Two strategies are available : either adding all the rasters
 at once, or a one by one approach where specific transformation can be specified and some rasters avoided.

```python
# adding a single raster
extractor.append('clc', nan=0, normalized=True, transform=some_user_defined_function)
```
or
```python
# adding all the raster at root_path
extractor.add_all(nan=0, normalized=True, transform=some_user_defined_function)
```
In addition, some rasters are preferably used through a one hot encoding representation 
thus increasing the depth of the environmental tensor. The global parameter ```raster_metadata```
enables to set some of these properties on a per raster basis.

If parameters are not set, default values are used. For instance, nan have a default value
on a per raster basis. If you want to change it, either modify the metadata or set the parameter.

Please check the ```environmental_raster_glc.py``` file for more details.

#### Using the Extractor
The extractor acts as an array. For instance, ```len(extractor)``` gives the number of availble rasters.
 Accessing to a specific vector or tensor is done in the following way, by giving latitude and longitude:

```python
env_tensor = extractor[43.61, 3.88]
# env_tensor is a numpy array
```
Attention the shape of ```env_tensor``` does not necessarily corresponds to ```len(extractor)``` as some variables using 
a one hot encoding representation actually correspond to a deeper representation.

The extractor also enables to plot a specific patch:

```python
extractor.plot((43.61, 3.88))
# accept an optional style parameter to modify the style temporarily
```
Resulting in images of the following type:
![Rasters Patchs](https://raw.githubusercontent.com/maximiliense/GLC/master/patchs.jpg)

The plot method accept a ```cancel_one_hot``` parameter which value is True by default thus representing a variable 
initially set to have a one hot encoding as a single patch. In the previous image, ```clc``` is set to have
a one hot encoding representation but is plotted as a single patch.

In addition, Land cover, altitude, near-IR and RGB are provided:
![Npy Patchs](https://raw.githubusercontent.com/maximiliense/GLC/master/patchs_2.jpeg)

### Command line use

Using the online extraction of patches is fast but requires a significant amount of memory to store all rasters. So, for those who would rather
export the patches on disk, an additional functionality is provided.

The patch corresponding to the csv dataset will be extracted using the following command:
``` 
python3.7 extract_offline.py rasters_directory dataset.csv destination_directory
````

The destination_directory will be created if it does not exist yet. 
Its content might be erased if two files have the same name.

The extractor code has been conceived for low memory usage but might be slower in that sense.

Notice that the patch will be exported in numpy format. R library ```RcppCNPy``` enables to read
numpy format.

Help command returns:
```
usage: environmental_raster_glc.py [-h] [--size SIZE] [--normalized NORM]
                                   rasters dataset destination

extract environmental patches to disk

positional arguments:
  rasters            the path to the raster directory
  dataset            the dataset in CSV format
  destination        The directory where the patches will be exported

optional arguments:
  -h, --help         show this help message and exit
  --size SIZE        size of the final patch (default : 64)
  --normalized NORM  true if patch normalized (False by default)
```

Notice that some rasters (```proxi_eau_fast```  in particular) require a lots of memory and can be removed from 
the extraction by using the exception variable (in the ```extract_offline.py``` file). 
