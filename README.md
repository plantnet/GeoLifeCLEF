# GeoLifeCLEF 2019
Automatically predicting the list of species that are the most likely to be observed at a given location is useful for 
many scenarios in biodiversity informatics. First of all, it could improve species identification processes and tools 
by reducing the list of candidate species that are observable at a given location (be they automated, semi-automated 
or based on classical field guides or flora). More generally, it could facilitate biodiversity inventories through the 
development of location-based recommendation services (typically on mobile phones) as well as the involvement of 
non-expert nature observers. Last but not least, it might serve educational purposes thanks to biodiversity discovery 
applications providing functionalities such as contextualized educational pathways.

The rest of this documents presents (1) the data, (2) the python code and the (3) R code.

## 1. Data
The data are composed of two parts: the environmental rasters and the actual dataset containing all the occurrences.
This section will describe both.
### Environmental Rasters
The rasters are available at the following 
[URL](http://otmedia.lirmm.fr/LifeCLEF/GeoLifeCLEF2019/). The following variables are
available:

| Name        | Description           | Nature  | Values |
| ------------- |:-------------| :-----:|-----:|
| CHBIO_1      | Annual Mean Temp. (mean of monthly) | quanti. |[-10.7,18.4]|
| CHBIO_2      | Max-temp - min-temp | quanti. |[7.8, 21.0]|
| CHBIO_3      | Isothermality (100*2/7) | quanti. |[41.1,60.0]|
| CHBIO_4      | Temp. seasonality (std.dev*100) | quanti. |[302.7, 777.8]|
| CHBIO_5      | Max Temp of warmest month | quanti. |[6.1,36.6]|
| CHBIO_6      | Min Temp of coldest month | quanti. |[-28.3,5.4]|
| CHBIO_7      | Temp. annual range | quanti. |[16.7,42.0]|
| CHBIO_8      | Mean temp. of wettest quarter | quanti. |[-14.2,23.0]|
| CHBIO_9      | Mean temp. of driest quarter | quanti. |[-17.7,26.5]|
| CHBIO_10      | Mean temp. of warmest quarter | quanti. |[-2.8, 26.5]|
| CHBIO_11      | Mean temp. of coldest quarter | quanti. |[-17.7, 11.8]|
| CHBIO_12      | Annual precipitations | quanti. |[318.3,2543.3]|
| CHBIO_13      | Precipitations of wettest month | quanti. |[43.0,285.5]| 
| CHBIO_14      | Precipitations of driest month | quanti. |[3.0,135.6]|
| CHBIO_15      | Precipitations seasonality (coef. of var.) | quanti. |[8.2,26.5]|
| CHBIO_16      | Precipitations of wettest quarter | quanti. |[121.6,855.6]|
| CHBIO_17      | Precipitations of driest quarter | quanti. |[19.8,421.3]|
| CHBIO_18      | Precipitations of warmest quarter | quanti. |[198,851.7]|
| CHBIO_19      | Precipitations of coldest quarter | quanti. |[60.5,520.4]|
| etp      | Potential evapo transpiration | quanti. |[133,1176]|
| alti      | Elevation | quanti. |[-188,4672]|
| awc_top      | Topsoil available water capacity | ordinal |{0,120,165,210}|
| bs_top      | Base saturation of the topsoil | ordinal |{35,62,85}|
| cec_top      | Topsoil cation exchange capacity | ordinal |{7,22,50}|
| crusting      | Soil crusting class | ordinal | [0,5] |
| dgh      | Depth to a gleyed horizon | ordinal |{20,60,140}|
| dimp      | Depth to an impermeable layer | ordinal |{60,100}|
| erodi      | Soil erodibility class | ordinal |[0,5]|
| oc_top      | Topsoil organic carbon content | ordinal |{1,2,4,8}|
| pd_top      | Topsoil packing density | ordinal |{1,2}|
| text      | Dominant surface textural class| ordinal |[0,5]|
| proxi_eau_fast      | <50 meters to fresh water| boolean |{0,1}|
| clc      | Ground occupation | categorial |[1,48]|

More details about each raster are available within the archive.
### Dataset of Occurrences 
The dataset itself is available at the following 
[URL](http://otmedia.lirmm.fr/LifeCLEF/GeoLifeCLEF2019/). The following variables are

The dataset is composed in multiple files:
- PL_complete.csv
- PL_trusted.csv
- noPlant.csv 
- GLC_2018.csv


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
![alt Patchs](https://raw.githubusercontent.com/maximiliense/GLC19/master/patchs.jpg)

The plot method accept a ```cancel_one_hot``` parameter which value is True by default thus representing a variable 
initially set to have a one hot encoding as a single patch. In the previous image, ```clc``` is set to have
a one hot encoding representation but is plotted as a single patch.

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

## 3. R
@Christophe