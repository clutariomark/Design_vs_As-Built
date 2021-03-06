# Design vs As-Built ( Japanese Standard )

## Description


This stand-alone python scripts ouputs As-Built Management Form PDF from LAS dataset and LandXML dataset.

As-Built Management Form PDF follows the rule set by the Ministry of Land, Infrastructure, Transport and Tourism in Japan.

Surveyors, constructors etc company in Japan needs to follow this form to show their design vs as-built of their building structure to their client.

LAS file (As-built) and LandXML (Design) file is fed into this code to create a raster image that indicates design vs as-built.

The result of design vs as-built is statistically calculated following the rule mentioned above and set to PDF output.


#### excerpt of design vs as-built japanese standard

https://www.pref.yamagata.jp/documents/18357/24_re-za-dekigata.pdf

  

## Demo

![dvsab](dvsab.png)

  
  

## Requirement

- ArcGIS PRO ( arcpy library )

- Windows server

  

## Usage

```

~\ArcGIS\Pro\bin\Python\Scripts>propy.bat "~\design_vs_as_built.py" inputLas inputLandXML type method area cellFill

```

### Parameters

| index | parameter |description|
|--|--|--|
| 1 | inputLas | Input path of LAS file. Above 10cm resolution is preferred. The coordinates must be aligned with the LandXML file. |
| 2 | inputLandXML | Input path of LandXML file. The coordinates must be aligned with the LAS file. |
| 3 | type | The type of the construction.   |
| 4 | method | The method of the construction.  |
| 5 | area | The area of the construction. |
| 6 | cellFill | How you thin LAZ format. Options:'MINIMUM', 'MAXIMUM', 'CLOSEST_TO_MEAN'  |


#### Construction type
| type | method | area |
|--|--|--|
| road | cut | ground, slope|
| | fill | crown, slope |
| riverbed | cut | ground, slope |
|  | fill | crown, slope |
| pavement | not ready | - |


### Rasterizing method

**LAS to Raster**

It thins the LAS dataset first, then create TIN for calculaion, finally rasterize.

**LandXML to Raster**

It create esri TIN surface and then rasterize.


## Install
1. prepare Windows server
2. install ArcGIS pro to the window server
3. login by your account
4. run the command line
