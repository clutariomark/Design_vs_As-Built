from Interpolator import RASTER
import arcpy
from utils import get_workspace

def validate_data_source(source_layer):

    if hasattr(source_layer, "dataSource"):
        source_layer = source_layer.dataSource
    else:
        source_layer = str(source_layer)
    return source_layer

class DesignVsAsBuilt:
    def __init__(self, aprx_path):
        self._aprx = arcpy.mp.ArcGISProject(aprx_path)
        self._layout = self._aprx.listLayouts("Design Vs As-Built")[0]
        self._mf = self._layout.listElements("MAPFRAME_ELEMENT")[0]

    @property
    def layouts(self):
        return {layout.name: layout for layout in self._aprx.listLayouts()}

    @property
    def maps(self):
        return {mp.name: mp for mp in self._aprx.listMaps()}

    def addlayer(self, layer):
        self._mf.map.addLayer(layer, "TOP")
        first_layer = self._mf.map.listLayers()[0]
        self._extent = self._mf.getLayerExtent(first_layer)
        self._mf.camera.setExtent(self._extent)
        
    def addDataFromPath (self, layer_path):
        first_layer = self._mf.map.addDataFromPath(layer_path)
        self._extent = self._mf.getLayerExtent(first_layer)
        self._mf.camera.setExtent(self._extent)

    def export(self, output):
        self._layout.exportToPDF(output)

    def setResults(self, result):
        element = self._layout.listElements("TEXT_ELEMENT")
        res_num = result["numbers"]
        res_det = result["determine"]
        for e in element:
            for name, val in res_num.items():
                if e.name == name:
                    e.text = val
            for name, val in res_det.items():
                if e.name == name:
                    if val:
                        e.text = "問題なし"
                    else:
                        e.text = "異常値あり"
                             
CONSTRUCTION = {
    "road": {
        "cut": {
            "ground": {
                "mean": [-0.05, 0.05],
                "minmax": [-0.15, 0.15],
                "jpname": "平場"
                },
            "slope": {
                "mean": [-0.07, 0.07],
                "minmax": [-0.16, 0.16],
                "jpname": "法面"
                },
            "jpname": "掘削工"
            },
        "fill" :{
            "crown":{
                "mean": [-0.05, 0.05],
                "minmax": [-0.15, 0.15],
                "jpname": "天端"
            },
            "slope": {
                "mean": [-0.08, 0.08],
                "minmax": [-0.19, 0.19],
                "jpname": "法面"
            },
            "jpname": "盛土工"
        },
        "jpname": "道路土工"
    },
    "riverbed": {
        "cut": {
            "ground": {
                "mean": [-0.05, 0.05],
                "minmax": [-0.15, 0.15],
                "jpname": "平場",        
                },
            "slope": {
                "mean": [-0.07, 0.07],
                "minmax": [-0.16, 0.16],
                "jpname": "法面",  
                },
            "jpname": "掘削工"
            },
        "fill" :{
            "crown":{
                "mean": [-0.05, 100.0],
                "minmax": [-0.15, 100.0],
                "jpname": "天端",
            },
            "slope": {
                "mean": [-0.05, 100.0],
                "minmax": [-0.17, 100.0],
                "jpname": "法面"
            },
            "jpname": "盛土工",
        },
        "jpname": "河川土工",
    },
    "pavement": {
        "surface_layer":{
            "small_scale":{
                "mean": [-0.003, 100.0],
                "max": [-0.020, 100.0],
                "min": [-0.020, 100.0],
                "jpname": "小規模",
            },
            "medium_scale":{
                "mean": [-0.002, 100.0],
                "max": [-0.017, 100.0],
                "min": [-0.017, 100.0],
                "jpname": "中規模",
            },
            "jpname": "表層",
        },
        "base_layer":{
            "small_scale":{
                "mean": [-0.004, 100.0],
                "max": [-0.024, 100.0],
                "min": [-0.024, 100.0],
                "jpname": "小規模",
            },
            "medium_scale":{
                "mean": [-0.003, 100.0],
                "max": [-0.024, 100.0],
                "min": [-0.024, 100.0],
                "jpname": "中規模",
            },
            "jpname": "基層",
        },
        "upper_roadbed":{
            "small_scale":{
                "mean": [-0.064, 100.0],
                "max": [-0.01, 100.0],
                "min": [-0.01, 100.0],
                "jpname": "小規模",
            },
            "medium_scale":{
                "mean": [-0.053, 100.0],
                "max": [-0.008, 100.0],
                "min": [-0.008, 100.0],
                "jpname": "中規模",
            },
            "jpname": "上層路盤",
            
        },
        "lower_roadbed":{
            "small_scale":{
                "mean": [-0.09, 0.09],
                "max": [-0.015, 0.04],
                "min": [-0.015, 0.04],
                "jpname": "小規模",
            },
            "medium_scale":{
                "mean": [-0.09, 0.09],
                "max": [-0.015, 0.04],
                "min": [-0.015, 0.04],
                "jpname": "中規模",
            },
            "jpname": "下層路盤",
        },
        "jpname": "舗装工",
    }
}



# def Eval_DvsAB(path, DvsAB_raster, standard):

#     per20, per50, per80 = 0.2, 0.5, 0.8
#     standard = float(standard)

#     remap = arcpy.sa.RemapRange([
#         [-999, -standard, 0],                  # as-built < design: ~ -100%
#         [-standard, -standard*per80, 1],       # on design: -100% ~ -80%
#         [-standard*per80, -standard*per50, 2], # on design: -80% ~ -50%
#         [-standard*per50, -standard*per20, 3], # on design: -50% ~ -20%
#         [-standard*per20, 0, 4],               # on design: -20% ~ 0%
#         [0, standard*per20, 5],                # on design:  0% ~ 20%
#         [standard*per20, standard*per50, 6],   # on design:  20% ~ 50%
#         [standard*per50, standard*per80, 7],   # on design:  50% ~ 80%
#         [standard*per80, standard, 8],         # on design:  80% ~ 100%
#         [standard, 999, 9]                     # as-built > design:  100% ~
#     ])
    
#     reclassified = arcpy.sa.Reclassify(DvsAB_raster, "Value", remap, "DATA")
    
#     reclassified.save(path)
#     arcpy.AddColormap_management(path, "#","./colormap.clr")
#     # arcpy.AddMessage("Added to Map and saved to:{}".format(path))

#     return reclassified