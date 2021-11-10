from configparser import Interpolation
from Interpolator import RASTER
import os
import arcpy.conversion
import math
from utils import debug
from utils import get_workspace
from design_vs_as_built_sidecar import *
import argparse

class DesignVsAsBuilt_Calc:
    def __init__(self, settings, workspace, i_lxml, i_las):
        self._settings = settings
        self._workspace = workspace
        self.i_lxml = i_lxml
        self.i_las = i_las
        self._lxml_ras = ""
        self._las_ras = ""
        self._tin_trian = ""
        self._design_vs_as_built = ""
        self._thresh = ""
        self._ws(workspace)
        self.landxml2ras()
        self.las2ras()
        self.diff_rasters()
        
    def _ws(self, ws):
        arcpy.env.workspace = ws
        
    def landxml2ras(self):
        tin_triangles = []
        tin_rasters = []
        arcpy.ddd.LandXMLToTin(self.i_lxml, 
                               self._workspace, 
                               "TINs")
       
        TINList = arcpy.ListDatasets("*", "Tin")
        suffix = 0
        for tin in TINList:
            out_tin_triangle = "design_triangle"+str(suffix)
            tinTriangle = arcpy.ddd.TinTriangle(tin, out_tin_triangle, 
                                                "DEGREE", 1.0)
            tin_triangles.append(tinTriangle)
            
            output_tin_raster = "design_raster"+str(suffix)+".tif"
            arcpy.ddd.TinRaster(tin,
                                output_tin_raster,
                                data_type="FLOAT",
                                method="LINEAR",
                                sample_distance="CELLSIZE 1",
                                z_factor=1
                                )
            if suffix==0:
                arcpy.env.snapRaster = output_tin_raster
                debug("arcpy.env.snapRaster", arcpy.env.snapRaster)
                # arcpy.AddMessage(f"arcpy.env.snapRaster {arcpy.env.snapRaster}")
            
            tin_rasters.append(output_tin_raster)
            
            suffix=suffix+1
        # merged features
        out_tin_triangles = arcpy.CreateUniqueName("tin_triangles_merged.shp", 
                                                   self._workspace)
        arcpy.management.Merge(tin_triangles, out_tin_triangles)
        self._tin_trian = out_tin_triangles
        
        # merged rasters
        design_raster_path = arcpy.CreateUniqueName("design_raster_merged.tif", 
                                                    self._workspace)
        out_tin_rasters = arcpy.ia.Merge(tin_rasters, "FIRST")
        out_tin_rasters.save(design_raster_path)
        self._lxml_ras = design_raster_path
           
    def las2ras(self):
        interpolation_type = f'TRIANGULATION LINEAR WINDOW_SIZE {self._settings["cell_assignment"]} 1'
        out_las_raster = os.path.join(self._workspace, "las_raster.tif")
        arcpy.conversion.LasDatasetToRaster(
            self.i_las,
            out_las_raster,
            value_field="ELEVATION",
            # interpolation_type="BINNING IDW LINEAR",
            interpolation_type=interpolation_type,
            data_type="FLOAT",
            sampling_type="CELLSIZE",
            sampling_value=1,
            z_factor=1)
        self._las_ras = out_las_raster
   
    def _extract_by_mask(self, area_type):
        if(area_type == "ground") or (area_type == "crown"):
            TIN_feature = os.path.join(workspace, "TIN_flat.shp")
            slop_deg='Slope_Deg <= 10'
            output_mask = arcpy.CreateUniqueName(
                "design_flat_mask.shp", workspace)
            extracted_dvab_ras = arcpy.CreateUniqueName(
                "extracted_flat_design_vs_as_built.tif", workspace)
        else:
            TIN_feature = os.path.join(workspace, "TIN_slope.shp")
            slop_deg='Slope_Deg > 10'
            output_mask = arcpy.CreateUniqueName(
                "design_slope_mask.shp", workspace)
            extracted_dvab_ras = arcpy.CreateUniqueName(
                "extracted_slope_design_vs_as_built.tif", workspace)
                  
        # extract the feature from merge TIN triangle based on given slope degree
        out_layer_tin_flat = TIN_feature
        input_feature = arcpy.management.MakeFeatureLayer(
            self._tin_trian, out_layer_tin_flat, slop_deg)

        # make a polygon around TIN triangle. it's not TIN triangle.
        reference_scale = "25000"
        spatial_reference = arcpy.SpatialReference(4326) 
        margin = "0 meters"
        method = "CONVEX_HULL"
        mask_for_non_placed_anno = "ONLY_PLACED"
        attributes = "ALL"
        arcpy.FeatureOutlineMasks_cartography(input_feature,
                                            output_mask,
                                            reference_scale,
                                            spatial_reference,
                                            margin, method,
                                            mask_for_non_placed_anno,
                                            attributes)
        extracted_dvab = arcpy.sa.ExtractByMask(self._design_vs_as_built, output_mask)
        extracted_dvab.save(extracted_dvab_ras)
        return extracted_dvab
        
    def diff_rasters(self):
        design_vs_as_built = arcpy.sa.Raster(self._las_ras) - arcpy.sa.Raster(self._lxml_ras)
        self._design_vs_as_built = arcpy.CreateUniqueName("design_vs_as_built.tif", workspace)
        design_vs_as_built.save(self._design_vs_as_built)        
        
    def _get_las_property(self):
        domain_area = arcpy.CreateUniqueName("domain_area.shp", self._workspace)
        arcpy.ddd.RasterDomain(self._design_vs_as_built, domain_area, "POLYGON")
        arcpy.ddd.ExtractLas(self.i_las, self._workspace, "MAXOF", domain_area)
        extracted_domain_las = os.path.join(workspace, os.path.basename(las_dataset))
        las_properties = arcpy.Describe(extracted_domain_las)
        
        return las_properties

    def calc_stats(self):
        # output: dvab extracted by area type
        dvab_part = self._extract_by_mask(self._settings["area_type"])
        pixel_block = dvab_part.blockSize
        elev_list = []
        for i in range(pixel_block[0]):
            for j in range(pixel_block[1]):
                if(math.isnan(dvab_part[i, j]) == False):
                    elev_list.append(dvab_part[i, j])
        count_point = self._get_las_property().pointCount
        mean_val = round(sum(elev_list)/len(elev_list), 12) #平均値
        max_val = round(max(elev_list), 12) #最大値差
        min_val= round(min(elev_list), 12) #最小値差
        stats = {
            "count": count_point,
            "mean": mean_val,
            "max": max_val,
            "min": min_val
            }
        
        thresh = self._determine_standard_value()
        reclassified = self._eval_dvab(reclass_raster_path, dvab_part, thresh)
        self._thresh = thresh

        arcpy.AddField_management(reclassified, "Area", "DOUBLE")
        cell_size = arcpy.Describe(reclassified).meanCellHeight
        arcpy.CalculateField_management(in_table=reclassified,
                                        field="Area",
                                        expression="!Count! * {}".format(cell_size ** 2))
        with arcpy.da.SearchCursor(reclass_raster_path, ["Value", "Area"]) as rows:
            data = {r[0]: r[1] for r in rows}
            debug("data", data)
            # arcpy.AddMessage("data: {}".format(data))

        # interpolation for missing statistic results
        coverage = [i for i in range(10)]
        data_keys = [key for (key, value) in data.items()]
        for n in coverage:
            if(n not in data_keys):
                data[n] = 0
        data_sorted = {}
        for key in sorted(data):
            data_sorted[key] = data[key]
        data=data_sorted
        
        return {"stats": stats, "data": data}
    
    def _determine_standard_value(self):
        cons_type1 =  self._settings["cons_type1"]
        cons_type2 =  self._settings["cons_type2"]
        area_type =  self._settings["area_type"]

        mean_thresh_low = CONSTRUCTION[cons_type1][cons_type2][area_type]["mean"][0]
        mean_thresh_up = CONSTRUCTION[cons_type1][cons_type2][area_type]["mean"][1]
        min = CONSTRUCTION[cons_type1][cons_type2][area_type]["minmax"][0]
        max = CONSTRUCTION[cons_type1][cons_type2][area_type]["minmax"][1]
        return {
            "mean_low": mean_thresh_low,
            "mean_up": mean_thresh_up,
            "min": min,
            "max": max}

    def _eval_dvab(self, path, DvsAB_raster, standard):
        
        per20, per50, per80 = 0.2, 0.5, 0.8
        remap = arcpy.sa.RemapRange([
            [-999, standard["min"], 0],                  # as-built < design: ~ -100%
            [standard["min"], standard["min"]*per80, 1],       # on design: -100% ~ -80%
            [standard["min"]*per80, standard["min"]*per50, 2], # on design: -80% ~ -50%
            [standard["min"]*per50, standard["min"]*per20, 3], # on design: -50% ~ -20%
            [standard["min"]*per20, 0, 4],               # on design: -20% ~ 0%
            [0, standard["max"]*per20, 5],                # on design:  0% ~ 20%
            [standard["max"]*per20, standard["max"]*per50, 6],   # on design:  20% ~ 50%
            [standard["max"]*per50, standard["max"]*per80, 7],   # on design:  50% ~ 80%
            [standard["max"]*per80, standard["max"], 8],         # on design:  80% ~ 100%
            [standard["max"], 999, 9]                     # as-built > design:  100% ~
        ])
        reclassified = arcpy.sa.Reclassify(DvsAB_raster, "Value", remap, "DATA")
        reclassified.save(path)
        debug("Added to Map and saved to:", path)
        colormap=r"C:\Users\kjohn\Documents\GitHub\Design_vs_As-Built\colormap.clr"
        debug("colormap:", colormap)
        arcpy.management.AddColormap(path, "#", colormap)
        return reclassified
        # customized_colormap=[
        #     [0,0,0,0],
        #     [1,255,0,0],
        #     [2,137,225,255],
        #     [3,255,140,102],
        #     [4,255, 207, 191],
        #     [5,0 ,255, 0],
        #     [6,50 ,205 ,0],
        #     [7,0 ,191 ,255],
        #     [8,0 ,0 ,255],
        #     [9,0 ,0 ,0],
        #     [10,255 ,255 ,0]]
        # recolored = arcpy.ia.Colormap(reclassified, colormap=customized_colormap)
        # recolored.save(path)
        # return recolored
        
    def determine_stats(self, stats):
        cons_type1 = self._settings["cons_type1"]
        cons_type2 = self._settings["cons_type2"]
        area_type = self._settings["area_type"]
        
        thresh = self._thresh
        stat = stats["stats"]
        data = stats["data"]
        eval_area = sum(data.values()) #評価面積
        rejected_points = data[0] + data[9]
        on_design_points = data[1] + data[2] + data[2] + data[3] + data[4] + data[5] + data[6] + data[7] + data[8]
        on_design_points_80 = data[2] + data[3] + data[4] + data[5] + data[6] + data[7]
        on_design_points_50 = data[2] + data[3] + data[4]
        percentage_80 = round(on_design_points_80 / on_design_points * 100, 1)
        percentage_50 = round(on_design_points_50 / on_design_points * 100, 1)
        percentage_rejected_points = round(rejected_points / (on_design_points + rejected_points) *100, 1)

        determine_mean_val = False
        determine_max_val = False
        determine_min_val = False
        determine_data_num_val = False
        determine_rejected_points = False
        determine_overall_txt = "不合格"

        if percentage_rejected_points <= 0.3:
            determine_rejected_points = True
        if (stat["mean"] <= thresh["mean_up"]) and (stat["mean"] >= thresh["mean_low"]):
            determine_mean_val = True
        if stat["max"] <= thresh["max"]:
            determine_max_val = True
        if stat["min"] >= thresh["min"]:
            determine_min_val = True
        if stat["count"] > eval_area:
            determine_data_num_val = True
        if determine_mean_val == True and determine_max_val == True and determine_min_val == True and determine_data_num_val == True:
            determine_overall_txt = "合格"

        return {
            "numbers": {
                "Standard1": abs(round(thresh["mean_low"]*1000,0)), 
                "Standard2": round(thresh["max"]*1000, 0),
                "Standard3": round(thresh["min"]*1000, 0),
                "Determine Overall": determine_overall_txt,
                "Type1": CONSTRUCTION[cons_type1]['jpname'],
                "Type2": CONSTRUCTION[cons_type1][cons_type2]['jpname'],
                "Measure Title 1": CONSTRUCTION[cons_type1][cons_type2][area_type]['jpname'],
                "Means Value": round(stat["mean"]*1000,2), 
                "Max Value": round(stat["max"]*1000,2), 
                "Min Value": round(stat["min"]*1000,2),
                "Data num Value": stat["count"], 
                "Area value": eval_area,
                "Min Data Num": eval_area, 
                "Rejected Value": rejected_points,
                "Percent 80 Data num": on_design_points_80, 
                "Percent 80 Value": str(percentage_80) + "%",
                "Percent 50 Data num": on_design_points_50,
                "Percent 50 Value": str(percentage_50) + "%"},
            "determine": {
                "Determine Means": determine_mean_val, 
                "Determine Max": determine_max_val, 
                "Determine Min": determine_min_val,
                "Determine Data num": determine_data_num_val,
                "Determine Rejected": determine_rejected_points
                }
            }

# =================== INPUT PARAMETERS ========================
# las_dataset = validate_data_source(arcpy.GetParameter(0))
# landxml = validate_data_source(arcpy.GetParameter(1))
# # Controller determines the construction types
# cons_type1 = arcpy.GetParameter(2)
# cons_type2 = arcpy.GetParameter(3)
# area_type = arcpy.GetParameter(4)
# cell_assgnmnt_las2ras = arcpy.GetParameter(5)# MAXIMUM, MINIMUM, CLOSEST_TO_MEAN

parser = argparse.ArgumentParser(description='Design vs As-Built')
parser.add_argument('lasfile', type=str, metavar='', help='input las file')
parser.add_argument('landxmlfile', type=str, metavar='', help='input landxml file')
parser.add_argument('cons_type', type=str, default='road', metavar='', help='input construction type')
parser.add_argument('cons_method', type=str, default='cut', metavar='', help='input construction method')
parser.add_argument('area_type', type=str, default='ground', metavar='', help='input area type')
parser.add_argument('las2ras_thin', type=str, default='MINIMUM', metavar='', help='input thinning method')
args = parser.parse_args()

las_dataset = args.lasfile
landxml = args.landxmlfile
cons_type1 = args.cons_type
cons_type2 = args.cons_method
area_type = args.area_type
cell_assignment = args.las2ras_thin
THIN_METHOD = ['MINIMUM', 'MAXIMUM', 'CLOSEST_TO_MEAN']
debug("las_dataset", las_dataset)
# arcpy.AddMessage(f"{las_dataset} {landxml} {cons_type1} {cons_type2} {area_type} {cell_assgnmnt_las2ras}")

workspace = get_workspace(os.path.splitext(os.path.basename(las_dataset))[0]+ "_" + str(area_type))
# if not os.path.exists(workspace):
#     os.mkdir(workspace)
arcpy.env.workspace = workspace
arcpy.env.overwriteOutput = True

reclass_raster_path = arcpy.CreateUniqueName("reclass.tif", workspace)
pdf_path = arcpy.CreateUniqueName(f"DesisgnvsAsbuilt_{area_type}.pdf", workspace)
# las2ras_path = arcpy.CreateUniqueName("las2ras.tif", workspace)
# hillshade_path = arcpy.CreateUniqueName("hillshade.tif", workspace)
# diff_raster_path = arcpy.CreateUniqueName("design_vs_as_built.tif", workspace)
# colormap_raster_path = arcpy.CreateUniqueName("colormap.tif  ", workspace)
# domain_raster_path = arcpy.CreateUniqueName("raster_domain.shp", workspace)
# TINs_path = arcpy.CreateUniqueName("TINs", workspace)
# design_raster_path = arcpy.CreateUniqueName("design_raster.tif", workspace)

# =================== PROCESS ========================================================
settings = {"cons_type1": cons_type1,
"cons_type2": cons_type2,
"area_type": area_type,
"cell_assignment": cell_assignment}
design_asbuilt_calc = DesignVsAsBuilt_Calc(
    settings, workspace, landxml, las_dataset)
stats = design_asbuilt_calc.calc_stats()
res = design_asbuilt_calc.determine_stats(stats)
# active_aprx = arcpy.mp.ArcGISProject("CURRENT")
# activeMap = active_aprx.activeMap
# arcpy.mp.ArcGISProject("CURRENT").activeMap.addDataFromPath(reclass_raster_path)
# lyr = activeMap.listLayers()[0]
# tempplate_path=r"HaulRoad.aprx"
# it needs to be full path
tempplate=r"C:\Users\kjohn\Documents\GitHub\Design_vs_As-Built\HaulRoad.aprx"
DesignVsAsBuilt_instance = DesignVsAsBuilt(tempplate)
# DesignVsAsBuilt_instance.addlayer(lyr)
DesignVsAsBuilt_instance.addDataFromPath(reclass_raster_path)
DesignVsAsBuilt_instance.setResults(res)
DesignVsAsBuilt_instance.export(pdf_path)
os.startfile(pdf_path)

