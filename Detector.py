import arcpy
import os
# import traceback
import numpy as np
import matplotlib.pyplot as plt
from Segment import Segment
from Interpolator import *
from config_manager import *
from centerlines import *
import xlwt
import math
from collections import namedtuple

IN_MEMORY       =  "in_memory"
ROAD_FC         =  "detected_roads.shp"
WINDROW_FC      =  "detected_windrows.shp"
SECTIONS_FC     =  "cross_sections.shp"
SLOPES_FC       =  "slopes.shp"
CLIP_POLY       =  "clip_polygon.shp"
GRADIENT_FC     =  "road_gradient.shp"
GRADIENT_3D     =  "road_gradient_3d.shp"
CROSSFALL_FC    =  "crossfall.shp"
STREAMS         =  "water_streams.shp"
SINKS           =  "sinks.shp"

Field = namedtuple("Field", 'name type')
road_fields = [Field(name="IDX", type="SHORT"),
               Field(name="Length", type="DOUBLE"),
               Field(name="Slope", type="DOUBLE"),
               Field(name="SlopeDeg", type="DOUBLE"),
               Field(name="Compliance", type="TEXT"),
               Field(name="Profile", type="TEXT")]

windrow_fields = [Field(name="Height", type="DOUBLE"),
                  Field(name="Compliance", type="TEXT")]

gradient_fields = [Field(name="Slope", type="DOUBLE"),
                   Field(name="SlopeDeg", type="DOUBLE")]

crossfall_fields = [Field(name="Slope", type="DOUBLE"),
                    Field(name="SlopeDeg", type="DOUBLE")]

config = Config()

class Detector:
    """ The detector object job it to detect and computes statistics
     for the roads segments, windrows, gradient and crossfall strings"""

    def __init__(self, center_lines: list, las_path: str, working_dir: str, type: str):

        arcpy.AddMessage(f"Initialize {self.__class__.__name__} object")
        self.las = las_path
        self.workspace = working_dir
        self.raster = os.path.join(IN_MEMORY, "las2ras")
        arcpy.AddMessage(f"Converting las to raster")
        arcpy.LasDatasetToRaster_conversion(in_las_dataset=self.las,
                                            out_raster=self.raster,
                                            value_field='ELEVATION',
                                            interpolation_type='BINNING AVERAGE NONE',
                                            data_type='FLOAT',
                                            sampling_type='CELLSIZE',
                                            sampling_value='0.5')

        # Filters & processing configurations
        self.min_area = center_lines[0].min_area
        self.windrow_min_distance_from_center = 0.1 if type == "ROM" else center_lines[0].width / 6.5
        self.windrow_max_distance_from_center = center_lines[0].windrow_max_distance_from_center
        self.gradient_spacing_distance = center_lines[0].gradient_spacing_distance
        self.spacing_distance = center_lines[0].cross_sections_spacing_distance
        self.max_spacing_distance_filter = center_lines[0].max_spacing_distance_filter
        self.cross_section_multiplier = center_lines[0].cross_section_multiplier
        self.road_min_length_filter = center_lines[0].road_min_length_filter
        self.road_max_length_filter = center_lines[0].road_max_length_filter
        self.windrow_min_height_filter = center_lines[0].windrow_min_height_filter
        self.windrow_max_height_filter = center_lines[0].windrow_max_height_filter
        self.generate_plots = bool(int(config.get("Windrow Heights", "generate_plots")))

        self.domain = self._getDomain(self.raster)
        self.center_lines = self._chopCenterLines(center_lines, buffer_distance=-5)
        self.interpolator = Interpolator(self.raster)
        self.spatial_reference = self.center_lines[0].geometry.spatialReference

        self.hydro_worker = None

        # Layers
        self.cross_sections_2d = None
        self.slopes = None
        self.cross_sections_3d = None
        self.detected_roads_2d = None
        self.detected_roads_3d = None
        self.detected_windrows = None

        gradient_fc = os.path.join(self.workspace, GRADIENT_3D)
        self.road_gradient = gradient_fc if os.path.exists(gradient_fc) else None
        self.road_crossfall = None
        self.road_streams = None
        self.road_sinks = None
        self.hydro_aoi = None
        self.hydro_raster = None

        self.windrow1_points = []
        self.windrow2_points = []

    def CreateCrossSections2D(self):

        arcpy.AddMessage("Creating cross sections in 2D")
        arcpy.CreateFeatureclass_management(out_path=self.workspace,
                                            out_name=SECTIONS_FC,
                                            geometry_type='POLYLINE',
                                            spatial_reference=self.spatial_reference)

        self.cross_sections_2d = os.path.join(self.workspace, SECTIONS_FC)

        for field in road_fields:
            arcpy.AddField_management(self.cross_sections_2d, field.name, field.type)

        idx = 0
        for centerline in self.center_lines:
            points = [pnt for pnt in centerline.geometry]
            segments = [Segment(origin, destination) for origin, destination in zip(points[0][:-1], points[0][1:])]

            distance = 0.0
            with arcpy.da.InsertCursor(self.cross_sections_2d, ["SHAPE@", "IDX"]) as iCur:
                for s in segments:
                    while distance <= s.Length:
                        cs = s.compute_perpendicular_segment(distance, (self.cross_section_multiplier * centerline.width))
                        iCur.insertRow((cs.toArcpyPolyline(self.spatial_reference), idx))
                        idx += 1
                        distance += self.spacing_distance
                    distance = distance - s.Length

        self._validateOutputFeatureClass(self.cross_sections_2d)
        return self.cross_sections_2d

    def DetectRoads(self):

        if not self.slopes:
            self.slopes = self._createSlopes()

        domain = None
        with arcpy.da.SearchCursor(self.domain, ["SHAPE@"]) as cur:
            for row in cur:
                domain = row[0]

        slopes = None
        with arcpy.da.SearchCursor(self.slopes, ["SHAPE@"]) as cur:
            for row in cur:
                slopes = row[0]

        arcpy.AddMessage("Computing symmetric difference")
        diff = domain.symmetricDifference(slopes)

        arcpy.CreateFeatureclass_management(self.workspace, CLIP_POLY,
                                            geometry_type='POLYGON',
                                            spatial_reference=self.spatial_reference)

        with arcpy.da.InsertCursor(os.path.join(self.workspace, CLIP_POLY), ["SHAPE@"]) as iCur:
            iCur.insertRow([diff])

        arcpy.AddMessage("Clipping the cross sections")
        self.detected_roads_2d = os.path.join(self.workspace, ROAD_FC)
        arcpy.Clip_analysis(in_features=self.cross_sections_2d,
                            clip_features=diff,
                            out_feature_class=self.detected_roads_2d)

        arcpy.AddMessage("Filtering the clipped cross sections")
        with arcpy.da.UpdateCursor(self.detected_roads_2d, ["SHAPE@"]) as cur:
            for row in cur:

                length_list = [arcpy.Polyline(g).length for g in row[0]]
                geom = arcpy.Polyline(row[0].getPart(length_list.index(max(length_list))))

                if geom.length < self.road_min_length_filter or geom.length > self.road_max_length_filter:
                    cur.deleteRow()
                    continue

                cur.updateRow([geom])

        self._validateOutputFeatureClass(self.detected_roads_2d)
        return self.detected_roads_2d

    def DetectWindrows(self):

        if self.detected_roads_2d is None:
            raise Exception("Cannot detect the windrows without 2D detected roads layer")

        if self.cross_sections_3d is None:
            self.cross_sections_3d = self.interpolator.Interpolate3D(self.cross_sections_2d, vertices_only=False)
            self._validateOutputFeatureClass(self.cross_sections_3d)

        if self.detected_roads_3d is None:
            self.detected_roads_3d = self._export3dRoads()

        with arcpy.da.UpdateCursor(self.detected_roads_3d, ["SHAPE@", "IDX", "Profile"]) as rCur:
            for road in rCur:
                roadId = road[1]
                where = f'"IDX" = {roadId}'
                with arcpy.da.SearchCursor(self.cross_sections_3d, ["SHAPE@", "IDX"], where_clause=where) as sCur:
                    for section in sCur:

                        # compute the barrier points
                        rp1, rp2 = road[0].firstPoint, road[0].lastPoint
                        bp1, bp2 = self._detectBarrierHeightsPointsWithinDistance(rp1=rp1, rp2=rp2, polyline=section[0])

                        if isinstance(bp1, type(None)) and isinstance(bp2, type(None)):
                            continue

                        elif isinstance(bp1, type(None)):
                            height2 = bp2.Z - rp2.Z
                            if self.windrow_min_height_filter < height2 < self.windrow_max_height_filter:
                                bp2.Z = height2
                                self.windrow2_points.append((bp2, roadId))

                        elif isinstance(bp2, type(None)):
                            height1 = bp1.Z - rp1.Z
                            if self.windrow_min_height_filter < height1 < self.windrow_max_height_filter:
                                bp1.Z = height1
                                self.windrow1_points.append((bp1, roadId))

                        else:
                            if self.generate_plots:
                                road[2] = self._saveProfile(section[0], bp1, bp2, rp1, rp2, roadId)

                            height1 = bp1.Z - rp1.Z
                            if self.windrow_min_height_filter < height1 < self.windrow_max_height_filter:
                                bp1.Z = height1
                                self.windrow1_points.append((bp1, roadId))

                            height2 = bp2.Z - rp2.Z
                            if self.windrow_min_height_filter < height2 < self.windrow_max_height_filter:
                                bp2.Z = height2
                                self.windrow2_points.append((bp2, roadId))

                        # Update the cursor with the updated list
                        rCur.updateRow(road)
                        break

            # Sorting by the original idx
            self.windrow1_points = sorted(self.windrow1_points, key=lambda x: x[1])
            self.windrow2_points = sorted(self.windrow2_points, key=lambda x: x[1])
        self.detected_windrows = self._exportWindrows()
        return self.detected_windrows

    def DetectWaterStreams(self, z_limit=0.2, filter_value=200):

        # if self.hydro_aoi is None:
        #     self.hydro_aoi = self._createHydrologyDomain(buffer_distance=5)
        #
        # if self.hydro_worker is None:
        #     self.hydro_worker = HydroWorker(domain=self.hydro_aoi, surface=self.raster, workspace=self.workspace)
        #
        # return self.hydro_worker.detect_water_streams(z_limit=0.2, filter_value=200)

        if self.hydro_aoi is None:
            self.hydro_aoi = self._createHydrologyDomain(buffer_distance=5)

        if self.hydro_raster is None:
            self.hydro_raster = arcpy.sa.ExtractByMask(self.raster, self.hydro_aoi)

        arcpy.AddMessage("Detecting water streams")
        fill_raster = arcpy.sa.Fill(in_surface_raster=self.hydro_raster, z_limit=z_limit)
        flow_direction = arcpy.sa.FlowDirection(in_surface_raster=fill_raster)
        flow_accumulation = arcpy.sa.FlowAccumulation(in_flow_direction_raster=flow_direction, data_type="FLOAT")
        above = arcpy.sa.Con(flow_accumulation, flow_accumulation, where_clause=f"value>{filter_value}")
        streams_raster = arcpy.sa.StreamLink(in_stream_raster=above, in_flow_direction_raster=flow_direction)

        return self._exportWaterStreams(streams_raster, flow_direction)

    def DetectSinks(self, z_limit=0.05, comparison_tolerance=0.0001, area_filter=2.0):

        # if self.hydro_aoi is None:
        #     self.hydro_aoi = self._createHydrologyDomain(buffer_distance=5)
        #
        # if self.hydro_worker is None:
        #     self.hydro_worker = HydroWorker(domain=self.hydro_aoi, surface=self.raster, workspace=self.workspace)
        #
        # return self.hydro_worker.detect_sinks(z_limit=0.05, comparison_tolerance=0.0001, area_filter=2.0)

        if self.hydro_aoi is None:
            self.hydro_aoi = self._createHydrologyDomain(buffer_distance=5)

        if self.hydro_raster is None:
            self.hydro_raster = arcpy.sa.ExtractByMask(self.raster, self.hydro_aoi)

        arcpy.AddMessage("Detecting sinks")
        fill_by_z = arcpy.sa.Fill(in_surface_raster=self.hydro_raster, z_limit=z_limit)
        fill_all = arcpy.sa.Fill(in_surface_raster=self.hydro_raster)
        substract_raster = arcpy.Minus_3d(in_raster_or_constant1=fill_all, in_raster_or_constant2=fill_by_z)
        sinks_raster = arcpy.sa.Con(substract_raster, "1", "0", where_clause=f"VALUE>{comparison_tolerance}")

        return self._exportSinks(sinks_raster)


    def Clean2dLayers(self):

        if self.cross_sections_3d is not None and self.cross_sections_2d is not None:
            arcpy.Delete_management(self.cross_sections_2d)
            self.cross_sections_2d = None

        if self.detected_roads_3d is not None and self.detected_roads_2d is not None:
            arcpy.Delete_management(self.detected_roads_2d)
            self.detected_roads_2d = None

    def RecalculateRoads3D(self, roads_2d):

        self.detected_roads_2d = roads_2d
        self._export3dRoads(recalculate=True)

    def ComputeRoadGradient(self):

        arcpy.AddMessage(f"Create segments along road within distance of {self.gradient_spacing_distance}")
        self.road_gradient = os.path.join(self.workspace, GRADIENT_FC)
        arcpy.CreateFeatureclass_management(self.workspace, GRADIENT_FC,
                                            geometry_type='POLYLINE',
                                            spatial_reference=self.spatial_reference)

        for field in gradient_fields:
            arcpy.AddField_management(self.road_gradient, field.name, field.type)

        for centerline in self.center_lines:
            points = [pnt for pnt in centerline.geometry]
            segments = [Segment(origin, destination) for origin, destination in zip(points[0][:-1], points[0][1:])]

            distance = 0
            with arcpy.da.InsertCursor(self.road_gradient, ["SHAPE@"]) as iCur:
                for s in segments:
                    while distance <= s.Length:
                        cs = Segment(s.compute_new_point_3D(distance, 0),
                                     s.compute_new_point_3D(distance + self.gradient_spacing_distance, 0))
                        iCur.insertRow([cs.toArcpyPolyline(self.spatial_reference)])
                        distance += self.gradient_spacing_distance
                    distance = distance - s.Length

        self.road_gradient = self.interpolator.Interpolate3D(self.road_gradient, vertices_only=True)

        with arcpy.da.UpdateCursor(self.road_gradient, ["SHAPE@"] + [f.name for f in gradient_fields]) as cursor:
            for row in cursor:
                p1 = row[0].firstPoint
                p2 = row[0].lastPoint
                dz = abs(p1.Z - p2.Z)
                slopeDegrees = np.degrees(np.arccos(row[0].length / row[0].length3D))

                row[1] = (dz / row[0].length) * 100.0
                row[2] = slopeDegrees

                # Update the cursor with the updated list
                cursor.updateRow(row)

        return self.road_gradient

    def ComputeRoadCrossFall(self):

        if self.detected_roads_2d is None and self.detected_roads_3d is None:
            raise Exception("Cannot compute cross fall without detected roads")

        if self.detected_roads_3d is None:
            self.detected_roads_3d = self.interpolator.Interplate3D(self.detected_roads_2d, vertices_only=True)
            self._validateOutputFeatureClass(self.detected_roads_3d)

        arcpy.AddMessage("Exporting the crossfall")
        self.road_crossfall = os.path.join(self.workspace, CROSSFALL_FC)

        interpolated_roads = self.interpolator.Interpolate3D(self.detected_roads_3d, vertices_only=False)

        arcpy.CreateFeatureclass_management(self.workspace, CROSSFALL_FC,
                                            geometry_type='POLYLINE',
                                            spatial_reference=self.spatial_reference)

        for field in crossfall_fields:
            arcpy.AddField_management(self.road_crossfall, field.name, field.type)

        with arcpy.da.InsertCursor(self.road_crossfall, ["SHAPE@"] + [f.name for f in crossfall_fields]) as iCur:
            with arcpy.da.SearchCursor(interpolated_roads, ["SHAPE@"]) as sCur:
                for row in sCur:
                    p1 = row[0].firstPoint
                    p2 = row[0].centroid
                    p3 = row[0].lastPoint

                    poly1 = arcpy.Polyline(arcpy.Array([p2, p1]))
                    poly2 = arcpy.Polyline(arcpy.Array([p2, p3]))

                    dz1 = abs(p1.Z - p2.Z)
                    dz2 = abs(p3.Z - p2.Z)

                    slope_deg1 = np.degrees(np.arccos(poly1.length / poly1.length3D))
                    slope_deg2 = np.degrees(np.arccos(poly2.length / poly2.length3D))

                    slope1 = (dz1 / poly1.length) * 100.0
                    slope2 = (dz2 / poly2.length) * 100.0

                    iCur.insertRow([poly1, slope1, slope_deg1])
                    iCur.insertRow([poly2, slope2, slope_deg2])

        arcpy.Delete_management(interpolated_roads)
        return self.road_crossfall

    def ComputeRoadWidthCompliance(self):

        if not os.path.exists(self.detected_roads_3d):
            raise FileExistsError(f"File {self.detected_roads_3d} does not exists")

        with arcpy.da.UpdateCursor(self.detected_roads_3d, ["Length", "Compliance", "SHAPE@"]) as cursor:
            for row in cursor:
                try:
                    distance = row[2].length3D
                    if distance >= self.center_lines[0].width:
                        row[1] = 'TRUE'
                    else:
                        row[1] = 'FALSE'

                    # Update the length
                    row[0] = distance

                    # Update the cursor with the updated list
                    cursor.updateRow(row)

                except AttributeError:
                    cursor.deleteRow()

    def ComputeBarrierHeightCompliance(self):

        if not os.path.exists(self.detected_windrows):
            raise FileExistsError(f"File {self.detected_windrows} does not exists")

        with arcpy.da.UpdateCursor(self.detected_windrows, ["Height", "Compliance"]) as cursor:
            for row in cursor:
                if row[0] >= self.center_lines[0].height:
                    row[1] = 'TRUE'
                else:
                    row[1] = 'FALSE'

                # Update the cursor with the updated list
                cursor.updateRow(row)

# ----------------------------------------------------------------------------------------------------------------------

    @property
    def DetectedRoads(self):
        return self.detected_roads_2d

    @property
    def DetectedRoads3D(self):
        return self.detected_roads_3d

    @property
    def DetectedWindrows(self):
        return self.detected_windrows

    @property
    def DetectedWaterStreams(self):
        return self.road_streams

    @property
    def DetectedSinks(self):
        return self.road_sinks

    @property
    def RoadCrossFall(self):
        return self.road_crossfall

    @property
    def RoadGradient(self):
        return self.road_gradient

# ----------------------------------------------------------------------------------------------------------------------

    def _chopCenterLines(self, centerlines, buffer_distance=0):
        """ Chop the center lines by the given polygon within given buffer distance"""

        arcpy.AddMessage("Buffering the LAS AOI")
        buffer = os.path.join(IN_MEMORY, "buffered_aoi")
        arcpy.Buffer_analysis(in_features=self.domain,
                              out_feature_class=buffer,
                              buffer_distance_or_field=buffer_distance)

        arcpy.AddMessage("Clipping the center lines")
        clipped = os.path.join(IN_MEMORY, "clipped_centerlines")
        arcpy.Clip_analysis(in_features=[cl.geometry for cl in centerlines],
                            clip_features=buffer,
                            out_feature_class=clipped)

        output = []
        for cl in centerlines:
            single_parts = os.path.join(IN_MEMORY, "single_parts")
            arcpy.MultipartToSinglepart_management(in_features=cl.geometry, out_feature_class=single_parts)

            with arcpy.da.SearchCursor(single_parts, ['SHAPE@']) as sCur:
                for l in sCur:
                    output.append(Centerline(road_name=cl.road_name,
                                             mission_id=cl.mission_id,
                                             width=cl.width,
                                             slope=cl.slope,
                                             height=cl.height,
                                             min_area=cl.min_area,
                                             type=cl.type,
                                             gradient_spacing_distance=cl.gradient_spacing_distance,
                                             cross_sections_spacing_distance=cl.cross_sections_spacing_distance,
                                             windrow_max_distance_from_center=cl.windrow_max_distance_from_center,
                                             max_spacing_distance_filter=cl.max_spacing_distance_filter,
                                             cross_section_multiplier=cl.cross_section_multiplier,
                                             road_min_length_filter=cl.road_min_length_filter,
                                             road_max_length_filter=cl.road_max_length_filter,
                                             windrow_min_height_filter=cl.windrow_min_height_filter,
                                             windrow_max_height_filter=cl.windrow_max_height_filter,
                                             pdf_dpi=cl.pdf_dpi,
                                             geometry=l[0]))

            arcpy.Delete_management(single_parts)

        # Delete all temp layers
        [arcpy.Delete_management(lyr) for lyr in [clipped, buffer]]
        return output

    def _getDomain(self, raster):

        return arcpy.RasterDomain_3d(in_raster=raster,
                                     out_feature_class=os.path.join(IN_MEMORY, "aoi"),
                                     out_geometry_type='POLYGON')

    def _createHydrologyDomain(self, buffer_distance):
        """ Create AOI polygon for hydrology analysis"""

        arcpy.AddMessage("Creating domain for Hydrology analysis")
        buffer = os.path.join(IN_MEMORY, "hydro_buffered_center_lines")
        arcpy.Buffer_analysis(in_features=[cl.geometry for cl in self.center_lines],
                              out_feature_class=buffer,
                              buffer_distance_or_field=max(cl.width for cl in self.center_lines) * 1.5)

        clipped = os.path.join(IN_MEMORY, "hydro_clipped_center_lines")
        arcpy.Clip_analysis(in_features=buffer,
                            clip_features=self.domain,
                            out_feature_class=clipped)

        buffer = os.path.join(IN_MEMORY, "hydro_minus_buffer")
        arcpy.Buffer_analysis(in_features=clipped,
                              out_feature_class=buffer,
                              buffer_distance_or_field=buffer_distance if buffer_distance < 0 else buffer_distance * -1)

        dissolve_aoi = os.path.join(IN_MEMORY, "hydro_aoi")
        arcpy.Dissolve_management(in_features=buffer, out_feature_class=dissolve_aoi)

        return dissolve_aoi

    def _createExcel(self, slope):

        book = xlwt.Workbook()
        sheet1 = book.add_sheet("Sheet1")

        slopeCode = -1
        for i in range(int(math.floor(slope))):

            row = sheet1.row(i)
            slopeCode = i + 1
            row.write(0, "{:.02f}".format(slopeCode))
            row.write(1, slopeCode)

        row = sheet1.row(slopeCode)
        row.write(0, "{:.02f}".format(slope))
        row.write(1, slopeCode + 1)

        row = sheet1.row(slopeCode + 1)
        row.write(0, "{:.02f}".format(90.0))
        row.write(1, slopeCode + 2)

        out = os.path.join(self.workspace, "slopes.xls")
        try:
            book.save(out)
        except Exception as e:
            arcpy.AddError("Could not saved excel due to {}".format(e))

        return out, slopeCode + 1

    def _createSlopes(self):

        # create excel table with thw given slope
        excel, code = self._createExcel(max(center_line.slope for center_line in self.center_lines))
        table = arcpy.ExcelToTable_conversion(excel)

        arcpy.AddMessage("Calculating slopes")
        slopes = os.path.join("in_memory", "temp_slopes")
        arcpy.SurfaceSlope_3d(self.las, slopes, units='DEGREE', class_breaks_table=table)

        arcpy.AddMessage("Filtering slopes")
        with arcpy.da.UpdateCursor(slopes, ["SHAPE@", "SlopeCode"]) as uCur:
            for row in uCur:
                if float(row[1]) <= code:
                    uCur.deleteRow()

                elif row[0].area <= self.min_area:
                    uCur.deleteRow()

        arcpy.AddMessage("Merging slopes")
        dissolved = os.path.join(IN_MEMORY, "dissolve_slopes")
        arcpy.Dissolve_management(slopes, dissolved, multi_part=True)

        arcpy.AddMessage("Extract by AOI")
        extract = os.path.join(IN_MEMORY, "buffer_for_slopes1")
        arcpy.Clip_analysis(in_features=dissolved,
                            clip_features=self.domain,
                            out_feature_class=extract)

        buffer = os.path.join(IN_MEMORY, "buffer_for_slopes2")
        arcpy.Buffer_analysis(in_features=[cl.geometry for cl in self.center_lines],
                              out_feature_class=buffer,
                              buffer_distance_or_field=int(1.5 * max([cl.width for cl in self.center_lines])))

        self.slopes = os.path.join(self.workspace, SLOPES_FC)
        arcpy.Clip_analysis(in_features=extract,
                            clip_features=buffer,
                            out_feature_class=self.slopes)
        self._validateOutputFeatureClass(self.slopes)
        return self.slopes

    def _export3dRoads(self, recalculate=False, vertices_only=True):

        if self.detected_roads_2d is None:
            raise Exception("Cannot export to 3d if 2d does not exist")

        if recalculate:
            try:
                arcpy.Delete_management(self.detected_roads_3d)
            except:
                pass

            self.detected_roads_3d = self.interpolator.Interpolate3D(self.detected_roads_2d, vertices_only=vertices_only)
            self._validateOutputFeatureClass(self.detected_roads_3d)

        arcpy.AddMessage("Exporting detected roads in 3D")

        fields = ["SHAPE@", "IDX", "Length", "Slope", "SlopeDeg"]
        with arcpy.da.UpdateCursor(self.detected_roads_3d, fields) as cursor:
            for row in cursor:
                p1 = row[0].firstPoint
                p2 = row[0].lastPoint

                dz = abs(p1.Z - p2.Z)
                row[2] = row[0].length3D
                slopeDegrees = np.degrees(np.arccos(row[0].length / row[0].length3D))
                row[3] = (dz / row[0].length) * 100.0
                row[4] = slopeDegrees

                # Update the cursor with the updated list
                cursor.updateRow(row)

        self._validateOutputFeatureClass(self.detected_roads_3d)
        return self.detected_roads_3d

    def _exportWindrows(self):

        if not len(self.windrow1_points) or not len(self.windrow2_points):
            raise Exception("Cannot export windrows while none windrow points had been detected!")

        arcpy.AddMessage("Exporting detected windrows")
        self.detected_windrows = os.path.join(self.workspace, WINDROW_FC)
        arcpy.CreateFeatureclass_management(self.workspace, WINDROW_FC,
                                            geometry_type='POLYLINE',
                                            spatial_reference=self.spatial_reference)

        for field in windrow_fields:
            arcpy.AddField_management(self.detected_windrows, field.name, field.type)

        buffer = os.path.join(self.workspace, "buffer_for_export_windorws.shp")
        arcpy.Buffer_analysis(in_features=[cl.geometry for cl in self.center_lines],
                              out_feature_class=buffer,
                              buffer_distance_or_field=self.windrow_min_distance_from_center)

        with arcpy.da.SearchCursor(buffer, ["SHAPE@"]) as sCur:
            for feature in sCur:
                domain = feature[0]
            with arcpy.da.InsertCursor(self.detected_windrows, ["SHAPE@"] + [f.name for f in windrow_fields]) as iCur:
                for windrow_points in [self.windrow1_points, self.windrow2_points]:
                    for p1, p2 in zip(windrow_points[:-1], windrow_points[1:]):
                        p1, p2 = p1[0], p2[0]
                        poly = arcpy.Polyline(arcpy.Array([p1, p2]))
                        if poly.length > self.max_spacing_distance_filter:
                            continue
                        if poly.within(domain):
                            continue
                        iCur.insertRow((poly, max(p1.Z, p2.Z), ''))
        self._validateOutputFeatureClass(self.detected_windrows)
        return self.detected_windrows

    def _exportWaterStreams(self, streams_raster, flow_direction):

        arcpy.AddMessage("Exporting the water streams")
        self.road_streams = os.path.join(self.workspace, STREAMS)
        arcpy.sa.StreamToFeature(in_stream_raster=streams_raster,
                                 in_flow_direction_raster=flow_direction,
                                 out_polyline_features=self.road_streams,
                                 simplify=True)

        self._validateOutputFeatureClass(self.road_streams)
        return self.road_streams

    def _exportSinks(self, raster):

        arcpy.AddMessage("Exporting the sinks")
        self.road_sinks = os.path.join(self.workspace, SINKS)
        arcpy.RasterToPolygon_conversion(in_raster=raster,
                                         out_polygon_features=self.road_sinks,
                                         simplify=True,
                                         create_multipart_features="SINGLE_OUTER_PART")

        with arcpy.da.UpdateCursor(self.road_sinks, ["SHAPE@", "gridcode"]) as uCur:
            for row in uCur:

                if row[1] == 0 or row[0].area <= self.min_area:
                    uCur.deleteRow()

        self._validateOutputFeatureClass(self.road_sinks)
        return self.road_sinks

    def _saveProfile(self, polyline, bp1, bp2, rp1, rp2, idx):

        my_dpi = 96
        plt.ioff()
        plt.figure(figsize=(800 / my_dpi, 600 / my_dpi), dpi=my_dpi)

        # the profile points
        y = [p.Y for p in polyline.getPart(0)]
        x = [p.X for p in polyline.getPart(0)]
        z = [p.Z for p in polyline.getPart(0)]

        # computing the polyline azimuth
        dy = y[-1] - y[0]
        dx = x[-1] - x[0]
        az = np.arctan2(dx, dy)

        # project the point along the line
        p = [yi / np.cos(az) for yi in y]
        start = min(p)
        p = [pi - start for pi in p]
        mid = (max(p) - min(p)) / 2.0
        plt.scatter(p, z, color='b', s=7.0, label='Profile point')

        # the barrier detected points
        p = [bp1.Y / np.cos(az) - start, bp2.Y / np.cos(az) - start]
        z = [bp1.Z, bp2.Z]
        Heights = [bp1.Z - rp1.Z, bp2.Z - rp2.Z]
        plt.scatter(p, z, color='r', s=30.0, label='Detected barrier point')
        for pi, zi, height in zip(p, z, Heights):
            plt.annotate('{0:.3f} m'.format(height),
                xy=(pi, zi), xytext=(0, 15),
                textcoords='offset points', va='bottom',
                bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5))

        # the road detected points
        p = [rp1.Y / np.cos(az) - start, rp2.Y / np.cos(az) - start]
        z = [rp1.Z, rp2.Z]
        plt.scatter(p, z, color='g', s=30.0, label='Detected road point')

        width = max(p) - min(p)
        plt.annotate('{0:.3f} m'.format(width),
            xy=(mid, min(z)), xytext=(0, -15),
            textcoords='offset points', va='bottom',
            bbox=dict(boxstyle='round,pad=0.5', fc='yellow', alpha=0.5))

        plt.title("Cross section #" + str(idx).zfill(5))
        plt.grid(True)

        # Set the axis aspect to exaggerate the Z
        plt.axes().set_aspect(4, 'datalim')
        plt.legend(loc='lower left')

        outputDir = os.path.join(self.workspace, "Plots")
        if not os.path.isdir(outputDir):
            os.mkdir(outputDir)

        path = os.path.join(outputDir, "Cross section #" + str(idx).zfill(5) + '.png')
        plt.savefig(path, format='png', dpi=my_dpi)
        plt.gcf().clear()
        plt.close()

        return path

    def _detectBarrierHeightsPointsWithinDistance(self, rp1: arcpy.Point, rp2: arcpy.Point, polyline: arcpy.Polyline):

        segments = self._splitPolylineIntoSegments3D(polyline)
        mid_point = polyline.trueCentroid
        if len(segments) < 3:
            return None, None

        # Split the segments into 2 groups
        segments1 = segments[:int(len(segments) / 2)]
        segments2 = segments[int(len(segments) / 2):]

        dis1 = Segment(segments1[int(len(segments1) / 2)].MidPoint, rp1).Length
        dis2 = Segment(segments2[int(len(segments2) / 2)].MidPoint, rp1).Length

        if dis2 < dis1:
            tmp = segments1
            segments1 = segments2
            segments2 = tmp

        barrier_point1 = None
        barrier_point2 = None

        # Initialize the points max
        maxZ1 = maxZ2 = -100.0

        for s in segments1:
            if s.Z > maxZ1 and Segment(mid_point, s.Destination).Length < self.windrow_max_distance_from_center:
                barrier_point1 = s.Origin if s.Origin.Z > s.Destination.Z else s.Destination
                maxZ1 = s.Z

        for s in segments2:
            if s.Z > maxZ2 and Segment(mid_point, s.Destination).Length < self.windrow_max_distance_from_center:
                barrier_point2 = s.Origin if s.Origin.Z > s.Destination.Z else s.Destination
                maxZ2 = s.Z

        return barrier_point1, barrier_point2

    def _splitPolylineIntoSegments3D(self, polyline: arcpy.Polyline):

        segments = []
        points = polyline.getPart(0)
        for p1, p2 in zip(points[:-1], points[1:]):
            segments.append(Segment(arcpy.Point(p1.X, p1.Y, p1.Z), arcpy.Point(p2.X, p2.Y, p2.Z)))
        return segments

    @staticmethod
    def _validateOutputFeatureClass(fc):
        """ validate that the given feature class has been created """

        feature_count = int(arcpy.GetCount_management(fc).getOutput(0))
        if feature_count == 0:
            arcpy.AddError(f"{fc} has no features.")
            sys.exit(0)
        else:
            arcpy.AddMessage(f"Created {os.path.basename(fc)} within "
                             f"{os.path.basename(os.path.dirname(fc))} "
                             f"({feature_count} features)")



