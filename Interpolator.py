import arcpy
import os
import sys
import time

LAS    = 'LasDataset'
RASTER = 'RasterDataset'

class Interpolator(object):

    def __init__(self, surface):
        
        if self._surfaceSupported(surface=surface):
            self.surface = surface
        else:
            raise ValueError(f"Surface {surface} is not supported")


    def Interpolate3D(self, inFC, vertices_only=False):

        outFC = os.path.splitext(inFC)[0] + "_3D.shp"
        surface_type = arcpy.Describe(self.surface).dataType

        arcpy.AddMessage("Interpolating into 3D")
        if surface_type == LAS:
            if vertices_only:
                arcpy.ddd.InterpolateShape(self.surface, inFC, outFC, method="BILINEAR", vertices_only=vertices_only)
            else:
                arcpy.ddd.InterpolateShape(self.surface, inFC, outFC, method="BILINEAR")

        elif surface_type == RASTER:
            arcpy.ddd.InterpolateShape(self.surface, inFC, outFC, method="BILINEAR", vertices_only=vertices_only)
        
        else:
            raise ValueError(f"Surface {self.surface} is not supported for interpolation")

        return outFC

    def _surfaceSupported(self, surface):
        
        surface_type = arcpy.Describe(surface).dataType
        
        if surface_type == LAS or surface_type == RASTER:
            return True
            
        return False

if __name__ == "__main__":
    arcpy.env.overwriteOutput = True

    inFC = "E:\GoogleDrive_CS\ScanX\RoadAnalytics\Development\HaulRoad\input\CenterLine.shp"
    surface = "E:\GoogleDrive_CS\ScanX\RoadAnalytics\Development\HaulRoad\input\LIDAR_HR.las"
    vertices_only = False

    interpolator = Interpolator(surface)

