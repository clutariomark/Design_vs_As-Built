import arcpy
import os
from pathlib import Path
from config_manager import Config
from utils import debug

path = Path(os.path.abspath(__file__))
resources_dir = os.path.join(path.parents[2], "resources")

class HaulRoadPDF:
    def __init__(self, aprx_path):
        self.config = Config()
        self._aprx_path = aprx_path
        self._aprx = arcpy.mp.ArcGISProject(self._aprx_path)

    @property
    def layouts(self):
        return {layout.name: layout for layout in self._aprx.listLayouts()}

    @property
    def maps(self):
        return {mp.name: mp for mp in self._aprx.listMaps()}

    @staticmethod
    def validate_data_source(source_layer):
        debug("Validating data source of {}".format(source_layer))
        if hasattr(source_layer, "dataSource"):
            debug("Source is a layer object")
            source_layer = source_layer.dataSource
        else:
            debug("Source is a real path")
            source_layer = str(source_layer)
        return source_layer
    
    def set_text(self, value, layout_name, element_name="Title", visible=True):

        try:
            layout = self.layouts[layout_name]
        except KeyError:
            raise ValueError("Layout {} was not found".format(layout_name))
        try:
            title_elm = [elm for elm in layout.listElements("TEXT_ELEMENT") if elm.name == element_name][0]
            title_elm.visible = visible
        except IndexError as e:
            raise ValueError("{} map element was not found".format(element_name))

        title_elm.text = value

    def zoom_to_layer(self, layout_name, layer_name):
        try:
            layout = self.layouts[layout_name]
        except KeyError:
            raise ValueError("Layout was not found")
        mf = layout.listElements("MAPFRAME_ELEMENT")[0]
        mapp = mf.map
        try:
            layer = [l for l in mapp.listLayers(layer_name) if not l.isBroken][0]
            extent = mf.getLayerExtent(layer)
            mf.camera.setExtent(extent)
        except IndexError:
            arcpy.AddError("Layer {} was not found".format(layer_name))

        allowed_scales = [value for value in range(500, 50000, 500)]
        debug("original scale: {}".format(mf.camera.scale))
        for scale in allowed_scales:
            if mf.camera.scale < scale:
                debug("rounding scale to {}".format(scale))
                mf.camera.scale = scale
                break

    def add_hillshade(self, map_name, raster_source, symbology_set):
        raster_source = self.validate_data_source(raster_source)

        try:
            the_map = self.maps[map_name]
        except KeyError:
            raise ValueError("Map named '{}' was not found in this project".format(map_name))

        analysis_layer = [lyr for lyr in the_map.listLayers() if lyr.isFeatureLayer][-1]

        hillshade_layer_file = arcpy.mp.LayerFile(
            os.path.join(resources_dir, "symbol", symbology_set, "hillshade.lyrx"))
        hillshade_layer = the_map.insertLayer(analysis_layer,  hillshade_layer_file, "AFTER")
        debug("Input hillshade source: {}".format(raster_source))
        debug("Old hillshade source: {}".format(hillshade_layer.dataSource))
        new_connection_properties = {
            'dataset': os.path.basename(raster_source),
            'connection_info': {
                'database': os.path.dirname(raster_source)
            }
        }
        hillshade_layer.updateConnectionProperties(hillshade_layer.connectionProperties, new_connection_properties)
        debug("New hillshade source: {}".format(hillshade_layer.dataSource))
        hillshade_layer.name = "Hillshade"

        return hillshade_layer.name

    def list_layers(self, layer_name=None, map_name="Map"):
        try:
            the_map = self.maps[map_name]
        except KeyError:
            raise ValueError("Map named '{}' was not found in this project".format(map_name))
        return the_map.listLayers(wildcard=layer_name)

    def change_data_source(self, map_name, layer_name, new_data_source):
        new_data_source = self.validate_data_source(new_data_source)
        try:
            layer = self.list_layers(layer_name, map_name)[0]
        except IndexError:
            raise ValueError("Layer {} was not found in {}".format(layer_name, map_name))

        new_connection_properties = {
            'dataset': os.path.basename(new_data_source),
            'connection_info': {
                'database': os.path.dirname(new_data_source)
            }
        }

        layer.updateConnectionProperties(layer.connectionProperties, new_connection_properties)
        debug("Data source is now: {}".format(layer.dataSource))
        if os.path.normpath(layer.dataSource) != os.path.normpath(new_data_source):
            arcpy.AddError("Failed to change data source from {} to {}".format(layer.dataSource, new_data_source))

    def apply_symbology(self, map_name, layer_name, lyrx, **kwargs):
        try:
            layer = self.list_layers(layer_name, map_name)[0]
        except IndexError:
            raise ValueError("Layer {} was not found in {}".format(layer_name, map_name))
        arcpy.ApplySymbologyFromLayer_management(layer, lyrx, **kwargs)

        if "road width" in layer.name.lower():
            # see https://community.esri.com/thread/216962-apply-lyrx-file-symbology-via-python
            debug("using BUG-000108497 workaround")

            try:
                the_map = self.maps[map_name]
            except KeyError:
                raise ValueError("Map named '{}' was not found in this project".format(map_name))
            
            added_layer = the_map.addLayer(arcpy.mp.LayerFile(lyrx))[0]

            dataset, database = os.path.basename(layer.dataSource), os.path.dirname(layer.dataSource)
            new_connection_properties = {
                'dataset': dataset,
                'connection_info': {
                    'database': database
                }
            }
            added_layer.updateConnectionProperties(added_layer.connectionProperties, new_connection_properties)
            the_map.removeLayer(layer)

            site = os.path.split(os.path.dirname(lyrx))[-1]
            for layout in self._aprx.listLayouts("Road Width*"):
                debug(f"updating graphic legend item Legend{site} position in {layout.name}")
                legend = layout.listElements("GRAPHIC_ELEMENT", f"Legend{site}")[0]

                position_y = 13.6204
                try:
                    debug("getting road width legend position from site config")
                    site_config = Config(site)
                    position_y = site_config.get("layout_settings", "width_legend_position_y")
                except:
                    debug("failed, using default value")
                finally:
                    debug(f"position Y: {position_y}")
                legend.elementPositionY = position_y

    def is_layout_fit(self, layout_name, map_name, layer_name):

        try:
            layout = self.layouts[layout_name]
        except KeyError:
            raise ValueError("Layout was not found")

        mf = layout.listElements("MAPFRAME_ELEMENT")[0]
        layer = self.list_layers(layer_name, map_name)[0]

        extent = mf.getLayerExtent(layer)
        is_layer_landscape = extent.width > extent.height
        is_layout_landscape = layout.pageWidth > layout.pageHeight
        if is_layer_landscape == is_layout_landscape:
            return True
        else:
            return False

    def save(self):
        self._aprx.save()

    def saveACopy(self, copy_path):
        self._aprx.saveACopy(copy_path)

    def export(self, output, layout_name="Layout", dpi=300):
        try:
            layout = self.layouts[layout_name]
        except KeyError:
            raise ValueError("Layout was not found")
        debug(f"exporting with resolution of {dpi} DPI")
        try:
            layout.exportToPDF(output, resolution=dpi)
        except OSError as e:
            arcpy.AddError("Failed to export PDF. Make sure you have permissions and the file is not already opened.")

    def remove_base_map(self, map_name):
        try:
            the_map = self.maps[map_name]
        except KeyError:
            raise ValueError("Map named '{}' was not found in this project".format(map_name))
        base_layer = [layer for layer in the_map.listLayers(wildcard='*') if layer.isBasemapLayer]
        if len(base_layer):
            base_layer = base_layer[0]
            base_layer.visible = False
            self._aprx.save()
