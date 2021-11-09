import datetime

import arcpy
from config_manager import Config
import os
import winreg


def file_from_onedrive(relative_path: str = "") -> str:
    """ Returns the path of a file in OneDrive local folder.

    :param relative_path: the file path relative to "Airo_GIS - Doc" folder (if empty, root folder is returned)
    :return: file path
    """
    try:
        with winreg.OpenKey(winreg.HKEY_CURRENT_USER,
                            r"Software\SyncEngines\Providers\OneDrive\a5b612b28ac94ab28b2d625adefbdc03+1") as regkey:
            for value in winreg.EnumValue(regkey, 0):
                if os.path.isdir(value):
                    onedrive = value
                    if not relative_path:
                        return onedrive
                    path = os.path.join(onedrive, relative_path)
                    if os.path.isfile(path) or os.path.isdir(path):
                        return path
                    raise FileNotFoundError("OneDrive folder was found, but the file does not exist")
        raise IOError(r"Failed to find ScanX GIS OneDrive directory")
    except (IOError, FileNotFoundError) as e:
        arcpy.AddError(e)


def extract_date(layer_name):
    if os.path.isfile(layer_name):
        layer_name = os.path.basename(layer_name).split(".")[0]
    if os.sep in layer_name:
        # in case the layer is in a group, the name is "[group]\[layer]", and only the last part is taken
        layer_name = layer_name.split(os.sep)[-1]
    if "." in layer_name:
        layer_name = layer_name.split(".")[0]
    separator = "_"
    for part in layer_name.split(separator):
        try:
            _ = int(part)  # only check if part is a sequence of numbers, no need for the actual number
            return datetime.datetime.strptime(part, "%Y%m%d")
        except ValueError:
            continue


def get_workspace(name: str = None, display: bool = False):
    config = Config()
    try:
        root = config.get("workspace", "workspace_root")
    except Exception:
        root = None
    if not root or not os.path.isdir(root):
        user_home_dir = os.path.expanduser("~")
        root = os.path.join(user_home_dir, "ScanX", "HaulRoad")
    if name is not None:
        ws = os.path.join(root, name)
    else:
        ws = root
    if not os.path.isdir(ws):
        os.makedirs(ws)
    if display:
        arcpy.AddMessage(f"Workspace: {ws}")
    else:
        _debug(f"Workspace: {ws}")
    return ws


def _debug(text):
    if os.getenv("DEV"):
        arcpy.AddMessage(f"DEBUG: {text} \n {var}")
        arcpy.AddMessage("*"*100)

def debug(text, var):
    arcpy.AddMessage("*"*100)
    arcpy.AddMessage(f"DEBUG: {str(text)} \n {var}")
    arcpy.AddMessage("*"*100)