import os
from pathlib import Path
from configparser import ConfigParser

class Config:
    def __init__(self, site=None):
        
        self._config = ConfigParser()
        path = Path(os.path.abspath(__file__))
        self.resources_dir = os.path.join(path.parents[2], "resources")

        if site is None:
            self._config_path = os.path.join(self.resources_dir, "HR.ini")
        else:
            self._config_path = os.path.join(self.resources_dir, "symbol", site, "site_config.ini")

        self._config.read(self._config_path)

    def get(self, section, option):
        return self._config.get(section, option)
    
    def __repr__(self):
        return "<Config: {}>".format(self._config_path)