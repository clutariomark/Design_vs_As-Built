from typing import Optional
from fastapi import FastAPI
from subprocess import *
import os

app = FastAPI()

@app.get("/dvab/")
def calc(las: Optional[str] = None, xml: Optional[str] = None ,c: Optional[str] = None, m: Optional[str] = None, a: Optional[str] = None, c_fill: Optional[str] = None):
    # subprocess.call(['C:\\Users\\kjohn\\Desktop\\Design_vs_AsBuilt_Beta\\test.bat {} {} {} {} {}'.format(las, xml, c, m, a)])
    # output = 'C:\\Users\\kjohn\\Desktop\\Design_vs_AsBuilt_Beta\\test.bat {}'.format(las)
    # subprocess.call([output])
    os.system(f"command.bat {las} {xml} {c} {m} {a} {c_fill}")
    return {'las': las,
            'xml': xml,
            'construction': c,
            'method': m,
            'area': a,
            'cell_fill': c_fill}

# "C:\Program Files\ArcGIS\Pro\bin\Python\Scripts\propy.bat" "C:\Users\kjohn\Documents\GitHub\Design_vs_As-Built\design_vs_as_built.py" 
# "C:\Users\kjohn\Desktop\Design_vs_AsBuilt_Beta\OGI_Infotec\dekigata.las" "C:\Users\kjohn\Desktop\Design_vs_AsBuilt_Beta\OGI_Infotec\design.xml" riverbed fill slope MINIMUM


@app.get("/")
def hello():
    return {'hello': 'world'}