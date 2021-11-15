from typing import Optional

from fastapi import FastAPI

app = FastAPI()


@app.get("/")
def read_root():
    return {"Hello": "World"}


# @app.get("/items/{item_id}")
# def read_item(item_id: int, q: Optional[str] = None):
#     return {"item_id": item_id, "q": q}


# e.g. http://127.0.0.1:8000/items/5?c=road&m=cut&a=ground
@app.get("/items/{item_id}")
def read_item(item_id: int, c: Optional[str] = None, m: Optional[str] = None, a: Optional[str] = None):
    return {"item_id": item_id, "construction": c, "method": m, "area": a}

