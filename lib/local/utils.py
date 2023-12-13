from typing import Any
import pandas as pd
import re

def pd_set_type(cols: list|str, t, df: pd.DataFrame):
    if isinstance(cols, str): cols = cols.split(', ')
    for col in cols: df[col] = df[col].astype(t)

def regex(r, s):
    for m in re.finditer(r, s):
        yield s[m.start():m.end()]
