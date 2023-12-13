# This file is part of LCAstar.
# 
# FabFos is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# FabFos is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of 
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with LCAstar. If not, see <https://www.gnu.org/licenses/>.

# copyright 2023 Niels W. Hanson, Kishori M. Konwar, Tony X. Liu, 
# Steven J. Hallam

import os, sys
import time
from datetime import datetime as dt
from pathlib import Path
import re
from typing import Iterable, Sequence

USER = "hallamlab" # github id
MODULE_ROOT = Path("/".join(os.path.realpath(__file__).split('/')[:-1]))
NAME = MODULE_ROOT.name.lower()
# ENTRY_POINTS = [NAME]
ENTRY_POINTS = []

def _get_version() -> str:
    with open(MODULE_ROOT.joinpath("version.txt")) as v:
        return v.readline()
VERSION = _get_version()

def regex(r, s):
    for m in re.finditer(r, s):
        yield s[m.start():m.end()]

def Batchify(lst: Sequence, size: int):
    for i in range(0, len(lst), size):
        yield lst[i:i + size]
        
class StdTime:
    FORMAT = '%Y-%m-%d_%H-%M-%S'

    @classmethod
    def Timestamp(cls, timestamp: dt|None = None):
        ts = dt.now() if timestamp is None else timestamp
        return f"{ts.strftime(StdTime.FORMAT)}"
    
    @classmethod
    def Parse(cls, timestamp: str|int):
        if isinstance(timestamp, str):
            return dt.strptime(timestamp, StdTime.FORMAT)
        else:
            return dt.fromtimestamp(timestamp/1000)
    
    @classmethod
    def CurrentTimeMillis(cls):
        return round(time.time() * 1000)
    
if __name__ == "__main__":
    sys.path = [str(p) for p in set([
        MODULE_ROOT.parents[1]
    ]+sys.path)]
    from setup import SHORT_SUMMARY
    if len(sys.argv)>1:
        k = sys.argv[1]
        meta = dict(
            USER = USER,
            NAME = NAME,
            ENTRY_POINTS = ENTRY_POINTS,
            VERSION = VERSION,
            SHORT_SUMMARY = SHORT_SUMMARY,
            MODULE_ROOT = MODULE_ROOT,
        )
        print(meta.get(k, ""))
