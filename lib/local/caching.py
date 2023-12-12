import os
import pickle
import gzip
import json
import sqlite3
from pathlib import Path
from io import BytesIO
from typing import Callable, TypeVar
from .constants import EXECUTION_DIR, WORKSPACE_ROOT

############################## pickling ##############################

CACHE = f'{EXECUTION_DIR}/cache'

_force_regenerate = False
def set_force_regenerate(b):
    global _force_regenerate
    _force_regenerate = b

def _get_paths(fname: str, alt=None):
    if isinstance(alt, str) and alt[-1] == '/': alt = alt[:-1]
    cache = CACHE if alt is None else f'{alt}/cache'
    fpath = f'{cache}/{fname}'
    return fpath, cache

def _ext_to_fpaths(fpath: str, compression=False):
    EXT = '.pkl.gz' if compression else '.pkl'
    fpath = fpath.replace(EXT, '')
    fpath += EXT
    fpath_str = fpath.replace(str(WORKSPACE_ROOT), "{WORKSPACE}") # for logging
    return fpath, fpath_str

def save_exists(name: str, alt_workspace=None):
    fpath_no_ext, cache = _get_paths(name, alt_workspace)
    fpath_comp, fpath_str_comp = _ext_to_fpaths(fpath_no_ext, compression=True)
    fpath, fpath_str = _ext_to_fpaths(fpath_no_ext, compression=True)
    return os.path.exists(fpath) or os.path.exists(fpath_comp)

def save(name, x, alt_workspace=None, compression_level=1, silent=False):
    fpath_no_ext, cache = _get_paths(name, alt_workspace)
    if not os.path.isdir(cache): os.system(f'mkdir {cache}')

    fpath, fpath_str = _ext_to_fpaths(fpath_no_ext, compression=compression_level>0)

    cmsg = 'compressing & ' if compression_level > 0 else ''
    if not silent: print(f'{cmsg}caching data to [{fpath_str}]')

    if compression_level == 0:
        with open(fpath, 'wb') as f:
            pickle.dump(x, f, protocol=pickle.HIGHEST_PROTOCOL)
    else:
        with gzip.open(fpath, "wb", compresslevel=compression_level) as f:
            pickle.dump(x, f, protocol=pickle.HIGHEST_PROTOCOL)


def load(name: str, alt_workspace=None, silent=False):
    fpath_no_ext, cache = _get_paths(name, alt_workspace)

    for c in [False, True]:
        fpath, fpath_str = _ext_to_fpaths(fpath_no_ext, compression=c)
        if not os.path.isfile(fpath): continue

        dcomp_msg = '& decompressing ' if c else ''
        if not silent: print(f'recovering {dcomp_msg}cached data from [{fpath_str}]')
        with gzip.open(fpath, "rb") if c else open(fpath, "rb") as f:
            return pickle.load(f)

    fpath, fpath_str = _ext_to_fpaths(fpath_no_ext)
    raise FileNotFoundError(f"{fpath_str} doesn't exist, nor can a compressed cache be found")

def cache(fname, regenerate, force_regenerate=None, compression_level=1):
    if force_regenerate is None: force_regenerate = _force_regenerate
    fpath_no_ext, cache = _get_paths(fname)
    fpath, fpath_str = _ext_to_fpaths(fpath_no_ext, compression=compression_level>0)

    if not force_regenerate and os.path.isfile(fpath):
        return load(fname)
    else:
        x = regenerate()
        save(fname, x, compression_level=compression_level)
        return x

############################## fn decorator ##############################

T = TypeVar('T')
def cache_fn_result(loader: Callable[..., T]) -> Callable[[], T]:
    data = None
    def getter(*args, **kargs):
        nonlocal data
        if data is None:
            data = loader(*args, **kargs)
        return data
    return getter

# #####################################################################################

class DictCache:
    EXT = ".db"
    def __init__(self, name: str, save_folder: Path|None=None, compression: int=9) -> None:
        if save_folder is None:
            save_folder = WORKSPACE_ROOT.joinpath(f"data/cache")
            if not save_folder.exists(): os.makedirs(save_folder, exist_ok=True)
        if not name.endswith(self.EXT): name += self.EXT
        # Connect to the SQLite database (or create it if it doesn't exist)
        self.conn = sqlite3.connect(save_folder.joinpath(name))

        # Create a table to store the compressed, cached JSON data
        self.conn.execute('''CREATE TABLE IF NOT EXISTS json_cache
                        (id TEXT PRIMARY KEY, data BLOB)''')
        
        self.compression = compression

    def save(self):
        # Commit the changes to the database
        self.conn.commit()

    def close(self):
        self.conn.commit()
        self.conn.close()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        self.save()
        self.close()

    def _get_iterator(self, fields: str):
        return self.conn.execute(f"SELECT {fields} FROM json_cache")

    def keys(self):
        for row in self._get_iterator("id"): yield row[0]

    def values(self):
        for row in self._get_iterator("data"): yield self._decompress(row[0])

    def items(self):
        for k, v in self._get_iterator("id, data"):
            yield k, self._decompress(v)

    def __iter__(self):
        return self.keys()

    def __contains__(self, key: str):
        return self.get(key) is not None

    # Define a function to cache JSON data (compressed with gzip)
    def __setitem__(self, key: str, data: dict):
        # Serialize the JSON data to a string
        json_data = json.dumps(data)
        
        # Compress the JSON data using gzip
        gzip_buffer = BytesIO()
        with gzip.GzipFile(mode='wb', fileobj=gzip_buffer, compresslevel=self.compression) as f:
            f.write(json_data.encode('utf-8'))
        compressed_data = gzip_buffer.getvalue()
        
        # Insert or replace the compressed data in the database
        self.conn.execute("INSERT OR REPLACE INTO json_cache (id, data) VALUES (?, ?)", (key, compressed_data))

    def _decompress(self, compressed_data):
        gzip_buffer = BytesIO(compressed_data)
        with gzip.GzipFile(mode='rb', fileobj=gzip_buffer) as f:
            return json.loads(f.read().decode('utf-8'))

    def get(self, key: str, default: dict|None=None) -> dict|None:
        # Query the database for the compressed JSON data
        cursor = self.conn.execute("SELECT data FROM json_cache WHERE id=?", (key,))
        
        # Get the first row of the result (or None if no rows are returned)
        row = cursor.fetchone()
        
        if row is not None:
            # Decompress the compressed JSON data
            compressed_data = row[0]
            # Deserialize the JSON data and return it
            return self._decompress(compressed_data)
        else:
            return default

    # Define a function to retrieve cached JSON data
    def __getitem__(self, key: str) -> dict:
        v = self.get(key)
        if v is None: raise KeyError(f"[{key}] not found")
        return v
