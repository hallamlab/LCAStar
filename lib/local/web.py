from pyexpat import ExpatError
import time
import requests
import xmltodict    
from typing import Any
import urllib.parse
from .constants import WORKSPACE_ROOT
from .caching import DictCache

# #####################################################################################################################
# utils

def chain_get(d: dict, path: list|str) -> None | list[Any]:
    if isinstance(path, str): path = path.split(', ')
    todo:list[tuple[dict, int]] = [(d, 0)]
    results = []
    while len(todo) > 0:
        curr, depth = todo.pop()
        if curr is None: continue
        if len(curr) == 0: continue
        if depth >= len(path):
            if isinstance(curr, list):
                results += curr
            else: 
                results.append(curr)
            continue
        
        if isinstance(curr, list):
            todo += [(x, depth) for x in curr]
        elif not isinstance(curr, dict):
            assert False, (curr, path[depth-1])
        else:
            todo.append((curr.get(path[depth], {}), depth+1))
    return results if len(results)>0 else None

# #####################################################################################################################
# ncbi

# https://www.ncbi.nlm.nih.gov/account/settings/
# https://www.ncbi.nlm.nih.gov/books/NBK25497/#_chapter2_Usage_Guidelines_and_Requiremen_
with open(WORKSPACE_ROOT.joinpath("secrets/ncbi_apikey")) as f:
    API_KEY = f.readline().replace("\n", "").strip()

def ncbi_get(action: str, db: str|None=None, params: list[tuple[str, str]]=list(), retry=False) -> tuple[Any, dict]:
    base_url = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
    for k, v in params:
        assert isinstance(k, str), (type(k), k)
        assert isinstance(v, str), (type(v), v)
    base = [
        f"api_key={API_KEY}",
    ]
    if db is not None:
        base.append(f"db={db}")
    
    url=f"{base_url}/{action}.fcgi?"+"&".join(base+[
        f"{k}={urllib.parse.quote(v)}" for k, v in params
    ])
    with DictCache("ncbi_requests") as request_cache:
        _s = url.index("api_key")
        _e = url.index("&", _s)+1
        ckey = url[:_s]+url[_e:]
        ckey = ckey.replace(base_url, "")
        _cached_r = request_cache.get(ckey)
        if retry or _cached_r is None:
            time.sleep(0.1)
            r = requests.get(url)
            try:
                d: dict = xmltodict.parse(r.text)
            except ExpatError as e:
                print(f"{e}")
                print(f"{r.text}")
                return "xml format", dict(html_text=r.text)
            if r.status_code != 200: return r.status_code, d
            request_cache[ckey] = dict(status_code = r.status_code, data=d)
        else:
            d: dict = _cached_r["data"]
        return 200, d

def ncbi_search(query: str, db: str, response_type: str="esummary", 
                search_params: list[tuple[str, str]]=list(), 
                response_params: list[tuple[str, str]]=list(), 
                silent: bool=False, retry: bool=False) -> tuple[str, dict]|tuple[None, list]:
    # db = "biosample", https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly

    _esearch = "esearch"
    code, d = ncbi_get(_esearch, db, params=[("term", query)]+search_params, retry=retry)
    if code != 200: return _esearch, d

    rdata = {} if not isinstance(d, dict) else d.get('eSearchResult', {})
    _ids = rdata.get('IdList', {})
    ids = None if _ids is None else _ids.get('Id')

    if ids is None: return "not found", d
    idcount = int(rdata.get('Count', 0))
    if idcount == 0: return "not found", d
    elif idcount == 1:
        assert isinstance(ids, str), (type(ids), ids)
        ids = [ids]
    ids = [id for id in ids if isinstance(id, str)]

    responses: list = []
    for i, id in enumerate(ids):
        assert id != ""
        if not silent: print(f"\rfetching result {i+1} of {len(ids)}", end="")
        status_code, data = ncbi_get(response_type, db, params=[("id", id)]+response_params, retry=retry)
        if status_code != 200: continue
        responses.append(data)
    if not silent and len(ids)>0: print()

    # if len(responses) == 0: 
    return None, responses

def ncbi_link(query: str, dba: str, dbb: str, response_type: str="esummary", 
                search_params: list[tuple[str, str]]=list(), 
                response_params: list[tuple[str, str]]=list(), 
                silent: bool=False, retry: bool=False):
    # example
    # dba = "biosample", https://www.ncbi.nlm.nih.gov/books/NBK25497/table/chapter2.T._entrez_unique_identifiers_ui/?report=objectonly
    # dbb = "bioproject"

    _esearch = "esearch"
    code, d = ncbi_get(_esearch, dba, params=[("term", query)]+search_params, retry=retry)
    if code != 200: return _esearch, d

    rdata = {} if not isinstance(d, dict) else d.get('eSearchResult', {})
    _ids = rdata.get('IdList', {})
    ids = None if _ids is None else _ids.get('Id')

    if ids is None: return "not found", d
    idcount = int(rdata.get('Count', 0))
    if idcount == 0: return "not found", d
    elif idcount == 1:
        assert isinstance(ids, str), (type(ids), ids)
        ids = [ids]
    ids = [id for id in ids if isinstance(id, str)]

    LK = f"{dba}_{dbb}"
    linked_ids = set()
    for i, id in enumerate(ids):
        assert id != ""
        if not silent: print(f"\rgetting link {i+1} of {len(ids)}", end="")
        status_code, data = ncbi_get("elink", db=dbb, params=[("id", id), ("dbfrom", dba), ("linkname", LK)], retry=retry)
        if status_code != 200: continue
        if not isinstance(data, dict): continue
        dbb_ids = chain_get(data, "eLinkResult, LinkSet, LinkSetDb, Link, Id")
        if dbb_ids is None: continue
        linked_ids.update(dbb_ids)

    if len(linked_ids) == 0: return "not found", d

    if response_type == "id":
        return None, linked_ids

    linked_ids = sorted(linked_ids)
    responses: list = []
    if not silent: print()
    for i, id in enumerate(linked_ids):
        assert id != ""
        if not silent: print(f"\rfetching result {i+1} of {len(linked_ids)}", end="")
        status_code, data = ncbi_get(response_type, dbb, params=[("id", id)]+response_params, retry=retry)
        if status_code != 200: continue
        responses.append(data)
    if not silent: print()
    return None, responses
