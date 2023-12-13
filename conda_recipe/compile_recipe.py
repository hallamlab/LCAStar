import os, sys
import stat
from pathlib import Path
import yaml

HERE = Path(os.path.realpath(__file__)).parent
sys.path = list(set([
    str(HERE.joinpath("../").absolute())
]+sys.path))

# import constants from setup.py
from setup import USER, NAME, VERSION, ENTRY_POINTS, SHORT_SUMMARY

# ======================================================
# parse dependencies
with open(HERE.joinpath(f"../envs/base.yml")) as y:
    raw_deps = yaml.safe_load(y)
def _parse_deps(level: list, compiled: str, depth: int):
    tabs_space = "  "*depth
    for item in level:
        # conda recipes can't have pip
        # instead, a few can be added into the template, but these will not be tracked!
        if not isinstance(item, str) or item in {"pip"}: continue
        if isinstance(item, str):
            compiled += f"{tabs_space}- {item}\n"
        else:
            k, v = list(item.items())[0]
            compiled += f"{tabs_space}- {k}:\n"
            compiled = _parse_deps(v, compiled, depth+1)
    compiled = compiled[:-1] # remove trailing \n
    return compiled
reqs = _parse_deps(raw_deps["dependencies"], "", 2)
python_dep = [d for d in raw_deps["dependencies"] if isinstance(d, str) and d.startswith("python=")]
if len(python_dep) < 1:
    python_dep = ["python=3.11"]
python_ver = _parse_deps(python_dep, "", 2)

# ======================================================
# entry points

entry_points = ""
for e in ENTRY_POINTS:
    tabs_space = "  "*2
    entry_points += f"{tabs_space}- {e}\n"
entry_points = entry_points[:-1] # remove trailing \n


# ======================================================
# path to tar archive of source code

dist_path = Path(os.path.abspath(HERE.joinpath("../dist")))
assert dist_path.exists(), f"did you forget to build the pip package first?"
tar_path = [dist_path.joinpath(f) for f in os.listdir(dist_path) if VERSION in f and ".tar.gz" in f][0]


# ======================================================
# generate recipe files

with open(HERE.joinpath("meta_template.yaml")) as f:
    template = "".join(f.readlines())
meta_values = {
    "USER": USER,
    "NAME": NAME,
    "SHORT_SUMMARY": SHORT_SUMMARY,
    "VERSION": VERSION,
    "ENTRY": entry_points,
    "REQUIREMENTS": reqs,
    "PYTHON": python_ver,
    "TAR": f"file://{tar_path}",
}
for k, v in meta_values.items():
    template = template.replace(f"<{k}>", v)
with open(HERE.joinpath("meta.yaml"), "w") as f:
    f.write(template)

build_file = HERE.joinpath("call_build.sh")
with open(build_file, "w") as f:
    channels = " ".join(f"-c {ch}" for ch in raw_deps["channels"])
    _here = 'HERE=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )'
    f.write(f"""\
        {_here}
        conda mambabuild {channels} --output-folder $HERE/../conda_build $HERE/
    """.replace("    ", ""))
st = os.stat(build_file)
os.chmod(build_file, st.st_mode | stat.S_IEXEC)
