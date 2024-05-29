from enum import Enum
from pathlib import Path
import os

EXECUTION_DIR = Path(os.getcwd())
WORKSPACE_ROOT = Path("/".join(os.path.realpath(__file__).split('/')[:-3]))
