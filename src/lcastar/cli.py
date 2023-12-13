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

import json
import os, sys
from pathlib import Path
import argparse
import inspect
import importlib

from .utils import NAME, USER, VERSION, ENTRY_POINTS, MODULE_ROOT, StdTime

CLI_ENTRY = ENTRY_POINTS[0]
    
class ArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        self.print_help(sys.stderr)
        self.exit(2, '\n%s: error: %s\n' % (self.prog, message))


class CommandLineInterface:
    def _get_fn_name(self):
        return inspect.stack()[1][3]

    def estimate(self, raw_args):
        parser = ArgumentParser(
            prog = f'{CLI_ENTRY} {self._get_fn_name()}',
        )

        print("CLI not implemented yet, please use as python module")

    def help(self, args=None):
        help = [
            f"{NAME} v{VERSION}",
            f"https://github.com/{USER}/{NAME}",
            f"",
            f"Syntax: {CLI_ENTRY} COMMAND [OPTIONS]",
            f"",
            f"Where COMMAND is one of:",
        ]+[f"- {k}" for k in COMMANDS]+[
            f"",
            f"for additional help, use:",
            f"{CLI_ENTRY} COMMAND -h/--help",
        ]
        help = "\n".join(help)
        print(help)
COMMANDS = {k:v for k, v in CommandLineInterface.__dict__.items() if k[0]!="_"}

def main():
    cli = CommandLineInterface()
    if len(sys.argv) <= 1:
        cli.help()
        return

    COMMANDS.get(# calls command function with args
        sys.argv[1], 
        CommandLineInterface.help # default
    )(cli, sys.argv[2:]) # cli is instance of "self"

if __name__ == "__main__":
    main()
