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

try:
    from .algorithm import LcaStar, Lineage, ResultNode, InitNcbi
except ImportError:
    print("Warning: could not import package, ignore if building package from git repo")
