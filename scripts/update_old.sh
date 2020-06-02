#!/bin/bash

# This file is part of ENSO.
# Copyright (C) 2020 Fabian Bohle
#
# ENSO is free software: you can redistribute it and/or modify it under
# the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ENSO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with ENSO. If not, see <https://www.gnu.org/licenses/>.


# The script is intended to update enso.json files prepared by enso.py below version 2.0
cp enso.json updated_enso.json
mv enso.json old_enso.json

sed -i 's/"nstruc":/"nconf":/g' updated_enso.json

# if single-point in gas phase (eg. no solvent, or additive solvation + gas phase single-point)
gas='False'
if [ "$gas" == 'True' ]; then
    sed -i's/"sp_part3": "calculated"/"sp_part3_gas": "calculated"/g' updated_enso.json
    sed -i 's/"sp_part3": "not_calculated"/"sp_part3_gas": "not_calculated"/g'  updated_enso.json
    sed -i 's/"energy_sp_part3":/"energy_sp_part3_gas":/g' updated_enso.json
fi

# if single-point in solution phase (NO additive solvation!!! )
solvent='True'
if [ "$solvent" == 'True' ]; then
    sed -i 's/"sp_part3": "calculated"/"sp_part3_solv": "calculated"/g' updated_enso.json
    sed -i 's/"sp_part3": "not_calculated"/"sp_part3_solv": "not_calculated"/g'  updated_enso.json
    sed -i 's/"energy_sp_part3":/"energy_sp_part3_solv":/g' updated_enso.json
fi

# if the rrho contribution was calculated with any GFNn-xTB method
rrho='xtb'
if [ "$rrho" == 'xtb' ]; then
    sed -i 's/"rrho": "calculated"/"rrho_xtb": "calculated"/g' updated_enso.json
    sed -i 's/"rrho": "not_calculated"/"rrho_xtb": "not_calculated"/g'  updated_enso.json
    sed -i 's/"energy_rrho":/"energy_rrho_xtb":/g' updated_enso.json
fi
cp updated_enso.json enso.json
