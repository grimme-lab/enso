#!/usr/bin/env python3

# This file is part of ENSO.
# Copyright (C) 2019 Fabian Bohle, Karola Schmitz
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

coding = "ISO-8859-1"
try:
    import argparse
except ImportError:
    raise ImportError(
        "ENSO requires the module argparse. Please install the module argparse."
    )
try:
    import sys
except ImportError:
    raise ImportError("ENSO requires the module sys. Please install the module sys.")
try:
    import os
except ImportError:
    raise ImportError("ENSO requires the module os. Please install the module sys.")
from os.path import expanduser

try:
    import shutil
except ImportError:
    raise ImportError(
        "ENSO requires the module shutil. Please install the module shutil."
    )
try:
    import subprocess
except ImportError:
    raise ImportError(
        "ENSO requires the module subprocess. Please install the module subprocess."
    )
try:
    from multiprocessing import JoinableQueue as Queue
except ImportError:
    raise ImportError(
        "ENSO requires the module multiprocessing. Please install the module "
        "multiprocessing."
    )
try:
    from threading import Thread
except ImportError:
    raise ImportError(
        "ENSO requires the module threading. Please install the module threading."
    )
try:
    import time
except ImportError:
    raise ImportError("ENSO requires the module time. Please install the module time.")
try:
    import csv
except ImportError:
    raise ImportError("ENSO requires the module csv. Please install the module csv.")
try:
    import math
except ImportError:
    raise ImportError("ENSO requires the module math. Please install the module math.")
try:
    import json
except ImportError:
    raise ImportError("ENSO requires the module json. Please install the module json.")
try:
    from collections import OrderedDict
except ImportError:
    raise ImportError(
        "ENSO requires the module OrderedDict. Please install the module" "OrderedDict."
    )
try:
    import traceback
except ImportError:
    raise ImportError("ENSO uses the module traceback.")


def cml(
    descr, solvents, func, func3, funcJ, funcS, gfnv, href, cref, fref, pref, siref
):
    """ Get args object from commandline interface.
        Needs argparse module."""

    parser = argparse.ArgumentParser(
        description=descr,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        usage=argparse.SUPPRESS,
    )

    group1 = parser.add_argument_group()
    group1.add_argument(
        "-checkinput",
        "--checkinput",
        dest="checkinput",
        action="store_true",
        required=False,
        help="check if input is correct and stops",
    )
    group1.add_argument(
        "-run",
        "--run",
        dest="run",
        action="store_true",
        required=False,
        help="necessary to run the full ENSO program",
    )
    group1.add_argument(
        "-write_ensorc",
        dest="writeensorc",
        default=False,
        action="store_true",
        required=False,
        help="Write new ensorc, which is placed into the current directory.",
    )

    group2 = parser.add_argument_group("System specific flags")
    group2.add_argument(
        "-chrg",
        "--charge",
        dest="chrg",
        action="store",
        required=False,
        help="charge of the investigated molecule.",
    )
    group2.add_argument(
        "-solv",
        "--solvent",
        dest="solv",
        choices=solvents,
        metavar="",
        action="store",
        required=False,
        help="solvent used in whole procedure, available solvents are: {}.".format(
            solvents
        ),
    )
    group2.add_argument(
        "-u",
        "--unpaired",
        dest="unpaired",
        action="store",
        required=False,
        type=int,
        help="integer number of unpaired electrons of the investigated molecule.",
    )
    group2.add_argument(
        "-T",
        "--temperature",
        dest="temperature",
        action="store",
        required=False,
        help="Temperature in Kelvin for thermostatistical evaluation.",
    )

    group3 = parser.add_argument_group("Programms")
    group3.add_argument(
        "-prog",
        "--program",
        choices=["tm", "orca"],
        dest="prog",
        required=False,
        help="either 'orca' or 'tm",
    )
    group3.add_argument(
        "-prog_rrho",
        "--programrrho",
        choices=["xtb", "prog"],
        dest="rrhoprog",
        required=False,
        help="program for RRHO contribution in part 2 and 3, either 'xtb' or 'prog'.",
    )
    group3.add_argument(
        "-prog3",
        "--programpart3",
        choices=["tm", "orca"],
        dest="prog3",
        required=False,
        help="program used for part3 either 'orca' or 'tm",
    )
    group3.add_argument(
        "-prog4",
        "--programpart4",
        choices=["tm", "orca"],
        dest="prog4",
        required=False,
        help="program for shieldings and couplings, either 'orca' or 'tm",
    )
    group3.add_argument(
        "-gfn",
        "--gfnversion",
        dest="gfnv",
        choices=gfnv,
        action="store",
        required=False,
        help="GFN-xTB version for calculating the RRHO contribution.",
    )
    group3.add_argument(
        "-ancopt",
        "--ancoptoptimization",
        choices=["on", "off"],
        dest="ancopt",
        required=False,
        help="if switched on, xtb will be used to run the optimization in part 1 and 2",
    )

    group4 = parser.add_argument_group("Part specific flags (sorted by part)")
    group4.add_argument(
        "-part1",
        "--part1",
        choices=["on", "off"],
        dest="part1",
        action="store",
        required=False,
        help=" switch part1 on or off",
    )
    group4.add_argument(
        "-thrpart1",
        "--thresholdpart1",
        dest="thr1",
        action="store",
        required=False,
        help="threshold for energetic sorting of part 1",
    )
    group4.add_argument(
        "-func",
        "--functional",
        dest="func",
        choices=func,
        action="store",
        required=False,
        help="functional for part 1 and 2",
    )
    group4.add_argument(
        "-sm",
        "--solventmodel",
        choices=["cpcm", "smd", "cosmo", "dcosmors"],
        dest="sm",
        action="store",
        required=False,
        help="solvent model for the single-point in part 1 and the "
        "optimizations in part 2. For TM, COSMO and DCOSMO-RS are "
        "available (default DCOSMO-RS) and for ORCA, CPCM and SMD "
        "(default SMD).",
    )
    group4.add_argument(
        "-part2",
        "--part2",
        choices=["on", "off"],
        dest="part2",
        action="store",
        required=False,
        help=" switch part2 on or off",
    )
    group4.add_argument(
        "-thrpart2",
        "--thresholdpart2",
        dest="thr2",
        action="store",
        required=False,
        help="threshold for energetic sorting of part 2",
    )
    group4.add_argument(
        "-smgsolv2",
        "--solventmodelgsolvpart2",
        choices=["sm", "cosmors", "gbsa_gsolv"],
        dest="gsolv2",
        action="store",
        required=False,
        help="solvent model for the Gsolv calculation in part 2. Either the solvent model of the "
        "optimizations (sm) or COSMO-RS can be used.",
    )
    group4.add_argument(
        "-part3",
        "--part3",
        choices=["on", "off"],
        dest="part3",
        action="store",
        required=False,
        help=" switch part3 on or off",
    )
    group4.add_argument(
        "-func3",
        "--functionalpart3",
        dest="func3",
        choices=func3,
        action="store",
        required=False,
        help="functional for part3. WARNING: With TM, only the functional PW6B95 is possible.",
    )
    group4.add_argument(
        "-basispart3",
        "--basispart3",
        dest="basis3",
        action="store",
        required=False,
        help="basis for sp in part 3",
    )
    group4.add_argument(
        "-sm3",
        "--solventmoldel3",
        choices=["smd", "dcosmors", "cosmors", "gbsa_gsolv"],
        dest="sm3",
        action="store",
        required=False,
        help="solvent model for part 3. For ORCA, SMD is available and for TM "
        "DCOSMO-RS. COSMO-RS can be used in combination with both "
        "programs.",
    )
    group4.add_argument(
        "-part4",
        "--part4",
        choices=["on", "off"],
        dest="part4",
        action="store",
        required=False,
        help=" switch part4 on or off",
    )
    group4.add_argument(
        "-J",
        "--couplings",
        choices=["on", "off"],
        dest="calcJ",
        action="store",
        required=False,
        help=" switch calculation of couplings on or off",
    )
    group4.add_argument(
        "-funcJ",
        "--functionalcoupling",
        dest="funcJ",
        choices=funcJ,
        action="store",
        required=False,
        help="functional for shielding calculation in part4. WARNING: With TM, "
        "only the functional TPSS is possible.",
    )
    group4.add_argument(
        "-basisJ",
        "--basisJ",
        dest="basisJ",
        action="store",
        required=False,
        help="basis for SSCC calculation",
    )
    group4.add_argument(
        "-S",
        "--shieldings",
        choices=["on", "off"],
        dest="calcS",
        action="store",
        required=False,
        help=" switch calculation of shieldings on or off",
    )
    group4.add_argument(
        "-funcS",
        "--functionalshielding",
        dest="funcS",
        choices=funcS,
        action="store",
        required=False,
        help="functional for shielding calculation in part4.",
    )
    group4.add_argument(
        "-basisS",
        "--basisS",
        dest="basisS",
        action="store",
        required=False,
        help="basis for shielding calculation",
    )
    group4.add_argument(
        "-sm4",
        "--solventmoldel4",
        choices=["cpcm", "smd", "cosmo"],
        dest="sm4",
        action="store",
        required=False,
        help="solvent model for part 4. With TM, only COSMO is possible, for "
        "ORCA, CPCM (default) or SMD can be chosen here.",
    )

    group5 = parser.add_argument_group("NMR specific flags")
    group5.add_argument(
        "-href",
        "-hydrogenreference",
        choices=href,
        dest="href",
        required=False,
        help="reference for 1H spectrum",
    )
    group5.add_argument(
        "-cref",
        "-carbonreference",
        choices=cref,
        dest="cref",
        required=False,
        help="reference for 13C spectrum",
    )
    group5.add_argument(
        "-fref",
        "-fluorinereference",
        choices=fref,
        dest="fref",
        required=False,
        help="reference for 19F spectrum",
    )
    group5.add_argument(
        "-pref",
        "-phosphorusreference",
        choices=pref,
        dest="pref",
        required=False,
        help="reference for 31P spectrum",
    )
    group5.add_argument(
        "-siref",
        "-siliconreference",
        choices=siref,
        dest="siref",
        required=False,
        help="reference for 29Si spectrum",
    )
    group5.add_argument(
        "-hactive",
        "-hydrogenactive",
        choices=["on", "off"],
        dest="hactive",
        required=False,
        help="calculate data for 1H NMR spectrum, either on or off",
    )
    group5.add_argument(
        "-cactive",
        "-carbonactive",
        choices=["on", "off"],
        dest="cactive",
        required=False,
        help="calculate data for 13C NMR spectrum, either on or off",
    )
    group5.add_argument(
        "-factive",
        "-fluorineactive",
        choices=["on", "off"],
        dest="factive",
        required=False,
        help="calculate data for 19F NMR spectrum, either on or off",
    )
    group5.add_argument(
        "-pactive",
        "-phosphorusactive",
        choices=["on", "off"],
        dest="pactive",
        required=False,
        help="calculate data for 31P NMR spectrum, either on or off",
    )
    group5.add_argument(
        "-siactive",
        "-siliconactive",
        choices=["on", "off"],
        dest="siactive",
        required=False,
        help="calculate data for 29Si NMR spectrum, either on or off",
    )
    group5.add_argument(
        "-mf",
        "--resonancefrequency",
        dest="mf",
        type=int,
        action="store",
        required=False,
        help="experimental resonance frequency in MHz",
    )

    group6 = parser.add_argument_group("Further options")
    group6.add_argument(
        "-backup",
        "--backup",
        dest="backup",
        action="store",
        required=False,
        help="recalculation with conformers slightly above the threshold of part 1",
    )
    group6.add_argument(
        "-boltzmann",
        "--boltzmann",
        dest="boltzmann",
        choices=["on", "off"],
        action="store",
        required=False,
        help="only boltzmann population in part 3 is calculated",
    )
    group6.add_argument(
        "-check",
        "--check",
        dest="check",
        choices=["on", "off"],
        action="store",
        required=False,
        help="cautious searching for errors and exiting if too many errors",
    )
    group6.add_argument(
        "-crestcheck",
        "--crestcheck",
        dest="crestcheck",
        choices=["on", "off"],
        action="store",
        required=False,
        help="Checking the DFT-optimized ensemble for identical structures or "
        "rotamers and sort them, be aware that this is done using the "
        "thresholds from CREST.",
    )
    group6.add_argument(
        "-nc",
        "--nconf",
        dest="nstruc",
        type=int,
        action="store",
        required=False,
        help=" # numbers of structures to investigate, default is all conformers.",
    )

    group7 = parser.add_argument_group("Options for parallel calculation")
    group7.add_argument(
        "-O",
        "--omp",
        dest="omp",
        type=int,
        action="store",
        help="sets omp for parallel calculation",
    )
    group7.add_argument(
        "-P",
        "--maxthreads",
        dest="maxthreads",
        type=int,
        action="store",
        help="sets maxthreads for parallel calculation",
    )
    group7.add_argument(
        "--debug",
        dest="debug",
        action="store_true",
        default=False,
        help=argparse.SUPPRESS,
    )
    args = parser.parse_args()
    return args


def mkdir_p(path):
    """ create mkdir -p like behaviour"""
    try:
        os.makedirs(path)
    except OSError as e:
        if os.path.isdir(path):
            pass
        else:
            raise e
    return


def print_block(strlist):
    """Print all elements of strlist in block mode
    e.g. within 80 characters then newline
    """
    length = 0
    try:
        maxlen = max([len(str(x)) for x in strlist])
    except:
        maxlen = 12
    for item in strlist:
        length += maxlen + 2
        if length <= 80:
            if not item == strlist[-1]:  # works only if item only once in list!
                print("{:>{digits}}, ".format(str(item), digits=maxlen), end="")
            else:
                print("{:>{digits}}".format(str(item), digits=maxlen), end="")
        else:
            print("{:>{digits}}".format(str(item), digits=maxlen))
            length = 0
    if length != 0:
        print("\n")
    return


def last_folders(path, number=1):
    """get last folder or last two folders of path, depending on number"""
    if number not in (1, 2):
        number = 1
    if number == 1:
        folder = os.path.basename(path)
    if number == 2:
        folder = os.path.join(
            os.path.basename(os.path.dirname(path)), os.path.basename(path)
        )
    return folder


def check_for_float(line):
    """ Go through line and check for float, return first float"""
    elements = line.split()
    value = None
    for element in elements:
        try:
            value = float(element)
            found = True
        except ValueError:
            found = False
            value = None
        if found:
            break
    return value


def read_json(nstruc, cwd, jsonfile, args):
    """Reading in jsonfile"""

    conf_data = [
        ("crude_opt", "not_calculated"),
        ("energy_crude_opt", None),
        ("backup_for_part2", False),
        ("consider_for_part2", True),
        ("opt", "not_calculated"),
        ("energy_opt", None),
        ("backup_for_part3", False),
        ("sp_part2", "not_calculated"),
        ("energy_sp_part2", None),
        ("consider_for_part3", False),
        ("sp_part3", "not_calculated"),
        ("energy_sp_part3", None),
        ("cosmo-rs", "not_calculated"),
        ("energy_cosmo-rs", None),
        ("gbsa_gsolv", "not_calculated"),
        ("energy_gbsa_gsolv", None),
        ("rrho", "not_calculated"),
        ("energy_rrho", None),
        ("symmetry", "C1"),
        ("consider_for_part4", False),
        ("1H_J", "not_calculated"),
        ("1H_S", "not_calculated"),
        ("13C_J", "not_calculated"),
        ("13C_S", "not_calculated"),
        ("19F_J", "not_calculated"),
        ("19F_S", "not_calculated"),
        ("29Si_J", "not_calculated"),
        ("29Si_S", "not_calculated"),
        ("31P_J", "not_calculated"),
        ("31P_S", "not_calculated"),
        ("removed_by_user", False),
    ]

    firstrun = True
    if os.path.isfile(os.path.join(cwd, jsonfile)):
        print("Reading file: {}".format(jsonfile))
        firstrun = False
        try:
            with open(
                os.path.join(cwd, jsonfile), "r", encoding=coding, newline=None
            ) as inp:
                json_dict = json.load(inp, object_pairs_hook=OrderedDict)
        except Exception as e:
            print("Your Jsonfile is corrupted!\n{}".format(e))
            sys.exit(1)
        # check if all args in flagskey are the same
        flags_unchangable = [
            "unpaired",
            "chrg",
            "solv",
            "prog",
            "gfnv",
            "ancopt",
            "func",
            "sm",
        ]  # must not be changed
        flags_changable = [
            "rrhoprog",
            "prog4",
            "func3",
            "funcJ",
            "funcS",
            "basis3",
            "basisJ",
            "basisS",
            "gsolv2",
            "sm3",
            "sm4",
        ]  # might be changed
        # all other flags of previous run are not affecting restart (e.g. OMP)
        flags_dict = vars(args)
        error_logical = False
        for i in flags_changable:
            if i not in json_dict["flags"]:
                print("ERROR: flag {} is missing in {}!".format(i, jsonfile))
                error_logical = True
        for i in flags_unchangable:
            if i not in json_dict["flags"]:
                print("ERROR: flag {} is missing in {}!".format(i, jsonfile))
                error_logical = True
            elif json_dict["flags"][i] != flags_dict[i]:
                print(
                    "ERROR: flag {} was changed from {} to {}!".format(
                        i, json_dict["flags"][i], flags_dict[i]
                    )
                )
                error_logical = True
        if error_logical:
            print(
                "One or multiple flags from the previous run were changed "
                "or are missing!\nGoing to exit."
            )
            sys.exit(1)
        if args.run:
            # add conformers if necessary
            for i in range(1, nstruc + 1):
                conf = "".join("CONF" + str(i))
                if conf not in [j for j in json_dict.keys()]:
                    json_dict[conf] = OrderedDict(conf_data.copy())
                    print("{} is added to the file enso.json.".format(conf))
        # check if all keys are set
        J_list = ["1H_J", "13C_J", "19F_J", "31P_J", "29Si_J"]
        S_list = ["1H_S", "13C_S", "19F_S", "31P_S", "29Si_S"]
        parts_keys = [
            "crude_opt",
            "opt",
            "sp_part2",
            "cosmo-rs",
            "gbsa_gsolv",
            "rrho",
            "sp_part3",
            "1H_J",
            "13C_J",
            "19F_J",
            "31P_J",
            "29Si_J",
            "1H_S",
            "13C_S",
            "19F_S",
            "31P_S",
            "29Si_S",
        ]
        parts_values = ["calculated", "failed", "not_calculated"]
        keys_numbers = [
            "energy_crude_opt",
            "energy_opt",
            "energy_sp_part2",
            "energy_sp_part3",
            "energy_cosmo-rs",
            "energy_gbsa_gsolv",
            "energy_rrho",
        ]
        keys_true_false = [
            "consider_for_part2",
            "consider_for_part3",
            "consider_for_part4",
            "backup_for_part2",
            "backup_for_part3",
        ]
        # removed_by_user  missing here!

        for item in json_dict.keys():
            if "CONF" in item:
                # add information which might be missing in older
                # enso.json files (<=version 1.2)
                if "29Si_S" not in json_dict[item]:
                    json_dict[item][ "29Si_S"] = "not_calculated"
                if "29Si_J" not in json_dict[item]:
                    json_dict[item][ "29Si_J"] = "not_calculated"
                if "removed_by_user" not in json_dict[item]:
                    json_dict[item]["removed_by_user"] = False
                if "gbsa_gsolv" not in json_dict[item]:
                    json_dict[item]["gbsa_gsolv"] = "not_calculated"
                if "energy_gbsa_gsolv" not in json_dict[item]:
                    json_dict[item]["energy_gbsa_gsolv"] = None
                for key in parts_keys:
                    if key not in json_dict[item]:
                        print(
                            "ERROR: information about {} for {} in the file "
                            "{} are missing!".format(key, item, jsonfile)
                        )
                        error_logical = True
                    elif json_dict[item][key] not in parts_values:
                        print(
                            "ERROR: information about {} for {} in the file "
                            "{} are not recognized!".format(key, item, jsonfile)
                        )
                        error_logical = True
                for key in keys_numbers:
                    if key not in json_dict[item]:
                        print(
                            "ERROR: information about {} for {} in the file "
                            "{} are missing!".format(key, item, jsonfile)
                        )
                        error_logical = True
                    elif json_dict[item][key] is None or isinstance(
                        json_dict[item][key], float
                    ):
                        pass
                    else:
                        print(
                            "ERROR: information about {} for {} in the file "
                            "{} are not recognized!".format(key, item, jsonfile)
                        )
                        error_logical = True
                for key in keys_true_false:
                    if key not in json_dict[item]:
                        print(
                            "ERROR: information about {} for {} in the file "
                            "{} are missing!".format(key, item, jsonfile)
                        )
                        error_logical = True
                    elif json_dict[item][key] not in (True, False):
                        print(
                            "ERROR: information about {} for {} in the file "
                            "{} are missing!".format(key, item, jsonfile)
                        )
                        error_logical = True
                if json_dict[item]["removed_by_user"] not in (True, False):
                    print(
                        "ERROR: information about {} for {} in the file "
                        "{} are not recognized!".format(
                            "removed_by_user", item, jsonfile
                        )
                    )
                    error_logical = True
                if "symmetry" not in json_dict[item]:
                    print(
                        "ERROR: information about {} for {} in the file {} "
                        "are not recognized!".format("symmetry", item, jsonfile)
                    )
                    error_logical = True
        if error_logical:
            print(
                "One or multiple errors were found in the file {}.\nGoing to "
                "exit.".format(jsonfile)
            )
            sys.exit(1)
        # modify json_dict if a flag was changed,
        # if gsolv2 is changed, no effect here
        if json_dict["flags"]["temperature"] != flags_dict["temperature"]:
            print(
                "WARNING: The temperature has been changed with respect to the "
                "previous run.\n       GRRHO is reset to not calculated for "
                "all conformers."
            )
            for i in json_dict.keys():
                if "CONF" in i:
                    json_dict[i]["rrho"] = "not_calculated"
                    json_dict[i]["energy_rrho"] = None
        if args.part2 == "on" or args.part3 == "on":
            if json_dict["flags"]["rrhoprog"] != flags_dict["rrhoprog"]:
                print(
                    "WARNING: The program for calculating GRRHO was changed "
                    "with respect to the previous run.\n         GRRHO is reset "
                    "to not_calculated for all conformers."
                )
                for i in json_dict.keys():
                    if "CONF" in i:
                        json_dict[i]["rrho"] = "not_calculated"
                        json_dict[i]["energy_rrho"] = None
        if args.part3 == "on":
            if (
                json_dict["flags"]["func3"] != flags_dict["func3"]
                or json_dict["flags"]["basis3"] != flags_dict["basis3"]
                or json_dict["flags"]["sm3"] != flags_dict["sm3"]
            ):
                tmpreset = False
                if json_dict["flags"]["func3"] != flags_dict["func3"]:
                    tmpreset = True
                if json_dict["flags"]["basis3"] != flags_dict["basis3"]:
                    tmpreset = True
                if json_dict["flags"]["sm3"] in (
                    "cosmors",
                    "gbsa_gsolv",
                ) and flags_dict["sm3"] in ("cosmors", "gbsa_gsolv"):
                    # dont recalculate SP if only the solent model is
                    # changed between additive corrections
                    pass
                else:
                    # if solvent model changed to eg cosmo dcosmors smd cpcm
                    tmpreset = True
                if tmpreset:
                    print(
                        "WARNING: The density functional / the basis set and/or"
                        " the solvent model for part 3 was changed \n"
                        "         with respect to the previous run.\n"
                        "         The single-point of part 3 is reset to "
                        "not_calculated for all conformers."
                    )
                    for i in json_dict.keys():
                        if "CONF" in i:
                            json_dict[i]["sp_part3"] = "not_calculated"
                            json_dict[i]["energy_sp_part3"] = None
        if args.part4 == "on":
            if (
                json_dict["flags"]["prog4"] != flags_dict["prog4"]
                or json_dict["flags"]["sm4"] != flags_dict["sm4"]
            ):
                print(
                    "WARNING: The program and/or the solvent model for part 4 "
                    "was changed with respect to the previous run. The "
                    "coupling and shielding calculation are reset to "
                    "not_calculated for all conformers."
                )
                for i in json_dict.keys():
                    if "CONF" in i:
                        for j in J_list:
                            json_dict[i][j] = "not_calculated"
                        for j in S_list:
                            json_dict[i][j] = "not_calculated"
            else:
                if (
                    json_dict["flags"]["funcJ"] != flags_dict["funcJ"]
                    or json_dict["flags"]["basisJ"] != flags_dict["basisJ"]
                ):
                    print(
                        "WARNING: The density functional and/or the basis set "
                        "for the coupling calculation was changed with respect "
                        "to the previous run.\n The coupling calculation is reset "
                        "to not_calculated for all conformers."
                    )
                    for i in json_dict.keys():
                        if "CONF" in i:
                            for j in J_list:
                                json_dict[i][j] = "not_calculated"
                if (
                    json_dict["flags"]["funcS"] != flags_dict["funcS"]
                    or json_dict["flags"]["basisS"] != flags_dict["basisS"]
                ):
                    print(
                        "WARNING: The density functional and/or the basis set "
                        "for the shielding calculation was changed with respect "
                        "to the previous run.\n The shielding calculation is "
                        "reset to not_calculated for all conformers."
                    )
                    for i in json_dict.keys():
                        if "CONF" in i:
                            for j in S_list:
                                json_dict[i][j] = "not_calculated"
        # copy json to json#1 and so on
        move_recursively(cwd, jsonfile)
    else:
        if args.run:
            print(
                "The file {} containing all results does not exist yet,\nand "
                "is written during this run.\n".format(jsonfile)
            )
        firstrun = True
        tmp = []
        tmp.append(("flags", OrderedDict()))
        for i in range(1, nstruc + 1):
            conf = "".join("CONF" + str(i))
            tmp.append((conf, OrderedDict(conf_data.copy())))
        json_dict = OrderedDict(tmp)
    json_dict["flags"] = OrderedDict(vars(args))
    return json_dict, firstrun


def write_json(instruction, json_dict, jsonfile):
    """Writing jsonfile"""
    if instruction == "save_and_exit":
        outfile = "".join("error_" + jsonfile)
    else:
        outfile = jsonfile
    print("Results are written to {}.".format(outfile))
    with open(outfile, "w") as out:
        json.dump(json_dict, out, indent=4, sort_keys=False)
    if instruction == "save_and_exit":
        print("\nGoing to exit.")
        sys.exit(1)


def conformersxyz2coord(conformersxyz, nat, directory, conflist, onlyenergy=False):
    """read crest_conformers.xyz and write coord into 
    designated folders, also get GFNx-xTB energies """
    with open(conformersxyz, "r", encoding=coding, newline=None) as xyzfile:
        stringfile_lines = xyzfile.readlines()
    if (len(conflist) * (nat + 2)) > len(stringfile_lines):
        print(
            "ERROR: Either the number of conformers ({}) or the number of "
            "atoms ({}) is wrong!".format(str(len(conflist)), str(nat))
        )
        write_json("save_and_exit", json_dict, jsonfile)
    gfne = [["", None] for _ in conflist]
    counter = 0
    for i in conflist:
        if "CONF" in str(i):
            i = int("".join(filter(str.isdigit, i)))
        gfne[counter][0] = "".join(("CONF", str(i)))
        gfne[counter][1] = check_for_float(stringfile_lines[(i - 1) * (nat + 2) + 1])
        if gfne[counter][1] is None:
            print(
                "Error in float conversion while reading file"
                " {}!".format(conformersxyz)
            )
        if not onlyenergy:
            atom = []
            x = []
            y = []
            z = []
            start = (i - 1) * (nat + 2) + 2
            end = i * (nat + 2)
            bohr2ang = 0.52917721067
            for line in stringfile_lines[start:end]:
                atom.append(str(line.split()[0].lower()))
                x.append(float(line.split()[1]) / bohr2ang)
                y.append(float(line.split()[2]) / bohr2ang)
                z.append(float(line.split()[3]) / bohr2ang)
            coordxyz = []
            for j in range(len(x)):
                coordxyz.append(
                    "{: 09.7f} {: 09.7f}  {: 09.7f}  {}".format(
                        x[j], y[j], z[j], atom[j]
                    )
                )
            if not os.path.isfile(os.path.join("CONF{}".format(i), directory, "coord")):
                print(
                    "Write new coord file in {}".format(
                        os.path.join("CONF{}".format(i), directory)
                    )
                )
                with open(
                    os.path.join("CONF{}".format(i), directory, "coord"),
                    "w",
                    newline=None,
                ) as coord:
                    coord.write("$coord\n")
                    for line in coordxyz:
                        coord.write(line + "\n")
                    coord.write("$end")
        counter = counter + 1
    return gfne


def crest_routine(results, func, crestcheck, json_dict):
    """check if two conformers are rotamers of each other,
    this check is always performed, but removing conformers depends on 
    the value of crestcheck"""
    error_logical = False
    cwd = os.getcwd()
    dirn = "conformer_rotamer_check"  ### directory name
    fn = "conformers.xyz"  ### file name
    dirp = os.path.join(cwd, dirn)  ### directory path
    fp = os.path.join(dirp, fn)  ### file path

    ### create new directory
    if not os.path.isdir(dirp):
        mkdir_p(dirp)
    ### delete file if it already exists
    if os.path.isfile(fp):
        os.remove(fp)

    ### sort conformers according to el. energy
    results.sort(key=lambda x: float(x.energy))
    ### cp coord file of lowest conformer in directory
    try:
        shutil.copy(
            os.path.join(results[0].name, func, "coord"), os.path.join(dirp, "coord")
        )
    except:
        print(
            "ERROR: while copying the coord file from {}! Probably, the "
            "corresponding directory or file does not exist.".format(
                os.path.join(results[0].name, func)
            )
        )
        error_logical = True

    ### write crest file in xyz
    with open(fp, "a", encoding=coding, newline=None) as inp:
        for i in results:
            conf_xyz, nat = coord2xyz(
                os.path.join(cwd, i.name, func)
            )  ### coordinates in xyz
            inp.write("  {}\n".format(nat))  ### number of atoms
            inp.write("{:20.8f}        !{}\n".format(i.energy, i.name))  ### energy
            for line in conf_xyz:
                inp.write(line + "\n")

    print(
        "\nChecking if conformers became rotamers of each other during "
        "the DFT-optimization.\nThe check is performed in the directory"
        " {}.".format(dirn)
    )

    ### call crest
    print("Calling CREST to identify rotamers.")
    with open(os.path.join(dirp, "crest.out"), "w", newline=None) as outputfile:
        subprocess.call(
            [crestpath, "-cregen", fn, "-enso"],
            shell=False,
            stdin=None,
            stderr=subprocess.STDOUT,
            universal_newlines=False,
            cwd=dirp,
            stdout=outputfile,
            env=environsettings,
        )
    time.sleep(0.05)
    ### read in crest results
    try:
        with open(
            os.path.join(dirp, "cregen.enso"), "r", encoding=coding, newline=None
        ) as inp:
            store = inp.readlines()
    except:
        print("ERROR: output file (cregen.enso) of CREST routine does not exist!")
        error_logical = True

    if error_logical:
        print("ERROR: CREST-CHECK can not be performed!")
    else:
        rotlist = []
        rotdict = {}
        if " ALL UNIQUE\n" in store:
            print("No conformers are identified as rotamers or identical.")
        elif " DUPLICATES FOUND\n" in store:
            if args.crestcheck:
                print(
                    "\nWARNING: The following conformers are identified as "
                    "rotamers or identical. They are sorted out."
                )
            else:
                print(
                    "\nWARNING: The following conformers are identified as "
                    "rotamers or identical.\nWARNING: They are NOT sorted out "
                    "since crestcheck is switched off."
                )
            try:
                length = max([len(str(j)) for i in store[1:] for j in i.split()]) + 4
            except ValueError:
                length = 4 + int(args.nstruc) + 1
            print(
                "{:{digits}} {:10}  {:5}<--> {:{digits}} {:10}  {:5}".format(
                    "CONFA", "E(A):", "G(A):", "CONFB", "E(B):", "G(B):", digits=length
                )
            )
            for line in store[1:]:
                confa = results[int(line.split()[0]) - 1]
                confb = results[int(line.split()[1]) - 1]
                rotlist.append(confa.name)
                rotdict[confa.name] = confb.name
                print(
                    "{:{digits}} {:>10.5f} {:>5.2f} <--> {:{digits}} {:>10.5f} {:>5.2f}".format(
                        confa.name,
                        confa.energy,
                        confa.rel_free_energy,
                        confb.name,
                        confb.energy,
                        confb.rel_free_energy,
                        digits=length,
                    )
                )

            if rotlist:
                for confa, confb in rotdict.items():
                    if json_dict[confa]["removed_by_user"]:
                        for i in list(results):
                            if i.name == confa:
                                json_dict[confa]["consider_for_part3"] = False
                                results.remove(i)
                    if json_dict[confb]["removed_by_user"]:
                        for i in list(results):
                            if i.name == confb:
                                json_dict[confb]["consider_for_part3"] = False
                                results.remove(i)

            if args.crestcheck and len(rotlist) > 0:
                for i in list(results):
                    if i.name in rotlist:
                        json_dict[i.name]["consider_for_part3"] = False
                        results.remove(i)
                        print("Removing {:{digits}}.".format(i.name, digits=length))
        else:
            print("ERROR: could not read CREST output (cregen.enso)!.")
            error_logical = True

    if error_logical:
        rotdict = {}

    return rotdict


def write_trj(results, outpath, optfolder):
    """Write trajectory to file"""
    try:
        with open(outpath, "a", encoding=coding, newline=None) as out:
            for i in results:
                conf_xyz, nat = coord2xyz(
                    os.path.join(cwd, i.name, optfolder)
                )  ### coordinates in xyz
                out.write("  {}\n".format(nat))  ### number of atoms
                out.write(
                    "E= {:20.8f}  G= {:20.8f}      !{}\n".format(
                        i.energy, i.free_energy, i.name
                    )
                )  ### energy
                for j in conf_xyz:
                    out.write(
                        " {:3} {:19.10f} {:19.10f} {:19.10f}\n".format(
                            j[0].upper(), float(j[1]), float(j[2]), float(j[3])
                        )
                    )
    except (FileExistsError, ValueError):
        print("Could not write trajectory: {}.".format(last_folders(outpath, 1)))
    return


def coord2xyz(path):
    """convert TURBOMOLE coord file to xyz"""
    time.sleep(0.1)
    bohr2ang = 0.52917721067
    with open(os.path.join(path, "coord"), "r", encoding=coding, newline=None) as f:
        coord = f.readlines()
    x = []
    y = []
    z = []
    atom = []
    for line in coord[1:]:
        if "$" in line:  # stop at $end ...
            break
        x.append(float(line.split()[0]) * bohr2ang)
        y.append(float(line.split()[1]) * bohr2ang)
        z.append(float(line.split()[2]) * bohr2ang)
        atom.append(str(line.split()[3].lower()))
    # natoms = int(len(coord[1:-1])) # unused
    coordxyz = []
    for i in range(len(x)):
        coordxyz.append(
            "{:3} {: 19.10f}  {: 19.10f}  {: 19.10f}".format(atom[i], x[i], y[i], z[i])
        )
    return coordxyz, len(coordxyz)


def RMSD_routine(workdir, nat):
    bohr2ang = 0.52917721067
    new = [[0.0, 0.0, 0.0] for _ in range(nat)]
    old = [[0.0, 0.0, 0.0] for _ in range(nat)]
    squared_difference_x = []
    squared_difference_y = []
    squared_difference_z = []
    with open(os.path.join(workdir, "coord"), "r", encoding=coding, newline=None) as f:
        coord = f.readlines()
        x = []
        y = []
        z = []
        for line in coord[1:]:
            if "$" in line:  # stop at $end ...
                break
            x.append(float(line.split()[0]) * bohr2ang)
            y.append(float(line.split()[1]) * bohr2ang)
            z.append(float(line.split()[2]) * bohr2ang)
    for i in range(len(x)):
        # old.append("{: 09.7f}  {: 09.7f}  {: 09.7f}".format(float(x[i]),float(y[i]),float(z[i])))
        old[i][0] = float(x[i])
        old[i][1] = float(y[i])
        old[i][2] = float(z[i])

    with open(
        os.path.join(workdir, "xtbopt.coord"), "r", encoding=coding, newline=None
    ) as f:
        coord = f.readlines()
        x = []
        y = []
        z = []
        for line in coord[1:]:
            if "$" in line:  # stop at $end ...
                break
            x.append(float(line.split()[0]) * bohr2ang)
            y.append(float(line.split()[1]) * bohr2ang)
            z.append(float(line.split()[2]) * bohr2ang)
    for i in range(len(x)):
        # old.append("{: 09.7f}  {: 09.7f}  {: 09.7f}".format(float(x[i]),float(y[i]),float(z[i])))
        new[i][0] = float(x[i])
        new[i][1] = float(y[i])
        new[i][2] = float(z[i])

    for i in range(len(old)):
        squared_difference_x.append(float((old[i][0] - new[i][0]) ** 2))
        squared_difference_y.append(float((old[i][1] - new[i][1]) ** 2))
        squared_difference_z.append(float((old[i][2] - new[i][2]) ** 2))

    sumdiff = math.fsum(
        squared_difference_x[i] + squared_difference_y[i] + squared_difference_z[i]
        for i in range(nat)
    )

    return math.sqrt(sumdiff / nat)


def RMSD_subprocess(workdir):
    rmsd = None
    with open(os.path.join(workdir, "rmsd"), "w", newline=None) as outputfile:
        callargs = ["rmsd", "coord", "xtbopt.coord"]
        subprocess.call(
            callargs,
            shell=False,
            stdin=None,
            stderr=subprocess.STDOUT,
            universal_newlines=False,
            cwd=workdir,
            stdout=outputfile,
            env=environsettings,
        )
    with open(os.path.join(workdir, "rmsd"), "r", encoding=coding, newline=None) as out:
        stor = out.readlines()
    for line in stor:
        if " rmsd (Angstroem) =" in line:
            try:
                rmsd = float(line.split()[3])
            except:
                rmsd = None
    return rmsd


def write_anmr_enso(cwd, results):
    """write anmr_enso"""
    try:
        length = max([len(i.name[4:]) for i in results])
        if length < 4:
            length = 4
        fmtenergy = max([len("{:.5f}".format(i.energy)) for i in results])
    except:
        length = 6
        fmtenergy = 7
    with open(os.path.join(cwd, "anmr_enso"), "w", newline="") as out:
        out.write(
            "{:5} {:{digits}} {:{digits}} {:6} {:{digits2}} {:7} {:7}\n".format(
                "ONOFF",
                "NMR",
                "CONF",
                "BW",
                "Energy",
                "Gsolv",
                "RRHO",
                digits=length,
                digits2=fmtenergy,
            )
        )
        for i in results:
            out.write(
                "{:<5} {:{digits}} {:{digits}} {:.4f} {:.5f} {:.5f} {:.5f}\n".format(
                    1,
                    i.name[4:],
                    i.name[4:],
                    i.new_bm_weight,
                    i.energy,
                    i.gsolv,
                    i.rrho,
                    digits=length,
                )
            )
    return


def writing_anmrrc(
    prog,
    func,
    sm,
    sm4,
    funcS,
    basisS,
    solv,
    mf,
    href,
    cref,
    fref,
    pref,
    siref,
    hactive,
    cactive,
    factive,
    pactive,
    siactive,
    calcJ,
    calcS,
    save_errors,
):
    h_tm_shieldings = {
        "TMS": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 32.0512048,
                    "acetone": 32.03971003333333,
                    "chcl3": 32.041133316666674,
                    "acetonitrile": 32.03617056666667,
                    "ch2cl2": 32.04777176666666,
                    "dmso": 32.039681316666666,
                    "h2o": 32.036860174999994,
                    "methanol": 32.04573335,
                    "thf": 32.04154705833333,
                    "toluene": 32.02829061666666,
                },
                "pbe0": {
                    "gas": 31.820450258333327,
                    "acetone": 31.801199816666667,
                    "chcl3": 31.807363400000003,
                    "acetonitrile": 31.797744033333334,
                    "ch2cl2": 31.815502166666665,
                    "dmso": 31.797286500000002,
                    "h2o": 31.801018416666665,
                    "methanol": 31.809920125,
                    "thf": 31.802681225,
                    "toluene": 31.790892416666665,
                },
                "pbeh-3c": {
                    "gas": 32.32369869999999,
                    "acetone": 32.30552229166667,
                    "chcl3": 32.30850654166667,
                    "acetonitrile": 32.3015773,
                    "ch2cl2": 32.31627083333333,
                    "dmso": 32.303862816666665,
                    "h2o": 32.30345545833333,
                    "methanol": 32.3130819,
                    "thf": 32.306951225,
                    "toluene": 32.29417180833333,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 32.099305599999994,
                    "acetone": 32.07685382499999,
                    "chcl3": 32.078372550000005,
                    "acetonitrile": 32.067920741666676,
                    "ch2cl2": 32.0876576,
                    "dmso": 32.07713496666667,
                    "h2o": 32.07222951666666,
                    "methanol": 32.085467083333334,
                    "thf": 32.07950451666667,
                    "toluene": 32.06162065,
                },
                "pbe0": {
                    "gas": 31.869211950000004,
                    "acetone": 31.83879448333333,
                    "chcl3": 31.845031441666663,
                    "acetonitrile": 31.829924375,
                    "ch2cl2": 31.855811533333338,
                    "dmso": 31.835178675000005,
                    "h2o": 31.83680665833334,
                    "methanol": 31.850090208333338,
                    "thf": 31.841073758333337,
                    "toluene": 31.824697675,
                },
                "pbeh-3c": {
                    "gas": 32.37107341666667,
                    "acetone": 32.341934458333334,
                    "chcl3": 32.34503841666666,
                    "acetonitrile": 32.332714675,
                    "ch2cl2": 32.35537393333334,
                    "dmso": 32.34058045833333,
                    "h2o": 32.338073200000004,
                    "methanol": 32.35207416666667,
                    "thf": 32.34418670833334,
                    "toluene": 32.32693729166667,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 31.86774000000001,
                    "acetone": 31.848927016666664,
                    "chcl3": 31.851003891666664,
                    "acetonitrile": 31.843538541666664,
                    "ch2cl2": 31.860415141666664,
                    "dmso": 31.849057266666673,
                    "h2o": 31.844762508333332,
                    "methanol": 31.857667625,
                    "thf": 31.851878716666665,
                    "toluene": 31.833541825,
                },
                "pbe0": {
                    "gas": 31.636587116666664,
                    "acetone": 31.60924136666667,
                    "chcl3": 31.616506625,
                    "acetonitrile": 31.604173191666664,
                    "ch2cl2": 31.62743169166667,
                    "dmso": 31.604975658333334,
                    "h2o": 31.607992624999994,
                    "methanol": 31.620864658333335,
                    "thf": 31.611675816666665,
                    "toluene": 31.59546233333333,
                },
                "pbeh-3c": {
                    "gas": 32.14311896666666,
                    "acetone": 32.11710325,
                    "chcl3": 32.12106585833333,
                    "acetonitrile": 32.11156126666667,
                    "ch2cl2": 32.1315459,
                    "dmso": 32.114840533333336,
                    "h2o": 32.11376850833333,
                    "methanol": 32.127508733333336,
                    "thf": 32.11950190833333,
                    "toluene": 32.1023676,
                },
            },
        }
    }
    h_orca_shieldings = {
        "TMS": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 32.17000000000001,
                    "acetone": 32.09433333333334,
                    "chcl3": 32.10649999999999,
                    "acetonitrile": 32.09366666666667,
                    "ch2cl2": 32.099,
                    "dmso": 32.09466666666666,
                    "h2o": 32.10341666666666,
                    "methanol": 32.09250000000001,
                    "thf": 32.10183333333333,
                    "toluene": 32.122833333333325,
                },
                "pbe0": {
                    "gas": 31.819000000000003,
                    "acetone": 31.732666666666663,
                    "chcl3": 31.747000000000003,
                    "acetonitrile": 31.73166666666667,
                    "ch2cl2": 31.738416666666666,
                    "dmso": 31.732666666666663,
                    "h2o": 31.741500000000002,
                    "methanol": 31.73066666666666,
                    "thf": 31.74116666666667,
                    "toluene": 31.765999999999995,
                },
                "dsd-blyp": {
                    "gas": 31.91416666666667,
                    "acetone": 31.83541666666667,
                    "chcl3": 31.84766666666667,
                    "acetonitrile": 31.834666666666667,
                    "ch2cl2": 31.839916666666667,
                    "dmso": 31.835583333333332,
                    "h2o": 31.844166666666666,
                    "methanol": 31.833166666666667,
                    "thf": 31.842583333333334,
                    "toluene": 31.86475,
                },
                "wb97x": {
                    "gas": 31.952,
                    "acetone": 31.867499999999996,
                    "chcl3": 31.880999999999997,
                    "acetonitrile": 31.866666666666664,
                    "ch2cl2": 31.872666666666664,
                    "dmso": 31.86758333333333,
                    "h2o": 31.876083333333337,
                    "methanol": 31.86533333333333,
                    "thf": 31.8755,
                    "toluene": 31.89966666666666,
                },
                "pbeh-3c": {
                    "gas": 32.324999999999996,
                    "acetone": 32.23866666666667,
                    "chcl3": 32.25299999999999,
                    "acetonitrile": 32.23783333333333,
                    "ch2cl2": 32.24466666666667,
                    "dmso": 32.23866666666667,
                    "h2o": 32.24733333333333,
                    "methanol": 32.23666666666667,
                    "thf": 32.24733333333333,
                    "toluene": 32.272,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 32.21800000000001,
                    "acetone": 32.140166666666666,
                    "chcl3": 32.152166666666666,
                    "acetonitrile": 32.140499999999996,
                    "ch2cl2": 32.145,
                    "dmso": 32.14183333333333,
                    "h2o": 32.175000000000004,
                    "methanol": 32.13766666666667,
                    "thf": 32.148,
                    "toluene": 32.168833333333325,
                },
                "pbe0": {
                    "gas": 31.868,
                    "acetone": 31.778999999999996,
                    "chcl3": 31.792583333333337,
                    "acetonitrile": 31.778666666666663,
                    "ch2cl2": 31.784333333333336,
                    "dmso": 31.78016666666667,
                    "h2o": 31.815166666666666,
                    "methanol": 31.77633333333333,
                    "thf": 31.787500000000005,
                    "toluene": 31.812,
                },
                "dsd-blyp": {
                    "gas": 31.962999999999997,
                    "acetone": 31.881250000000005,
                    "chcl3": 31.89325,
                    "acetonitrile": 31.881583333333335,
                    "ch2cl2": 31.886000000000006,
                    "dmso": 31.882583333333333,
                    "h2o": 31.916833333333333,
                    "methanol": 31.878500000000003,
                    "thf": 31.889,
                    "toluene": 31.910750000000004,
                },
                "wb97x": {
                    "gas": 32.00091666666666,
                    "acetone": 31.913416666666663,
                    "chcl3": 31.9265,
                    "acetonitrile": 31.9135,
                    "ch2cl2": 31.918499999999995,
                    "dmso": 31.914666666666665,
                    "h2o": 31.94883333333333,
                    "methanol": 31.910666666666668,
                    "thf": 31.921500000000005,
                    "toluene": 31.94516666666667,
                },
                "pbeh-3c": {
                    "gas": 32.373,
                    "acetone": 32.28366666666667,
                    "chcl3": 32.29716666666666,
                    "acetonitrile": 32.28333333333333,
                    "ch2cl2": 32.288666666666664,
                    "dmso": 32.284499999999994,
                    "h2o": 32.317166666666665,
                    "methanol": 32.28066666666667,
                    "thf": 32.29183333333334,
                    "toluene": 32.31616666666667,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 31.97300000000001,
                    "acetone": 31.898,
                    "chcl3": 31.909500000000005,
                    "acetonitrile": 31.897833333333338,
                    "ch2cl2": 31.902666666666665,
                    "dmso": 31.898999999999997,
                    "h2o": 31.910666666666668,
                    "methanol": 31.89566666666667,
                    "thf": 31.90516666666667,
                    "toluene": 31.925,
                },
                "pbe0": {
                    "gas": 31.625,
                    "acetone": 31.537166666666668,
                    "chcl3": 31.550499999999996,
                    "acetonitrile": 31.536666666666665,
                    "ch2cl2": 31.542500000000004,
                    "dmso": 31.537666666666667,
                    "h2o": 31.549500000000005,
                    "methanol": 31.53458333333334,
                    "thf": 31.545499999999993,
                    "toluene": 31.569,
                },
                "dsd-blyp": {
                    "gas": 31.718000000000004,
                    "acetone": 31.639666666666667,
                    "chcl3": 31.651416666666663,
                    "acetonitrile": 31.639499999999998,
                    "ch2cl2": 31.644083333333338,
                    "dmso": 31.640416666666667,
                    "h2o": 31.65216666666667,
                    "methanol": 31.636916666666664,
                    "thf": 31.64683333333333,
                    "toluene": 31.667833333333334,
                },
                "wb97x": {
                    "gas": 31.757,
                    "acetone": 31.672250000000002,
                    "chcl3": 31.68516666666667,
                    "acetonitrile": 31.67166666666667,
                    "ch2cl2": 31.6775,
                    "dmso": 31.67266666666666,
                    "h2o": 31.68466666666666,
                    "methanol": 31.66966666666667,
                    "thf": 31.680166666666665,
                    "toluene": 31.703,
                },
                "pbeh-3c": {
                    "gas": 32.13400000000001,
                    "acetone": 32.047333333333334,
                    "chcl3": 32.06066666666667,
                    "acetonitrile": 32.04666666666666,
                    "ch2cl2": 32.05266666666666,
                    "dmso": 32.047666666666665,
                    "h2o": 32.059,
                    "methanol": 32.044666666666664,
                    "thf": 32.05566666666666,
                    "toluene": 32.079,
                },
            },
        }
    }
    c_tm_shieldings = {
        "TMS": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 186.6465687,
                    "acetone": 187.27903107499998,
                    "chcl3": 187.238498325,
                    "acetonitrile": 187.372512775,
                    "ch2cl2": 187.0771589,
                    "dmso": 187.243299225,
                    "h2o": 187.37157565,
                    "methanol": 187.10988087500002,
                    "thf": 187.19458635,
                    "toluene": 187.48276635,
                },
                "pbe0": {
                    "gas": 188.859355325,
                    "acetone": 189.6196798,
                    "chcl3": 189.4971041,
                    "acetonitrile": 189.698041075,
                    "ch2cl2": 189.318608125,
                    "dmso": 189.68253637499998,
                    "h2o": 189.65553119999998,
                    "methanol": 189.409198575,
                    "thf": 189.55889105,
                    "toluene": 189.776394325,
                },
                "pbeh-3c": {
                    "gas": 198.41611147499998,
                    "acetone": 199.13367970000002,
                    "chcl3": 199.054179875,
                    "acetonitrile": 199.250248325,
                    "ch2cl2": 198.845265825,
                    "dmso": 199.185056825,
                    "h2o": 199.2289907,
                    "methanol": 198.917945675,
                    "thf": 199.076003325,
                    "toluene": 199.3931504,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 186.97419324999998,
                    "acetone": 187.496073025,
                    "chcl3": 187.45393565,
                    "acetonitrile": 187.554538075,
                    "ch2cl2": 187.31238564999998,
                    "dmso": 187.469466275,
                    "h2o": 187.57139320000002,
                    "methanol": 187.344972675,
                    "thf": 187.42200885,
                    "toluene": 187.671731225,
                },
                "pbe0": {
                    "gas": 189.169130675,
                    "acetone": 189.816064175,
                    "chcl3": 189.69082477499998,
                    "acetonitrile": 189.860330875,
                    "ch2cl2": 189.532330975,
                    "dmso": 189.88587445000002,
                    "h2o": 189.8368566,
                    "methanol": 189.62332455,
                    "thf": 189.76569125,
                    "toluene": 189.94371412499999,
                },
                "pbeh-3c": {
                    "gas": 198.7168509,
                    "acetone": 199.3308802,
                    "chcl3": 199.25125382500002,
                    "acetonitrile": 199.41320919999998,
                    "ch2cl2": 199.06108425,
                    "dmso": 199.390014125,
                    "h2o": 199.41478467500002,
                    "methanol": 199.13192775,
                    "thf": 199.28161922500001,
                    "toluene": 199.562540575,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 185.410099625,
                    "acetone": 185.99193982499997,
                    "chcl3": 185.949648475,
                    "acetonitrile": 186.0799505,
                    "ch2cl2": 185.80363820000002,
                    "dmso": 185.97415155,
                    "h2o": 186.07484635,
                    "methanol": 185.839592875,
                    "thf": 185.9190184,
                    "toluene": 186.17204557500003,
                },
                "pbe0": {
                    "gas": 187.626469575,
                    "acetone": 188.34549135,
                    "chcl3": 188.212218325,
                    "acetonitrile": 188.413268225,
                    "ch2cl2": 188.04820440000003,
                    "dmso": 188.42875420000001,
                    "h2o": 188.3724699,
                    "methanol": 188.14698049999998,
                    "thf": 188.2963985,
                    "toluene": 188.46803717499998,
                },
                "pbeh-3c": {
                    "gas": 197.27823677499998,
                    "acetone": 197.953274625,
                    "chcl3": 197.871683925,
                    "acetonitrile": 198.0615831,
                    "ch2cl2": 197.6764831,
                    "dmso": 198.014841225,
                    "h2o": 198.048432475,
                    "methanol": 197.75143105,
                    "thf": 197.905333025,
                    "toluene": 198.186480775,
                },
            },
        }
    }
    c_orca_shieldings = {
        "TMS": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 188.604,
                    "acetone": 189.7395,
                    "chcl3": 189.5435,
                    "acetonitrile": 189.77,
                    "ch2cl2": 189.6625,
                    "dmso": 189.8015,
                    "h2o": 189.8495,
                    "methanol": 189.77,
                    "thf": 189.647,
                    "toluene": 189.30400000000003,
                },
                "pbe0": {
                    "gas": 188.867,
                    "acetone": 190.265,
                    "chcl3": 190.02224999999999,
                    "acetonitrile": 190.298,
                    "ch2cl2": 190.16649999999998,
                    "dmso": 190.33175,
                    "h2o": 190.38799999999998,
                    "methanol": 190.29875,
                    "thf": 190.1445,
                    "toluene": 189.73375,
                },
                "dsd-blyp": {
                    "gas": 191.37099999999998,
                    "acetone": 192.606,
                    "chcl3": 192.385,
                    "acetonitrile": 192.63599999999997,
                    "ch2cl2": 192.51575000000003,
                    "dmso": 192.66625000000002,
                    "h2o": 192.7205,
                    "methanol": 192.63524999999998,
                    "thf": 192.4955,
                    "toluene": 192.12275,
                },
                "wb97x": {
                    "gas": 190.36075,
                    "acetone": 191.689,
                    "chcl3": 191.453,
                    "acetonitrile": 191.72175000000001,
                    "ch2cl2": 191.5935,
                    "dmso": 191.753,
                    "h2o": 191.8085,
                    "methanol": 191.72150000000002,
                    "thf": 191.57150000000001,
                    "toluene": 191.17225,
                },
                "pbeh-3c": {
                    "gas": 198.458,
                    "acetone": 199.905,
                    "chcl3": 199.649,
                    "acetonitrile": 199.94,
                    "ch2cl2": 199.8025,
                    "dmso": 199.9715,
                    "h2o": 200.0265,
                    "methanol": 199.93900000000002,
                    "thf": 199.7775,
                    "toluene": 199.3395,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 188.908,
                    "acetone": 190.0265,
                    "chcl3": 189.83749999999998,
                    "acetonitrile": 190.062,
                    "ch2cl2": 189.954,
                    "dmso": 190.103,
                    "h2o": 190.07774999999998,
                    "methanol": 190.0595,
                    "thf": 189.9445,
                    "toluene": 189.614,
                },
                "pbe0": {
                    "gas": 189.18025,
                    "acetone": 190.57025000000002,
                    "chcl3": 190.33075,
                    "acetonitrile": 190.60525,
                    "ch2cl2": 190.47,
                    "dmso": 190.65175,
                    "h2o": 190.59925000000004,
                    "methanol": 190.60775,
                    "thf": 190.456,
                    "toluene": 190.058,
                },
                "dsd-blyp": {
                    "gas": 191.66199999999998,
                    "acetone": 192.88025,
                    "chcl3": 192.66174999999998,
                    "acetonitrile": 192.915,
                    "ch2cl2": 192.79025,
                    "dmso": 192.95425,
                    "h2o": 192.91275000000002,
                    "methanol": 192.91250000000002,
                    "thf": 192.77625,
                    "toluene": 192.4135,
                },
                "wb97x": {
                    "gas": 190.65525,
                    "acetone": 191.97199999999998,
                    "chcl3": 191.73825,
                    "acetonitrile": 192.00725,
                    "ch2cl2": 191.875,
                    "dmso": 192.04950000000002,
                    "h2o": 191.99675000000002,
                    "methanol": 192.007,
                    "thf": 191.86025,
                    "toluene": 191.47125,
                },
                "pbeh-3c": {
                    "gas": 198.752,
                    "acetone": 200.196,
                    "chcl3": 199.9445,
                    "acetonitrile": 200.23250000000002,
                    "ch2cl2": 200.0925,
                    "dmso": 200.277,
                    "h2o": 200.15925,
                    "methanol": 200.23350000000002,
                    "thf": 200.075,
                    "toluene": 199.65050000000002,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 187.22,
                    "acetone": 188.442,
                    "chcl3": 188.214,
                    "acetonitrile": 188.4745,
                    "ch2cl2": 188.351,
                    "dmso": 188.5115,
                    "h2o": 188.58350000000002,
                    "methanol": 188.473,
                    "thf": 188.33950000000002,
                    "toluene": 187.965,
                },
                "pbe0": {
                    "gas": 187.5725,
                    "acetone": 188.99225,
                    "chcl3": 188.73424999999997,
                    "acetonitrile": 189.0295,
                    "ch2cl2": 188.8875,
                    "dmso": 189.06875,
                    "h2o": 189.14175,
                    "methanol": 189.0275,
                    "thf": 188.8665,
                    "toluene": 188.4305,
                },
                "dsd-blyp": {
                    "gas": 190.06825,
                    "acetone": 191.39,
                    "chcl3": 191.15425,
                    "acetonitrile": 191.42600000000002,
                    "ch2cl2": 191.29475000000002,
                    "dmso": 191.461,
                    "h2o": 191.53225,
                    "methanol": 191.4225,
                    "thf": 191.27499999999998,
                    "toluene": 190.87675000000002,
                },
                "wb97x": {
                    "gas": 189.04575,
                    "acetone": 190.45225000000002,
                    "chcl3": 190.20074999999997,
                    "acetonitrile": 190.4885,
                    "ch2cl2": 190.35025000000002,
                    "dmso": 190.52525,
                    "h2o": 190.5975,
                    "methanol": 190.4855,
                    "thf": 190.32899999999998,
                    "toluene": 189.904,
                },
                "pbeh-3c": {
                    "gas": 197.184,
                    "acetone": 198.7195,
                    "chcl3": 198.449,
                    "acetonitrile": 198.75799999999998,
                    "ch2cl2": 198.611,
                    "dmso": 198.7955,
                    "h2o": 198.8655,
                    "methanol": 198.755,
                    "thf": 198.587,
                    "toluene": 198.1245,
                },
            },
        }
    }
    f_tm_shieldings = {
        "CFCl3": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 163.5665883,
                    "acetone": 165.9168679,
                    "chcl3": 165.043061,
                    "acetonitrile": 166.377831,
                    "ch2cl2": 164.776383,
                    "dmso": 166.1839641,
                    "h2o": 166.880495,
                    "methanol": 165.4364879,
                    "thf": 165.7384153,
                    "toluene": 165.7812123,
                },
                "pbe0": {
                    "gas": 179.4820255,
                    "acetone": 181.9743764,
                    "chcl3": 181.1338758,
                    "acetonitrile": 182.4438224,
                    "ch2cl2": 180.8751895,
                    "dmso": 182.2224636,
                    "h2o": 182.9958356,
                    "methanol": 181.5031528,
                    "thf": 181.7669891,
                    "toluene": 181.7963177,
                },
                "pbeh-3c": {
                    "gas": 225.045234,
                    "acetone": 226.6335916,
                    "chcl3": 226.0133192,
                    "acetonitrile": 226.9371636,
                    "ch2cl2": 225.8300352,
                    "dmso": 226.8061873,
                    "h2o": 227.4000142,
                    "methanol": 226.3012569,
                    "thf": 226.5247654,
                    "toluene": 226.555523,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 150.4514566,
                    "acetone": 151.5612999,
                    "chcl3": 150.5819485,
                    "acetonitrile": 151.9884593,
                    "ch2cl2": 150.2953968,
                    "dmso": 151.8818575,
                    "h2o": 151.6179136,
                    "methanol": 151.0439011,
                    "thf": 151.4207377,
                    "toluene": 151.4686522,
                },
                "pbe0": {
                    "gas": 167.7783433,
                    "acetone": 169.09491,
                    "chcl3": 168.1354478,
                    "acetonitrile": 169.5416871,
                    "ch2cl2": 167.8558489,
                    "dmso": 169.3950732,
                    "h2o": 169.2178304,
                    "methanol": 168.5860848,
                    "thf": 168.9136991,
                    "toluene": 168.9347931,
                },
                "pbeh-3c": {
                    "gas": 213.6651892,
                    "acetone": 214.1284506,
                    "chcl3": 213.4293417,
                    "acetonitrile": 214.4297108,
                    "ch2cl2": 213.2298905,
                    "dmso": 214.366451,
                    "h2o": 214.1162368,
                    "methanol": 213.76845,
                    "thf": 214.0512078,
                    "toluene": 214.0924969,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 146.4091676,
                    "acetone": 148.7113398,
                    "chcl3": 147.7715256,
                    "acetonitrile": 149.1854535,
                    "ch2cl2": 147.4708159,
                    "dmso": 148.9781692,
                    "h2o": 148.8407317,
                    "methanol": 148.1815132,
                    "thf": 148.5140784,
                    "toluene": 148.6001306,
                },
                "pbe0": {
                    "gas": 163.4654205,
                    "acetone": 165.9356023,
                    "chcl3": 165.0269644,
                    "acetonitrile": 166.4188044,
                    "ch2cl2": 164.7336009,
                    "dmso": 166.1830401,
                    "h2o": 166.0858984,
                    "methanol": 165.4145633,
                    "thf": 165.7038144,
                    "toluene": 165.7726604,
                },
                "pbeh-3c": {
                    "gas": 209.8752809,
                    "acetone": 211.4025693,
                    "chcl3": 210.7286529,
                    "acetonitrile": 211.7120494,
                    "ch2cl2": 210.5166504,
                    "dmso": 211.5990015,
                    "h2o": 211.4250312,
                    "methanol": 211.0321396,
                    "thf": 211.2798891,
                    "toluene": 211.3499046,
                },
            },
        }
    }
    f_orca_shieldings = {
        "CFCl3": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 166.028,
                    "acetone": 167.858,
                    "chcl3": 167.569,
                    "acetonitrile": 167.92,
                    "ch2cl2": 167.732,
                    "dmso": 167.992,
                    "h2o": 168.239,
                    "methanol": 167.889,
                    "thf": 167.737,
                    "toluene": 167.278,
                },
                "pbe0": {
                    "gas": 178.99,
                    "acetone": 181.086,
                    "chcl3": 180.741,
                    "acetonitrile": 181.154,
                    "ch2cl2": 180.939,
                    "dmso": 181.224,
                    "h2o": 181.464,
                    "methanol": 181.123,
                    "thf": 180.934,
                    "toluene": 180.377,
                },
                "dsd-blyp": {
                    "gas": 225.542,
                    "acetone": 227.877,
                    "chcl3": 227.478,
                    "acetonitrile": 227.949,
                    "ch2cl2": 227.712,
                    "dmso": 228.007,
                    "h2o": 228.213,
                    "methanol": 227.919,
                    "thf": 227.691,
                    "toluene": 227.033,
                },
                "wb97x": {
                    "gas": 193.433,
                    "acetone": 195.381,
                    "chcl3": 195.059,
                    "acetonitrile": 195.445,
                    "ch2cl2": 195.245,
                    "dmso": 195.508,
                    "h2o": 195.733,
                    "methanol": 195.415,
                    "thf": 195.239,
                    "toluene": 194.719,
                },
                "pbeh-3c": {
                    "gas": 224.834,
                    "acetone": 226.308,
                    "chcl3": 226.076,
                    "acetonitrile": 226.36,
                    "ch2cl2": 226.207,
                    "dmso": 226.424,
                    "h2o": 226.639,
                    "methanol": 226.333,
                    "thf": 226.215,
                    "toluene": 225.843,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 153.325,
                    "acetone": 153.259,
                    "chcl3": 152.987,
                    "acetonitrile": 153.326,
                    "ch2cl2": 153.137,
                    "dmso": 153.425,
                    "h2o": 153.729,
                    "methanol": 153.292,
                    "thf": 153.16,
                    "toluene": 152.731,
                },
                "pbe0": {
                    "gas": 167.245,
                    "acetone": 167.447,
                    "chcl3": 167.121,
                    "acetonitrile": 167.52,
                    "ch2cl2": 167.31,
                    "dmso": 167.626,
                    "h2o": 167.92,
                    "methanol": 167.486,
                    "thf": 167.322,
                    "toluene": 166.785,
                },
                "dsd-blyp": {
                    "gas": 216.287,
                    "acetone": 217.144,
                    "chcl3": 216.726,
                    "acetonitrile": 217.223,
                    "ch2cl2": 216.969,
                    "dmso": 217.304,
                    "h2o": 217.555,
                    "methanol": 217.19,
                    "thf": 216.957,
                    "toluene": 216.272,
                },
                "wb97x": {
                    "gas": 182.767,
                    "acetone": 182.921,
                    "chcl3": 182.602,
                    "acetonitrile": 182.99,
                    "ch2cl2": 182.783,
                    "dmso": 183.077,
                    "h2o": 183.351,
                    "methanol": 182.957,
                    "thf": 182.792,
                    "toluene": 182.279,
                },
                "pbeh-3c": {
                    "gas": 213.421,
                    "acetone": 213.215,
                    "chcl3": 212.997,
                    "acetonitrile": 213.271,
                    "ch2cl2": 213.116,
                    "dmso": 213.36,
                    "h2o": 213.627,
                    "methanol": 213.241,
                    "thf": 213.14,
                    "toluene": 212.796,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 148.387,
                    "acetone": 149.573,
                    "chcl3": 149.247,
                    "acetonitrile": 149.647,
                    "ch2cl2": 149.43,
                    "dmso": 149.748,
                    "h2o": 150.066,
                    "methanol": 149.609,
                    "thf": 149.446,
                    "toluene": 148.927,
                },
                "pbe0": {
                    "gas": 162.075,
                    "acetone": 163.638,
                    "chcl3": 163.239,
                    "acetonitrile": 163.71,
                    "ch2cl2": 163.472,
                    "dmso": 163.807,
                    "h2o": 164.125,
                    "methanol": 163.671,
                    "thf": 163.476,
                    "toluene": 162.835,
                },
                "dsd-blyp": {
                    "gas": 211.635,
                    "acetone": 213.66,
                    "chcl3": 213.199,
                    "acetonitrile": 213.746,
                    "ch2cl2": 213.469,
                    "dmso": 213.828,
                    "h2o": 214.092,
                    "methanol": 213.71,
                    "thf": 213.451,
                    "toluene": 212.692,
                },
                "wb97x": {
                    "gas": 177.986,
                    "acetone": 179.452,
                    "chcl3": 179.093,
                    "acetonitrile": 179.528,
                    "ch2cl2": 179.299,
                    "dmso": 179.616,
                    "h2o": 179.902,
                    "methanol": 179.491,
                    "thf": 179.302,
                    "toluene": 178.721,
                },
                "pbeh-3c": {
                    "gas": 208.73,
                    "acetone": 209.687,
                    "chcl3": 209.429,
                    "acetonitrile": 209.749,
                    "ch2cl2": 209.573,
                    "dmso": 209.825,
                    "h2o": 210.102,
                    "methanol": 209.716,
                    "thf": 209.592,
                    "toluene": 209.176,
                },
            },
        }
    }
    p_tm_shieldings = {
        "PH3": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 560.9783608,
                    "acetone": 559.5567974,
                    "chcl3": 555.7297268,
                    "acetonitrile": 558.7420853,
                    "ch2cl2": 555.9207578,
                    "dmso": 559.0317956,
                    "h2o": 551.9868157,
                    "methanol": 557.7229598,
                    "thf": 559.4070044,
                    "toluene": 558.9538264,
                },
                "pbe0": {
                    "gas": 573.7889709,
                    "acetone": 572.6807308,
                    "chcl3": 568.6200619,
                    "acetonitrile": 572.0156003,
                    "ch2cl2": 568.6775273,
                    "dmso": 572.2984368,
                    "h2o": 564.8512663,
                    "methanol": 570.6948985,
                    "thf": 572.4491708,
                    "toluene": 572.2945282,
                },
                "pbeh-3c": {
                    "gas": 622.6149401,
                    "acetone": 624.221383,
                    "chcl3": 622.2460822,
                    "acetonitrile": 624.0839458,
                    "ch2cl2": 622.3660073,
                    "dmso": 623.8685076,
                    "h2o": 622.54767,
                    "methanol": 623.1569748,
                    "thf": 623.7253948,
                    "toluene": 623.2733775,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 559.5296772,
                    "acetone": 557.5438599,
                    "chcl3": 553.7653249,
                    "acetonitrile": 556.735552,
                    "ch2cl2": 554.1613395,
                    "dmso": 557.010476,
                    "h2o": 550.1185847,
                    "methanol": 555.82703,
                    "thf": 557.2207586,
                    "toluene": 556.8427805,
                },
                "pbe0": {
                    "gas": 572.4232552,
                    "acetone": 570.7398164,
                    "chcl3": 566.7271447,
                    "acetonitrile": 570.0779914,
                    "ch2cl2": 566.9826221,
                    "dmso": 570.3456887,
                    "h2o": 563.05667,
                    "methanol": 568.8622417,
                    "thf": 570.3305746,
                    "toluene": 570.2507738,
                },
                "pbeh-3c": {
                    "gas": 621.2286124,
                    "acetone": 622.356702,
                    "chcl3": 620.3365742,
                    "acetonitrile": 622.2263079,
                    "ch2cl2": 620.6570087,
                    "dmso": 621.9912341,
                    "h2o": 620.7021951,
                    "methanol": 621.3567408,
                    "thf": 621.7091401,
                    "toluene": 621.3088355,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 558.1589032,
                    "acetone": 556.5475548,
                    "chcl3": 553.3273579,
                    "acetonitrile": 555.6559443,
                    "ch2cl2": 553.600544,
                    "dmso": 556.0983125,
                    "h2o": 548.970911,
                    "methanol": 555.4535832,
                    "thf": 556.3191064,
                    "toluene": 555.9299261,
                },
                "pbe0": {
                    "gas": 571.012794,
                    "acetone": 569.7250563,
                    "chcl3": 566.2936179,
                    "acetonitrile": 568.9923465,
                    "ch2cl2": 566.4237381,
                    "dmso": 569.4236946,
                    "h2o": 561.898531,
                    "methanol": 568.4989088,
                    "thf": 569.4140377,
                    "toluene": 569.3191735,
                },
                "pbeh-3c": {
                    "gas": 620.0674752,
                    "acetone": 621.5116584,
                    "chcl3": 619.9397925,
                    "acetonitrile": 621.2898165,
                    "ch2cl2": 620.15928,
                    "dmso": 621.2154327,
                    "h2o": 619.7280828,
                    "methanol": 621.0126668,
                    "thf": 620.9449236,
                    "toluene": 620.5363442,
                },
            },
        },
        "TMP": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 281.6302978,
                    "acetone": 265.4354914,
                    "chcl3": 257.5409613,
                    "acetonitrile": 263.2430698,
                    "ch2cl2": 257.0543221,
                    "dmso": 262.8752182,
                    "h2o": 242.4838211,
                    "methanol": 245.6431135,
                    "thf": 266.7188352,
                    "toluene": 269.0597797,
                },
                "pbe0": {
                    "gas": 277.8252556,
                    "acetone": 261.5502528,
                    "chcl3": 254.1109855,
                    "acetonitrile": 259.5059377,
                    "ch2cl2": 253.6358478,
                    "dmso": 258.7821425,
                    "h2o": 239.5329333,
                    "methanol": 242.1687948,
                    "thf": 262.8378646,
                    "toluene": 265.4050199,
                },
                "pbeh-3c": {
                    "gas": 390.6073841,
                    "acetone": 378.6668397,
                    "chcl3": 373.2000393,
                    "acetonitrile": 377.1343123,
                    "ch2cl2": 372.9163524,
                    "dmso": 376.6203422,
                    "h2o": 362.7163813,
                    "methanol": 364.8220379,
                    "thf": 379.5051748,
                    "toluene": 381.2789752,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 276.8654211,
                    "acetone": 259.8829696,
                    "chcl3": 251.5648819,
                    "acetonitrile": 257.7225804,
                    "ch2cl2": 251.0880934,
                    "dmso": 256.90761,
                    "h2o": 234.4800595,
                    "methanol": 237.4630709,
                    "thf": 261.291204,
                    "toluene": 263.9827571,
                },
                "pbe0": {
                    "gas": 273.0911933,
                    "acetone": 256.1507446,
                    "chcl3": 248.2072561,
                    "acetonitrile": 254.0571117,
                    "ch2cl2": 247.7513367,
                    "dmso": 253.0100842,
                    "h2o": 231.7425518,
                    "methanol": 234.1695454,
                    "thf": 257.4644157,
                    "toluene": 260.3717755,
                },
                "pbeh-3c": {
                    "gas": 386.2437698,
                    "acetone": 373.8145109,
                    "chcl3": 368.1719462,
                    "acetonitrile": 372.350904,
                    "ch2cl2": 367.8934403,
                    "dmso": 371.4995766,
                    "h2o": 355.9965281,
                    "methanol": 358.0517851,
                    "thf": 374.7716841,
                    "toluene": 376.8283779,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 278.0447826,
                    "acetone": 261.4382678,
                    "chcl3": 253.5317417,
                    "acetonitrile": 259.5831076,
                    "ch2cl2": 253.0735218,
                    "dmso": 258.8205488,
                    "h2o": 236.9938311,
                    "methanol": 240.0596152,
                    "thf": 262.646474,
                    "toluene": 265.5482099,
                },
                "pbe0": {
                    "gas": 274.1582231,
                    "acetone": 257.5976215,
                    "chcl3": 250.0455696,
                    "acetonitrile": 255.8739799,
                    "ch2cl2": 249.6032437,
                    "dmso": 254.7109046,
                    "h2o": 234.1066151,
                    "methanol": 236.6658834,
                    "thf": 258.6914971,
                    "toluene": 261.8410368,
                },
                "pbeh-3c": {
                    "gas": 387.4697022,
                    "acetone": 375.2569197,
                    "chcl3": 369.9533245,
                    "acetonitrile": 374.0256406,
                    "ch2cl2": 369.6688695,
                    "dmso": 373.1520781,
                    "h2o": 358.1827766,
                    "methanol": 360.3168296,
                    "thf": 376.0015788,
                    "toluene": 378.3153047,
                },
            },
        },
    }
    p_orca_shieldings = {
        "PH3": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 578.49,
                    "acetone": 577.53,
                    "chcl3": 577.773,
                    "acetonitrile": 577.631,
                    "ch2cl2": 577.63,
                    "dmso": 577.688,
                    "h2o": 577.764,
                    "methanol": 577.506,
                    "thf": 577.671,
                    "toluene": 577.946,
                },
                "pbe0": {
                    "gas": 573.639,
                    "acetone": 573.637,
                    "chcl3": 573.71,
                    "acetonitrile": 573.764,
                    "ch2cl2": 573.67,
                    "dmso": 573.829,
                    "h2o": 573.914,
                    "methanol": 573.632,
                    "thf": 573.688,
                    "toluene": 573.665,
                },
                "dsd-blyp": {
                    "gas": 569.431,
                    "acetone": 567.575,
                    "chcl3": 567.994,
                    "acetonitrile": 567.65,
                    "ch2cl2": 567.746,
                    "dmso": 567.695,
                    "h2o": 567.745,
                    "methanol": 567.531,
                    "thf": 567.809,
                    "toluene": 568.372,
                },
                "wb97x": {
                    "gas": 568.27,
                    "acetone": 568.185,
                    "chcl3": 568.261,
                    "acetonitrile": 568.31,
                    "ch2cl2": 568.218,
                    "dmso": 568.375,
                    "h2o": 568.459,
                    "methanol": 568.18,
                    "thf": 568.236,
                    "toluene": 568.231,
                },
                "pbeh-3c": {
                    "gas": 622.505,
                    "acetone": 626.377,
                    "chcl3": 625.536,
                    "acetonitrile": 626.609,
                    "ch2cl2": 626.042,
                    "dmso": 626.709,
                    "h2o": 626.85,
                    "methanol": 626.48,
                    "thf": 625.933,
                    "toluene": 624.513,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 574.673,
                    "acetone": 575.587,
                    "chcl3": 575.672,
                    "acetonitrile": 575.6,
                    "ch2cl2": 575.619,
                    "dmso": 575.662,
                    "h2o": 575.948,
                    "methanol": 575.57,
                    "thf": 575.668,
                    "toluene": 575.8,
                },
                "pbe0": {
                    "gas": 569.721,
                    "acetone": 571.667,
                    "chcl3": 571.577,
                    "acetonitrile": 571.703,
                    "ch2cl2": 571.631,
                    "dmso": 571.774,
                    "h2o": 572.075,
                    "methanol": 571.67,
                    "thf": 571.656,
                    "toluene": 571.48,
                },
                "dsd-blyp": {
                    "gas": 565.936,
                    "acetone": 565.88,
                    "chcl3": 566.179,
                    "acetonitrile": 565.866,
                    "ch2cl2": 566.012,
                    "dmso": 565.915,
                    "h2o": 566.166,
                    "methanol": 565.843,
                    "thf": 566.084,
                    "toluene": 566.506,
                },
                "wb97x": {
                    "gas": 564.429,
                    "acetone": 566.244,
                    "chcl3": 566.161,
                    "acetonitrile": 566.279,
                    "ch2cl2": 566.206,
                    "dmso": 566.349,
                    "h2o": 566.646,
                    "methanol": 566.247,
                    "thf": 566.233,
                    "toluene": 566.086,
                },
                "pbeh-3c": {
                    "gas": 618.99,
                    "acetone": 624.483,
                    "chcl3": 623.499,
                    "acetonitrile": 624.639,
                    "ch2cl2": 624.087,
                    "dmso": 624.744,
                    "h2o": 625.072,
                    "methanol": 624.593,
                    "thf": 623.983,
                    "toluene": 622.448,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 574.839,
                    "acetone": 574.09,
                    "chcl3": 574.267,
                    "acetonitrile": 574.11,
                    "ch2cl2": 574.167,
                    "dmso": 574.166,
                    "h2o": 574.435,
                    "methanol": 574.084,
                    "thf": 574.22,
                    "toluene": 574.478,
                },
                "pbe0": {
                    "gas": 569.911,
                    "acetone": 570.088,
                    "chcl3": 570.127,
                    "acetonitrile": 570.133,
                    "ch2cl2": 570.135,
                    "dmso": 570.198,
                    "h2o": 570.482,
                    "methanol": 570.103,
                    "thf": 570.164,
                    "toluene": 570.119,
                },
                "dsd-blyp": {
                    "gas": 566.08,
                    "acetone": 564.411,
                    "chcl3": 564.793,
                    "acetonitrile": 564.406,
                    "ch2cl2": 564.583,
                    "dmso": 564.448,
                    "h2o": 564.684,
                    "methanol": 564.385,
                    "thf": 564.658,
                    "toluene": 565.213,
                },
                "wb97x": {
                    "gas": 564.63,
                    "acetone": 564.706,
                    "chcl3": 564.726,
                    "acetonitrile": 564.75,
                    "ch2cl2": 564.72,
                    "dmso": 564.813,
                    "h2o": 565.093,
                    "methanol": 564.721,
                    "thf": 564.752,
                    "toluene": 564.742,
                },
                "pbeh-3c": {
                    "gas": 619.182,
                    "acetone": 623.189,
                    "chcl3": 622.29,
                    "acetonitrile": 623.352,
                    "ch2cl2": 622.833,
                    "dmso": 623.451,
                    "h2o": 623.764,
                    "methanol": 623.308,
                    "thf": 622.734,
                    "toluene": 621.304,
                },
            },
        },
        "TMP": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 291.33,
                    "acetone": 276.264,
                    "chcl3": 277.254,
                    "acetonitrile": 275.207,
                    "ch2cl2": 276.171,
                    "dmso": 276.988,
                    "h2o": 262.671,
                    "methanol": 263.366,
                    "thf": 278.685,
                    "toluene": 283.761,
                },
                "pbe0": {
                    "gas": 277.761,
                    "acetone": 262.673,
                    "chcl3": 263.634,
                    "acetonitrile": 261.631,
                    "ch2cl2": 262.58,
                    "dmso": 263.406,
                    "h2o": 249.27,
                    "methanol": 249.931,
                    "thf": 265.061,
                    "toluene": 270.123,
                },
                "dsd-blyp": {
                    "gas": 299.195,
                    "acetone": 286.35,
                    "chcl3": 287.213,
                    "acetonitrile": 285.469,
                    "ch2cl2": 286.302,
                    "dmso": 286.997,
                    "h2o": 274.843,
                    "methanol": 275.42,
                    "thf": 288.362,
                    "toluene": 292.724,
                },
                "wb97x": {
                    "gas": 277.52,
                    "acetone": 262.317,
                    "chcl3": 263.295,
                    "acetonitrile": 261.26,
                    "ch2cl2": 262.227,
                    "dmso": 263.036,
                    "h2o": 248.805,
                    "methanol": 249.485,
                    "thf": 264.716,
                    "toluene": 269.816,
                },
                "pbeh-3c": {
                    "gas": 390.602,
                    "acetone": 379.7,
                    "chcl3": 380.279,
                    "acetonitrile": 378.978,
                    "ch2cl2": 379.593,
                    "dmso": 380.317,
                    "h2o": 368.831,
                    "methanol": 369.216,
                    "thf": 381.391,
                    "toluene": 384.986,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 286.404,
                    "acetone": 270.748,
                    "chcl3": 271.725,
                    "acetonitrile": 269.462,
                    "ch2cl2": 270.524,
                    "dmso": 271.355,
                    "h2o": 256.342,
                    "methanol": 257.122,
                    "thf": 273.469,
                    "toluene": 278.676,
                },
                "pbe0": {
                    "gas": 272.706,
                    "acetone": 257.164,
                    "chcl3": 258.119,
                    "acetonitrile": 255.895,
                    "ch2cl2": 256.94,
                    "dmso": 257.797,
                    "h2o": 242.92,
                    "methanol": 243.667,
                    "thf": 259.855,
                    "toluene": 264.973,
                },
                "dsd-blyp": {
                    "gas": 294.405,
                    "acetone": 281.158,
                    "chcl3": 282.018,
                    "acetonitrile": 280.073,
                    "ch2cl2": 280.993,
                    "dmso": 281.703,
                    "h2o": 269.086,
                    "methanol": 269.737,
                    "thf": 283.464,
                    "toluene": 287.882,
                },
                "wb97x": {
                    "gas": 272.595,
                    "acetone": 256.861,
                    "chcl3": 257.836,
                    "acetonitrile": 255.578,
                    "ch2cl2": 256.643,
                    "dmso": 257.483,
                    "h2o": 242.627,
                    "methanol": 243.389,
                    "thf": 259.577,
                    "toluene": 264.773,
                },
                "pbeh-3c": {
                    "gas": 385.991,
                    "acetone": 374.828,
                    "chcl3": 375.394,
                    "acetonitrile": 373.92,
                    "ch2cl2": 374.61,
                    "dmso": 375.349,
                    "h2o": 363.431,
                    "methanol": 363.874,
                    "thf": 376.762,
                    "toluene": 380.401,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 286.331,
                    "acetone": 271.022,
                    "chcl3": 271.947,
                    "acetonitrile": 269.751,
                    "ch2cl2": 270.768,
                    "dmso": 271.616,
                    "h2o": 256.882,
                    "methanol": 257.6,
                    "thf": 273.659,
                    "toluene": 278.687,
                },
                "pbe0": {
                    "gas": 272.619,
                    "acetone": 257.298,
                    "chcl3": 258.198,
                    "acetonitrile": 256.053,
                    "ch2cl2": 257.051,
                    "dmso": 257.926,
                    "h2o": 243.408,
                    "methanol": 244.095,
                    "thf": 259.935,
                    "toluene": 264.977,
                },
                "dsd-blyp": {
                    "gas": 294.334,
                    "acetone": 281.319,
                    "chcl3": 282.131,
                    "acetonitrile": 280.265,
                    "ch2cl2": 281.144,
                    "dmso": 281.852,
                    "h2o": 269.472,
                    "methanol": 270.068,
                    "thf": 283.556,
                    "toluene": 287.875,
                },
                "wb97x": {
                    "gas": 272.586,
                    "acetone": 257.148,
                    "chcl3": 258.069,
                    "acetonitrile": 255.901,
                    "ch2cl2": 256.919,
                    "dmso": 257.755,
                    "h2o": 243.195,
                    "methanol": 243.894,
                    "thf": 259.785,
                    "toluene": 264.863,
                },
                "pbeh-3c": {
                    "gas": 385.897,
                    "acetone": 374.881,
                    "chcl3": 375.407,
                    "acetonitrile": 373.999,
                    "ch2cl2": 374.652,
                    "dmso": 375.391,
                    "h2o": 363.697,
                    "methanol": 364.097,
                    "thf": 376.757,
                    "toluene": 380.319,
                },
            },
        },
    }
    si_tm_shieldings = {
        "TMS": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 334.2579542,
                    "acetone": 334.1639413,
                    "chcl3": 334.1459912,
                    "acetonitrile": 334.1644763,
                    "ch2cl2": 334.143167,
                    "dmso": 334.2355086,
                    "h2o": 334.1700712,
                    "methanol": 334.1638302,
                    "thf": 334.1765686,
                    "toluene": 334.1672644,
                },
                "pbe0": {
                    "gas": 332.1432161,
                    "acetone": 332.0806043,
                    "chcl3": 332.027555,
                    "acetonitrile": 332.070525,
                    "ch2cl2": 332.0181509,
                    "dmso": 332.1389588,
                    "h2o": 332.0768365,
                    "methanol": 332.082777,
                    "thf": 332.0989747,
                    "toluene": 332.0655251,
                },
                "pbeh-3c": {
                    "gas": 425.4500968,
                    "acetone": 425.4194168,
                    "chcl3": 425.3783658,
                    "acetonitrile": 425.4187809,
                    "ch2cl2": 425.3492293,
                    "dmso": 425.4302912,
                    "h2o": 425.4004059,
                    "methanol": 425.3865089,
                    "thf": 425.4157351,
                    "toluene": 425.4555181,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 334.5698984,
                    "acetone": 334.0803779,
                    "chcl3": 334.1093328,
                    "acetonitrile": 334.0665281,
                    "ch2cl2": 334.1280337,
                    "dmso": 334.1272572,
                    "h2o": 334.0495564,
                    "methanol": 334.1137413,
                    "thf": 334.1251606,
                    "toluene": 334.1235476,
                },
                "pbe0": {
                    "gas": 332.3546979,
                    "acetone": 331.9058869,
                    "chcl3": 331.8955148,
                    "acetonitrile": 331.8800833,
                    "ch2cl2": 331.9140658,
                    "dmso": 331.948424,
                    "h2o": 331.8617288,
                    "methanol": 331.9375391,
                    "thf": 331.9562723,
                    "toluene": 331.9253075,
                },
                "pbeh-3c": {
                    "gas": 426.0062656,
                    "acetone": 425.7811084,
                    "chcl3": 425.7602588,
                    "acetonitrile": 425.745999,
                    "ch2cl2": 425.7473718,
                    "dmso": 425.779427,
                    "h2o": 425.7365851,
                    "methanol": 425.7713265,
                    "thf": 425.7964293,
                    "toluene": 425.8200844,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 333.7779314,
                    "acetone": 333.3511708,
                    "chcl3": 333.3794838,
                    "acetonitrile": 333.3298692,
                    "ch2cl2": 333.3946486,
                    "dmso": 333.3881767,
                    "h2o": 333.3406562,
                    "methanol": 333.3784136,
                    "thf": 333.3860666,
                    "toluene": 333.3885135,
                },
                "pbe0": {
                    "gas": 331.5820841,
                    "acetone": 331.1904714,
                    "chcl3": 331.1839521,
                    "acetonitrile": 331.1565218,
                    "ch2cl2": 331.1982524,
                    "dmso": 331.2347884,
                    "h2o": 331.1670301,
                    "methanol": 331.2231923,
                    "thf": 331.2383692,
                    "toluene": 331.2108329,
                },
                "pbeh-3c": {
                    "gas": 425.0726297,
                    "acetone": 424.9009564,
                    "chcl3": 424.8706079,
                    "acetonitrile": 424.8831877,
                    "ch2cl2": 424.8554965,
                    "dmso": 424.9143792,
                    "h2o": 424.8579037,
                    "methanol": 424.8851226,
                    "thf": 424.9146175,
                    "toluene": 424.9330242,
                },
            },
        }
    }
    si_orca_shieldings = {
        "TMS": {
            "pbeh-3c": {
                "tpss": {
                    "gas": 344.281,
                    "acetone": 344.239,
                    "chcl3": 344.311,
                    "acetonitrile": 344.198,
                    "ch2cl2": 344.231,
                    "dmso": 344.292,
                    "h2o": 344.228,
                    "methanol": 344.291,
                    "thf": 344.283,
                    "toluene": 344.452,
                },
                "pbe0": {
                    "gas": 332.181,
                    "acetone": 332.067,
                    "chcl3": 332.162,
                    "acetonitrile": 332.033,
                    "ch2cl2": 332.082,
                    "dmso": 332.122,
                    "h2o": 332.048,
                    "methanol": 332.122,
                    "thf": 332.134,
                    "toluene": 332.298,
                },
                "dsd-blyp": {
                    "gas": 357.874,
                    "acetone": 357.762,
                    "chcl3": 357.864,
                    "acetonitrile": 357.726,
                    "ch2cl2": 357.783,
                    "dmso": 357.798,
                    "h2o": 357.715,
                    "methanol": 357.809,
                    "thf": 357.826,
                    "toluene": 358.001,
                },
                "wb97x": {
                    "gas": 335.739,
                    "acetone": 335.641,
                    "chcl3": 335.74,
                    "acetonitrile": 335.606,
                    "ch2cl2": 335.659,
                    "dmso": 335.687,
                    "h2o": 335.608,
                    "methanol": 335.692,
                    "thf": 335.707,
                    "toluene": 335.879,
                },
                "pbeh-3c": {
                    "gas": 425.385,
                    "acetone": 425.52,
                    "chcl3": 425.527,
                    "acetonitrile": 425.511,
                    "ch2cl2": 425.508,
                    "dmso": 425.578,
                    "h2o": 425.566,
                    "methanol": 425.557,
                    "thf": 425.54,
                    "toluene": 425.556,
                },
            },
            "b97-3c": {
                "tpss": {
                    "gas": 344.503,
                    "acetone": 344.558,
                    "chcl3": 344.676,
                    "acetonitrile": 344.487,
                    "ch2cl2": 344.537,
                    "dmso": 344.67,
                    "h2o": 344.542,
                    "methanol": 344.662,
                    "thf": 344.637,
                    "toluene": 344.919,
                },
                "pbe0": {
                    "gas": 332.338,
                    "acetone": 332.293,
                    "chcl3": 332.442,
                    "acetonitrile": 332.236,
                    "ch2cl2": 332.31,
                    "dmso": 332.4,
                    "h2o": 332.288,
                    "methanol": 332.392,
                    "thf": 332.403,
                    "toluene": 332.676,
                },
                "dsd-blyp": {
                    "gas": 357.729,
                    "acetone": 357.628,
                    "chcl3": 357.774,
                    "acetonitrile": 357.578,
                    "ch2cl2": 357.655,
                    "dmso": 357.692,
                    "h2o": 357.632,
                    "methanol": 357.703,
                    "thf": 357.725,
                    "toluene": 357.985,
                },
                "wb97x": {
                    "gas": 335.744,
                    "acetone": 335.688,
                    "chcl3": 335.837,
                    "acetonitrile": 335.633,
                    "ch2cl2": 335.71,
                    "dmso": 335.774,
                    "h2o": 335.704,
                    "methanol": 335.776,
                    "thf": 335.792,
                    "toluene": 336.064,
                },
                "pbeh-3c": {
                    "gas": 425.911,
                    "acetone": 426.14,
                    "chcl3": 426.185,
                    "acetonitrile": 426.113,
                    "ch2cl2": 426.124,
                    "dmso": 426.254,
                    "h2o": 426.162,
                    "methanol": 426.22,
                    "thf": 426.196,
                    "toluene": 426.294,
                },
            },
            "tpss": {
                "tpss": {
                    "gas": 343.24,
                    "acetone": 343.388,
                    "chcl3": 343.506,
                    "acetonitrile": 343.343,
                    "ch2cl2": 343.385,
                    "dmso": 343.48,
                    "h2o": 343.378,
                    "methanol": 343.47,
                    "thf": 343.449,
                    "toluene": 343.647,
                },
                "pbe0": {
                    "gas": 331.055,
                    "acetone": 331.217,
                    "chcl3": 331.313,
                    "acetonitrile": 331.175,
                    "ch2cl2": 331.224,
                    "dmso": 331.303,
                    "h2o": 331.205,
                    "methanol": 331.296,
                    "thf": 331.293,
                    "toluene": 331.461,
                },
                "dsd-blyp": {
                    "gas": 357.099,
                    "acetone": 357.125,
                    "chcl3": 357.231,
                    "acetonitrile": 357.081,
                    "ch2cl2": 357.141,
                    "dmso": 357.179,
                    "h2o": 357.075,
                    "methanol": 357.188,
                    "thf": 357.195,
                    "toluene": 357.379,
                },
                "wb97x": {
                    "gas": 334.802,
                    "acetone": 334.886,
                    "chcl3": 334.987,
                    "acetonitrile": 334.842,
                    "ch2cl2": 334.897,
                    "dmso": 334.957,
                    "h2o": 334.855,
                    "methanol": 334.958,
                    "thf": 334.959,
                    "toluene": 335.134,
                },
                "pbeh-3c": {
                    "gas": 424.346,
                    "acetone": 424.653,
                    "chcl3": 424.66,
                    "acetonitrile": 424.64,
                    "ch2cl2": 424.633,
                    "dmso": 424.74,
                    "h2o": 424.718,
                    "methanol": 424.709,
                    "thf": 424.681,
                    "toluene": 424.701,
                },
            },
        }
    }
    # shieldings and solvent models
    if solv is None:
        solv = "gas"
    if prog == "tm":
        # print('NMR data: func {}, funcS {}, href {}, solv {}'.format(str(func), str(funcS), str(href), str(solv)))
        try:
            hshielding = "{:4.3f}".format(h_tm_shieldings[href][func][funcS][solv])
            cshielding = "{:4.3f}".format(c_tm_shieldings[cref][func][funcS][solv])
            fshielding = "{:4.3f}".format(f_tm_shieldings[fref][func][funcS][solv])
            pshielding = "{:4.3f}".format(p_tm_shieldings[pref][func][funcS][solv])
            sishielding = "{:4.3f}".format(si_tm_shieldings[siref][func][funcS][solv])
        except KeyError:
            hshielding = 0
            cshielding = 0
            fshielding = 0
            pshielding = 0
            sishielding = 0
            print(
                "ERROR! The reference absolute shielding constant could not be"
                " found!\n You have to edit the file .anmrrc by hand!"
            )
            save_errors.append(
                "ERROR! The reference absolute shielding constant could not be"
                " found!\n You have to edit the file .anmrrc by hand!"
            )
        if sm == "cosmo":
            print(
                "WARNING: The geometry optimization of the reference molecule "
                "was calculated with DCOSMO-RS instead of COSMO as solvent "
                "model (sm)!"
            )
        if basisS != "def2-TZVP":
            print(
                "WARNING: The reference shielding was calculated with the "
                "basis def2-TZVP (basisS)!"
            )
        sm = "DCOSMO-RS"
        sm4 = "DCOSMO-RS"
        basisS = "def2-TZVP"
    elif prog == "orca":
        try:
            hshielding = "{:4.3f}".format(h_orca_shieldings[href][func][funcS][solv])
            cshielding = "{:4.3f}".format(c_orca_shieldings[cref][func][funcS][solv])
            fshielding = "{:4.3f}".format(f_orca_shieldings[fref][func][funcS][solv])
            pshielding = "{:4.3f}".format(p_orca_shieldings[pref][func][funcS][solv])
            sishielding = "{:4.3f}".format(si_orca_shieldings[siref][func][funcS][solv])
        except KeyError:
            hshielding = 0
            cshielding = 0
            fshielding = 0
            pshielding = 0
            sishielding = 0
            print(
                "ERROR! The reference absolute shielding constant could not be"
                " found!\n You have to edit the file .anmrrc by hand!"
            )
            save_errors.append(
                "ERROR! The reference absolute shielding constant could not be"
                " found!\n You have to edit the file .anmrrc by hand!"
            )
        if sm == "cpcm":
            print(
                "WARNING: The geometry optimization of the reference molecule "
                "was calculated with SMD instead of CPCM as solvent model (sm)!"
            )
        if sm4 == "cpcm":
            print(
                "WARNING: The reference shielding was calculated with SMD "
                "instead of CPCM as solvent model (sm)!"
            )
        if basisS != "def2-TZVP":
            print(
                "WARNING: The reference shielding was calculated with the "
                "basis def2-TZVP (basisS)!"
            )
        sm = "SMD"
        sm4 = "SMD"
        basisS = "def2-TZVP"
    # basis set for optimization
    if func == "pbeh-3c":
        basisfunc = "def2-mSVP"
    elif func == "b97-3c":
        basisfunc = "def2-mTZVP"
    elif func == "tpss":
        basisfunc = "def2-TZVP"
    else:
        basisfunc = "unknown"
    # write .anmrrc
    if args.prog4 == "tm":
        prog = "TM"
    elif args.prog4 == "orca":
        prog = "ORCA"
    ha = 0
    ca = 0
    fa = 0
    pa = 0
    sia = 0
    if hactive == "on":
        ha = 1
    if cactive == "on":
        ca = 1
    if factive == "on":
        fa = 1
    if pactive == "on":
        pa = 1
    if siactive == "on":
        sia = 1
    with open(os.path.join(cwd, ".anmrrc"), "w", newline=None) as arc:
        arc.write("7 8 XH acid atoms\n")
        if args.mf is not None:
            arc.write(
                "ENSO qm= {} mf= {} lw= 1.0  J= {} S= {}\n".format(
                    str(prog).upper(), str(mf), calcJ, calcS
                )
            )
        else:
            arc.write("ENSO qm= {} lw= 1.2\n".format(str(prog).upper()))
        try:
            length = max(
                [
                    len(i)
                    for i in [
                        hshielding,
                        cshielding,
                        fshielding,
                        pshielding,
                        sishielding,
                    ]
                ]
            )
        except:
            length = 6
        arc.write(
            "{}[{}] {}[{}]/{}//{}[{}]/{}\n".format(
                href, solv, funcS, sm4, basisS, func, sm, basisfunc
            )
        )
        arc.write(
            "1  {:{digits}}    0.0    {}\n".format(hshielding, ha, digits=length)
        )  # hydrogen
        arc.write(
            "6  {:{digits}}    0.0    {}\n".format(cshielding, ca, digits=length)
        )  # carbon
        arc.write(
            "9  {:{digits}}    0.0    {}\n".format(fshielding, fa, digits=length)
        )  # fluorine
        arc.write(
            "14 {:{digits}}    0.0    {}\n".format(sishielding, sia, digits=length)
        )  # silicon
        arc.write(
            "15 {:{digits}}    0.0    {}\n".format(pshielding, pa, digits=length)
        )  # phosphorus
    return


def splitting(item):
    """Used in move recursively"""
    try:
        return int(item.rsplit(".", 1)[1])
    except ValueError:
        return 0


def move_recursively(path, filename):
    """Check if file or file.x exists and move them to file.x+1
       ignores e.g. file.save"""
    files = [
        f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))
    ]  # list of all files in directory
    newfiles = []  # list of all files in directory that contain filename and '.'
    for item in files:
        if filename + "." in item:
            newfiles.append(item)
    newfiles.sort(key=splitting, reverse=True)
    for item in newfiles:
        try:
            data = item.rsplit(".", 1)  # splits only at last '.'
            int(data[1])
        except ValueError:
            continue
        tmp_from = os.path.join(path, item)
        newfilename = str(data[0]) + "." + str(int(data[1]) + 1)
        tmp_to = os.path.join(path, newfilename)
        print("Backing up {} to {}.".format(item, newfilename))
        shutil.move(tmp_from, tmp_to)

    if filename in files:
        print("Backing up {} to {}.".format(filename, filename + ".1"))
        shutil.move(os.path.join(path, filename), os.path.join(path, filename + ".1"))
    return


class qm_job:
    name = ""
    xtb_energy = None
    rel_xtb_energy = None
    energy = None
    rel_energy = None
    sp3_energy = None
    success = False
    workdir = ""
    func = ""
    solv = None
    sm = ""
    chrg = 0
    jobtype = ""
    full = True
    cycles = 0
    rrho = None
    gsolv = 0.0
    free_energy = None
    rel_free_energy = None
    bm_weight = None
    new_bm_weight = None
    NMR = False
    boltzmann = False
    basis = None
    gfnv = None
    unpaired = 0
    spenergy = None
    hactive = False
    cactive = False
    pactive = False
    siactive = False
    factive = False
    symmetry = "C1"
    temperature = 298.15
    nat = 0
    progsettings = {
        "tempprogpath": "",
        "xtbpath": "",
        "orca_old": "",
        "omp": 1,
        "cosmorssetup": None,
    }

    def execute(self):
        pass

    def _sp(self):
        pass

    def _opt(self):
        pass

    def _genericoutput(self):
        """ READ shielding and coupling constants and write them to plain output"""
        pass

    def _gbsa_rs(self, environsettings):
        """ Calculate GBSA-RS"""
        if not self.boltzmann:
            tmp_gas = 0
            tmp_solv = 0
            print("Running GBSA-RS calculation in " + last_folders(self.workdir, 2))
            if os.path.isfile(os.path.join(self.workdir, "xtbrestart")):
                os.remove(os.path.join(self.workdir, "xtbrestart"))
            with open(
                os.path.join(self.workdir, "gas.out"), "w", newline=None
            ) as outputfile:
                callargs = [
                    self.progsettings["xtbpath"],
                    "coord",
                    "-" + str(self.gfnv),
                    "-sp",
                    "-chrg",
                    str(self.chrg),
                    "--norestart",
                ]
                returncode = subprocess.call(
                    callargs,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                    env=environsettings,
                )
            if returncode is not 0:
                self.gsolv = None
                self.success = False
                print(
                    "ERROR: GFN-xTB error in {:18}".format(
                        last_folders(self.workdir, 2)
                    ),
                    file=sys.stderr,
                )
                return
            with open(
                os.path.join(self.workdir, "solv.out"), "w", newline=None
            ) as outputfile:
                callargs = [
                    self.progsettings["xtbpath"],
                    "coord",
                    "-" + str(self.gfnv),
                    "-sp",
                    "-gbsa",
                    self.solv,
                    "-chrg",
                    str(self.chrg),
                    "--norestart",
                ]
                returncode = subprocess.call(
                    callargs,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                    env=environsettings,
                )
            if returncode is not 0:
                self.gsolv = None
                self.success = False
                print(
                    "ERROR: GFN-xTB error in {:18}".format(
                        last_folders(self.workdir, 2)
                    ),
                    file=sys.stderr,
                )
                return
            time.sleep(0.05)

            if os.path.isfile(os.path.join(self.workdir, "gas.out")):
                with open(
                    os.path.join(self.workdir, "gas.out"),
                    "r",
                    encoding=coding,
                    newline=None,
                ) as inp:
                    store = inp.readlines()
                    for line in store:
                        if "| TOTAL ENERGY" in line:
                            try:
                                tmp_gas = float(line.split()[3])
                                self.success = True
                            except:
                                print(
                                    "Error while converting gas phase "
                                    "single-point in: {}".format(
                                        last_folders(self.workdir, 2)
                                    ),
                                    file=sys.stderr,
                                )
                                tmp_gas = None
                                self.success = False
                            break
                        if (
                            "external code error" in line
                            or "|grad| > 500, something is totally wrong!" in line
                            or "abnormal termination of xtb" in line
                        ):
                            print(
                                "ERROR: GFN-xTB error in {:18}".format(
                                    last_folders(self.workdir, 2)
                                ),
                                file=sys.stderr,
                            )
                            self.gsolv = None
                            self.success = False
                            break
            else:
                print(
                    "WARNING: File {} doesn't exist!".format(
                        os.path.join(self.workdir, "gas.out")
                    )
                )
                self.gsolv = False
                self.success = None
            if os.path.isfile(os.path.join(self.workdir, "solv.out")):
                with open(
                    os.path.join(self.workdir, "solv.out"),
                    "r",
                    encoding=coding,
                    newline=None,
                ) as inp:
                    store = inp.readlines()
                    for line in store:
                        if "| TOTAL ENERGY" in line:
                            try:
                                tmp_solv = float(line.split()[3])
                                self.success = True
                            except:
                                print(
                                    "Error while converting solvation "
                                    "single-point in: {}".format(
                                        last_folders(self.workdir, 2)
                                    ),
                                    file=sys.stderr,
                                )
                                tmp_solv = None
                                self.success = False
                            break
                        if (
                            "external code error" in line
                            or "|grad| > 500, something is totally wrong!" in line
                            or "abnormal termination of xtb" in line
                        ):
                            print(
                                "ERROR: GFN-xTB error in {:18}".format(
                                    last_folders(self.workdir, 2)
                                ),
                                file=sys.stderr,
                            )
                            self.gsolv = None
                            self.success = False
                            break
            else:
                print(
                    "WARNING: File {} doesn't exist!".format(
                        os.path.join(self.workdir, "gas.out")
                    )
                )
                self.gsolv = None
                self.success = False
            if not self.success is False:
                if tmp_solv is None or tmp_gas is None:
                    self.gsolv = None
                    self.success = False
                else:
                    self.gsolv = tmp_solv - tmp_gas
                    self.success = True
            return

    def _xtbrrho(self, environsettings):
        """
        RRHO contribution with GFN-XTB, available both to ORCA and TM
        """
        if not self.boltzmann:
            print("Running xtb RRHO in " + last_folders(self.workdir, 2))
            if os.path.isfile(os.path.join(self.workdir, "xtbrestart")):
                os.remove(os.path.join(self.workdir, "xtbrestart"))
            elif os.path.isfile(os.path.join(self.workdir, "xcontrol-inp")):
                os.remove(os.path.join(self.workdir, "xcontrol-inp"))
            with open(
                os.path.join(self.workdir, "xcontrol-inp"), "w", newline=None
            ) as xcout:
                xcout.write("$thermo\n")
                xcout.write("    temp={}\n".format(self.temperature))
                xcout.write("$end")
            time.sleep(0.05)
            with open(
                os.path.join(self.workdir, "ohess.out"), "w", newline=None
            ) as outputfile:
                if self.solv:
                    callargs = [
                        self.progsettings["xtbpath"],
                        "coord",
                        "-" + str(self.gfnv),
                        "-ohess",
                        "-gbsa",
                        self.solv,
                        "-chrg",
                        str(self.chrg),
                        "-enso",
                        "--norestart",
                        "-I",
                        "xcontrol-inp",
                    ]
                else:
                    callargs = [
                        self.progsettings["xtbpath"],
                        "coord",
                        "-" + str(self.gfnv),
                        "-ohess",
                        "-chrg",
                        str(self.chrg),
                        "-enso",
                        "--norestart",
                        "-I",
                        "xcontrol-inp",
                    ]
                returncode = subprocess.call(
                    callargs,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                    env=environsettings,
                )
            time.sleep(0.05)
            # check if converged:
            if returncode is not 0:
                self.rrho = None
                self.success = False
                print(
                    "ERROR: GFN-xTB ohess error in {:18}".format(
                        last_folders(self.workdir, 2)
                    ),
                    file=sys.stderr,
                )
                return
        if os.path.isfile(os.path.join(self.workdir, "xtb_enso.json")):
            with open(
                os.path.join(self.workdir, "xtb_enso.json"),
                "r",
                encoding=coding,
                newline=None,
            ) as f:
                data = json.load(f)
            if "ZPVE" not in data:
                data["ZPVE"] = 0.0
            if "G(T)" in data:
                if float(self.temperature) == 0:
                    self.rrho = data["ZPVE"]
                    self.success = True
                else:
                    self.rrho = data["G(T)"]
                    self.success = True
                if "point group" in data:
                    self.symmetry = data["point group"]
            else:
                print(
                    "Error while converting rrho in: {}".format(
                        last_folders(self.workdir, 2)
                    ),
                    file=sys.stderr,
                )
                self.rrho = None
                self.success = False
        else:
            print(
                "WARNING: File {} doesn't exist!".format(
                    os.path.join(self.workdir, "xtb_enso.json")
                )
            )
            self.rrho = None
            self.success = False
        return

    def _TMSP(self, environsettings):
        """Turbomole single-point calculation"""
        if not self.boltzmann:
            with open(
                os.path.join(self.workdir, "ridft.out"), "w", newline=None
            ) as outputfile:
                subprocess.call(
                    ["ridft"],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                    env=environsettings,
                )
        time.sleep(0.02)
        # check if scf is converged:
        if os.path.isfile(os.path.join(self.workdir, "ridft.out")):
            with open(
                os.path.join(self.workdir, "ridft.out"),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                stor = inp.readlines()
                if " ENERGY CONVERGED !\n" not in stor:
                    print(
                        "ERROR: scf in {:18} not converged!".format(
                            last_folders(self.workdir, 2)
                        ),
                        file=sys.stderr,
                    )
                    self.success = False
                    self.spenergy = None
                    return 1
        else:
            print(
                "WARNING: File {} doesn't exist!".format(
                    os.path.join(self.workdir, "ridft.out")
                )
            )
            self.success = False
            self.spenergy = None
            return 1
        if os.path.isfile(os.path.join(self.workdir, "energy")):
            with open(
                os.path.join(self.workdir, "energy"), "r", encoding=coding, newline=None
            ) as energy:
                storage = energy.readlines()
            try:
                self.spenergy = float(storage[-2].split()[1])
                self.success = True
            except ValueError:
                print(
                    "ERROR while converting energy in: {:18}".format(
                        last_folders(self.workdir, 2)
                    ),
                    file=sys.stderr,
                )
        else:
            self.spenergy = None
            self.success = False
        return

    def _polyfit(self, n, m, X, Y):
        """ translation of fortran code from cosmothermrd.f90 into python code"""
        try:
            n1 = n + 1
            m1 = m + 1
            m2 = m + 2
            Xc = [None] * m2
            for k in range(0, m2, 1):
                Xc[k] = 0.0
                for i in range(0, n1, 1):
                    Xc[k] = Xc[k] + X[i] ** (k + 1)
            Yc = math.fsum(Y)
            Yx = [None] * m
            for k in range(0, m, 1):
                Yx[k] = 0.0
                for i in range(0, n1, 1):
                    Yx[k] = Yx[k] + Y[i] * X[i] ** (k + 1)
            C = [[None] * m1 for _ in range(m1)]
            for i in range(0, m1, 1):
                for j in range(0, m1, 1):
                    ij = i + j - 1
                    if i == 0 and j == 0:
                        C[0][0] = n1
                    elif ij > len(Xc) - 1:
                        C[i][j] = 0.0
                    else:
                        C[i][j] = Xc[ij]
            B = [None] * m1
            B[0] = Yc
            for i in range(1, m1, 1):
                B[i] = Yx[i - 1]
            for k in range(0, m, 1):
                for i in range(k + 1, m1, 1):
                    B[i] = B[i] - C[i][k] / C[k][k] * B[k]
                    for j in range(k + 1, m1, 1):
                        C[i][j] = C[i][j] - C[i][k] / C[k][k] * C[k][j]
            A = [None] * m1
            A[m1 - 1] = B[m1 - 1] / C[m1 - 1][m1 - 1]
            for i in range(m - 1, -1, -1):
                s = 0.0
                for k in range(i + 1, m1, 1):
                    s = s + C[i][k] * A[k]
                A[i] = (B[i] - s) / C[i][i]

            a0 = []
            for item in A:
                a0.append(item)
        except ZeroDivisionError:
            print("Error in Gsolv!")
            a0 = False
        return a0

    def _solv_complete(self, environsettings):
        """complete COSMO-RS within the ENSO script"""
        if not self.boltzmann:
            print(
                "Running COSMO-RS calculation in {:18}".format(
                    last_folders(self.workdir, 2)
                )
            )
            self.spenergy = None
            # run two single-points:
            cosmo_cefine = {
                "normal": [
                    "cefine",
                    "-chrg",
                    str(self.chrg),
                    "-uhf",
                    str(self.unpaired),
                    "-func",
                    "b-p",
                    "-bas",
                    "def-TZVP",
                    "-noopt",
                    "-grid",
                    " m3",
                    "-sym",
                    "c1",
                    "-novdw",
                ],
                "fine": [
                    "cefine",
                    "-chrg",
                    str(self.chrg),
                    "-uhf",
                    str(self.unpaired),
                    "-func",
                    "b-p",
                    "-bas",
                    "def2-TZVPD",
                    "-noopt",
                    "-grid",
                    " m3",
                    "-sym",
                    "c1",
                    "-novdw",
                ],
                "19-fine": [
                    "cefine",
                    "-chrg",
                    str(self.chrg),
                    "-uhf",
                    str(self.unpaired),
                    "-func",
                    "b-p",
                    "-bas",
                    "def2-TZVPD",
                    "-noopt",
                    "-grid",
                    "m3",
                    "-sym",
                    "c1",
                ],
            }

            old_workdir = self.workdir  # starting workdir
            self.workdir = os.path.join(self.workdir, "COSMO")  # COSMO folder
            mkdir_p(self.workdir)
            shutil.copy(
                os.path.join(old_workdir, "coord"), os.path.join(self.workdir, "coord")
            )
            # parametrization and version:
            if "FINE" in self.progsettings["cosmorssetup"].split()[2]:
                fine = True
                param = "fine"
                if self.progsettings["cosmothermversion"] == 19:
                    param = "19-fine"
            else:
                fine = False
                param = "normal"

            # cefine
            for k in range(2):
                s = subprocess.check_output(
                    cosmo_cefine[param],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                )
                time.sleep(0.15)
                output = s.decode("utf-8").splitlines()
                # checkoutput for errors
                for line in output:
                    if "define ended abnormally" in line:
                        self.success = False
                        return 1
                # check if wrong functional was written by cefine
                with open(
                    os.path.join(self.workdir, "control"),
                    "r",
                    encoding=coding,
                    newline=None,
                ) as control:
                    checkup = control.readlines()
                for line in checkup:
                    if "functional" in line:
                        if "b-p" not in line:
                            print(
                                "Wrong functional in control file"
                                " in {}".format(last_folders(self.workdir, 2))
                            )
                            self.success = False
                        else:
                            break
            # end cefine
            # running single-point in gas phase
            self._TMSP(environsettings)
            if not self.success:
                print("Error in COSMO-RS calculation!")
                return 1
            with open(
                os.path.join(self.workdir, "out.energy"), "w", newline=None
            ) as out:
                out.write(str(self.spenergy) + "\n")
            self.spenergy = None
            # running single-point in ideal conductor!
            with open(
                os.path.join(self.workdir, "control"),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                tmp = inp.readlines()
            with open(os.path.join(self.workdir, "control"), "w", newline=None) as out:
                for line in tmp[:-1]:
                    out.write(line + "\n")
                if not fine:
                    out.write("$cosmo \n")
                    out.write(" epsilon=infinity \n")
                    out.write("$cosmo_out file=out.cosmo \n")
                    out.write("$end \n")
                else:
                    # fine
                    out.write("$cosmo \n")
                    out.write(" epsilon=infinity \n")
                    out.write(" use_contcav \n")
                    out.write(" cavity closed \n")
                    out.write("$cosmo_out file=out.cosmo \n")
                    out.write("$cosmo_isorad \n")
                    out.write("$end \n")
            self._TMSP(environsettings)
            if not self.success:
                print("Error in COSMO-RS calculation during single-point calculation!")
                return 1
            # info from .ensorc # replacement for cosmothermrc
            # fdir=/software/cluster/COSMOthermX16/COSMOtherm/DATABASE-COSMO/BP-TZVP-COSMO autoc
            cosmors_solv = {
                "acetone": "f = propanone.cosmo ",
                "h2o": "f = h2o.cosmo ",
                "chcl3": "f = chcl3.cosmo ",
                "ch2cl2": "f = ch2cl2.cosmo ",
                "acetonitrile": "f = acetonitrile_c.cosmo",
                "dmso": "f = dimethylsulfoxide.cosmo ",
                "methanol": "f = methanol.cosmo ",
                "thf": "f = thf.cosmo ",
                "toluene": "f = toluene_c0.cosmo ",
            }
            if fine:
                solv_data = os.path.join(
                    os.path.split(
                        self.progsettings["cosmorssetup"].split()[5].strip('"')
                    )[0],
                    "DATABASE-COSMO/BP-TZVPD-FINE",
                )
            else:
                solv_data = os.path.join(
                    os.path.split(
                        self.progsettings["cosmorssetup"].split()[5].strip('"')
                    )[0],
                    "DATABASE-COSMO/BP-TZVP-COSMO",
                )
            # test = ['ctd = BP_TZVP_C30_1601.ctd cdir = "/software/cluster/COSMOthermX16/COSMOtherm/CTDATA-FILES"']

            henry = [
                "henry  xh={ 1.0 0.0     }  tc=-50.0 Gsolv",
                "henry  xh={ 1.0 0.0     }  tc=-10.0 Gsolv",
                "henry  xh={ 1.0 0.0     }  tc=0.0  Gsolv",
                "henry  xh={ 1.0 0.0     }  tc=10.0 Gsolv",
                "henry  xh={ 1.0 0.0     }  tc=20.0 Gsolv",
                "henry  xh={ 1.0 0.0     }  tc=25.0 Gsolv",
                "henry  xh={ 1.0 0.0     }  tc=30.0 Gsolv",
                "henry  xh={ 1.0 0.0     }  tc=40.0 Gsolv",
                "henry  xh={ 1.0 0.0     }  tc=50.0 Gsolv",
                "henry  xh={ 1.0 0.0     }  tc=60.0 Gsolv",
            ]
            with open(
                os.path.join(self.workdir, "cosmotherm.inp"), "w", newline=None
            ) as out:
                out.write(self.progsettings["cosmorssetup"] + "\n")
                # write from ensorc
                out.write("EFILE VPFILE \n")
                if self.progsettings["cosmothermversion"] > 16:
                    pass
                else:  # cosmothermX16
                    out.write("\n")  # needs empty line!
                out.write(cosmors_solv[self.solv] + "fdir=" + solv_data + " autoc \n")
                out.write("f = out.cosmo \n")
                for line in henry:
                    out.write(line + "\n")
            # running COSMOtherm
            with open(
                os.path.join(self.workdir, "cosmotherm.out"), "w", newline=None
            ) as outputfile:
                subprocess.call(
                    ["cosmotherm", "cosmotherm.inp"],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                    env=environsettings,
                )
            time.sleep(0.1)
            # get T and Gsolv for version > cosmothermX16
            try:
                with open(
                    os.path.join(self.workdir, "cosmotherm.tab"),
                    "r",
                    encoding=coding,
                    newline=None,
                ) as inp:
                    stor = inp.readlines()
                T = []
                gsolv = []
                for line in stor:
                    if "T=" in line:
                        T.append(float(line.split()[5]))
                    elif " out " in line:
                        gsolv.append(float(line.split()[5]))
            except (FileNotFoundError, ValueError):
                print(
                    "ERROR: cosmotherm.tab was not written, this error can be due to a missing licensefile "
                    "information, or wrong path to the COSMO-RS Database."
                )
                self.gsolv = None
                self.success = False
                return 1
            # cosmothermrd
            if os.stat(os.path.join(self.workdir, "cosmotherm.tab")).st_size == 0:
                print(
                    "ERROR: cosmotherm.tab was not written, this error can be due to a missing licensefile "
                    "information, or wrong path to the COSMO-RS Database."
                )
                self.gsolv = None
                self.success = False
                return 1
            z = []
            for i in range(len(T)):
                z.append(float(gsolv[i]) / float(T[i]))
            a0 = self._polyfit(len(T) - 1, 4, T, gsolv)
            if not a0:
                self.gsolv = None
                self.success = False
                return 1
            temp = float(self.temperature)
            gsolv_out = (
                a0[0]
                + a0[1] * temp
                + a0[2] * temp ** 2
                + a0[3] * temp ** 3
                + a0[4] * temp ** 4
            )
            dh = a0[0]
            ssolv = -(
                a0[1]
                + 2.0 * a0[2] * temp
                + 3.0 * a0[3] * temp ** 2
                + 4.0 * a0[4] * temp ** 3
            )
            hsolv = gsolv_out + temp * ssolv
            a0 = self._polyfit(len(T) - 1, 4, T, z)
            hsolv2 = (
                -temp
                * temp
                * (
                    a0[1]
                    + 2.0 * a0[2] * temp
                    + 3.0 * a0[3] * temp ** 2
                    + 4.0 * a0[4] * temp ** 3
                )
            )
            ssolv2 = (hsolv2 - gsolv_out) / temp

            ## volumework:
            R = 1.987203585e-03  # kcal/(mol*K)
            videal = (
                24.789561955 / 298.15
            )  # molar volume for ideal gas at 298.15 K 100.0 kPa
            volwork = R * temp * math.log(videal * temp)

            self.workdir = old_workdir
            with open(
                os.path.join(self.workdir, "cosmors.out"), "w", newline=None
            ) as out:
                out.write(
                    "This is cosmothermrd (python version in ENSO) (SG,FB,SAW, 06/18)\n"
                )
                out.write("final thermochemical solvation properties in kcal/mol\n")
                out.write(
                    "derived from {} different temperatures\n".format(str(len(T)))
                )
                out.write("big differences between the two values for S and H\n")
                out.write("indicate numerical problems (change cosmotherm T range)\n")
                out.write(
                    "----------------------------------------------------------\n"
                )
                out.write(" Hsolv(0 K,extrapolated)  = {:10.3f}\n".format(dh))
                out.write(
                    " Ssolv({} K)= {:10.5f} {:10.5f}\n".format(temp, ssolv, ssolv2)
                )
                out.write(
                    " Hsolv({} K)= {:10.3f} {:10.3f}\n".format(temp, hsolv, hsolv2)
                )
                out.write(" Gsolv({} K)= {:10.3f}\n".format(temp, gsolv_out))
                out.write(" VWork({} K)= {:10.3f}\n".format(temp, volwork))
                out.write(
                    " Gsolv+VWork({} K)= {:10.3f}\n".format(temp, (gsolv_out + volwork))
                )
            time.sleep(0.05)
            self.gsolv = gsolv_out / 627.50947428
            self.success = True
        else:  # read only output if boltzmann
            if os.path.isfile(os.path.join(self.workdir, "cosmors.out")):
                with open(
                    os.path.join(self.workdir, "cosmors.out"),
                    "r",
                    encoding=coding,
                    newline=None,
                ) as inp:
                    stor = inp.readlines()
                    for line in stor:
                        if " Gsolv(" in line:
                            try:
                                self.gsolv = (
                                    float(line.split()[2]) / 627.50947428
                                )  # with volume work from cosmothermrd
                                self.success = True
                            except:
                                self.success = False
                                self.gsolv = None
                                print(
                                    "\nERROR: could not get Gsolv from COSMO-RS in {:18}!".format(
                                        last_folders(self.workdir, 2)
                                    ),
                                    file=sys.stderr,
                                )
                                return 1
                    if math.isnan(self.gsolv):
                        self.success = False
                    if not self.success:
                        self.gsolv = None
                        print(
                            "\nERROR: COSMO-RS in {:18} not converged!".format(
                                last_folders(self.workdir, 2)
                            ),
                            file=sys.stderr,
                        )
            else:
                print(
                    "WARNING: {} doesn't exist!".format(
                        os.path.join(self.workdir, "cosmors.out")
                    )
                )
                self.success = False
                self.gsolv = None
                return 1
        return
        # END QM_JOB ######################################


class tm_job(qm_job):
    def cefine(self, environsettings):
        """Do cefine for func"""
        removegf = False
        if self.basis == "def2-QZVP(-gf)":
            self.basis = "def2-QZVP"
            removegf = True
        cef_calls = {
            "b97-3c": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "b973c",
                "-bas",
                "def2-mTZVP",
                "-noopt",
                "-grid",
                " m4",
                "-scfconv",
                "6",
                "-sym",
                "c1",
                "-novdw",
            ],
            "pbeh-3c": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "pbeh-3c",
                "-bas",
                "def2-mSVP",
                "-noopt",
                "-grid",
                " m4",
                "-scfconv",
                "6",
                "-sym",
                "c1",
            ],
            "tpss": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "tpss",
                "-fpol",
                "-bas",
                "def2-TZVP",
                "-noopt",
                "-grid",
                " m4",
                "-scfconv",
                "6",
                "-sym",
                "c1",
                "-d3",
            ],
            "pw6b95": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "pw6b95",
                "-bas",
                str(self.basis),
                "-noopt",
                "-grid",
                " m5",
                "-scfconv",
                "7",
                "-sym",
                "c1",
                "-d3",
                "-ri",
            ],
            "pbe0": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "pbe0",
                "-bas",
                str(self.basis),
                "-noopt",
                "-grid",
                " m4",
                "-scfconv",
                "6",
                "-sym",
                "c1",
                "-d3",
            ],
        }
        cef_rrhocalls = {
            "b97-3c": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "b973c",
                "-bas",
                "def2-mTZVP",
                "-noopt",
                "-grid",
                " m5",
                "-scfconv",
                "7",
                "-sym",
                "c1",
                "-novdw",
            ],
            "pbeh-3c": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "pbeh-3c",
                "-bas",
                "def2-mSVP",
                "-noopt",
                "-grid",
                " m5",
                "-scfconv",
                "7",
                "-sym",
                "c1",
            ],
            "tpss": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "tpss",
                "-fpol",
                "-bas",
                "def2-TZVP",
                "-noopt",
                "-grid",
                " m5",
                "-scfconv",
                "7",
                "-sym",
                "c1",
                "-d3",
            ],
        }
        cef_callsnmr = {
            "tpss": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "tpss",
                "-fpol",
                "-bas",
                str(self.basis),
                "-noopt",
                "-grid",
                " 4",
                "-scfconv",
                "7",
                "-sym",
                "c1",
                "-d3",
            ],
            "pbe0": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "pbe0",
                "-fpol",
                "-bas",
                str(self.basis),
                "-noopt",
                "-grid",
                " 4",
                "-scfconv",
                "7",
                "-sym",
                "c1",
                "-d3",
            ],
            "pbeh-3c": [
                "cefine",
                "-chrg",
                str(self.chrg),
                "-func",
                "pbeh-3c",
                "-fpol",
                "-bas",
                "def2-mSVP",
                "-noopt",
                "-grid",
                " 4",
                "-scfconv",
                "7",
                "-sym",
                "c1",
                "-d3",
            ],
        }
        # remove -fg functions from def2-QZVP basis set
        if removegf:
            cef_calls[self.func] = cef_calls[self.func] + ["-gf"]
        for k in range(2):
            if self.jobtype == "rrhotm":
                # use higher grid
                s = subprocess.check_output(
                    cef_rrhocalls[self.func],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                )
            else:
                if self.NMR:
                    s = subprocess.check_output(
                        cef_callsnmr[self.func],
                        shell=False,
                        stdin=None,
                        stderr=subprocess.STDOUT,
                        universal_newlines=False,
                        cwd=self.workdir,
                    )
                elif self.unpaired > 0:
                    s = subprocess.check_output(
                        cef_calls[self.func] + ["-uhf", str(self.unpaired)],
                        shell=False,
                        stdin=None,
                        stderr=subprocess.STDOUT,
                        universal_newlines=False,
                        cwd=self.workdir,
                    )
                elif self.unpaired == 0:
                    s = subprocess.check_output(
                        cef_calls[self.func],
                        shell=False,
                        stdin=None,
                        stderr=subprocess.STDOUT,
                        universal_newlines=False,
                        cwd=self.workdir,
                    )

            time.sleep(0.15)
            output = s.decode("utf-8").splitlines()
            # checkoutput for errors
            for line in output:
                if "define ended abnormally" in line:
                    self.success = False
                    return 1
            # check if wrong functional was written by cefine
            with open(
                os.path.join(self.workdir, "control"),
                "r",
                encoding=coding,
                newline=None,
            ) as control:
                checkup = control.readlines()
            for line in checkup:
                if "functional" in line:
                    if self.func == "b97-3c":
                        testfunc = "b973c"
                    else:
                        testfunc = self.func
                    if testfunc not in line:
                        print(
                            "Wrong functional in control file"
                            " in {}".format(last_folders(self.workdir, 2))
                        )
                        self.success = False
                    else:
                        self.success = True
                        break
        # continue with modifications to control
        solvent_dcosmors = {
            "acetone": [
                "$cosmo",
                " epsilon= 20.7",
                " cavity closed",
                " use_contcav",
                "$dcosmo_rs file=propanone_25.pot",
            ],
            "chcl3": [
                "$cosmo",
                " epsilon= 4.8",
                " cavity closed",
                " use_contcav",
                "$dcosmo_rs file=chcl3_25.pot",
            ],
            "acetonitrile": [
                "$cosmo",
                " epsilon= 36.6",
                " cavity closed",
                " use_contcav",
                "$dcosmo_rs file=acetonitrile_25.pot",
            ],
            "ch2cl2": [
                "$cosmo",
                " epsilon= 9.1",
                " cavity closed",
                " use_contcav",
                "$dcosmo_rs file=chcl3_25.pot",
            ],
            "dmso": [
                "$cosmo",
                " epsilon= 47.2",
                " cavity closed",
                " use_contcav",
                "$dcosmo_rs file=dimethylsulfoxide_25.pot",
            ],
            "h2o": [
                "$cosmo",
                " epsilon= 80.1",
                " cavity closed",
                " use_contcav",
                "$dcosmo_rs file=h2o_25.pot",
            ],
            "methanol": [
                "$cosmo",
                " epsilon= 32.7",
                " cavity closed",
                " use_contcav",
                "$dcosmo_rs file=methanol_25.pot",
            ],
            "thf": [
                "$cosmo",
                " epsilon= 7.6",
                " cavity closed",
                " use_contcav",
                "$dcosmo_rs file=thf_25.pot",
            ],
            "toluene": [
                "$cosmo",
                " epsilon= 2.4",
                " cavity closed",
                " use_contcav",
                "$dcosmo_rs file=toluene_25.pot",
            ],
        }
        solvent_cosmo = {
            "acetone": ["$cosmo", " epsilon= 20.7", " cavity closed", " use_contcav"],
            "chcl3": ["$cosmo", " epsilon= 4.8", " cavity closed", " use_contcav"],
            "acetonitrile": [
                "$cosmo",
                " epsilon= 36.6",
                " cavity closed",
                " use_contcav",
            ],
            "ch2cl2": ["$cosmo", " epsilon= 9.1", " cavity closed", " use_contcav"],
            "dmso": ["$cosmo", " epsilon= 47.2", " cavity closed", " use_contcav"],
            "h2o": ["$cosmo", " epsilon= 80.1", " cavity closed", " use_contcav"],
            "methanol": ["$cosmo", " epsilon= 32.7", " cavity closed", " use_contcav"],
            "thf": ["$cosmo", " epsilon= 7.6", " cavity closed", " use_contcav"],
            "toluene": ["$cosmo", " epsilon= 2.4", " cavity closed", " use_contcav"],
        }
        # add 3 body correction (atm) for dispersion in b97-3c
        if self.func == "b97-3c":
            with open(
                os.path.join(self.workdir, "control"),
                "r",
                encoding=coding,
                newline=None,
            ) as control:
                tmp = control.readlines()
            with open(
                os.path.join(self.workdir, "control"), "w", newline=None
            ) as newcontrol:
                for line in tmp[:-1]:  # works because of novdw
                    newcontrol.write(line)
                newcontrol.write("$disp3 -bj -abc\n")
                newcontrol.write("$end")
        # add solvent part and NMR part to control file if necessary
        if self.solv and self.sm not in ["gas", None]:
            with open(
                os.path.join(self.workdir, "control"),
                "r",
                encoding=coding,
                newline=None,
            ) as control:
                tmp = control.readlines()
            with open(
                os.path.join(self.workdir, "control"), "w", newline=None
            ) as newcontrol:
                for line in tmp[:-1]:  # change thresholds if DCOSMO-RS is used
                    newcontrol.write(line)
                if self.sm == "dcosmors":
                    for line in solvent_dcosmors[self.solv]:
                        newcontrol.write(line + " \n")
                if self.sm == "cosmo":
                    for line in solvent_cosmo[self.solv]:
                        newcontrol.write(line + " \n")
                newcontrol.write("$end")
        # add NMR part if cefine for part4
        if self.NMR:
            numactive = 0
            if self.hactive:
                numactive += 1
            if self.cactive:
                numactive += 1
            if self.factive:
                numactive += 1
            if self.pactive:
                numactive += 1
            if self.siactive:
                numactive += 1
            with open(
                os.path.join(self.workdir, "control"),
                "r",
                encoding=coding,
                newline=None,
            ) as control:
                tmp = control.readlines()
            with open(
                os.path.join(self.workdir, "control"), "w", newline=None
            ) as newcontrol:
                for line in tmp:
                    if "rpacor" in line:
                        tmp[tmp.index(line)] = "$rpacor 4000 \n"
                for line in tmp[:-1]:
                    newcontrol.write(line)
                newcontrol.write("$ncoupling\n")
                newcontrol.write(" fconly\n")
                # newcontrol.write(' sdonly\n')
                # newcontrol.write(' psoonly\n')
                newcontrol.write(" thr=0.1\n")
                if numactive == 1:
                    if self.hactive or self.pactive:
                        newcontrol.write('$nucsel "h" "f" "p" \n')
                        newcontrol.write('$nucsel2 "h" "f" "p" \n')
                    elif self.cactive:
                        newcontrol.write('$nucsel "c" "f"  \n')
                        newcontrol.write('$nucsel2 "c" "f"  \n')
                    elif self.factive:
                        newcontrol.write('$nucsel "h" "c" "f" "p"\n')
                        newcontrol.write('$nucsel2 "h" "c" "f" "p"\n')
                else:
                    newcontrol.write('$nucsel "h" "c" "f" "p" "si"\n')
                    newcontrol.write('$nucsel2 "h" "c" "f" "p" "si"\n')
                newcontrol.write("$rpaconv 8\n")
                newcontrol.write("$end")
        time.sleep(0.15)
        return 0

    def _sp(self, environsettings):
        """Turbomole single-point calculation"""
        if not self.boltzmann:
            print("Running single-point in {:18}".format(last_folders(self.workdir, 2)))
            with open(
                os.path.join(self.workdir, "ridft.out"), "w", newline=None
            ) as outputfile:
                subprocess.call(
                    ["ridft"],
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                    env=environsettings,
                )
        time.sleep(0.02)
        # check if scf is converged:
        if os.path.isfile(os.path.join(self.workdir, "ridft.out")):
            with open(
                os.path.join(self.workdir, "ridft.out"),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                stor = inp.readlines()
                if " ENERGY CONVERGED !\n" not in stor:
                    print(
                        "ERROR: scf in {:18} not converged!".format(
                            last_folders(self.workdir, 2)
                        ),
                        file=sys.stderr,
                    )
                    self.success = False
                    self.energy = None
                    return 1
        else:
            print(
                "WARNING: {} doesn't exist!".format(
                    os.path.join(self.workdir, "ridft.out")
                )
            )
            self.success = False
            self.energy = None
            return 1
        if os.path.isfile(os.path.join(self.workdir, "energy")):
            with open(
                os.path.join(self.workdir, "energy"), "r", encoding=coding, newline=None
            ) as energy:
                storage = energy.readlines()
            try:
                self.energy = float(storage[-2].split()[1])
                self.success = True
            except ValueError:
                print(
                    "ERROR while converting energy in: {:18}".format(
                        last_folders(self.workdir, 2)
                    ),
                    file=sys.stderr,
                )
        else:
            self.energy = None
            self.success = False
        return

    def _xtbopt(self, environsettings):
        """Turbomole optimization using ANCOPT implemented in GFN-xTB"""
        if self.full:
            output = "opt-part2.out"
        else:
            output = "opt-part1.out"
        if not self.boltzmann:
            print("Running optimization in {:18}".format(last_folders(self.workdir, 2)))
            if os.path.isfile(os.path.join(self.workdir, "xtbrestart")):
                os.remove(os.path.join(self.workdir, "xtbrestart"))
            if self.full:
                if self.sm == "dcosmors" and self.solv:
                    callargs = [
                        self.progsettings["xtbpath"],
                        "coord",
                        "-opt",
                        "lax",
                        "-tm",
                    ]
                else:  # gas phase
                    callargs = [self.progsettings["xtbpath"], "coord", "-opt", "-tm"]
            else:  # only crude
                callargs = [
                    self.progsettings["xtbpath"],
                    "coord",
                    "-opt",
                    "crude",
                    "-tm",
                ]

            with open(
                os.path.join(self.workdir, output), "w", newline=None
            ) as outputfile:
                returncode = subprocess.call(
                    callargs,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                    env=environsettings,
                )
            if returncode is not 0:
                self.energy = None
                self.success = False
                print(
                    "ERROR: optimization in {:18} not converged".format(
                        last_folders(self.workdir, 2)
                    ),
                    file=sys.stderr,
                )
                return
            time.sleep(0.02)
        # check if converged:
        if os.path.isfile(os.path.join(self.workdir, output)):
            with open(
                os.path.join(self.workdir, output), "r", encoding=coding, newline=None
            ) as inp:
                stor = inp.readlines()
                for line in stor:
                    if (
                        "external code error" in line
                        or "|grad| > 500, something is totally wrong!" in line
                        or "abnormal termination of xtb" in line
                    ):
                        print(
                            "ERROR: optimization in {:18} not converged".format(
                                last_folders(self.workdir, 2)
                            ),
                            file=sys.stderr,
                        )
                        self.success = False
                        self.energy = None
                        return 1
                    if "   *** GEOMETRY OPTIMIZATION CONVERGED AFTER " in line:
                        self.cycles = int(line.split()[5])
        else:
            print(
                "WARNING: {} doesn't exist!".format(os.path.join(self.workdir, output))
            )
            self.success = False
            self.energy = None
            return 1
        if os.path.isfile(os.path.join(self.workdir, "energy")):
            with open(
                os.path.join(self.workdir, "energy"), "r", encoding=coding, newline=None
            ) as energy:
                storage = energy.readlines()
            try:
                self.energy = float(storage[-2].split()[1])
                self.success = True
            except ValueError:
                print(
                    "ERROR while converting energy in {:18}".format(
                        last_folders(self.workdir, 2)
                    ),
                    file=sys.stderr,
                )
        else:
            self.energy = None
            self.success = False
        return

    def _opt(self, environsettings):
        """Turbomole optimization using JOBEX, using adapted thresholds!"""

        # part2 = full optimization at either: normal (Econv/Eh = 5 × 10⁻⁶; Gconv/Eh·α⁻¹ = 1 × 10⁻³) or
        #                     lax if dcosmors is used (Econv/Eh = 2 × 10⁻⁵; Gconv/Eh·α⁻¹ = 2 × 10⁻³)
        # part1 = crude optimization : crude =        (Econv/Eh = 5 × 10⁻⁴; Gconv/Eh·α⁻¹ = 1 × 10⁻²)

        normal = {
            "threchange": "   threchange  5.0d-6 \n",
            "thrrmsgrad": "   thrrmsgrad  1.0d-3 \n",
            "thrmaxdispl": "   thrmaxdispl  1.0d-1 \n",
            "thrrmsdispl": "   thrrmsdispl  1.0d-1 \n",
        }
        lax = {
            "threchange": "   threchange  2.0d-5 \n",
            "thrrmsgrad": "   thrrmsgrad  2.0d-3 \n",
            "thrmaxdispl": "   thrmaxdispl  1.0d-1 \n",
            "thrrmsdispl": "   thrrmsdispl  1.0d-1 \n",
        }
        crude = {
            "threchange": "   threchange  5.0d-4 \n",
            "thrrmsgrad": "   thrrmsgrad  1.0d-2 \n",
            "thrmaxdispl": "   thrmaxdispl  1.0d-1 \n",
            "thrrmsdispl": "   thrrmsdispl  1.0d-1 \n",
        }

        if self.full:
            output = "jobex-part2.out"
            callargs = ["jobex", "-c", "200"]
            if self.sm == "dcosmors":
                thresholds = lax
            else:
                thresholds = normal
        else:
            output = "jobex-part1.out"
            callargs = ["jobex", "-c", "100"]
            thresholds = crude
        if not self.boltzmann:
            print("Running optimization in {:18}".format(last_folders(self.workdir, 2)))
            with open(
                os.path.join(self.workdir, "control"),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                stor = inp.readlines()
            with open(
                os.path.join(self.workdir, "control"), "w", newline=None
            ) as newcontrol:
                for line in stor:
                    if any([key in line for key in thresholds.keys()]):
                        for key in thresholds.keys():
                            if key in line:
                                newcontrol.write(thresholds[key])
                    else:
                        newcontrol.write(line + " \n")
            with open(
                os.path.join(self.workdir, output), "w", newline=None
            ) as outputfile:
                subprocess.call(
                    callargs,
                    shell=False,
                    stdin=None,
                    stderr=None,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                    env=environsettings,
                )

            time.sleep(0.02)
            # check if scf is converged:
        if os.path.isfile(os.path.join(self.workdir, "job.last")) and os.path.isfile(
            os.path.join(self.workdir, "energy")
        ):
            with open(
                os.path.join(self.workdir, "job.last"),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                stor = inp.readlines()
                for line in stor:
                    if "                 |  total energy      =" in line:
                        self.energy = float(line.split()[4])
                    if "CONVERGENCY CRITERIA FULFILLED IN CYCLE" in line:
                        self.success = True
            self.cycles = (
                sum(1 for line in open(os.path.join(self.workdir, "energy"))) - 2
            )
        else:
            print(
                "WARNING: {} or {} doesn't exist!".format(
                    os.path.join(self.workdir, "job.last"),
                    os.path.join(self.workdir, "energy"),
                )
            )
            self.success = False
            self.energy = None
            self.cycles = 0
            return 1
        if not self.full:
            self.success = True
        if self.energy is None and self.success:
            self.success = False
            print(
                "ERROR: optimization in {:18} not converged!".format(
                    last_folders(self.workdir, 2)
                ),
                file=sys.stderr,
            )
            return 1
        if not self.success:
            self.energy = None
            print(
                "ERROR: optimization in {:18} not converged!".format(
                    last_folders(self.workdir, 2)
                ),
                file=sys.stderr,
            )
            return 1
        return

    def _rrho(self, environsettings):
        """Turbomole RRHO correction in the gas phase"""
        thermo_dict = {
            "b97-3c": ["100", self.temperature, "1.0"],
            "tpss": ["100", self.temperature, "1.0"],
            "pbeh-3c": ["100", self.temperature, "0.95"],
        }
        opt = False  # parameter to check whether reoptimization is required
        if self.boltzmann:  # only for boltzmann evaluation
            time.sleep(0.01)
            # check if aoforce.out exists
            if os.path.isfile(os.path.join(self.workdir, "aoforce.out")):
                with open(
                    os.path.join(self.workdir, "thermo.out"),
                    "r",
                    encoding=coding,
                    newline=None,
                ) as inp:
                    stor = inp.readlines()
                    for line in stor:
                        if "G(T)           " in line:
                            try:
                                self.rrho = float(line.split()[1])  # a.u.
                                self.success = True
                            except:
                                print(
                                    "Error while converting energy in: {}".format(
                                        last_folders(self.workdir, 2)
                                    ),
                                    file=sys.stderr,
                                )
                                self.rrho = None
                                self.success = False
                            break
            else:
                print(
                    "ERROR: could not read aoforce.out in {}!".format(
                        last_folders(self.workdir, 2)
                    )
                )
                self.rrho = None
                self.success = False
                ## endif boltzmann
        elif self.solv and not self.boltzmann:
            self.solv = None
            opt = True
        if not self.boltzmann:  # only gas phase frequencies implemented!
            self.cefine(environsettings)
            # if optimization was carried out in solution,
            # need new optimization in gas phase
            if opt:
                self.full = True
                self._opt(environsettings)
            else:  # if optimization was carried out in gas phase,
                # coord is sufficient without reoptimization
                self._sp(environsettings)
                # rdgrad
                print("Running rdgrad in {}".format(last_folders(self.workdir, 2)))
                with open(
                    os.path.join(self.workdir, "rdgrad.out"), "w", newline=None
                ) as outputfile:
                    subprocess.call(
                        ["rdgrad"],
                        shell=False,
                        stdin=None,
                        stderr=None,
                        universal_newlines=False,
                        cwd=self.workdir,
                        stdout=outputfile,
                        env=environsettings,
                    )
                time.sleep(0.02)
                with open(
                    os.path.join(self.workdir, "rdgrad.out"),
                    "r",
                    encoding=coding,
                    newline=None,
                ) as inp:
                    stor = inp.readlines()
                    if (
                        "     --- calculation of the energy gradient finished ---\n"
                        not in stor
                    ):
                        print(
                            "ERROR: rdgrad calculation in {:18} not converged!".format(
                                last_folders(self.workdir, 2)
                            ),
                            file=sys.stderr,
                        )
                        self.success = False
                        return 1
            # AOFORCE
            print("Running aoforce in {}".format(last_folders(self.workdir, 2)))
            with open(
                os.path.join(self.workdir, "aoforce.out"), "w", newline=None
            ) as outputfile:
                subprocess.call(
                    ["aoforce", "-smpcpus", str(self.progsettings["omp"])],
                    shell=False,
                    stdin=None,
                    stderr=None,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                    env=environsettings,
                )
            time.sleep(0.05)
            # run thermo
            if os.path.isfile(os.path.join(self.workdir, "aoforce.out")):
                with open(
                    os.path.join(self.workdir, "thermo.out"), "w", newline=None
                ) as outputfile:
                    sthr, temp, scalefactor = thermo_dict[self.func]
                    subprocess.call(
                        ["thermo", sthr, temp, scalefactor],
                        shell=False,
                        stdin=None,
                        stderr=None,
                        universal_newlines=False,
                        cwd=self.workdir,
                        stdout=outputfile,
                        env=environsettings,
                    )
                with open(
                    os.path.join(self.workdir, "thermo.out"),
                    "r",
                    encoding=coding,
                    newline=None,
                ) as inp:
                    stor = inp.readlines()
                for line in stor:
                    if "G(T)           " in line:
                        try:
                            self.rrho = float(line.split()[1])  # a.u.
                            self.success = True
                        except:
                            print(
                                "Error while converting energy in: {}".format(
                                    last_folders(self.workdir, 2)
                                ),
                                file=sys.stderr,
                            )
                            self.rrho = None
                            self.success = False
                        break
            else:  # Aoforce not found
                print("ERROR: aoforce output not found!")
                self.rrho = None
                self.success = False
        return

    def _nmrJ(self, environsettings):
        """ TM NMR coupling calculation"""
        print(
            "Running couplings calculation in {}".format(last_folders(self.workdir, 2))
        )
        # print('Using: {}'.format(escfpath))
        with open(
            os.path.join(self.workdir, "escf.out"), "w", newline=None
        ) as outputfile:
            subprocess.call(
                [
                    self.progsettings["tempprogpath"],
                    "-smpcpus",
                    str(self.progsettings["omp"]),
                ],
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.workdir,
                stdout=outputfile,
                # env=environsettings,
            )
        time.sleep(0.02)
        # check for convergence
        with open(
            os.path.join(self.workdir, "escf.out"), "r", encoding=coding, newline=None
        ) as inp:
            stor = inp.readlines()
            if "   ****  escf : all done  ****\n" in stor:
                self.success = True
            else:
                print(
                    "ERROR: coupling calculation failed in {:18}".format(
                        last_folders(self.workdir, 1)
                    ),
                    file=sys.stderr,
                )
                self.success = False
        return

    def _nmrS(self, environsettings):
        """TM NMR shielding calculation"""
        print(
            "Running shielding calculation in {}".format(last_folders(self.workdir, 2))
        )
        # print('Using: {}'.format(mpshiftpath))
        with open(
            os.path.join(self.workdir, "mpshift.out"), "w", newline=None
        ) as outputfile:
            subprocess.call(
                [
                    self.progsettings["tempprogpath"],
                    "-smpcpus",
                    str(self.progsettings["omp"]),
                ],
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.workdir,
                stdout=outputfile,
                env=environsettings,
            )
            time.sleep(0.02)
            # check if shift calculation is converged:
            with open(
                os.path.join(self.workdir, "mpshift.out"),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                stor = inp.readlines()
                if (
                    "   ***  nmr shielding constants written onto general input/output file!  ***\n"
                    in stor
                ):
                    self.success = True
                else:
                    print(
                        "ERROR: shielding calculation failed in {:18}".format(
                            last_folders(self.workdir, 1)
                        ),
                        file=sys.stderr,
                    )
                    self.success = False
        return

    def _genericoutput(self):
        """ READ shielding and coupling constants and write them to plain output"""
        fnameshield = "mpshift.out"
        atom = []
        sigma = []
        try:
            with open(
                os.path.join(self.workdir, fnameshield),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                data = inp.readlines()
            for line in data:
                if ">>>>> DFT MAGNETIC SHIELDINGS <<<<<" in line:
                    start = data.index(line)
            for line in data[start:]:
                if "ATOM" in line:
                    splitted = line.split()
                    atom.append(int(splitted[2]))
                    sigma.append(float(splitted[4]))

        except FileNotFoundError:
            print(
                "Missing file: {} in {}".format(
                    fnameshield, last_folders(self.workdir, 2)
                )
            )
            self.success = False
        self.success = True
        fnamecoupl = "escf.out"
        atom1 = []
        atom2 = []
        jab = []
        try:
            with open(
                os.path.join(self.workdir, fnamecoupl),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                data = inp.readlines()
            for line in data:
                if "Nuclear coupling constants in Hz" in line:
                    start = int(data.index(line)) + 3
                if "-----------------------------------" in line:
                    end = int(data.index(line))
            for line in data[start:end]:
                if len(line.split()) == 6:
                    splitted = line.split()
                    atom1.append(int(splitted[1]))
                    atom2.append(int(splitted[4].split(":")[0]))
                    jab.append(float(splitted[5]))
        except FileNotFoundError:
            print(
                "Missing file: {} in {}".format(
                    fnamecoupl, last_folders(self.workdir, 2)
                )
            )
            self.success = False
        self.success = True

        with open(os.path.join(self.workdir, "nmrprop.dat"), "w", newline=None) as out:
            for i in range(len(atom)):
                out.write("{:{digits}} {}\n".format(atom[i], sigma[i], digits=4))
            for i in range(1, self.nat - len(atom)):
                out.write("\n")
            for i in range(len(atom1)):
                out.write(
                    "{:{digits}} {:{digits}}   {}\n".format(
                        atom1[i], atom2[i], jab[i], digits=4
                    )
                )
        time.sleep(0.02)
        return

    def execute(self):
        """Choose what to execute for the jobtype"""
        if self.jobtype == "prep":
            if self.boltzmann:
                self.success = True
                pass
            else:
                self.cefine(environsettings)
                print("{}  ".format(self.name))
        elif self.jobtype == "sp":
            self._sp(environsettings)
        elif self.jobtype == "xtbopt":
            self._xtbopt(environsettings)
        elif self.jobtype == "opt":
            self._opt(environsettings)
        elif self.jobtype == "rrhoxtb":
            self._xtbrrho(environsettings)
        elif self.jobtype == "rrhotm":
            self._rrho(environsettings)
        elif self.jobtype == "solv":
            self._solv_complete(environsettings)
        elif self.jobtype == "gbsa_gsolv":
            self._gbsa_rs(environsettings)
        elif self.jobtype == "nmrJ":
            self._nmrJ(environsettings)
        elif self.jobtype == "nmrS":
            self._nmrS(environsettings)
        elif self.jobtype == "genericout":
            self._genericoutput()


class orca_job(qm_job):
    coord = []
    solvent_smd_new = {
        "acetone": ["%cpcm", " smd     true", ' smdsolvent "acetone"', "end"],
        "chcl3": ["%cpcm", " smd     true", ' smdsolvent "chloroform"', "end"],
        "acetonitrile": ["%cpcm", " smd     true", ' smdsolvent "ACETONITRILE"', "end"],
        "ch2cl2": ["%cpcm", " smd     true", ' smdsolvent "ch2cl2"', "end"],
        "dmso": ["%cpcm", " smd     true", ' smdsolvent "dmso"', "end"],
        "h2o": ["%cpcm", " smd     true", ' smdsolvent "water"', "end"],
        "methanol": ["%cpcm", " smd     true", ' smdsolvent "methanol"', "end"],
        "thf": ["%cpcm", " smd     true", ' smdsolvent "thf"', "end"],
        "toluene": ["%cpcm", " smd     true", ' smdsolvent "toluene"', "end"],
    }

    solvent_smd_old = {
        "acetone": ["%cpcm", " smd     true", ' solvent "acetone"', "end"],
        "chcl3": ["%cpcm", " smd     true", ' solvent "chloroform"', "end"],
        "acetonitrile": ["%cpcm", " smd     true", ' solvent "ACETONITRILE"', "end"],
        "ch2cl2": ["%cpcm", " smd     true", ' solvent "ch2cl2"', "end"],
        "dmso": ["%cpcm", " smd     true", ' solvent "dmso"', "end"],
        "h2o": ["%cpcm", " smd     true", ' solvent "water"', "end"],
        "methanol": ["%cpcm", " smd     true", ' solvent "methanol"', "end"],
        "thf": ["%cpcm", " smd     true", ' solvent "thf"', "end"],
        "toluene": ["%cpcm", " smd     true", ' solvent "toluene"', "end"],
    }

    solvent_cpcm = {
        "acetone": ["!CPCM(Acetone)"],
        "chcl3": ["!CPCM(Chloroform)"],
        "acetonitrile": ["!CPCM(Acetonitrile)"],
        "ch2cl2": ["!CPCM(CH2Cl2)"],
        "dmso": ["!CPCM(DMSO)"],
        "h2o": ["!CPCM(Water)"],
        "methanol": ["!CPCM(Methanol)"],
        "thf": ["!CPCM(THF)"],
        "toluene": ["!CPCM(Toluene)"],
    }

    def xyz2coord(self):
        """convert file inp.xyz to TURBOMOLE coord file"""
        time.sleep(0.1)
        bohr2ang = 0.52917721067
        with open(
            os.path.join(self.workdir, "inp.xyz"), "r", encoding=coding, newline=None
        ) as f:
            xyz = f.readlines()
            atom = []
            x = []
            y = []
            z = []
            # natoms = xyz[0]     # unused variable
            for line in xyz[2:]:
                atom.append(str(line.split()[0].lower()))
                x.append(float(line.split()[1]) / bohr2ang)
                y.append(float(line.split()[2]) / bohr2ang)
                z.append(float(line.split()[3]) / bohr2ang)
            coordxyz = []
            for i in range(len(x)):
                coordxyz.append(
                    "{: 09.7f} {: 09.7f}  {: 09.7f}  {}".format(
                        x[i], y[i], z[i], atom[i]
                    )
                )
            with open(os.path.join(self.workdir, "coord"), "w", newline=None) as coord:
                coord.write("$coord\n")
                for line in coordxyz:
                    coord.write(line + "\n")
                coord.write("$end")
        return

    def _sp(self):
        """ORCA input generation and single-point calculation"""
        if not self.boltzmann:
            # convert coord to xyz
            self.coord, self.nat = coord2xyz(self.workdir)
            # input generation
            osp_calls = {
                "b97-3c": [
                    "%MaxCore 8000",
                    "! def2-mTZVP b97-3c grid4",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pbeh-3c": [
                    "%MaxCore 8000",
                    "! def2-mSVP pbeh-3c grid3",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "tpss": [
                    "%MaxCore 8000",
                    "! def2-TZVP(-f) tpss  grid4 d3bj",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pw6b95": [
                    "%MaxCore 8000",
                    "! "
                    + str(self.basis)
                    + " pw6b95 grid5 rijcosx def2/j d3bj loosescf nososcf gridx4",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pbe0": [
                    "%MaxCore 8000",
                    "! "
                    + str(self.basis)
                    + " pbe0 grid5 rijcosx def2/j d3bj loosescf nososcf gridx4",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "dsd-blyp": [
                    "%MaxCore 8000",
                    "! def2-TZVPP def2-TZVPP/c ri-dsd-blyp frozencore grid5 rijk def2/jk d3bj",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "wb97x": [
                    "%MaxCore 8000",
                    "! "
                    + str(self.basis)
                    + " wb97x-d3 grid5 rijcosx def2/j loosescf nososcf gridx4",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
            }
            with open(os.path.join(self.workdir, "inp"), "w", newline=None) as inp:
                # functionals
                for line in osp_calls[self.func]:
                    inp.write(line + "\n")
                # nprocs
                if int(self.progsettings["omp"]) >= 1:
                    inp.write("%pal\n")
                    inp.write("    nprocs {}\n".format(self.progsettings["omp"]))
                    inp.write("end\n")
                # add solvent correction to input if solvent is specified and using smd
                if self.solv and self.sm == "smd":
                    if self.progsettings["orca_old"]:
                        for line in self.solvent_smd_old[self.solv]:
                            inp.write(line + "\n")
                    else:
                        for line in self.solvent_smd_new[self.solv]:
                            inp.write(line + "\n")
                if self.solv and self.sm == "cpcm":
                    for line in self.solvent_cpcm[self.solv]:
                        inp.write(line + "\n")
                # unpaired, charge, and coordinates
                inp.write(
                    "*xyz {} {}\n".format(str(self.chrg), str(self.unpaired + 1))
                )  # unpaired == number of unpaired electrons
                for line in self.coord:
                    inp.write(line + "\n")
                inp.write("*")
            # Done writing input!
            time.sleep(0.02)
            print("Running single-point in {}".format(last_folders(self.workdir, 2)))
            # start SP calculation
            with open(
                os.path.join(self.workdir, "sp.out"), "w", newline=None
            ) as outputfile:
                call = [os.path.join(self.progsettings["tempprogpath"], "orca"), "inp"]
                subprocess.call(
                    call,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                )
            time.sleep(0.1)
        # check if scf is converged:
        self.success = False
        if os.path.isfile(os.path.join(self.workdir, "sp.out")):
            with open(
                os.path.join(self.workdir, "sp.out"), "r", encoding=coding, newline=None
            ) as inp:
                stor = inp.readlines()
                fspe = "FINAL SINGLE POINT ENERGY"
                for line in stor:
                    if fspe in line:
                        self.energy = float(line.split()[4])
                    if "ORCA TERMINATED NORMALLY" in line:
                        self.success = True
        else:
            self.energy = None
            self.success = False
            print(
                "WARNING: {} doesn't exist!".format(
                    os.path.join(self.workdir, "sp.out")
                )
            )
        if self.energy is None and self.success:
            self.success = False
            print(
                "ERROR: scf in {} not converged!".format(last_folders(self.workdir, 2)),
                file=sys.stderr,
            )
        if not self.success:
            self.energy = None
            print(
                "ERROR: in single point calculation {}!".format(
                    last_folders(self.workdir, 2)
                ),
                file=sys.stderr,
            )
        return

    def _xtbopt(self, environsettings):
        """
        ORCA input generation and geometry optimization using ANCOPT
        implemented within xtb, generates inp.xyz, inp (orca-input) 
        and adds information to coord (xtb can then tell which file 
        orca has to use).
        """
        if self.full:
            output = "opt-part2.out"
        else:
            output = "opt-part1.out"
        if not self.boltzmann:
            # convert coord to xyz, write inp.xyz
            self.coord, self.nat = coord2xyz(self.workdir)
            with open(os.path.join(self.workdir, "inp.xyz"), "w", newline=None) as inp:
                inp.write("{}\n\n".format(self.nat))
                for line in self.coord:
                    inp.write(line + "\n")
            # add input file to coord
            with open(os.path.join(self.workdir, "coord"), "r", newline=None) as coord:
                tmp = coord.readlines()
            with open(
                os.path.join(self.workdir, "coord"), "w", newline=None
            ) as newcoord:
                for line in tmp[:-1]:
                    newcoord.write(line)
                newcoord.write("$external\n")
                newcoord.write("   orca input file= inp\n")
                newcoord.write(
                    "   orca bin= {} \n".format(
                        os.path.join(self.progsettings["tempprogpath"], "orca")
                    )
                )
                newcoord.write("$end")
            # input generation
            ogo_calls = {
                "b97-3c": [
                    "%MaxCore 8000",
                    "! def2-mTZVP b97-3c grid4",
                    "!ENGRAD",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pbeh-3c": [
                    "%MaxCore 8000",
                    "! def2-mSVP pbeh-3c grid4",
                    "!ENGRAD",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "tpss": [
                    "%MaxCore 8000",
                    "! def2-TZVP(-f) tpss grid4",
                    "!ENGRAD",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pbe0": [
                    "%MaxCore 8000",
                    "! def2-TZVP(-f) pbe0 grid4",
                    "!ENGRAD",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
            }
            with open(os.path.join(self.workdir, "inp"), "w", newline=None) as inp:
                # functionals
                for line in ogo_calls[self.func]:
                    inp.write(line + "\n")
                # limitation of iterations
                # if not self.full:
                #    inp.write('%geom\n    TolE 5e-4\n    TolRMSG 1e-2\n    TolMaxG ??\n    TolRMSD ??\n    TolMaxD 1.0\nend\n')
                # nprocROR: in the optimization of
                if int(self.progsettings["omp"]) >= 1:
                    inp.write("%pal\n")
                    inp.write("    nprocs {}\n".format(self.progsettings["omp"]))
                    inp.write("end\n")
                # add solvent correction to input if solvent is specified and using smd
                if self.solv and self.sm == "smd":
                    if self.progsettings["orca_old"]:
                        for line in self.solvent_smd_old[self.solv]:
                            inp.write(line + "\n")
                    else:
                        for line in self.solvent_smd_new[self.solv]:
                            inp.write(line + "\n")
                if self.solv and self.sm == "cpcm":
                    for line in self.solvent_cpcm[self.solv]:
                        inp.write(line + "\n")
                # unpaired, charge, and coordinates
                inp.write(
                    "* xyzfile {} {} inp.xyz\n".format(
                        str(self.chrg), str(self.unpaired + 1)
                    )
                )
            # Done writing input!
            time.sleep(0.02)
            print("Running optimization in {:18}".format(last_folders(self.workdir, 2)))
            if self.full:
                if self.sm == "smd" and self.solv:
                    callargs = [
                        self.progsettings["xtbpath"],
                        "coord",
                        "--opt",
                        "lax",
                        "--orca",
                    ]
                else:  # gas phase
                    callargs = [
                        self.progsettings["xtbpath"],
                        "coord",
                        "--opt",
                        "--orca",
                    ]
            else:  # only crude
                callargs = [
                    self.progsettings["xtbpath"],
                    "coord",
                    "--opt",
                    "crude",
                    "--orca",
                ]
            with open(
                os.path.join(self.workdir, output), "w", newline=None
            ) as outputfile:
                subprocess.call(
                    callargs,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                )
            time.sleep(0.1)
        # check if optimization finished correctly:
        if os.path.isfile(os.path.join(self.workdir, output)):
            with open(
                os.path.join(self.workdir, output), "r", encoding=coding, newline=None
            ) as inp:
                stor = inp.readlines()
            for line in stor:
                if (
                    "external code error" in line
                    or "|grad| > 500, something is totally wrong!" in line
                    or "abnormal termination of xtb" in line
                ):
                    print(
                        "ERROR: optimization in {:18} not converged".format(
                            last_folders(self.workdir, 2)
                        ),
                        file=sys.stderr,
                    )
                    self.success = False
                    self.energy = None
                    return 1
                if "   *** GEOMETRY OPTIMIZATION CONVERGED AFTER " in line:
                    self.cycles = int(line.split()[5])
                if "ORCA TERMINATED NORMALLY" in line:
                    self.success = True
                if "FINAL SINGLE POINT ENERGY" in line:
                    self.energy = float(line.split()[4])
        else:
            self.success = False
            self.energy = None
            self.cycles = 0
            print(
                "WARNING: {} doesn't exist!".format(os.path.join(self.workdir, output))
            )
        if self.energy is None and self.success:
            self.success = False
            print(
                "ERROR: scf in {:18} not converged!".format(
                    last_folders(self.workdir, 2)
                ),
                file=sys.stderr,
            )
        ## if energy but not success
        if not self.success:
            self.energy = None
            print(
                "ERROR: in the optimization of {:18}".format(
                    last_folders(self.workdir, 2)
                ),
                file=sys.stderr,
            )
        # convert optimized xyz to coord file
        self.xyz2coord()
        return

    def _opt(self):
        """ORCA input generation and geometry optimization"""
        if self.full:
            output = "opt-part2.out"
        else:
            output = "opt-part1.out"
        if not self.boltzmann:
            # convert coord to xyz
            self.coord, self.nat = coord2xyz(self.workdir)
            # input generation
            ogo_calls = {
                "b97-3c": [
                    "%MaxCore 8000",
                    "! def2-mTZVP b97-3c Opt grid4 LOOSEOPT",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pbeh-3c": [
                    "%MaxCore 8000",
                    "! def2-mSVP pbeh-3c Opt grid4 LOOSEOPT",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "tpss": [
                    "%MaxCore 8000",
                    "! def2-TZVP(-f) tpss Opt grid4 LOOSEOPT",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pbe0": [
                    "%MaxCore 8000",
                    "! def2-TZVP(-f) pbe0 Opt grid4 LOOSEOPT",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
            }

            thresholds = {
                "normal": [
                    "%geom",
                    "  TolE  5e-6",
                    "  TolRMSG  1e-3",
                    "  TolMaxG  2e-3",
                    "  TolRMSD  1e-1",
                    "  TolMaxD 1e-1",
                    "end",
                ],
                "lax": [
                    "%geom",
                    "  TolE  2e-5",
                    "  TolRMSG  2e-3",
                    "  TolMaxG  2e-3",
                    "  TolRMSD  1e-1",
                    "  TolMaxD 1e-1",
                    "end",
                ],
                "crude": [
                    "%geom",
                    "  TolE  5e-4",
                    "  TolRMSG  1e-2",
                    "  TolMaxG  5e-3",
                    "  TolRMSD  1e-1",
                    "  TolMaxD 1e-1",
                    "end",
                ],
            }

            with open(os.path.join(self.workdir, "inp"), "w", newline=None) as inp:
                # functionals
                for line in ogo_calls[self.func]:
                    inp.write(line + "\n")
                # set thresholds:
                if self.full:
                    if self.sm == "smd":
                        for line in thresholds["lax"]:
                            inp.write(line + "\n")
                    else:
                        for line in thresholds["normal"]:
                            inp.write(line + "\n")
                else:
                    for line in thresholds["crude"]:
                        inp.write(line + "\n")
                if int(self.progsettings["omp"]) >= 1:
                    inp.write("%pal\n")
                    inp.write("    nprocs {}\n".format(self.progsettings["omp"]))
                    inp.write("end\n")
                # add solvent correction to input if solvent is specified and using smd
                if self.solv and self.sm == "smd":
                    if self.progsettings["orca_old"]:
                        for line in self.solvent_smd_old[self.solv]:
                            inp.write(line + "\n")
                    else:
                        for line in self.solvent_smd_new[self.solv]:
                            inp.write(line + "\n")
                if self.solv and self.sm == "cpcm":
                    for line in self.solvent_cpcm[self.solv]:
                        inp.write(line + "\n")
                # unpaired, charge, and coordinates
                inp.write("*xyz {} {}\n".format(str(self.chrg), str(self.unpaired + 1)))
                for line in self.coord:
                    inp.write(line + "\n")
                inp.write("*")
            # Done writing input!
            time.sleep(0.02)
            # start geometry optimization
            print("Running optimization in {:18}".format(last_folders(self.workdir, 2)))
            with open(
                os.path.join(self.workdir, output), "w", newline=None
            ) as outputfile:
                call = [os.path.join(self.progsettings["tempprogpath"], "orca"), "inp"]
                subprocess.call(
                    call,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                )
            time.sleep(0.1)
        # check if optimization finished correctly:
        self.success = False
        if os.path.isfile(os.path.join(self.workdir, output)):
            with open(
                os.path.join(self.workdir, output), "r", encoding=coding, newline=None
            ) as inp:
                stor = inp.readlines()
                fspe = "FINAL SINGLE POINT ENERGY"
                for line in stor:
                    if "GEOMETRY OPTIMIZATION CYCLE" in line:
                        self.cycles += 1
                    if fspe in line:
                        self.energy = float(line.split()[4])
                    if "ORCA TERMINATED NORMALLY" in line:
                        self.success = True
        else:
            self.success = False
            self.energy = None
            self.cycles = 0
            print(
                "WARNING: {} doesn't exist!".format(os.path.join(self.workdir, output))
            )
        if self.energy is None and self.success:
            self.success = False
            print(
                "ERROR: scf in {:18} not converged!".format(
                    last_folders(self.workdir, 2)
                ),
                file=sys.stderr,
            )
        ## if energy but not success
        if not self.success:
            self.energy = None
            print(
                "ERROR: in the optimization of {:18}".format(
                    last_folders(self.workdir, 2)
                ),
                file=sys.stderr,
            )
        # convert optimized xyz to coord file
        self.xyz2coord()
        return

    def _rrho(self):
        """ ORCA RRHO calculation to free energy in gas phase"""
        if not self.boltzmann:
            self.coord, self.nat = coord2xyz(self.workdir)
            # input generation
            ofreq_calls = {
                "b97-3c": [
                    "%MaxCore 8000",
                    "! def2-mTZVP b97-3c grid5 NumFreq",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pbeh-3c": [
                    "%MaxCore 8000",
                    "! def2-mSVP pbeh-3c grid5 NumFreq",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "tpss": [
                    "%MaxCore 8000",
                    "! def2-TZVP(-f) tpss grid5 NumFreq",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
            }

            ooptfreq_calls = {
                "b97-3c": [
                    "%MaxCore 8000",
                    "! def2-mTZVP b97-3c Opt grid4 NumFreq",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pbeh-3c": [
                    "%MaxCore 8000",
                    "! def2-mSVP pbeh-3c Opt grid4 NumFreq",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "tpss": [
                    "%MaxCore 8000",
                    "! def2-TZVP(-f) tpss Opt grid4 NumFreq",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
            }
            with open(os.path.join(self.workdir, "inp"), "w", newline=None) as inp:
                # if the optimization was performed in solution, the structure has to be optimized again in the gas phase
                if self.solv:
                    for line in ooptfreq_calls[self.func]:
                        inp.write(line + "\n")
                else:
                    for line in ofreq_calls[self.func]:
                        inp.write(line + "\n")
                # nprocs
                if int(self.progsettings["omp"]) >= 1:
                    inp.write("%pal\n")
                    inp.write("    nprocs {}\n".format(self.progsettings["omp"]))
                    inp.write("end\n")
                # unpaired, charge, and coordinates
                inp.write(
                    "*xyz {} {}\n".format(str(self.chrg), str(self.unpaired + 1))
                )  # int(args.unpaired) + 1
                for line in self.coord:
                    inp.write(line + "\n")
                inp.write("*")
            # Done writing input!
            time.sleep(0.02)
            # start frequency calculation
            print(
                "Running ORCA frequency calculation in {:18}".format(
                    last_folders(self.workdir, 2)
                )
            )
            with open(
                os.path.join(self.workdir, "freq.out"), "w", newline=None
            ) as outputfile:
                call = [os.path.join(self.progsettings["tempprogpath"], "orca"), "inp"]
                subprocess.call(
                    call,
                    shell=False,
                    stdin=None,
                    stderr=subprocess.STDOUT,
                    universal_newlines=False,
                    cwd=self.workdir,
                    stdout=outputfile,
                )
            time.sleep(0.1)
        # check if scf is converged:
        if os.path.isfile(os.path.join(self.workdir, "freq.out")):
            with open(
                os.path.join(self.workdir, "freq.out"),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                stor = inp.readlines()
                for line in stor:
                    if (
                        "G-E(el)                           ..." in line
                    ):  # tested with ORCA4.1
                        self.rrho = float(line.split()[2])
                    if "ORCA TERMINATED NORMALLY" in line:
                        self.success = True
        else:
            self.success = False
            self.rrho = None
            print(
                "WARNING: {} doesn't exist!".format(
                    os.path.join(self.workdir, "freq.out")
                )
            )
        if self.rrho is None and self.success:
            self.success = False
            print(
                "ERROR: RRHO calculation in {:18} failed!".format(
                    last_folders(self.workdir, 2)
                ),
                file=sys.stderr,
            )
        if not self.success:
            self.rrho = None
            print(
                "ERROR: in RRHO calculation of {:18}!".format(
                    last_folders(self.workdir, 2)
                ),
                file=sys.stderr,
            )
        return

    def _nmrJ(self):
        """ORCA NMR coupling calculation"""
        self.coord, self.nat = coord2xyz(self.workdir)
        # generate input   # double hybrids not implemented
        nmrj_calls = {
            "tpss": [
                "%MaxCore 8000",
                "! tpss " + str(self.basis) + " grid5 def2/j  nofinalgrid nososcf ",
                "!     smallprint printgap noloewdin",
                "%output",
                "       print[P_BondOrder_M] 1",
                "       print[P_Mayer] 1",
                "       print[P_basis] 2",
                "end",
            ],
            "pbe0": [
                "%MaxCore 8000",
                "! pbe0 "
                + str(self.basis)
                + " grid5 rijk def2/jk nofinalgrid nososcf ",
                "!     smallprint printgap noloewdin",
                "%output",
                "       print[P_BondOrder_M] 1",
                "       print[P_Mayer] 1",
                "       print[P_basis] 2",
                "end",
            ],
            "pbeh-3c": [
                "%MaxCore 8000",
                "! PBEh-3c grid5 rijk def2/jk nofinalgrid nososcf",
                "!     smallprint printgap noloewdin",
                "%output",
                "       print[P_BondOrder_M] 1",
                "       print[P_Mayer] 1",
                "       print[P_basis] 2",
                "end",
            ],
        }

        # 'dsd-blyp': ['%MaxCore 8000', '! dsd-blyp def2-TZVPP grid5 rijk def2/jk def2-TZVPP/C nofinalgrid nososcf',
        #              '!     smallprint printgap noloewdin', '%output', '       print[P_BondOrder_M] 1',
        #              '       print[P_Mayer] 1', '       print[P_basis] 2', 'end', '%mp2', '    RI true',
        #              '    density relaxed', 'end']
        numactive = 0
        if self.hactive:
            numactive += 1
        if self.cactive:
            numactive += 1
        if self.factive:
            numactive += 1
        if self.pactive:
            numactive += 1
        if self.siactive:
            numactive += 1
        with open(os.path.join(self.workdir, "inpJ"), "w", newline=None) as inp:
            # functional
            for line in nmrj_calls[self.func]:
                inp.write(line + "\n")
            # additional basis functions
            if self.basis == "pcJ-0":
                inp.write("%basis\n")
                inp.write("NewGTO 1\n")
                inp.write('"pcJ-0"\n')
                inp.write("p 1\n")
                inp.write("1 1.0 1.0\n")
                inp.write("end\n")
                inp.write("NewGTO 6\n")
                inp.write('"pcJ-0"\n')
                inp.write("d 1\n")
                inp.write("1 0.8 1.0\n")
                inp.write("end\n")
                inp.write("NewGTO 7\n")
                inp.write('"pcJ-0"\n')
                inp.write("d 1\n")
                inp.write("1 0.9 1.0\n")
                inp.write("end\n")
                inp.write("NewGTO 8\n")
                inp.write('"pcJ-0"\n')
                inp.write("d 1\n")
                inp.write("1 1.0 1.0\n")
                inp.write("end\n")
                inp.write("end\n")
            # nprocs
            if int(self.progsettings["omp"]) >= 1:
                inp.write("%pal\n")
                inp.write("    nprocs {}\n".format(self.progsettings["omp"]))
                inp.write("end\n")
            # add solvent correction to input if solvent is specified and using smd
            if self.solv and self.sm == "smd":
                if self.progsettings["orca_old"]:
                    for line in self.solvent_smd_old[self.solv]:
                        inp.write(line + "\n")
                else:
                    for line in self.solvent_smd_new[self.solv]:
                        inp.write(line + "\n")
            if self.solv and self.sm == "cpcm":
                for line in self.solvent_cpcm[self.solv]:
                    inp.write(line + "\n")
            # unpaired,charge, and coordinates
            inp.write("*xyz {} {}\n".format(str(self.chrg), str(self.unpaired + 1)))
            for line in self.coord:
                inp.write(line + "\n")
            inp.write("*\n")
            # couplings
            inp.write("%eprnmr\n")
            if numactive == 1:
                if self.hactive or self.pactive:
                    inp.write(" Nuclei = all H { ssfc }\n")
                    inp.write(" Nuclei = all F { ssfc }\n")
                    inp.write(" Nuclei = all P { ssfc }\n")
                elif self.cactive:
                    inp.write(" Nuclei = all C { ssfc }\n")
                    inp.write(" Nuclei = all F { ssfc }\n")
                elif self.factive:
                    inp.write(" Nuclei = all H { ssfc }\n")
                    inp.write(" Nuclei = all C { ssfc }\n")
                    inp.write(" Nuclei = all F { ssfc }\n")
                    inp.write(" Nuclei = all P { ssfc }\n")
            else:
                inp.write(" Nuclei = all H { ssfc }\n")
                inp.write(" Nuclei = all C { ssfc }\n")
                inp.write(" Nuclei = all F { ssfc }\n")
                inp.write(" Nuclei = all P { ssfc }\n")
                inp.write(" Nuclei = all Si { ssfc }\n")
            # if self.hactive or self.factive:
            #    if self.cactive or self.pactive or self.factive:
            #        inp.write(' Nuclei = all H { ssfc }\n')
            #        inp.write(' Nuclei = all C { ssfc }\n')
            #        inp.write(' Nuclei = all F { ssfc }\n')
            #        inp.write(' Nuclei = all P { ssfc }\n')
            #    else:
            #        inp.write(' Nuclei = all H { ssfc }\n')
            #        inp.write(' Nuclei = all F { ssfc }\n')
            # else:
            #    inp.write(' Nuclei = all C { ssfc }\n')
            #    inp.write(' Nuclei = all F { ssfc }\n')
            #    inp.write(' Nuclei = all P { ssfc }\n')
            inp.write(" SpinSpinRThresh 8.0\n")
            inp.write("end\n")
        # Done input!
        # start coupling calculation
        print(
            "Running coupling calculation in {}".format(last_folders(self.workdir, 2))
        )
        with open(
            os.path.join(self.workdir, "orcaJ.out"), "w", newline=None
        ) as outputfile:
            call = [os.path.join(self.progsettings["tempprogpath"], "orca"), "inpJ"]
            subprocess.call(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.workdir,
                stdout=outputfile,
            )
        time.sleep(0.1)
        # check if calculation was successfull:
        with open(
            os.path.join(self.workdir, "orcaJ.out"), "r", encoding=coding, newline=None
        ) as inp:
            store = inp.readlines()
            self.success = False
            for line in store:
                if "ORCA TERMINATED NORMALLY" in line:
                    self.success = True
        if not self.success:
            print(
                "ERROR: coupling calculation in {:18} failed!".format(
                    last_folders(self.workdir, 1)
                ),
                file=sys.stderr,
            )
        return

    def _nmrS(self):
        """ORCA NMR shielding calculation"""
        self.coord, self.nat = coord2xyz(self.workdir)
        nmrs_calls = {
            "tpss": [
                "%MaxCore 8000",
                "! tpss " + str(self.basis) + " grid5 def2/j nofinalgrid nososcf",
                "!     smallprint printgap noloewdin",
                "%output",
                "       print[P_BondOrder_M] 1",
                "       print[P_Mayer] 1",
                "       print[P_basis] 2",
                "end",
            ],
            "pbe0": [
                "%MaxCore 8000",
                "! pbe0 "
                + str(self.basis)
                + " grid5 rijcosx def2/j nofinalgrid nososcf gridx6",
                "!     smallprint printgap noloewdin",
                "%output",
                "       print[P_BondOrder_M] 1",
                "       print[P_Mayer] 1",
                "       print[P_basis] 2",
                "end",
            ],
            "dsd-blyp": [
                "%MaxCore 8000",
                "! dsd-blyp def2-TZVPP grid5 rijcosx def2/j def2-TZVPP/C nofinalgrid nososcf gridx7 nofinalgridx",
                "!     nofrozencore smallprint printgap noloewdin",
                "%output",
                "       print[P_BondOrder_M] 1",
                "       print[P_Mayer] 1",
                "       print[P_basis] 2",
                "end",
                "%mp2",
                "    RI true",
                "    density relaxed",
                "end",
            ],
            "wb97x": [
                "%MaxCore 8000",
                "! wb97x-d3 "
                + str(self.basis)
                + " grid5 rijcosx def2/j nofinalgrid nososcf gridx6",
                "!     smallprint printgap noloewdin",
                "%output",
                "       print[P_BondOrder_M] 1",
                "       print[P_Mayer] 1",
                "       print[P_basis] 2",
                "end",
            ],
            "pbeh-3c": [
                "%MaxCore 8000",
                "! PBEh-3c grid5 def2/j nofinalgrid nososcf",
                "!     smallprint printgap noloewdin",
                "%output",
                "       print[P_BondOrder_M] 1",
                "       print[P_Mayer] 1",
                "       print[P_basis] 2",
                "end",
            ],
        }

        with open(os.path.join(self.workdir, "inpS"), "w", newline=None) as inp:
            # functional
            try:
                for line in nmrs_calls[self.func]:
                    inp.write(line + "\n")
            except KeyError:
                print("Error: Functional was not found!")
                self.success = False
                return
            # nprocs
            if int(self.progsettings["omp"]) >= 1:
                inp.write("%pal\n")
                inp.write("    nprocs {}\n".format(self.progsettings["omp"]))
                inp.write("end\n")
            # add solvent correction to input if solvent is specified and using smd
            if self.solv and self.sm == "smd":
                if self.progsettings["orca_old"]:
                    for line in self.solvent_smd_old[self.solv]:
                        inp.write(line + "\n")
                else:
                    for line in self.solvent_smd_new[self.solv]:
                        inp.write(line + "\n")
            if self.solv and self.sm == "cpcm":
                for line in self.solvent_cpcm[self.solv]:
                    inp.write(line + "\n")
            # unpaired, charge, and coordinates
            inp.write("*xyz {} {}\n".format(str(self.chrg), str(self.unpaired + 1)))
            for line in self.coord:
                inp.write(line + "\n")
            inp.write("*\n")
            # shielding
            inp.write("%eprnmr\n")
            if self.hactive:
                inp.write(" Nuclei = all H { shift }\n")
            if self.cactive:
                inp.write(" Nuclei = all C { shift }\n")
            if self.factive:
                inp.write(" Nuclei = all F { shift }\n")
            if self.pactive:
                inp.write(" Nuclei = all P { shift }\n")
            if self.siactive:
                inp.write(" Nuclei = all Si { shift }\n")
            inp.write(" origin giao\n")
            inp.write(" giao_2el giao_2el_same_as_scf\n")
            inp.write(" giao_1el giao_1el_analytic\n")
            inp.write("end\n")
        # Done input!
        # shielding calculation
        print(
            "Running shielding calculation in {:18}".format(
                last_folders(self.workdir, 2)
            )
        )
        with open(
            os.path.join(self.workdir, "orcaS.out"), "w", newline=None
        ) as outputfile:
            call = [os.path.join(self.progsettings["tempprogpath"], "orca"), "inpS"]
            subprocess.call(
                call,
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.workdir,
                stdout=outputfile,
            )
        time.sleep(0.1)
        # check if calculation was successfull:
        with open(
            os.path.join(self.workdir, "orcaS.out"), "r", encoding=coding, newline=None
        ) as inp:
            store = inp.readlines()
            self.success = False
            for line in store:
                if "ORCA TERMINATED NORMALLY" in line:
                    self.success = True
        if not self.success:
            print(
                "ERROR: shielding calculation in {:18} failed!".format(
                    last_folders(self.workdir, 1)
                ),
                file=sys.stderr,
            )
        return

    def _genericoutput(self):
        """ORCA READ shielding and coupling constants and write them to plain output"""
        fnameshield = "orcaS.out"
        atom = []
        sigma = []
        try:
            with open(
                os.path.join(self.workdir, fnameshield),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                data = inp.readlines()
            for line in data:
                if "CHEMICAL SHIELDING SUMMARY (ppm)" in line:
                    start = data.index(line)
            for line in data[(start + 6) :]:
                splitted = line.split()
                if len(splitted) == 4:
                    atom.append(int(splitted[0]) + 1)
                    sigma.append(float(splitted[2]))
                else:
                    break
        except FileNotFoundError:
            print(
                "Missing file: {} in {}".format(
                    fnameshield, last_folders(self.workdir, 2)
                )
            )
            self.success = False
        self.success = True
        fnamecoupl = "orcaJ.out"
        atom1 = []
        atom2 = []
        jab = []
        try:
            with open(
                os.path.join(self.workdir, fnamecoupl),
                "r",
                encoding=coding,
                newline=None,
            ) as inp:
                data = inp.readlines()
            for line in data:
                if "NMR SPIN-SPIN COUPLING CONSTANTS" in line:
                    start = int(data.index(line)) + 6
                if " ****ORCA TERMINATED NORMALLY****" in line:
                    end = int(data.index(line))

            for line in data[start:end]:
                if "NUCLEUS" in line:
                    tmpsplitted = line.split()
                    atom1.append(int(tmpsplitted[4]) + 1)
                    atom2.append(int(tmpsplitted[9]) + 1)
                elif "Total" in line and "iso= " in line:
                    splitted = line.split()
                    jab.append(float(splitted[5]))
                else:
                    pass
        except FileNotFoundError:
            print(
                "Missing file: {} in {}".format(
                    fnamecoupl, last_folders(self.workdir, 2)
                )
            )
            self.success = False
        self.success = True

        with open(os.path.join(self.workdir, "nmrprop.dat"), "w", newline=None) as out:
            for i in range(len(atom)):
                out.write("{:{digits}} {}\n".format(atom[i], sigma[i], digits=4))
            for i in range(1, self.nat - len(atom)):
                out.write("\n")
            for i in range(len(atom1)):
                out.write(
                    "{:{digits}} {:{digits}}   {}\n".format(
                        atom1[i], atom2[i], jab[i], digits=4
                    )
                )
        time.sleep(0.02)
        return

    def execute(self):
        """Choose what to execute for the jobtype"""
        if self.jobtype == "prep":
            self.success = True
            pass
        elif self.jobtype == "sp":
            self._sp()
        elif self.jobtype == "opt":
            self._opt()
        elif self.jobtype == "xtbopt":
            self._xtbopt(environsettings)
        elif self.jobtype == "rrhoxtb":
            self._xtbrrho(environsettings)
        elif self.jobtype == "rrhoorca":
            self._rrho()
        elif self.jobtype == "solv":
            self._solv_complete(environsettings)
        elif self.jobtype == "gbsa_gsolv":
            self._gbsa_rs(environsettings)
        elif self.jobtype == "nmrJ":
            self._nmrJ()
        elif self.jobtype == "nmrS":
            self._nmrS()
        elif self.jobtype == "genericout":
            self._genericoutput()
        return


def execute_data(q, resultq):
    """code that the worker has to execute """
    while True:
        if q.empty():
            # print('Queue empty')
            break
        task = q.get()
        # print('Working in dir: {} with {}'.format(task.workdir, task.jobtype))
        task.execute()
        resultq.put(task)
        time.sleep(0.02)
        q.task_done()
    return


def run_in_parallel(q, resultq, job, maxthreads, loopover, instructdict, foldername=""):
    """Run jobs in parallel
    q = queue to put assemble tasks
    resultq = queue to retrieve results
    job = information which kind of job is to be performed tm_job , orca_job
    loopover either nested list with [[path/CONFX, path/CONFX/func]...] or
    loopover is list of qm_class objects
    instrucdict example : {'jobtype': 'prep', 'chrg': args.chrg}
    foldername is for existing objects to change the workdir
    results = list of qm_class objects with results from calculations
    """

    if instructdict.get("jobtype", None) is None:
        raise KeyError("jobtype is missing in instructdict!")

    cwd = os.getcwd()
    tmp_len = []
    if all(isinstance(x, qm_job) for x in loopover):
        # for already existing qm_job objects
        for item in loopover:
            if isinstance(item, tm_job) and job == orca_job:
                item.__class__ = job
                # print('AFTER change: name: {}, type: {}'.format(item.name, type(item)))
            elif isinstance(item, orca_job) and job == tm_job:
                item.__class__ = job
                # print('AFTER change: name: {}, type: {}'.format(item.name, type(item)))
            else:
                pass
            # name is already known
            item.workdir = os.path.normpath(
                os.path.join(cwd, os.path.join(item.name, foldername))
            )
            for instruction in instructdict:
                # update instructions
                setattr(item, instruction, instructdict[instruction])
            tmp_len.append(last_folders(item.workdir, 2))
            q.put(item)
    else:
        # for new creation of qm_job objects
        for item in loopover:
            task = job()
            task.name = os.path.basename(item[0])
            task.workdir = os.path.normpath(
                os.path.join(cwd, os.path.join(item[0], foldername))
            )
            for instruction in instructdict:
                setattr(task, instruction, instructdict[instruction])
            tmp_len.append(last_folders(task.workdir, 2))
            q.put(task)
    njobs = q.qsize()
    if instructdict.get("boltzmann", False):
        print(
            "\nReading data from {} conformers calculated in previous run.".format(
                njobs
            )
        )
    else:
        if instructdict["jobtype"] is "prep":
            print("\nPreparing {} calculations.".format(njobs))
        elif instructdict["jobtype"] in ("xtbopt", "opt"):
            if instructdict["full"]:
                print("\nStarting {} optimizations.".format(njobs))
            else:
                print("\nStarting {} crude optimization calculations.".format(njobs))
        elif instructdict["jobtype"] is "sp":
            print("\nStarting {} single-point calculations.".format(njobs))
        elif instructdict["jobtype"] is "solv":
            print("\nStarting {} Gsolv calculations.".format(njobs))
        elif instructdict["jobtype"] in ("rrhoxtb", "rrhoorca", "rrhotm"):
            print("\nStarting {} GRRHO calculations".format(njobs))
        elif instructdict["jobtype"] is "nmrJ":
            print("\nStarting {} Coupling calculations".format(njobs))
        elif instructdict["jobtype"] is "nmrS":
            print("\nStarting {} Shielding calculations".format(njobs))
        elif instructdict["jobtype"] is "gbsa_gsolv":
            print("\nStarting {} GBSA-Gsolv calculations".format(njobs))
        elif instructdict["jobtype"] is "genericout":
            print("\nWriting generic output!")
    time.sleep(0.5)

    # start working in parallel
    for i in range(maxthreads):
        worker = Thread(target=execute_data, args=(q, resultq))
        worker.setDaemon(True)
        worker.start()
    q.join()

    if not instructdict.get("boltzmann", False):
        print("Tasks completed!\n")
    else:
        print("Reading data from in previous run completed!\n")

    # Get results
    results = []
    ### get formatting information:
    try:
        maxworkdirlen = max([len(i) for i in tmp_len])
    except Exception as e:
        print(e)
        maxworkdirlen = 20
    ### end formatting information
    while not resultq.empty():
        results.append(resultq.get())
        if instructdict["jobtype"] is "prep":
            if results[-1].success:
                print(
                    "Preparation in {:{digits}} was successful".format(
                        last_folders(results[-1].workdir, 2), digits=maxworkdirlen
                    )
                )
            else:
                print(
                    "Preparation in {:{digits}} FAILED".format(
                        last_folders(results[-1].workdir, 2), digits=maxworkdirlen
                    )
                )
        elif instructdict["jobtype"] in ("xtbopt", "opt"):
            if instructdict["full"]:
                if results[-1].energy is not None and results[-1].success:
                    print(
                        "Finished optimization for {:{digits}} after {:>3} cycles: {:.6f} ".format(
                            last_folders(results[-1].workdir, 2),
                            str(results[-1].cycles),
                            results[-1].energy,
                            digits=maxworkdirlen,
                        )
                    )
                else:
                    print(
                        "Optimization FAILED for {:{digits}} after {:>3} cycles: {} ".format(
                            last_folders(results[-1].workdir, 2),
                            str(results[-1].cycles),
                            str(results[-1].energy),
                            digits=maxworkdirlen,
                        )
                    )
            else:  # crude optimization
                if results[-1].energy is not None and results[-1].success:
                    print(
                        "Finished crude optimization for {:{digits}} after {:>3} cycles: {:.6f} ".format(
                            last_folders(results[-1].workdir, 2),
                            str(results[-1].cycles),
                            results[-1].energy,
                            digits=maxworkdirlen,
                        )
                    )
                else:
                    print(
                        "Crude optimization FAILED for {:{digits}} after {:>3} cycles: {} ".format(
                            last_folders(results[-1].workdir, 2),
                            str(results[-1].cycles),
                            str(results[-1].energy),
                            digits=maxworkdirlen,
                        )
                    )
        elif instructdict["jobtype"] is "sp":
            if results[-1].energy is not None:
                print(
                    "Finished single-point calculation for {:{digits}}: {:.6f} ".format(
                        last_folders(results[-1].workdir, 2),
                        results[-1].energy,
                        digits=maxworkdirlen,
                    )
                )
            else:
                print(
                    "Single-point calculation FAILED for {:{digits}}: {} ".format(
                        last_folders(results[-1].workdir, 2),
                        str(results[-1].energy),
                        digits=maxworkdirlen,
                    )
                )
        elif instructdict["jobtype"] is "solv":
            if results[-1].gsolv is not None:
                print(
                    "Finished Gsolv for {:{digits}}: {:.6f}".format(
                        last_folders(results[-1].workdir, 2),
                        results[-1].gsolv,
                        digits=maxworkdirlen,
                    )
                )
            else:
                print(
                    "Gsolv FAILED for {:{digits}}: {}".format(
                        last_folders(results[-1].workdir, 2),
                        str(results[-1].gsolv),
                        digits=maxworkdirlen,
                    )
                )
        elif instructdict["jobtype"] in ("rrhoxtb", "rrhoorca", "rrhotm"):
            if results[-1].rrho is not None:
                print(
                    "Finished RRHO for {:{digits}}:  {:.6f} in sym: {}".format(
                        last_folders(results[-1].workdir, 2),
                        results[-1].rrho,
                        results[-1].symmetry,
                        digits=maxworkdirlen,
                    )
                )
            else:
                print(
                    "FAILED RRHO for {:{digits}}:  {}".format(
                        last_folders(results[-1].workdir, 2),
                        str(results[-1].rrho),
                        digits=maxworkdirlen,
                    )
                )
        elif instructdict["jobtype"] is "nmrJ":
            if results[-1].success:
                print(
                    "NMR-J calculation in {:{digits}} was successful.".format(
                        last_folders(results[-1].workdir, 2), digits=maxworkdirlen
                    )
                )
            else:
                print(
                    " NMR-J calculation FAILED for {:{digits}}!".format(
                        last_folders(results[-1].workdir, 2), digits=maxworkdirlen
                    )
                )
        elif instructdict["jobtype"] is "nmrS":
            if results[-1].success:
                print(
                    "NMR-S calculation in {:{digits}} was successful.".format(
                        last_folders(results[-1].workdir, 2), digits=maxworkdirlen
                    )
                )
            else:
                print(
                    "NMR-S calculation FAILED for {:{digits}}!".format(
                        last_folders(results[-1].workdir, 2), digits=maxworkdirlen
                    )
                )
        elif instructdict["jobtype"] is "gbsa_gsolv":
            if results[-1].gsolv is not None:
                print(
                    "Finished Gsolv-GBSA-Gsolv for {:{digits}}: {:.6f}".format(
                        last_folders(results[-1].workdir, 2),
                        results[-1].gsolv,
                        digits=maxworkdirlen,
                    )
                )
            else:
                print(
                    "Gsolv-GBSA-Gsolv FAILED for {:{digits}}: {}".format(
                        last_folders(results[-1].workdir, 2),
                        str(results[-1].gsolv),
                        digits=maxworkdirlen,
                    )
                )
        time.sleep(0.01)  # sleep is important don't remove it!!!

    # sort results by name
    results.sort(key=lambda x: int(x.name[4:]))  # (CONFX)
    return results


def check_tasks(results, args, thresh=0.25):
    """ Check if too many tasks failed and exit if so!"""
    # Check if preparation failed too often:
    counter = 0
    exit_log = False
    fail_rate = None
    for item in results:
        if not item.success:
            counter += 1
    fail_rate = float(counter) / float(len(results)) * 100
    if float(counter) / float(len(results)) >= thresh and args.check:
        exit_log = True
    return exit_log, fail_rate


class handle_input:
    """ Read ensorc, write ensorc, write flags.dat"""

    #add read_json write_json at some point, because they should be connected

    known_keys = (
        "reference for 1H",
        "reference for 13C",
        "reference for 19F",
        "reference for 31P",
        "reference for 29Si",
        "1H active",
        "13C active",
        "19F active",
        "31P active",
        "29Si active",
        "resonance frequency",
        "nconf",
        "charge",
        "unpaired",
        "solvent",
        "prog",
        "ancopt",
        "prog_rrho",
        "gfn_version",
        "prog3",
        "prog4",
        "part1",
        "part2",
        "part3",
        "part4",
        "boltzmann",
        "backup",
        "func",
        "func3",
        "basis3",
        "funcJ",
        "basisJ",
        "funcS",
        "basisS",
        "couplings",
        "shieldings",
        "part1_threshold",
        "part2_threshold",
        "sm",
        "sm3",
        "sm4",
        "check",
        "crestcheck",
        "maxthreads",
        "omp",
        "smgsolv2",
        "temperature",
    )
    key_args_dict = {
        "nconf": "nstruc",
        "charge": "chrg",
        "unpaired": "unpaired",
        "solvent": "solv",
        "prog": "prog",
        "ancopt": "ancopt",
        "prog_rrho": "rrhoprog",
        "gfn_version": "gfnv",
        "temperature": "temperature",
        "prog3": "prog3",
        "prog4": "prog4",
        "part1": "part1",
        "part2": "part2",
        "part3": "part3",
        "part4": "part4",
        "boltzmann": "boltzmann",
        "backup": "backup",
        "func": "func",
        "func3": "func3",
        "basis3": "basis3",
        "couplings": "calcJ",
        "funcJ": "funcJ",
        "basisJ": "basisJ",
        "shieldings": "calcS",
        "funcS": "funcS",
        "basisS": "basisS",
        "part1_threshold": "thr1",
        "part2_threshold": "thr2",
        "sm": "sm",
        "smgsolv2": "gsolv2",
        "sm3": "sm3",
        "sm4": "sm4",
        "check": "check",
        "crestcheck": "crestcheck",
        "maxthreads": "maxthreads",
        "omp": "omp",
        "1H active": "hactive",
        "13C active": "cactive",
        "19F active": "factive",
        "31P active": "pactive",
        "29Si active": "siactive",
        "resonance frequency": "mf",
        "reference for 1H": "href",
        "reference for 13C": "cref",
        "reference for 31P": "pref",
        "reference for 19F": "fref",
        "reference for 29Si": "siref",
    }
    # später für unittest!
    # print(list(set(known_keys) - set(key_args_dict.keys())))

    enso_internal_defaults = {
        "reference for 1H": "TMS",
        "reference for 13C": "TMS",
        "reference for 19F": "CFCl3",
        "reference for 31P": "TMP",
        "reference for 29Si": "TMS",
        "1H active": "on",
        "13C active": "off",
        "19F active": "off",
        "31P active": "off",
        "29Si active": "off",
        "resonance frequency": "None",
        "nconf": "all",
        "charge": "0",
        "unpaired": "0",
        "solvent": "gas",
        "prog": "None",
        "ancopt": "on",
        "prog_rrho": "xtb",
        "gfn_version": "gfn2",
        "temperature": "298.15",
        "prog3": "prog",
        "prog4": "prog",
        "part1": "on",
        "part2": "on",
        "part3": "on",
        "part4": "on",
        "boltzmann": "off",
        "backup": "off",
        "func": "pbeh-3c",
        "func3": "pw6b95",
        "basis3": "def2-TZVPP",
        "funcJ": "pbe0",
        "basisJ": "def2-TZVP",
        "funcS": "pbe0",
        "basisS": "def2-TZVP",
        "couplings": "on",
        "shieldings": "on",
        "part1_threshold": "4.0",
        "part2_threshold": "2.0",
        "sm": "default",
        "smgsolv2": "sm",
        "sm3": "default",
        "sm4": "default",
        "check": "on",
        "crestcheck": "off",
        "maxthreads": "1",
        "omp": "1",
    }
    # später für unittest
    # print(list(set(known_keys) - set(enso_internal_defaults.keys())))

    knownbasissets3 = (
        "SVP",
        "SV(P)",
        "TZVP",
        "TZVPP",
        "QZVP",
        "QZVPP",
        "def2-SV(P)",
        "def2-SVP",
        "def2-TZVP",
        "def2-TZVPP",
        "def-SVP",
        "def-SV(P)",
        "def2-QZVP",
        "DZ",
        "QZV",
        "cc-pVDZ",
        "cc-pVTZ",
        "cc-pVQZ",
        "cc-pV5Z",
        "aug-cc-pVDZ",
        "aug-cc-pVTZ",
        "aug-cc-pVQZ",
        "aug-cc-pV5Z",
        "def2-QZVPP",
    )
    knownbasissetsJ = (
        "SVP",
        "SV(P)",
        "TZVP",
        "TZVPP",
        "QZVP",
        "QZVPP",
        "def2-SV(P)",
        "def2-SVP",
        "def2-TZVP",
        "def2-TZVPP",
        "def-SVP",
        "def-SV(P)",
        "def2-QZVP",
        "DZ",
        "QZV",
        "cc-pVDZ",
        "cc-pVTZ",
        "cc-pVQZ",
        "cc-pV5Z",
        "aug-cc-pVDZ",
        "aug-cc-pVTZ",
        "aug-cc-pVQZ",
        "aug-cc-pV5Z",
        "def2-QZVPP",
        "pcJ-0",
        "pcJ-1",
        "pcJ-2",
    )
    knownbasissetsS = (
        "SVP",
        "SV(P)",
        "TZVP",
        "TZVPP",
        "QZVP",
        "QZVPP",
        "def2-SV(P)",
        "def2-SVP",
        "def2-TZVP",
        "def2-TZVPP",
        "def-SVP",
        "def-SV(P)",
        "def2-QZVP",
        "DZ",
        "QZV",
        "cc-pVDZ",
        "cc-pVTZ",
        "cc-pVQZ",
        "cc-pV5Z",
        "aug-cc-pVDZ",
        "aug-cc-pVTZ",
        "aug-cc-pVQZ",
        "aug-cc-pV5Z",
        "def2-QZVPP",
        "pcSseg-0",
        "pcSseg-1",
        "pcSseg-2",
        "pcSseg-3",
        "x2c-SVPall-s",
        "x2c-TZVPall-s",
    )

    def __init__(
        self, solvents, gfnv, func, func3, funcJ, funcS, href, cref, fref, pref, siref
    ):
        self.configdata = {}
        self.value_options = {
            "part1": ("on", "off"),
            "nconf": ("all", "number e.g. 10"),
            "charge": "number e.g. 0",
            "unpaired": "number e.g. 0",
            "solvent": ("gas", solvents),
            "prog": ("tm", "orca"),
            "part2": ("on", "off"),
            "part3": ("on", "off"),
            "part4": ("on", "off"),
            "prog3": ("tm", "orca", "prog"),
            "prog4": ("tm", "orca", "prog"),
            "ancopt": ("on", "off"),
            "prog_rrho": ("xtb", "prog", "off"),
            "gfn_version": gfnv,
            "temperature": "temperature in K e.g. 298.15",
            "boltzmann": ("on", "off"),
            "backup": ("on", "off"),
            "func": func,
            "func3": func3,
            "couplings": ("on", "off"),
            "basis3": self.knownbasissets3,
            "funcJ": funcJ,
            "basisJ": self.knownbasissetsJ,
            "funcS": funcS,
            "basisS": self.knownbasissetsS,
            "reference for 1H": href,
            "reference for 13C": cref,
            "reference for 19F": fref,
            "reference for 31P": pref,
            "reference for 29Si": siref,
            "1H active": ("on", "off"),
            "13C active": ("on", "off"),
            "19F active": ("on", "off"),
            "31P active": ("on", "off"),
            "29Si active": ("on", "off"),
            "resonance frequency": "MHz number of your experimental spectrometer",
            "shieldings": ("on", "off"),
            "part1_threshold": "number e.g. 4.0",
            "part2_threshold": "number e.g. 2.0",
            "sm": ("cosmo", "dcosmors", "cpcm", "smd"),
            "sm3": ("cosmors", "smd", "gbsa_gsolv"),
            "sm4": ("cosmo", "cpcm", "smd"),
            "check": ("on", "off"),
            "crestcheck": ("on", "off"),
            "maxthreads": "number e.g. 2",
            "omp": "number e.g. 4",
            "smgsolv2": ("sm", "cosmors", "gbsa_gsolv"),
        }

    def _decomment(self, csvfile):
        """ remove any comments from file before parsing with csv.DictReader"""
        for row in csvfile:
            raw = row.split("#")[0].strip()
            raw2 = raw.split("$")[0].strip()
            if raw2:
                yield raw2
        return

    def _read_config(self, path, test):
        """ Read from config data from file (here flags.dat or .ensorc"""
        configdata = {}
        with open(path, "r") as csvfile:
            # skip header:
            while True:
                line = csvfile.readline()
                if line.startswith(test) or line.startswith("$"):
                    break
                else:
                    pass
            reader = csv.DictReader(
                self._decomment(csvfile),
                fieldnames=("key", "value"),
                skipinitialspace=True,
                delimiter=":",
            )
            for row in reader:
                if "end" in row["key"]:
                    break
                else:
                    configdata[row["key"]] = row["value"]
        self.configdata = configdata
        return

    def read_ensorc_data(self, args, path, startreading):
        """ check configdata for typos"""
        print("Reading user set defaults (ensorc).")
        self._read_config(path, startreading)
        error_logical = False
        # check keys:
        if "end" in self.configdata:
            del self.configdata["end"]
        readinkeys = []
        for item in self.configdata:
            if item not in self.known_keys:
                print("WARNING: {} is not a known keyword in ensorc.".format(item))
                error_logical = True
            else:
                readinkeys.append(item)
        diff = list(set(self.known_keys) - set(readinkeys))
        if diff:
            print(
                "WARNING: These keywords were not found in your ensorc \n"
                "         and therefore default values are taken "
                "for: {}".format(diff)
            )
            for item in diff:
                self.configdata[item] = self.enso_internal_defaults[item]
                print(
                    'WARNING: for keyword: "{}" the default: "{}" has been set!'.format(
                        item, self.enso_internal_defaults[item]
                    )
                )
        # choices:
        for key in self.configdata:
            if self.configdata[key] == "":
                self.configdata[key] = None
        # if error_logical:
        # print('You can choose from: {}'.format(self.known_keys))
        # print('An error has occured!')

        # ensorc value data is not checked --> done later on
        # --> cml >> ensorc
        for key, value in self.key_args_dict.items():
            if key in self.configdata and not getattr(args, value):
                try:
                    setattr(args, value, self.configdata[key])
                except:
                    pass
        return

    def read_program_paths(self, path):
        """Get absolute paths of programs employed in enso"""
        # read program paths
        print("Reading absolute paths of programs employed in ENSO.")
        dbpath = None
        cosmothermversion = None
        cosmorssetup = None
        orcapath = None
        crestpath = None
        xtbpath = None
        mpshiftpath = None
        escfpath = None
        orca_old = None
        orcaversion = None

        with open(path, "r") as inp:
            stor = inp.readlines()
        for line in stor:
            if "ctd =" in line:
                try:
                    cosmorssetup = str(line.rstrip(os.linesep))
                except:
                    print(
                        "WARNING: could not read settings for COSMO-RS from" " .ensorc!"
                    )
                try:
                    dbpath = os.path.join(
                        os.path.split(cosmorssetup.split()[5].strip('"'))[0],
                        "DATABASE-COSMO/BP-TZVP-COSMO",
                    )
                    os.path.isdir(dbpath)
                except:
                    print(
                        "WARNING: could not read settings for COSMO-RS from "
                        ".ensorc!\nMost probably there is a user "
                        "input error."
                    )
            if "cosmothermversion:" in line:
                try:
                    cosmothermversion = int(line.split()[1])
                except:
                    print(
                        "WARNING: cosmothermversion could not be read! This "
                        "is necessary to prepare the cosmotherm.inp! "
                    )
            if "ORCA:" in line:
                try:
                    orcapath = str(line.split()[1])
                except:
                    print("WARNING: could not read path for ORCA from " ".ensorc!.")
            if "ORCA version:" in line:
                try:
                    tmp = line.split()[2]
                    tmp = tmp.split(".")
                    tmp.insert(1, ".")
                    tmp = "".join(tmp)
                    if float(tmp) < 4.1:
                        orca_old = True
                    elif float(tmp) >= 4.1:
                        orca_old = False
                    orcaversion = tmp
                except:
                    print("WARNING: could not read ORCA version from " ".ensorc!")
            if "GFN-xTB:" in line:
                try:
                    xtbpath = str(line.split()[1])
                except:
                    print("WARNING: could not read path for GFNn-xTB from .ensorc!")
                    if shutil.which("xtb") is not None:
                        xtbpath = shutil.which("xtb")
                        print("Going to use {} instead.".format(xtbpath))
            if "CREST:" in line:
                try:
                    crestpath = str(line.split()[1])
                except:
                    print("WARNING: could not read path for CREST from .ensorc!")
                    if shutil.which("crest") is not None:
                        crestpath = shutil.which("crest")
                        print("Going to use {} instead.".format(crestpath))
            if "mpshift:" in line:
                try:
                    mpshiftpath = str(line.split()[1])
                except:
                    print("ẂARNING: could not read path for mpshift from " ".ensorc!")
            if "escf:" in line:
                try:
                    escfpath = str(line.split()[1])
                except:
                    print("WARNING: could not read path for escf from " ".ensorc!")
        return (
            dbpath,
            cosmothermversion,
            cosmorssetup,
            orcapath,
            crestpath,
            xtbpath,
            mpshiftpath,
            escfpath,
            orca_old,
            orcaversion,
        )

    def write_ensorc(self):
        """ write an ensorc-new file into your current directory"""

        # name of file
        # remember to update paths in ensorc and adjust personal settings
        # move ensorc to /home/$USER/
        try:
            with open("ensorc-new", "w", newline=None) as outdata:
                outdata.write(".ENSORC\n")
                outdata.write("\n")
                outdata.write("ORCA: /path/excluding/binary/\n")
                outdata.write("ORCA version: 4.2.1\n")
                outdata.write("GFN-xTB: /path/including/binary/xtb-binary\n")
                outdata.write("CREST: /path/including/binary/crest-binary\n")
                outdata.write("mpshift: /path/including/binary/mpshift-binary\n")
                outdata.write("escf: /path/including/binary/xtb-binary\n")
                outdata.write("\n")
                outdata.write("#COSMO-RS\n")
                outdata.write(
                    "ctd = BP_TZVP_C30_1601.ctd cdir = "
                    '"/software/cluster/COSMOthermX16/COSMOtherm/CTDATA-FILES" ldir = '
                    '"/software/cluster/COSMOthermX16/COSMOtherm/CTDATA-FILES"\n'
                )
                outdata.write("cosmothermversion: 16\n")
                outdata.write("\n")
                outdata.write("$ENDPROGRAMS\n\n")
                outdata.write("#NMR data\n")
                for item in self.enso_internal_defaults:
                    keylen = len(item)
                    outdata.write(
                        "{}: {:{digits}} # {} \n".format(
                            item,
                            self.enso_internal_defaults[item],
                            self.value_options.get(item, "possibilities"),
                            digits=40 - keylen,
                        )
                    )
        except:
            pass
        return

    def write_flags(
        self,
        args,
        solvents,
        func,
        func3,
        funcJ,
        funcS,
        gfnv,
        href,
        cref,
        fref,
        pref,
        siref,
    ):
        """ Write file "flags.dat" with restart options"""
        if args.boltzmann == "on":
            args.part1 = "off"
            args.part2 = "off"
            args.part3 = "on"
            args.part4 = "off"
        if args.solv and args.solv != "gas":  # if solvent model is used
            if args.part1 == "on" or args.part2 == "on":
                if args.prog == "tm":
                    if args.sm == "dcosmors" or args.sm == "cosmo":
                        pass
                    elif args.sm is None or args.sm == "default":
                        args.sm = "dcosmors"
                        print(
                            "WARNING: DCOSMO-RS is used as default solvent model "
                            "with TM."
                        )
                    else:
                        args.sm = "dcosmors"
                        print(
                            """WARNING: solvent model for part 1 and 2 is not recognized! Options for TM are COSMO 
                               (cosmo) or DCOSMO-RS (dcosmors). Solvent model is set to TM default DCOSMO-RS."""
                        )
                elif args.prog == "orca":
                    if args.sm == "cpcm" or args.sm == "smd":
                        pass
                    elif args.sm is None or args.sm == "default":
                        args.sm = "smd"
                        print(
                            "WARNING: SMD is used as default solvent model with ORCA."
                        )
                    else:
                        args.sm = "smd"
                        print(
                            """WARNING: solvent model for part 1 and 2 is not recognized! Options for ORCA are CPCM
                              (cpcm) or SMD (smd). solvent model is set to ORCA default SMD."""
                        )
        with open(os.path.join(os.getcwd(), "flags.dat"), "w", newline=None) as inp:
            inp.write("FLAGS\n")
            if args.nstruc:
                inp.write(
                    "nconf: {:{digits}} # all or integer between 0 and total number of conformers\n".format(
                        str(args.nstruc), digits=40 - len("nconf")
                    )
                )
            else:
                inp.write(
                    "nconf: {:{digits}} # all or integer between 0 and total number of conformers\n".format(
                        "all", digits=40 - len("nconf")
                    )
                )
            inp.write(
                "charge: {:{digits}} # integer\n".format(
                    str(args.chrg), digits=40 - len("charge")
                )
            )
            inp.write(
                "unpaired: {:{digits}} # integer\n".format(
                    str(args.unpaired), digits=40 - len("unpaired")
                )
            )
            if args.solv:
                inp.write(
                    "solvent: {:{digits}} # {} \n".format(
                        str(args.solv), ", ".join(solvents), digits=40 - len("solvent")
                    )
                )
            else:
                inp.write(
                    "solvent: {:{digits}} # {} \n".format(
                        "gas", ", ".join(solvents), digits=40 - len("solvent")
                    )
                )
            inp.write(
                "prog: {:{digits}} # tm, orca\n".format(
                    str(args.prog), digits=40 - len("prog")
                )
            )
            inp.write(
                "ancopt: {:{digits}} # on, off\n".format(
                    str(args.ancopt), digits=40 - len("ancopt")
                )
            )
            if args.rrhoprog:
                inp.write(
                    "prog_rrho: {:{digits}} # xtb, prog\n".format(
                        str(args.rrhoprog), digits=40 - len("prog_rrho")
                    )
                )
            else:
                inp.write(
                    "prog_rrho: {:{digits}} # xtb, prog\n".format(
                        "prog", digits=40 - len("prog_rrho")
                    )
                )
            inp.write(
                "gfn_version: {:{digits}} # {}\n".format(
                    str(args.gfnv), ", ".join(gfnv), digits=40 - len("gfn_version")
                )
            )
            if args.temperature:
                inp.write(
                    "temperature: {:{digits}} # in Kelvin\n".format(
                        str(args.temperature), digits=40 - len("temperature")
                    )
                )
            else:
                inp.write(
                    "temperature: {:{digits}} # in Kelvin\n".format(
                        self.enso_internal_defaults["temperature"],
                        digits=40 - len("temperature"),
                    )
                )
            if args.prog3:
                inp.write(
                    "prog3: {:{digits}} # tm, orca, prog\n".format(
                        str(args.prog3), digits=40 - len("prog3")
                    )
                )
            else:
                inp.write(
                    "prog3: {:{digits}} # tm, orca, prog\n".format(
                        "prog", digits=40 - len("prog3")
                    )
                )
            if args.prog4:
                inp.write(
                    "prog4: {:{digits}} # prog, tm, orca\n".format(
                        str(args.prog4), digits=40 - len("prog4")
                    )
                )
            else:
                inp.write(
                    "prog4: {:{digits}} # tm, orca\n".format(
                        "prog", digits=40 - len("prog4")
                    )
                )
            inp.write(
                "part1: {:{digits}} # on, off\n".format(
                    str(args.part1), digits=40 - len("part1")
                )
            )
            inp.write(
                "part2: {:{digits}} # on, off\n".format(
                    str(args.part2), digits=40 - len("part2")
                )
            )
            inp.write(
                "part3: {:{digits}} # on, off\n".format(
                    str(args.part3), digits=40 - len("part3")
                )
            )
            inp.write(
                "part4: {:{digits}} # on, off\n".format(
                    str(args.part4), digits=40 - len("part4")
                )
            )
            if args.boltzmann == "on":
                inp.write(
                    "boltzmann: {:{digits}} # on, off\n".format(
                        "on", digits=40 - len("boltzmann")
                    )
                )
            elif args.boltzmann == "off" or not args.boltzmann:
                inp.write(
                    "boltzmann: {:{digits}} # on, off\n".format(
                        "off", digits=40 - len("boltzmann")
                    )
                )
            if args.backup == "on":
                inp.write(
                    "backup: {:{digits}} # on, off\n".format(
                        "on", digits=40 - len("backup")
                    )
                )
            elif args.backup == "off" or not args.backup:
                inp.write(
                    "backup: {:{digits}} # on, off\n".format(
                        "off", digits=40 - len("backup")
                    )
                )
            inp.write(
                "func: {:{digits}} # {} \n".format(
                    str(args.func), ", ".join(func), digits=40 - len("func")
                )
            )
            inp.write(
                "func3: {:{digits}} # {} \n".format(
                    str(args.func3), ", ".join(func3), digits=40 - len("func3")
                )
            )
            inp.write(
                "basis3: {:{digits}} #  \n".format(
                    str(args.basis3), digits=40 - len("basis3")
                )
            )
            inp.write(
                "couplings: {:{digits}} # on, off\n".format(
                    str(args.calcJ), digits=40 - len("couplings")
                )
            )
            inp.write(
                "funcJ: {:{digits}} # {} \n".format(
                    str(args.funcJ), ", ".join(funcJ), digits=40 - len("funcJ")
                )
            )
            inp.write(
                "basisJ: {:{digits}} #  \n".format(
                    str(args.basisJ), digits=40 - len("basisJ")
                )
            )
            inp.write(
                "shieldings: {:{digits}} # on, off\n".format(
                    str(args.calcS), digits=40 - len("shieldings")
                )
            )
            inp.write(
                "funcS: {:{digits}} # {} \n".format(
                    str(args.funcS), ", ".join(funcS), digits=40 - len("funcS")
                )
            )
            inp.write(
                "basisS: {:{digits}} # \n".format(
                    str(args.basisS), digits=40 - len("basisS")
                )
            )
            inp.write(
                "part1_threshold: {:{digits}} # integer or real number\n".format(
                    str(args.thr1), digits=40 - len("part1_threshold")
                )
            )
            inp.write(
                "part2_threshold: {:{digits}} # integer or real number\n".format(
                    str(args.thr2), digits=40 - len("part2_threshold")
                )
            )
            if args.sm:
                inp.write(
                    "sm: {:{digits}} # cosmo, dcosmors, cpcm, smd\n".format(
                        str(args.sm), digits=40 - len("sm")
                    )
                )
            else:
                inp.write(
                    "sm: {:{digits}} # cosmo, dcosmors, cpcm, smd\n".format(
                        "default", digits=40 - len("sm")
                    )
                )
            if args.gsolv2:
                inp.write(
                    "smgsolv2: {:{digits}} # sm, cosmors\n".format(
                        str(args.gsolv2), digits=40 - len("smgsolv2")
                    )
                )
            else:
                inp.write(
                    "smgsolv2: {:{digits}} # sm, cosmors\n".format(
                        "sm", digits=40 - len("smgsolv2")
                    )
                )
            if args.sm3:
                inp.write(
                    "sm3: {:{digits}} # dcosmors, cosmors, smd\n".format(
                        str(args.sm3), digits=40 - len("sm3")
                    )
                )
            else:
                inp.write(
                    "sm3: {:{digits}} # cosmors, smd\n".format(
                        "default", digits=40 - len("sm3")
                    )
                )
            if args.sm4:
                inp.write(
                    "sm4: {:{digits}} # cosmo, cpcm, smd\n".format(
                        str(args.sm4), digits=40 - len("sm4")
                    )
                )
            else:
                inp.write(
                    "sm4: {:{digits}} # cosmo, cpcm, smd\n".format(
                        "default", digits=40 - len("sm4")
                    )
                )
            if args.check == "on":
                inp.write(
                    "check: {:{digits}} # on, off\n".format(
                        "on", digits=40 - len("check")
                    )
                )
            elif args.check == "off":
                inp.write(
                    "check: {:{digits}} # on, off\n".format(
                        "off", digits=40 - len("check")
                    )
                )
            if args.crestcheck == "on":
                inp.write(
                    "crestcheck: {:{digits}} # on, off\n".format(
                        "on", digits=40 - len("crestcheck")
                    )
                )
            elif args.crestcheck == "off":
                inp.write(
                    "crestcheck: {:{digits}} # on, off\n".format(
                        "off", digits=40 - len("crestcheck")
                    )
                )
            inp.write(
                "maxthreads: {:{digits}} # integer larger than 0\n".format(
                    str(args.maxthreads), digits=40 - len("maxthreads")
                )
            )
            inp.write(
                "omp: {:{digits}} # integer larger than 0\n".format(
                    str(args.omp), digits=40 - len("omp")
                )
            )
            inp.write(
                "reference for 1H: {:{digits}} # {}\n".format(
                    str(args.href), ", ".join(href), digits=40 - len("reference for 1H")
                )
            )
            inp.write(
                "reference for 13C: {:{digits}} # {}\n".format(
                    str(args.cref),
                    ", ".join(cref),
                    digits=40 - len("reference for 13C"),
                )
            )
            inp.write(
                "reference for 19F: {:{digits}} # {}\n".format(
                    str(args.fref),
                    ", ".join(fref),
                    digits=40 - len("reference for 19F"),
                )
            )
            inp.write(
                "reference for 31P: {:{digits}} # {}\n".format(
                    str(args.pref),
                    ", ".join(pref),
                    digits=40 - len("reference for 31P"),
                )
            )
            inp.write(
                "reference for 29Si: {:{digits}} # {}\n".format(
                    str(args.siref),
                    ", ".join(siref),
                    digits=40 - len("reference for 29Si"),
                )
            )
            inp.write(
                "1H active: {:{digits}} # on, off\n".format(
                    str(args.hactive), digits=40 - len("1H active")
                )
            )
            inp.write(
                "13C active: {:{digits}} # on, off\n".format(
                    str(args.cactive), digits=40 - len("13C active")
                )
            )
            inp.write(
                "19F active: {:{digits}} # on, off\n".format(
                    str(args.factive), digits=40 - len("19F active")
                )
            )
            inp.write(
                "31P active: {:{digits}} # on, off\n".format(
                    str(args.pactive), digits=40 - len("31P active")
                )
            )
            inp.write(
                "29Si active: {:{digits}} # on, off\n".format(
                    str(args.siactive), digits=40 - len("29Si active")
                )
            )
            inp.write(
                "resonance frequency: {:{digits}} # integer\n".format(
                    str(args.mf), digits=40 - len("resonance frequency")
                )
            )
            inp.write("end\n")
        return

    def read_flags(
        self,
        args,
        solvents,
        gfnv,
        func,
        func3,
        funcJ,
        funcS,
        href,
        cref,
        fref,
        pref,
        siref,
        path,
        startread,
    ):
        """Read and check flags.dat"""
        print("\nReading information from flags.dat.")
        self._read_config(path, startread)
        error_logical = (
            False
        )  # is set true if error occurs and program has to be aborted
        if "end" in self.configdata:
            del self.configdata["end"]
        readinkeys = []
        for item in self.configdata:
            if item not in self.known_keys:
                print("WARNING: {} is not a known keyword in flags.dat.".format(item))
            else:
                readinkeys.append(item)
        diff = list(set(self.known_keys) - set(readinkeys))
        if diff:
            print(
                "WARNING: These keywords were not found in the configuration "
                "file flags.dat\n         and therefore default "
                "values are taken for:"
            )
            for item in diff:
                print("         {}".format(item))
                self.configdata[item] = self.enso_internal_defaults[item]
        if "part1" in self.configdata:
            if self.configdata["part1"] in ["on", "off"]:
                args.part1 = self.configdata["part1"]
            else:
                print(
                    "\nERROR: Part 1 is not recognized! Options are {}.".format(
                        self.value_options["part1"]
                    )
                )
                error_logical = True
        if "part2" in self.configdata:
            if self.configdata["part2"] in ["on", "off"]:
                args.part2 = self.configdata["part2"]
            else:
                print(
                    "\nERROR: Part 2 is not recognized! Options are {}.".format(
                        self.value_options["part2"]
                    )
                )
                error_logical = True
        if "part3" in self.configdata:
            if self.configdata["part3"] in ["on", "off"]:
                args.part3 = self.configdata["part3"]
            else:
                print(
                    "\nERROR: Part 3 is not recognized! Options are {}.".format(
                        self.value_options["part3"]
                    )
                )
                error_logical = True
        if "part4" in self.configdata:
            if self.configdata["part4"] in ["on", "off"]:
                args.part4 = self.configdata["part4"]
            else:
                print(
                    "\nERROR: Part 4 is not recognized! Options are {}.".format(
                        self.value_options["part4"]
                    )
                )
                error_logical = True
        if "nconf" in self.configdata:
            if self.configdata["nconf"] == "all":
                args.nstruc = None
            else:
                try:
                    args.nstruc = int(self.configdata["nconf"])
                except ValueError:
                    if args.part1 == "on":
                        print(
                            "\nERROR: Can not convert number of conformers "
                            "(nconf)! Number of conformers has to be "
                            "an integer or set to all."
                        )
                        error_logical = True
        if "charge" in self.configdata:
            try:
                args.chrg = int(self.configdata["charge"])
            except ValueError:
                print(
                    "\nERROR: Can not read the charge information! Charge has to"
                    " be an integer."
                )
                error_logical = True
        if "unpaired" in self.configdata:
            try:
                args.unpaired = int(self.configdata["unpaired"])
            except ValueError:
                print(
                    "\nERROR: Can not convert number of unpaired electrons! Number "
                    "of unpaired electrons has to be an integer."
                )
                error_logical = True
        if "solvent" in self.configdata:
            if self.configdata["solvent"] == "gas":
                args.solv = None
            elif self.configdata["solvent"] in solvents:
                args.solv = self.configdata["solvent"]
            else:
                print(
                    '\nERROR: Solvent "{}" is not implemented! Options are {}.'.format(
                        self.configdata["solvent"], self.value_options["solvent"]
                    )
                )
                error_logical = True
        if "prog" in self.configdata:
            if self.configdata["prog"] in ["tm", "orca"]:
                args.prog = self.configdata["prog"]
            else:
                print(
                    "\nERROR: Prog is not recognized! Options for the main "
                    "program are ORCA (orca) or TURBOMOLE (tm)."
                )
                error_logical = True
        if "prog3" in self.configdata:
            if self.configdata["prog3"] in ["tm", "orca"]:
                args.prog3 = self.configdata["prog3"]
            elif self.configdata["prog3"] == "prog":
                args.prog3 = args.prog
            else:
                print(
                    "\nERROR: Prog3 is not recognized! Options for the main "
                    "program are ORCA (orca) or TURBOMOLE (tm)."
                )
                error_logical = True
        if "prog4" in self.configdata:
            if self.configdata["prog4"] in ["tm", "orca"]:
                args.prog4 = self.configdata["prog4"]
            elif self.configdata["prog4"] == "prog":
                args.prog4 = args.prog
            else:
                if args.part4 == "on":
                    print(
                        "\nERROR: Program for part 4 is not recognized! Options "
                        "are ORCA (orca), TURBOMOLE (tm), "
                        "or the main program (prog)."
                    )
                    error_logical = True
                else:
                    print(
                        "\nWARNING: Program for part 4 is not recognized! Options "
                        "are ORCA (orca), TURBOMOLE (tm), "
                        "or the main program (prog)."
                    )
        if "ancopt" in self.configdata:
            if self.configdata["ancopt"] in ["on", "off"]:
                args.ancopt = self.configdata["ancopt"]
            else:
                print(
                    "\nWARNING: The use of the ANCOPT-optimizer implemented in "
                    "GFN-xTB is not recognized! Options are on or off."
                )
                error_logical = True
        if "prog_rrho" in self.configdata:
            if self.configdata["prog_rrho"] == "xtb":
                args.rrhoprog = self.configdata["prog_rrho"]
            elif self.configdata["prog_rrho"] == "off":
                args.rrhoprog = self.configdata["prog_rrho"]
                print(
                    "WARNING: Thermostatistical contribution to free energy will"
                    "not be calculated!"
                )
            elif (
                self.configdata["prog_rrho"] in ["prog", "orca"] and args.prog == "orca"
            ):
                args.rrhoprog = "orca"
            elif self.configdata["prog_rrho"] in ["prog", "tm"] and args.prog == "tm":
                if shutil.which("thermo") is not None:
                    print("Found thermo, therefore tm hessian calculation allowed!")
                    args.rrhoprog = "tm"
                else:
                    args.rrhoprog = "xtb"
                    print(
                        "WARNING: Currently are only GFN-xTB hessians possible "
                        "and no TM hessians"
                    )
            else:
                print(
                    "\nERROR: Program for RRHO contribution in part 2 and 3 "
                    "(rrhoprog) is not recognized!\n       Options are "
                    "GFN-xTB (xtb) or the main program (prog)."
                )
                error_logical = True
        if "gfn_version" in self.configdata:
            if self.configdata["gfn_version"] in gfnv:
                args.gfnv = self.configdata["gfn_version"]
            else:
                print(
                    "\nWARNING: GFN-xTB version is not recognized! Options are "
                    "{}.\nGFN2-xTB is used as default.".format(", ".join(gfnv))
                )
                args.gfnv = "gfn2"
        if "temperature" in self.configdata:
            try:
                float(self.configdata["temperature"])
                args.temperature = self.configdata["temperature"]
            except ValueError:
                print(
                    "\nWARNING: Temperature {} could not be converted to a "
                    "float! Using 298.15 K instead as default!".format(
                        self.configdata["temperature"]
                    )
                )
                args.temperature = "298.15"
        if "boltzmann" in self.configdata:
            if self.configdata["boltzmann"] == "off":
                args.boltzmann = False
            elif self.configdata["boltzmann"] == "on":
                args.boltzmann = True
                args.part1 = "off"
                args.part2 = "off"
                args.part3 = "on"
                args.part4 = "off"
            else:
                print(
                    "\nERROR: Boltzmann is not recognized! Options are off for "
                    "calculating everything or on for "
                    "calculating only the boltzmann population in part3."
                )
                error_logical = True
        if "backup" in self.configdata:
            if self.configdata["backup"] == "off":
                args.backup = False
            elif self.configdata["backup"] == "on":
                args.backup = True
                args.part1 = "off"
            else:
                print(
                    "\nERROR: Backup is not recognized! Options are off for the "
                    "normal procedure or on for "
                    "reevaluating the conformers sorted out."
                )
                error_logical = True
        if "func" in self.configdata:
            if self.configdata["func"] in func:
                args.func = self.configdata["func"]
            else:
                print(
                    "\nERROR: Chosen functional for part 1 and 2 (func) is not "
                    "implemented! Options are {}.".format(", ".join(func))
                )
                error_logical = True
        if "func3" in self.configdata:
            if self.configdata["func3"] in func3:
                args.func3 = self.configdata["func3"]
                if args.prog3 == "tm" and args.func3 != "pw6b95":
                    args.func3 = "pw6b95"
                    print(
                        "\nWARNING: Only PW6B95 is implemented as functional for"
                        " part 3 for TM!"
                    )
            else:
                print(
                    "\nERROR: Chosen functional for part 3 (func3) is not "
                    "implemented! Options are {}.".format(", ".join(func3))
                )
                error_logical = True
        if "basis3" in self.configdata:
            try:
                args.basis3 = str(self.configdata["basis3"])
                if args.basis3 not in self.knownbasissets3:
                    print(
                        "WARNING! Basis for part3: {} is used but could not be "
                        "checked for correct input!".format(args.basis3)
                    )
            except:
                print("\nERROR: Can not read in basis set for part 3 (basis3)!")
                error_logical = True
        if "couplings" in self.configdata:
            if self.configdata["couplings"] in ["on", "off"]:
                args.calcJ = self.configdata["couplings"]
            else:
                print(
                    '\nERROR: Couplings is not recognized! Options are "on" or "off".'
                )
                error_logical = True
        if "basisJ" in self.configdata:
            try:
                args.basisJ = str(self.configdata["basisJ"])
            except:
                print(
                    "\nERROR: Can not read in basis set for coupling calculation (basisJ)!"
                )
                error_logical = True
            if args.basisJ == "default":
                args.basisJ = "def2-TZVP"
            if args.basisJ not in self.knownbasissetsJ:
                print(
                    "WARNING! Basis for coupling calculation: {} is used but "
                    "could not be checked for correct "
                    "input!".format(args.basisJ)
                )
        if "funcJ" in self.configdata:
            if self.configdata["funcJ"] in funcJ:
                args.funcJ = self.configdata["funcJ"]
            else:
                print(
                    "\nERROR: Chosen functional for coupling calculation (funcJ)"
                    " is not implemented! Options are {}.".format(", ".join(funcJ))
                )
                error_logical = True
            if self.configdata["funcJ"] == "pbeh-3c":
                args.basisJ = "def2-mSVP"
        if "shieldings" in self.configdata:
            if self.configdata["shieldings"] in ["on", "off"]:
                args.calcS = self.configdata["shieldings"]
            else:
                print(
                    '\nERROR: Shieldings is not recognized! Options are "on" or "off".'
                )
                error_logical = True
        if "basisS" in self.configdata:
            try:
                args.basisS = str(self.configdata["basisS"])
            except:
                print(
                    "\nERROR: Can not read in basis set for shielding calculation (basisS)!"
                )
                error_logical = True
            if args.basisS == "default":
                args.basisS = "def2-TZVP"
            if args.basisS not in self.knownbasissetsS:
                print(
                    "WARNING! Basis for shielding calculation: {} is used but could not be checked for correct "
                    "input!".format(args.basisS)
                )
        if "funcS" in self.configdata:
            if self.configdata["funcS"] in funcS:
                args.funcS = self.configdata["funcS"]
                if args.funcS == "dsd-blyp" and args.prog4 == "tm":
                    print(
                        "\nError: DSD-BLYP for shielding calculations is only "
                        "possible with ORCA!"
                    )
                    error_logical = True
            else:
                print(
                    "\nERROR: Chosen functional for shielding calculation (funcS)"
                    " is not implemented! Options are {}.".format(", ".join(funcS))
                )
                error_logical = True
                if args.debug:
                    args.funcS = self.configdata["funcS"]
            if self.configdata["funcS"] == "pbeh-3c":
                args.basisS = "def2-mSVP"
        if "part1_threshold" in self.configdata:
            try:
                args.thr1 = float(self.configdata["part1_threshold"])
            except ValueError:
                print(
                    "\nERROR: Threshold for part 1 could not be converted! Threshold has to be number."
                )
                error_logical = True
        if "part2_threshold" in self.configdata:
            try:
                args.thr2 = float(self.configdata["part2_threshold"])
            except ValueError:
                print(
                    "\nERROR: Threshold for part 2 could not be converted! Threshold has to be number."
                )
                error_logical = True
        if "sm" in self.configdata:
            if args.prog == "tm" and self.configdata["sm"] in ("cosmo", "dcosmors"):
                args.sm = self.configdata["sm"]
            elif args.prog == "tm" and self.configdata["sm"] == "cpcm":
                args.sm = "cosmo"
                print(
                    "WARNING: CPCM is not available with TM! COSMO is used instead as solvent model for part 1 and 2."
                )
            elif args.prog == "tm" and self.configdata["sm"] not in (
                "cosmo",
                "dcosmors",
                "default",
            ):
                args.sm = "dcosmors"
                print(
                    "WARNING: Solvent model for part1 and 2 is not recognized! DCOSMO-RS is used as default with TM."
                )
            elif args.prog == "orca" and self.configdata["sm"] in ("cpcm", "smd"):
                args.sm = self.configdata["sm"]
            elif args.prog == "orca" and self.configdata["sm"] == "cosmo":
                args.sm = "cpcm"
                print(
                    "WARNING: COSMO is not available with ORCA! CPCM is used instead as solvent model for part 1 and 2."
                )
            elif args.prog == "orca" and self.configdata["sm"] not in (
                "cpcm",
                "smd",
                "default",
            ):
                args.sm = "smd"
                print(
                    "WARNING: Solvent model for part1 and 2 is not recognized! SMD is used as default with ORCA."
                )
            elif self.configdata["sm"] == "default":
                if args.prog == "tm":
                    args.sm = "dcosmors"
                elif args.prog == "orca":
                    args.sm = "smd"
        if "smgsolv2" in self.configdata:
            if self.configdata["smgsolv2"] in ("cosmors", "sm", "gbsa_gsolv"):
                args.gsolv2 = self.configdata["smgsolv2"]
            else:
                args.gsolv2 = args.sm
                print(
                    "WARNING: Solvent model for gsolv contribution of part 2 is "
                    "not recognized! The same solvent "
                    "model as for the optimizations in part 2 is used: {}.".format(
                        args.gsolv2
                    )
                )
        if "sm3" in self.configdata:
            if args.solv is None:
                args.sm3 = None
            elif args.prog3 == "orca":
                if self.configdata["sm3"] in ("cosmors", "smd", "gbsa_gsolv"):
                    args.sm3 = self.configdata["sm3"]
                elif self.configdata["sm3"] == "default":
                    args.sm3 = "smd"
                else:
                    args.sm3 = "smd"
                    print(
                        "WARNING: Solvent model for part 3 is not recognized! Options are SMD (smd) or COSMO-RS "
                        "(cosmors) for ORCA. Solvent model for part 3 is set to SMD."
                    )
            elif args.prog3 == "tm":
                if self.configdata["sm3"] in ("cosmors", "dcosmors", "gbsa_gsolv"):
                    args.sm3 = self.configdata["sm3"]
                elif self.configdata["sm3"] == "default":
                    args.sm3 = "dcosmors"
                else:
                    args.sm3 = "dcosmors"
                    print(
                        "WARNING: solvent model for part 3 is not recognized! Options are DCOSMO-RS (dcosmors) "
                        "or COSMO-RS (cosmors) for TM. Solvent model for part 3 is set to DCOSMO-RS."
                    )
        if "sm4" in self.configdata:
            if args.solv is None:
                args.sm4 = None
            elif args.prog4 == "tm":
                if self.configdata["sm4"] in ("dcosmors", "cosmo"):
                    args.sm4 = self.configdata["sm4"]
                elif self.configdata["sm4"] == "default":
                    args.sm4 = "dcosmors"
                else:
                    args.sm4 = "dcosmors"
                    print(
                        "WARNING: solvent model for part 4 is not recognized! Options are DCOSMO-RS (dcosmors) "
                        "or COSMO for TM. Solvent model for part 4 is set to DCOSMO-RS."
                    )
            elif args.prog4 == "orca":
                if self.configdata["sm4"] in ("smd", "cpcm"):
                    args.sm4 = self.configdata["sm4"]
                elif self.configdata["sm4"] == "default":
                    args.sm4 = "smd"
                else:
                    args.sm4 = "smd"
                    print(
                        "WARNING: Solvent model for part 4 is not recognized! Options are SMD (smd) or CPCM "
                        "for ORCA. Solvent model for part 4 is set to SMD."
                    )
        if "check" in self.configdata:
            if self.configdata["check"] == "on":
                args.check = True
            elif self.configdata["check"] == "off":
                args.check = False
            else:
                print(
                    "WARNING: Chosen option for cautious searching for errors and exiting if too many errors occur ("
                    "check) is not recognized! Options are on or off. Going to switch it off."
                )
                args.check = False
        if "crestcheck" in self.configdata:
            if self.configdata["crestcheck"] == "on":
                args.crestcheck = True
            elif self.configdata["crestcheck"] == "off":
                args.crestcheck = False
            else:
                print(
                    "WARNING: Chosen option for DFT-ensemble check using CREST is not recognized! Options are on or "
                    "off. Going to switch it off."
                )
                args.crestcheck = False
        if "maxthreads" in self.configdata:
            try:
                args.maxthreads = int(self.configdata["maxthreads"])
            except ValueError:
                print(
                    "\nERROR: Can not convert maxthreads! Maxthreads has to be an integer."
                )
                error_logical = True
        if "omp" in self.configdata:
            try:
                args.omp = int(self.configdata["omp"])
            except ValueError:
                print("\nERROR: Can not convert OMP! OMP has to be an integer.")
                error_logical = True
        spectrumlist = []
        if "reference for 1H" in self.configdata:
            if self.configdata["reference for 1H"] in href:
                args.href = self.configdata["reference for 1H"]
            else:
                print(
                    "\nWARNING: Reference for 1H is not implemented! Choices are {}. TMS is used.".format(
                        ", ".join(href)
                    )
                )
                args.href = "TMS"
        if "reference for 13C" in self.configdata:
            if self.configdata["reference for 13C"] in cref:
                args.cref = self.configdata["reference for 13C"]
            else:
                print(
                    "\nWARNING: Reference for 13C is not implemented! Choices are {}. TMS is used.".format(
                        ", ".join(cref)
                    )
                )
                args.cref = "TMS"
        if "reference for 19F" in self.configdata:
            if self.configdata["reference for 19F"] in fref:
                args.fref = self.configdata["reference for 19F"]
            elif self.configdata["reference for 19F"] == "default":
                args.fref = "CFCl3"
            else:
                print(
                    "\nWARNING: Reference for 19F is not implemented! Choices are {}. CFCl3 is used.".format(
                        ", ".join(fref)
                    )
                )
                args.fref = "CFCl3"
        if "reference for 31P" in self.configdata:
            if self.configdata["reference for 31P"] in pref:
                args.pref = self.configdata["reference for 31P"]
            elif self.configdata["reference for 31P"] == "default":
                args.pref = "TMP"
            else:
                print(
                    "\nWARNING: Reference for 31P is not implemented! Choices are {}. TMP is used.".format(
                        ", ".join(pref)
                    )
                )
                args.pref = "TMP"
        if "reference for 29Si" in self.configdata:
            if self.configdata["reference for 29Si"] in siref:
                args.siref = self.configdata["reference for 29Si"]
            else:
                print(
                    "\nWARNING: Reference for 29Si is not implemented! Choices are {}. TMS is used.".format(
                        ", ".join(siref)
                    )
                )
                args.siref = "TMS"
        if "1H active" in self.configdata:
            if self.configdata["1H active"] == "on":
                args.hactive = self.configdata["1H active"]
                spectrumlist.append("1H")
            elif self.configdata["1H active"] == "off":
                args.hactive = self.configdata["1H active"]
            else:
                args.hactive = "on"
                spectrumlist.append("1H")
                print(
                    "\nWARNING: Can not read 1H active! Data for 1H NMR spectrum is calculated."
                )
        if "13C active" in self.configdata:
            if self.configdata["13C active"] == "on":
                args.cactive = self.configdata["13C active"]
                spectrumlist.append("13C")
            elif self.configdata["13C active"] == "off":
                args.cactive = self.configdata["13C active"]
            else:
                args.cactive = "off"
                print(
                    "\nWARNING: Can not read 13C active! Data for 13C NMR spectrum is not calculated."
                )
        if "19F active" in self.configdata:
            if self.configdata["19F active"] == "on":
                args.factive = self.configdata["19F active"]
                spectrumlist.append("19F")
            elif self.configdata["19F active"] == "off":
                args.factive = self.configdata["19F active"]
            else:
                args.factive = "off"
                print(
                    "\nWARNING: Can not read 19F active! Data for 19F NMR spectrum is not calculated."
                )
        if "31P active" in self.configdata:
            if self.configdata["31P active"] == "on":
                args.pactive = self.configdata["31P active"]
                spectrumlist.append("31P")
            elif self.configdata["31P active"] == "off":
                args.pactive = self.configdata["31P active"]
            else:
                args.pactive = "off"
                print(
                    "\nWARNING: Can not read 31P active! Data for 31P NMR spectrum is not calculated."
                )
        if "29Si active" in self.configdata:
            if self.configdata["29Si active"] == "on":
                args.siactive = self.configdata["29Si active"]
                spectrumlist.append("29Si")
            elif self.configdata["29Si active"] == "off":
                args.siactive = self.configdata["29Si active"]
            else:
                args.siactive = "off"
                print(
                    "\nWARNING: Can not read 31P active! Data for 31P NMR spectrum is not calculated."
                )
        if "resonance frequency" in self.configdata:
            try:
                args.mf = int(self.configdata["resonance frequency"])
            except ValueError:
                args.mf = None
                print(
                    "\nWARNING: Can not convert resonance frequency! Flag -mf has to be set while executing the "
                    "anmr program."
                )
        if error_logical and not args.debug:
            print("Going to exit due to input errors!")
            sys.exit(1)
        return spectrumlist


def check_for_folder(conflist, functional, debug):
    """Check if folders exist (of conformers calculated in previous run) """
    error_logical = False
    for i in conflist:
        tmp_dir = os.path.join(i, args.func)
        if not os.path.exists(tmp_dir):
            print(
                "ERROR: directory of conformer {} does not exist, although it was calculated before!".format(
                    i
                )
            )
            error_logical = True
    if error_logical and not debug:
        print("One or multiple directories are missing.")
        write_json("save_and_exit", json_dict, jsonfile)
    return


def new_folders(conflist, functional, save_errors):
    """ """
    directories = []
    for conf in list(conflist):
        if "CONF" in str(conf):
            tmp_conf = conf
        else:
            tmp_conf = "CONF" + str(conf)
        tmp_dir = [tmp_conf, os.path.join(tmp_conf, args.func)]
        directories.append(tmp_dir)
        try:
            mkdir_p(tmp_dir[1])
        except:
            if not os.path.isdir(tmp_dir[1]):
                print("ERROR: Could not create folder for {}!".format(tmp_dir[0]))
                print("ERROR: Removing {}!".format(tmp_dir[0]))
                save_errors.append(
                    "{} was removed, because IO failed!".format(tmp_dir[0])
                )
                directories.remove(tmp_dir)
                conflist.remove(conf)
    print("Constructed all {} CONF-folders.".format(len(directories)))
    return directories


def rrho_part23(
    args,
    q,
    resultq,
    job,
    in_part2,
    in_part3,
    results,
    json_dict,
    jsonfile,
    save_errors,
    digilen,
    cwd,
):
    """calculate RRHO for part2 or part3"""

    if in_part2:
        pass
    elif in_part3:
        pass

    print("\nRRHO calculation:")
    if args.solv and args.rrhoprog == "orca":
        print(
            "WARNING: ORCA RRHO calculation with solvation correction is currently not "
            "implemented! \nFor a RRHO contribution in solution, please use GFNn-xTB."
            "\nRRHO is calculated in the gas phase.\n"
        )
    elif args.solv and args.rrhoprog == "tm":
        print(
            "WARNING: TM RRHO calculation with solvation correction is currently not "
            "implemented! \nFor a RRHO contribution in solution, please use GFNn-xTB."
            "\nRRHO is calculated in the gas phase.\n"
        )
    tmp_results = []
    for conf in list(results):
        if json_dict[conf.name]["rrho"] == "calculated":
            # this conformer was already calculated successfully,
            # append rrho energy and move it to tmp_results
            conf.rrho = json_dict[conf.name]["energy_rrho"]
            conf.symmetry = json_dict[conf.name]["symmetry"]
            tmp_results.append(results.pop(results.index(conf)))
        elif json_dict[conf.name]["rrho"] == "failed":
            # this conformer was already calculated,
            # but the calculation failed, remove it
            print(
                "The RRHO calculation of part 2 failed for {} in the previous run."
                " The conformer is sorted out.".format(conf.name)
            )
            results.remove(conf)
        elif json_dict[conf.name]["rrho"] == "not_calculated":
            # this conformer has to be calculated now
            pass
    # check if there is at least one conformer
    if len(tmp_results) == 0 and len(results) == 0:
        print("ERROR: There are no conformers left!")
        write_json("save_and_exit", json_dict, jsonfile)
    print(
        "number of further considered conformers: {:{digits}} {}".format(
            "",
            str(len(results) + len(tmp_results)),
            digits=digilen - len("number of further considered conformers"),
        )
    )
    if len(tmp_results) > 0:
        print(
            "The RRHO calculation was already carried out for {} conformers:".format(
                str(len(tmp_results))
            )
        )
        print_block([i.name for i in tmp_results])
        if len(results) > 0:
            print(
                "The RRHO calculation is carried out now for {} conformers:".format(
                    str(len(results))
                )
            )
            print_block([i.name for i in results])
        else:
            print("No conformers are considered additionally.")
    else:
        print("Considered conformers:")
        print_block([i.name for i in results])

    # calculate RRHO for all conformers that remain in list results
    if len(results) > 0:
        if args.rrhoprog == "xtb":
            # set omp num threads for GFN-xTB calculation
            # os.environ["OMP_NUM_THREADS"] = "{:d}".format(args.omp)
            environsettings["OMP_NUM_THREADS"] = "{:d}".format(args.omp)
        # create new folder and copy optimized coord file in CONFx/rrho
        print("Setting up new directories.")
        for conf in list(results):
            tmp_from = os.path.join(cwd, conf.name, args.func)
            tmp_to = os.path.join(cwd, conf.name, "rrho")
            if not os.path.isdir(tmp_to):
                mkdir_p(tmp_to)
            try:
                shutil.copy(
                    os.path.join(tmp_from, "coord"), os.path.join(tmp_to, "coord")
                )
            except FileNotFoundError:
                if not os.path.isfile(os.path.join(tmp_from, "coord")):
                    print(
                        "ERROR: while copying the coord file from {}! The corresponding "
                        "file does not exist.".format(tmp_from)
                    )
                elif not os.path.isdir(tmp_to):
                    print("ERROR: Could not create folder {}!".format(tmp_to))
                print("ERROR: Removing conformer {}!".format(conf.name))
                json_dict[conf.name]["rrho"] = "failed"
                json_dict[conf.name]["energy_rrho"] = None
                if in_part2:
                    json_dict[conf.name]["consider_for_part3"] = False
                    json_dict[conf.name]["backup_part3"] = False
                elif in_part3:
                    json_dict[conf.name]["consider_for_part4"] = False
                save_errors.append(
                    "Conformer {} was removed, because IO failed!".format(conf.name)
                )
                results.remove(conf)
        print("Constructed all new folders.")

        # check if at least one conformer is left
        if len(results) == 0:
            print("ERROR: No conformers left!")
            write_json("save_and_exit", json_dict, jsonfile)

        instructrrho = {
            "jobtype": "rrhoxtb",
            "chrg": args.chrg,
            "unpaired": args.unpaired,
            "func": args.func,
            "solv": args.solv,
            "boltzmann": args.boltzmann,
            "temperature": args.temperature,
            "progsettings": {"omp": args.omp, "tempprogpath": ""},
        }
        if args.rrhoprog == "xtb":
            instructrrho["jobtype"] = "rrhoxtb"
            instructrrho["func"] = "GFN-xTB"
            instructrrho["gfnv"] = args.gfnv
            instructrrho["progsettings"]["tempprogpath"] = ""
            instructrrho["progsettings"]["xtbpath"] = xtbpath
        elif args.rrhoprog == "orca":
            instructrrho["jobtype"] = "rrhoorca"
            instructrrho["func"] = args.func
            instructrrho["progsettings"]["tempprogpath"] = orcapath
        elif args.rrhoprog == "tm":
            instructrrho["jobtype"] = "rrhotm"
            instructrrho["func"] = args.func
            instructrrho["progsettings"]["tempprogpath"] = ""

        results = run_in_parallel(
            q, resultq, job, maxthreads, results, instructrrho, "rrho"
        )
        exit_log, fail_rate = check_tasks(results, args)

        # sort out conformers with crashed optimizations
        for conf in list(results):
            if not conf.success:
                print(
                    "\nERROR: A problem has occurred in the RRHO calculation of {}!"
                    " The conformer is removed.\n".format(conf.name)
                )
                json_dict[conf.name]["rrho"] = "failed"
                json_dict[conf.name]["energy_rrho"] = None
                if in_part2:
                    json_dict[conf.name]["consider_for_part3"] = False
                    json_dict[conf.name]["backup_part3"] = False
                elif in_part3:
                    json_dict[conf.name]["consider_for_part4"] = False
                save_errors.append(
                    "Conformer {} was removed, because the RRHO calculation "
                    "failed!".format(conf.name)
                )
                results.remove(conf)

        for conf in results:
            json_dict[conf.name]["rrho"] = "calculated"
            json_dict[conf.name]["energy_rrho"] = conf.rrho
            json_dict[conf.name]["symmetry"] = conf.symmetry

        if exit_log:
            print(
                "\nERROR: too many RRHO calculations failed ({:.2f} %)!".format(
                    fail_rate
                )
            )
            write_json("save_and_exit", json_dict, jsonfile)

        if len(results) != 0 and args.rrhoprog == "xtb":
            # check rmsd between GFN-xTB RRHO and DFT structure
            rmsd = []
            for item in results:
                if item.rrho is not None:
                    workdir = str(os.path.join(cwd, os.path.join(item.name, "rrho")))
                    rmsd.append([item.name, RMSD_routine(workdir, nat)])
            for item in list(rmsd):
                if not isinstance(item[1], float):
                    rmsd.remove(item)
                    print("ERROR: could not calculate RMSD in {}".format(item[0]))
                elif item[1] < 0.4:
                    rmsd.remove(item)
            if len(rmsd) >= 1:
                print(
                    "\nWARNING: The RMSD between the DFT and the GFN-xTB structure "
                    "used for the RRHO\n contribution is larger than 0.4 Angstrom for "
                    "the following conformers:"
                )
                print("#CONF    RMSD\n")
                for line in rmsd:
                    if line[1] > 0.4:
                        print("{:8} {:.3f}\n".format(line[0], line[1]))
                        save_errors.append(
                            "WARNING: Large RMSD between DFT and "
                            "GFNn-xTB geometry for {}.".format(line[0])
                        )
    # adding conformers calculated before to results
    try:
        length = max([len(str(i.name)) for i in tmp_results])
    except:
        length = 4 + int(args.nstruc)
    for item in tmp_results:
        print(
            "RRHO of {:>{digits}}: {:.7f} in sym: {}".format(
                item.name, float(item.rrho), item.symmetry, digits=length
            )
        )
        results.append(item)
    results.sort(key=lambda x: int(x.name[4:]))

    # check if at least one conformer is left
    if len(results) == 0:
        print("ERROR: No conformers left!")
        write_json("save_and_exit", json_dict, jsonfile)

    return results, save_errors


def additive_gsolv(
    args,
    q,
    resultq,
    job,
    in_part2,
    in_part3,
    results,
    json_dict,
    jsonfile,
    save_errors,
    digilen,
):
    """calculate GBSA-RS or COSMO-RS for part2 or part3"""

    if in_part2:
        if args.gsolv2 == "cosmors":
            js_smodel = "cosmo-rs"
            js_sm_energy = "energy_cosmo-rs"
            folder = "gsolv"
            sm_capital = "COSMO-RS"
        elif args.gsolv2 == "gbsa_gsolv":
            js_smodel = "gbsa_gsolv"
            js_sm_energy = "energy_gbsa_gsolv"
            folder = "gbsa_gsolv"
            sm_capital = "GBSA-Gsolv"
        else:
            print("ERROR: If this occured, we are on the wrong track.")
        part = "part2"
    elif in_part3:
        if args.sm3 == "cosmors":
            js_smodel = "cosmo-rs"
            js_sm_energy = "energy_cosmo-rs"
            folder = "gsolv"
            sm_capital = "COSMO-RS"
        elif args.sm3 == "gbsa_gsolv":
            js_smodel = "gbsa_gsolv"
            js_sm_energy = "energy_gbsa_gsolv"
            folder = "gbsa_gsolv"
            sm_capital = "GBSA-Gsolv"
        else:
            print("ERROR: If this occured, we are on the wrong track.")
        part = "part3"

    print("\nCalculating {} contribution to free energy.".format(sm_capital))
    tmp_results = []
    for i in list(results):
        if json_dict[i.name][js_smodel] == "calculated":
            # this conformer was already calculated successfully, it is moved to tmp_results
            i.gsolv = json_dict[i.name][js_sm_energy]
            tmp_results.append(results.pop(results.index(i)))
        elif json_dict[i.name][js_smodel] == "failed":
            # this conformer was already calculated but the calculation failed, it is removed
            print(
                "The {} calculation of {} failed for {} in the previous run. The conformer is "
                "sorted out.".format(sm_capital, part, i.name)
            )
            results.remove(i)
        elif json_dict[i.name][js_smodel] == "not_calculated":
            # this conformer has to be calculated now
            pass
    print(
        "number of further considered conformers: {:{digits}} {}".format(
            "",
            str(len(results) + len(tmp_results)),
            digits=digilen - len("number of further considered conformers"),
        )
    )
    if len(tmp_results) > 0:
        print(
            "The {} calculation was already carried out for {} conformers:".format(
                sm_capital, str(len(tmp_results))
            )
        )
        print_block([i.name for i in tmp_results])
        if len(results) > 0:
            print(
                "The {} calculation is carried out now for {} conformers:".format(
                    sm_capital, str(len(results))
                )
            )
            print_block([i.name for i in results])
        else:
            print("No conformers are considered additionally.")
    else:
        print("Considered conformers:")
        print_block([i.name for i in results])

    if len(results) > 0:
        # calculate solvation contribution for all conformers that remain in list results
        print("Setting up new directories.")
        for conf in list(results):
            tmp_from = os.path.join(cwd, conf.name, args.func)
            tmp_to = os.path.join(cwd, conf.name, folder)
            if not os.path.isdir(tmp_to):
                mkdir_p(tmp_to)
            try:
                shutil.copy(
                    os.path.join(tmp_from, "coord"), os.path.join(tmp_to, "coord")
                )
            except FileNotFoundError:
                if not os.path.isfile(os.path.join(tmp_from, "coord")):
                    print(
                        "ERROR: while copying the coord file from {}! The corresponding "
                        "file does not exist.".format(tmp_from)
                    )
                elif not os.path.isdir(tmp_to):
                    print("ERROR: Could not create folder {}!".format(tmp_to))
                print("ERROR: Removing conformer {}!".format(conf.name))
                json_dict[conf.name][js_smodel] = "failed"
                json_dict[conf.name][js_sm_energy] = None
                if in_part2:
                    json_dict[conf.name]["consider_for_part3"] = False
                    json_dict[conf.name]["backup_part3"] = False
                if in_part3:
                    json_dict[conf.name]["consider_for_part4"] = False
                save_errors.append(
                    "Conformer {} was removed, because IO failed!".format(conf.name)
                )
                results.remove(conf)
        print("Constructed all new folders.")
        if len(results) == 0:
            print("ERROR: No conformers left!")
            write_json("save_and_exit", json_dict, jsonfile)
        if js_smodel == "cosmo-rs":
            instructsolv = {
                "jobtype": "solv",
                "chrg": args.chrg,
                "unpaired": args.unpaired,
                "solv": args.solv,
                "boltzmann": args.boltzmann,
                "temperature": args.temperature,
                "progsettings": {
                    "omp": args.omp,
                    "tempprogpath": "",
                    "cosmorssetup": cosmorssetup,
                    "cosmothermversion": cosmothermversion,
                },
            }
        elif js_smodel == "gbsa_gsolv":
            instructsolv = {
                "jobtype": "gbsa_gsolv",
                "chrg": args.chrg,
                "unpaired": args.unpaired,
                "solv": args.solv,
                "gfnv": args.gfnv,
                "boltzmann": args.boltzmann,
                "temperature": args.temperature,
                "progsettings": {
                    "omp": args.omp,
                    "tempprogpath": xtbpath,
                    "xtbpath": xtbpath,
                },
            }
        results = run_in_parallel(
            q, resultq, job, maxthreads, results, instructsolv, folder
        )
        exit_log, fail_rate = check_tasks(results, args)

        for i in list(results):
            if not i.success:
                print(
                    "\nERROR: A problem has occurred in the {} calculation of {}! The conformer is "
                    "removed.\n".format(sm_capital, i.name)
                )
                json_dict[i.name][js_smodel] = "failed"
                json_dict[i.name][js_sm_energy] = None
                if in_part2:
                    json_dict[i.name]["consider_for_part3"] = False
                    json_dict[i.name]["backup_part3"] = False
                if in_part3:
                    json_dict[i.name]["consider_for_part4"] = False
                save_errors.append(
                    "Error in {} calculation for conformer {}".format(
                        sm_capital, i.name
                    )
                )
                results.remove(i)
        # write energies to json_dict
        for i in results:
            json_dict[i.name][js_smodel] = "calculated"
            json_dict[i.name][js_sm_energy] = i.gsolv
        if exit_log:
            print(
                "\nERROR: too many {} calculations failed ({:.2f} %)!".format(
                    sm_capital, fail_rate
                )
            )
            write_json("save_and_exit", json_dict, jsonfile)
    # adding conformers calculated before to results
    try:
        length = max([len(str(i.name)) for i in tmp_results])
    except ValueError:
        length = 4 + int(args.nstruc)
    for item in tmp_results:
        print(
            "Gsolv of {:>{digits}}: {:.7f}".format(item.name, item.gsolv, digits=length)
        )
        results.append(item)
    results.sort(key=lambda x: int(x.name[4:]))
    # check if at least one conformer is left:
    if len(results) == 0:
        print("ERROR: No conformers left!")
        write_json("save_and_exit", json_dict, jsonfile)

    return results


def correct_anmr_enso(cwd, removelist):
    """remove conformers with failed calculations from anmr_enso,
    removelist contains the name of conformers to remove"""
    with open(
        os.path.join(cwd, "anmr_enso"), "r", encoding=coding, newline=None
    ) as inp:
        data = inp.readlines()
    for name in removelist:
        for line in data[1:]:
            if int(line.split()[1]) == int(name[4:]):
                index = data.index(line)
                newline = "0" + line.lstrip()[1:]
        data[index] = newline
    with open(os.path.join(cwd, "anmr_enso"), "w", newline=None) as out:
        for line in data:
            out.write(line)


def enso_startup(cwd):
    """
    1)print header
    2) read cml input,
    3) read or create ensorc
    4) read or write flags.dat control file 
    5) check for crest_conformers.xyz
    6) check program settings
    7) read or write enso.json
    """

    # all implemented functionals and solvents
    impsolvents = (
        "acetone",
        "chcl3",
        "acetonitrile",
        "ch2cl2",
        "dmso",
        "h2o",
        "methanol",
        "thf",
        "toluene",
        "gas",
    )
    impfunc = ("pbeh-3c", "b97-3c", "tpss")
    impfunc3 = ("pw6b95", "wb97x", "dsd-blyp")
    impfuncJ = ("tpss", "pbe0", "pbeh-3c")
    impfuncS = ("tpss", "pbe0", "dsd-blyp", "pbeh-3c")
    impgfnv = ("gfn1", "gfn2")
    imphref = ("TMS",)
    impcref = ("TMS",)
    impfref = ("CFCl3",)
    imppref = ("TMP", "PH3")
    impsiref = ("TMS",)
    orca_old = False

    descr = """
     __________________________________________________
    |                                                  |
    |                                                  |
    |                       ENSO -                     |
    |          energetic sorting of CREST CRE          |
    |          for automated NMR calculations          |
    |             University of Bonn, MCTC             |
    |                    July 2018                     |
    |                   version 1.300                  |
    |  F. Bohle, K. Schmitz, J. Pisarek and S. Grimme  |
    |                                                  |
    |__________________________________________________|"""

    args = cml(
        descr,
        impsolvents,
        impfunc,
        impfunc3,
        impfuncJ,
        impfuncS,
        impgfnv,
        imphref,
        impcref,
        impfref,
        imppref,
        impsiref,
    )
    print(descr + "\n")  # print header

    # read .ensorc
    if os.path.isfile(os.path.join(cwd, ".ensorc")):
        ensorc = os.path.join(cwd, ".ensorc")
    else:
        home = expanduser("~")
        ensorc = os.path.join(home, ".ensorc")
        if not os.path.isfile(ensorc):
            ensorc = ""

    orcapath = ""
    orcaversion = ""
    xtbpath = "xtb"
    crestpath = "crest"
    mpshiftpath = ""
    escfpath = ""
    cosmorssetup = ""
    # hactive = "off"
    # cactive = "off"
    # factive = "off"
    # pactive = "off"
    # siactive = "off"
    spectrumlist = []

    if args.writeensorc:
        ensorcsettings = handle_input(
            impsolvents,
            impgfnv,
            impfunc,
            impfunc3,
            impfuncJ,
            impfuncS,
            imphref,
            impcref,
            impfref,
            imppref,
            impsiref,
        )
        ensorcsettings.write_ensorc()
        print(
            "A new ensorc was written into the current directory!\n"
            "You have to adjust the settings to your needs and it is"
            " mandatory to correctly set the program paths!"
        )
        print("All done!")
        sys.exit(0)
    elif ensorc:
        print(".ensorc is taken from: {}".format(str(ensorc)))
        ensorcsettings = handle_input(
            impsolvents,
            impgfnv,
            impfunc,
            impfunc3,
            impfuncJ,
            impfuncS,
            imphref,
            impcref,
            impfref,
            imppref,
            impsiref,
        )
        (
            dbpath,
            cosmothermversion,
            cosmorssetup,
            orcapath,
            crestpath,
            xtbpath,
            mpshiftpath,
            escfpath,
            orca_old,
            orcaversion,
        ) = ensorcsettings.read_program_paths(ensorc)
        print("\nThe following pathways were read in:")
        print("    ORCA:         {}".format(orcapath))
        print("    ORCA Version: {}".format(orcaversion))
        # print('    TURBOMOLE:  {}'.format(tmpath))
        print("    GFN-xTB:      {}".format(xtbpath))
        print("    CREST:        {}".format(crestpath))
        print("    mpshift:      {}".format(mpshiftpath))
        print("    escf:         {}".format(escfpath))
        try:
            tmp = cosmorssetup.split()
            if len(tmp) == 9:
                print("    Set up of COSMO-RS:")
                print("        {}".format(" ".join(tmp[0:3])))
                print("        {}".format(" ".join(tmp[3:6])))
                print("        {}".format(" ".join(tmp[6:9])))
            else:
                print("    Set up of COSMO-RS: {}".format(str(cosmorssetup)))
        except:
            print("    Set up of COSMO-RS: {}".format(str(cosmorssetup)))
        print("    Using {}\n    as path to the COSMO-RS DATABASE.\n".format(dbpath))
        ensorcsettings.read_ensorc_data(args, ensorc, "#NMR data")
    else:
        print("ERROR: could not find the file .ensorc!")
        print(
            "File has to be either in the user's home directory,"
            "\n or in the current working directory.\n Going to exit!"
        )
        sys.exit(1)

    # Write flags.dat or run program
    if not args.run and not args.checkinput:
        print("\nWriting data to flags.dat.")
        ensorcsettings.write_flags(
            args,
            impsolvents,
            impfunc,
            impfunc3,
            impfuncJ,
            impfuncS,
            impgfnv,
            imphref,
            impcref,
            impfref,
            imppref,
            impsiref,
        )
    elif (
        args.run or args.checkinput
    ):  # if restart is 'on', flags are read from file flags.dat
        flagsdata = handle_input(
            impsolvents,
            impgfnv,
            impfunc,
            impfunc3,
            impfuncJ,
            impfuncS,
            imphref,
            impcref,
            impfref,
            imppref,
            impsiref,
        )
        spectrumlist = flagsdata.read_flags(
            args,
            impsolvents,
            impgfnv,
            impfunc,
            impfunc3,
            impfuncJ,
            impfuncS,
            imphref,
            impcref,
            impfref,
            imppref,
            impsiref,
            "flags.dat",
            "FLAGS",
        )

    # check whether crest_conformers.xyz file is available
    if os.path.isfile(os.path.join(cwd, "crest_conformers.xyz")):
        conformersxyz = os.path.join(cwd, "crest_conformers.xyz")
        print("\nUsing conformers from file crest_conformers.xyz.")
    else:
        print("\nERROR: Could not find crest_conformers.xyz file!" "\nGoing to exit.")
        sys.exit(1)

    # get nstruc and nat from conformersxyz
    try:
        with open(conformersxyz, "r", encoding=coding, newline=None) as xyzfile:
            stringfile_lines = xyzfile.readlines()
    except:
        print("\nERROR: Can not read file {}!" "\nGoing to exit.".format(conformersxyz))
        sys.exit(1)
    try:
        nat = int(stringfile_lines[0])
    except ValueError:
        print(
            "Could not get the number of atoms from file {}, something "
            "is wrong!".format(conformersxyz)
        )
        sys.exit(1)
    try:
        nelements = int(int(len(stringfile_lines)) / (int(stringfile_lines[0]) + 2))
    except ValueError:
        print("Could not determine the number of conformers from the crest ensemble!")
        sys.exit(1)

    error_logical = False
    # check if nstruc has been set otherwise read from conformersxyz
    if not args.nstruc or args.nstruc == "all":
        args.nstruc = nelements
    if args.nstruc > nelements:
        print(
            "WARNING: the chosen number of conformers is larger than the number"
            " of conformers of the crest_conformers.xyz file! \nOnly the "
            "conformers in the crest_conformers.xyz file are used."
        )
        args.nstruc = nelements
    if crestpath is None:
        print("\nERROR: path for CREST is not correct!")
        error_logical = True
    elif shutil.which(crestpath) is None:
        print("\nERROR: path for CREST is not correct!")
        error_logical = True
    if (
        args.prog == "tm"
        or args.prog4 == "tm"
        or args.gsolv2 == "cosmors"
        or args.sm3 == "cosmors"
    ):
        if shutil.which("cefine") is not None:
            print("Using cefine from {}".format(shutil.which("cefine")))
        else:
            print("\nERROR: cefine has not been found!")
            error_logical = True
    if args.prog == "orca" or args.prog3 == "orca" or args.prog4 == "orca":
        if orcapath is None:
            print("\nERROR: path for ORCA is not correct!")
            error_logical = True
        elif shutil.which(os.path.join(orcapath, "orca")) is None:
            print("\nERROR: path for ORCA is not correct!")
            error_logical = True
        if orca_old is None:
            print("\nERROR: ORCA version was not found!")
            error_logical = True
    if args.rrhoprog == "xtb" or args.ancopt == "on":
        if xtbpath is None:
            print("\nERROR: path for xTB is not correct!")
            error_logical = True
        elif shutil.which(xtbpath) is None:
            print("\nERROR: path for xTB is not correct!")
            error_logical = True
    if args.part4 == "on" and args.calcJ == "on":
        if args.prog4 == "tm" or (args.prog == "tm" and args.prog4 == "prog"):
            if escfpath is None:
                print("\nERROR: path for escf is not correct!")
                error_logical = True
            elif shutil.which(escfpath) is None:
                print("\nERROR: path for escf is not correct!")
                error_logical = True
    if args.part4 == "on" and args.calcS == "on":
        if args.prog4 == "tm" or (args.prog == "tm" and args.prog4 == "prog"):
            if mpshiftpath is None:
                print("\nERROR: path for mpshift is not correct!")
                error_logical = True
            elif shutil.which(mpshiftpath) is None:
                print("\nERROR: path for mpshift is not correct!")
                error_logical = True
    if args.gsolv2 == "cosmors" or args.sm3 == "cosmors":
        if cosmorssetup is None:
            print("\nERROR: Set up for COSMO-RS has to be written to .ensorc!")
            error_logical = True
        if cosmothermversion is None:
            print("\nERROR: Version of COSMO-RS has to be written to .ensorc!")
            error_logical = True
        home = expanduser("~")
        if shutil.which("cosmotherm") is not None:
            print("Using COSMOtherm from {}".format(shutil.which("cosmotherm")))
        else:
            print("\nERROR: COSMOtherm has not been found!")
            error_logical = True
    environsettings = os.environ
    if (
        args.prog == "tm"
        or args.prog3 == "tm"
        or args.prog4 == "tm"
        or args.sm3 == "cosmors"
        or args.gsolv2 == "cosmors"
    ):
        # preparation of parallel calculation with TM
        try:
            if environsettings["PARA_ARCH"] == "SMP":
                try:
                    environsettings["PARNODES"] = str(args.omp)
                    print(
                        "PARNODES for TM or COSMO-RS calculation was set "
                        "to {}".format(environsettings["PARNODES"])
                    )
                except:
                    print("\nERROR: PARNODES can not be changed!")
                    error_logical = True
            else:
                print(
                    "\nERROR: PARA_ARCH has to be set to SMP for parallel TM "
                    "calculations!"
                )
                if args.run:
                    error_logical = True
        except:
            print(
                "\nERROR: PARA_ARCH has to be set to SMP and PARNODES have to "
                "be set\n       for parallel TM calculations!."
            )
            if args.run:
                error_logical = True

    if args.run or args.checkinput:
        if args.part4 == "on" and args.calcJ == "off" and args.calcS == "off":
            args.part4 = "off"
            print(
                "WARNING: Neither calculating coupling nor shielding "
                "constants is activated! Part 4 is not executed."
            )
        elif len(spectrumlist) == 0:
            if args.part4 == "on":
                print(
                    "WARNING: No type of NMR spectrum is activated in the "
                    ".ensorc! Part 4 is not executed."
                )
                args.part4 = "off"
            else:
                print(
                    "WARNING: No type of NMR spectrum is activated in " "the .ensorc!"
                )
        else:
            args.fref = "CFCl3"
            args.pref = "TMP"

    if error_logical and not args.debug:
        print("\nERROR: ENSO can not continue due to input errors!" "\nGoing to exit!")
        sys.exit(1)

    # print parameter setting
    print("\n-----------------------------------------------------------")
    print(" PARAMETERS")
    print("-----------------------------------------------------------\n")

    digilen = 60
    print(
        "number of atoms in system: {:{digits}} {}".format(
            "", nat, digits=digilen - len("number of atoms in system")
        )
    )
    print(
        "number of conformers: {:{digits}} {}".format(
            "", args.nstruc, digits=digilen - len("number of conformers")
        )
    )
    print(
        "charge: {:{digits}} {}".format("", args.chrg, digits=digilen - len("charge"))
    )
    print(
        "unpaired: {:{digits}} {}".format(
            "", args.unpaired, digits=digilen - len("unpaired")
        )
    )
    if args.solv:
        print(
            "solvent: {:{digits}} {}".format(
                "", args.solv, digits=digilen - len("solvent")
            )
        )
    else:
        print(
            "solvent: {:{digits}} {}".format("", "gas", digits=digilen - len("solvent"))
        )
    print(
        "program for part1 and part2: {:{digits}} {}".format(
            "", args.prog, digits=digilen - len("program for part1 and part2")
        )
    )
    print(
        "program for part 3: {:{digits}} {}".format(
            "", args.prog3, digits=digilen - len("program for part 3")
        )
    )
    if args.prog4:
        print(
            "program for part 4: {:{digits}} {}".format(
                "", args.prog4, digits=digilen - len("program for part 4")
            )
        )
    else:
        print(
            "program for part 4: {:{digits}} {}".format(
                "", "main program", digits=digilen - len("program for part 4")
            )
        )
    print(
        "using ANCOPT implemented in GFN-xTB for the optimization: {:{digits}} {}".format(
            "",
            str(args.ancopt),
            digits=digilen
            - len("using ANCOPT implemented in GFN-xTB for the optimization"),
        )
    )
    print(
        "program for RRHO in part 2 and part 3: {:{digits}} {}".format(
            "",
            args.rrhoprog,
            digits=digilen - len("program for RRHO in part 2 and part 3"),
        )
    )
    print(
        "temperature: {:{digits}} {}".format(
            "", args.temperature, digits=digilen - len("temperature")
        )
    )
    if args.rrhoprog == "xtb":
        print(
            "GFN-xTB version for RRHO in part 2 and part 3: {:{digits}} {}".format(
                "",
                args.gfnv,
                digits=digilen - len("GFN-xTB version for RRHO in part 2 and part 3"),
            )
        )
    print(
        "part 1: {:{digits}} {}".format("", args.part1, digits=digilen - len("part 1"))
    )
    print(
        "part 2: {:{digits}} {}".format("", args.part2, digits=digilen - len("part 2"))
    )
    print(
        "part 3: {:{digits}} {}".format("", args.part3, digits=digilen - len("part 3"))
    )
    print(
        "part 4: {:{digits}} {}".format("", args.part4, digits=digilen - len("part 4"))
    )
    if args.boltzmann:
        print(
            "only boltzmann population: {:{digits}} {}".format(
                "", "on", digits=digilen - len("only boltzmann population")
            )
        )
    else:
        print(
            "only boltzmann population: {:{digits}} {}".format(
                "", "off", digits=digilen - len("only boltzmann population")
            )
        )
    if args.backup:
        print(
            "calculate backup conformers: {:{digits}} {}".format(
                "", "on", digits=digilen - len("calculate backup conformers")
            )
        )
    else:
        print(
            "calculate backup conformers: {:{digits}} {}".format(
                "", "off", digits=digilen - len("calculate backup conformers")
            )
        )
    print(
        "functional for part 1 and 2: {:{digits}} {}".format(
            "", args.func, digits=digilen - len("functional for part 1 and 2")
        )
    )
    print(
        "functional for part 3: {:{digits}} {}".format(
            "", args.func3, digits=digilen - len("functional for part 3")
        )
    )
    print(
        "basis set for part 3: {:{digits}} {}".format(
            "", str(args.basis3), digits=digilen - len("basis set for part 3")
        )
    )
    print(
        "calculate couplings: {:{digits}} {}".format(
            "", str(args.calcJ), digits=digilen - len("calculate couplings")
        )
    )
    print(
        "functional for coupling calculation: {:{digits}} {}".format(
            "", args.funcJ, digits=digilen - len("functional for coupling calculation")
        )
    )
    print(
        "basis set for coupling calculation: {:{digits}} {}".format(
            "", args.basisJ, digits=digilen - len("basis set for coupling calculation")
        )
    )
    print(
        "calculate shieldings: {:{digits}} {}".format(
            "", str(args.calcS), digits=digilen - len("calculate shieldings")
        )
    )
    print(
        "functional for shielding calculation: {:{digits}} {}".format(
            "", args.funcS, digits=digilen - len("functional for shielding calculation")
        )
    )
    print(
        "basis set for shielding calculation: {:{digits}} {}".format(
            "", args.basisS, digits=digilen - len("basis set for shielding calculation")
        )
    )
    print(
        "threshold for part 1: {:{digits}} {} kcal/mol".format(
            "", args.thr1, digits=digilen - len("threshold for part 1")
        )
    )
    print(
        "threshold for part 2: {:{digits}} {} kcal/mol".format(
            "", args.thr2, digits=digilen - len("threshold for part 2")
        )
    )
    if args.sm:
        print(
            "solvent model for part 1 and part 2: {:{digits}} {}".format(
                "", args.sm, digits=digilen - len("solvent model for part 1 and part 2")
            )
        )
    else:
        print(
            "solvent model for part 1 and part 2: {:{digits}} {}".format(
                "",
                "default",
                digits=digilen - len("solvent model for part 1 and part 2"),
            )
        )
    if args.gsolv2:
        print(
            "solvent model for Gsolv contribution of part 2: {:{digits}} {}".format(
                "",
                args.gsolv2,
                digits=digilen - len("solvent model for Gsolv contribution of part 2"),
            )
        )
    else:
        print(
            "solvent model for Gsolv contribution of part 2: {:{digits}} {}".format(
                "",
                "default",
                digits=digilen - len("solvent model for Gsolv contribution of part 2"),
            )
        )
    if args.sm3:
        print(
            "solvent model for part 3: {:{digits}} {}".format(
                "", args.sm3, digits=digilen - len("solvent model for part 3")
            )
        )
    else:
        print(
            "solvent model for part 3: {:{digits}} {}".format(
                "", "default", digits=digilen - len("solvent model for part 3")
            )
        )
    if args.sm4:
        print(
            "solvent model for part 4: {:{digits}} {}".format(
                "", args.sm4, digits=digilen - len("solvent model for part 4")
            )
        )
    else:
        print(
            "solvent model for part 4: {:{digits}} {}".format(
                "", "default", digits=digilen - len("solvent model for part 4")
            )
        )
    if args.check:
        print(
            "cautious checking for error and failed calculations: {:{digits}} {}".format(
                "",
                "on",
                digits=digilen
                - len("cautious checking for error and failed calculations"),
            )
        )
    else:
        print(
            "cautious checking for error and failed calculations: {:{digits}} {}".format(
                "",
                "off",
                digits=digilen
                - len("cautious checking for error and failed calculations"),
            )
        )
    if args.crestcheck:
        print(
            "checking the DFT-ensemble using CREST: {:{digits}} {}".format(
                "", "on", digits=digilen - len("checking the DFT-ensemble using CREST")
            )
        )
    else:
        print(
            "checking the DFT-ensemble using CREST: {:{digits}} {}".format(
                "", "off", digits=digilen - len("checking the DFT-ensemble using CREST")
            )
        )
    print(
        "maxthreads: {:{digits}} {}".format(
            "", args.maxthreads, digits=digilen - len("maxthreads")
        )
    )
    print("omp: {:{digits}} {}".format("", args.omp, digits=digilen - len("omp")))
    if len(spectrumlist) > 0:
        print(
            "calculating spectra for: {:{digits}} {}".format(
                "",
                (", ").join(spectrumlist),
                digits=digilen - len("calculating spectra for"),
            )
        )
        print(
            "reference for 1H: {:{digits}} {}".format(
                "", str(args.href), digits=digilen - len("reference for 1H")
            )
        )
        print(
            "reference for 13C: {:{digits}} {}".format(
                "", str(args.cref), digits=digilen - len("reference for 13C")
            )
        )
        print(
            "reference for 19F: {:{digits}} {}".format(
                "", str(args.fref), digits=digilen - len("reference for 19F")
            )
        )
        print(
            "reference for 31P: {:{digits}} {}".format(
                "", str(args.pref), digits=digilen - len("reference for 31P")
            )
        )
        print(
            "reference for 29Si: {:{digits}} {}".format(
                "", str(args.siref), digits=digilen - len("reference for 29Si")
            )
        )
        print(
            "resonance frequency: {:{digits}} {}".format(
                "", str(args.mf), digits=digilen - len("resonance frequency")
            )
        )
    print("\nEND of parameters\n")

    maxthreads = args.maxthreads

    if not args.run and not args.checkinput:
        print("\nFURTHER USER STEPS:\n")
        print("    1) carefully check/adjust the file flags.dat\n")
        print(
            "  ( 2) optional but recommended: check your input\n"
            "       by executing the ENSO program with the flag "
            "-checkinput )\n"
        )
        print("    3) execute the ENSO program with the flag -run\n")
        print("\nGoing to exit.\n")
        sys.exit(0)

    jsonfile = "enso.json"
    json_dict, firstrun = read_json(args.nstruc, cwd, jsonfile, args)
    if args.checkinput:
        write_json("save", json_dict, jsonfile)
        print(
            "Input check is finished. The ENSO program can be executed with "
            "the flag -run.\n"
        )
        sys.exit(0)

    return (
        args,
        jsonfile,
        json_dict,
        firstrun,
        conformersxyz,
        nat,
        maxthreads,
        xtbpath,
        environsettings,
        dbpath,
        cosmothermversion,
        cosmorssetup,
        orcapath,
        orca_old,
        crestpath,
        mpshiftpath,
        escfpath,
        spectrumlist,
    )


def part1(
    args, jsonfile, json_dict, conformersxyz, nat, maxthreads, xtbpath, environsettings
):
    """ Run crude optimization"""
    save_errors = (
        []
    )  # list to store all relevant errors and print them bundled for user convenience
    print("-----------------------------------------------------------")
    print(" PART 1 - Crude optimization for initial evaluation")
    print("-----------------------------------------------------------\n")
    if args.part1 == "on":
        # print flags for part 1
        digilen = 60
        print(
            "program: {:{digits}} {}".format(
                "", args.prog, digits=digilen - len("program")
            )
        )
        print(
            "functional for part1 and 2: {:{digits}} {}".format(
                "", args.func, digits=digilen - len("functional for part1 and 2")
            )
        )
        print(
            "using ANCOPT algorithm for optimization: {:{digits}} {}".format(
                "",
                args.ancopt,
                digits=digilen - len("using ANCOPT algorithm for optimization"),
            )
        )
        if args.solv:
            print(
                "solvent model: {:{digits}} {}".format(
                    "", args.sm, digits=digilen - len("solvent model")
                )
            )
        print(
            "threshold: {:{digits}} {} kcal/mol".format(
                "", args.thr1, digits=digilen - len("threshold")
            )
        )
        print(
            "starting number of considered conformers: {:{digits}} {}\n".format(
                "",
                args.nstruc,
                digits=digilen - len("starting number of considered conformers"),
            )
        )

        if args.prog == "tm":
            job = tm_job
        elif args.prog == "orca":
            job = orca_job
        else:
            print("ERROR jobtype (ORCA,TM) is not set correctly! Going to exit!")
            sys.exit(1)
        # list of numbers of all considered conformers, needed for converting
        # conformersxyz in coord files
        all_confs_list = []
        new_confs_list = []  # list of numbers of considered conformers
        old_confs_list = []  # calculated before
        for i in range(1, args.nstruc + 1):
            conf = "".join("CONF" + str(i))
            if json_dict[conf]["crude_opt"] == "calculated":
                # this conformer was calculated successfully before
                all_confs_list.append(i)
                old_confs_list.append(i)
            elif json_dict[conf]["crude_opt"] == "failed":
                # this conformer was calculated before but the crude optimization failed
                print(
                    "The calculation of part 1 failed for {} in the previous run."
                    " The conformer is sorted out.".format(conf)
                )
            elif json_dict[conf]["crude_opt"] == "not_calculated":
                # this conformer was not calculated before, so it has be calculated now
                new_confs_list.append(i)
                all_confs_list.append(i)
            if json_dict[conf]["removed_by_user"]:
                print("CONF{} is removed as requested by the user!".format(i))
                try:
                    if i in new_confs_list:
                        new_confs_list.remove(i)
                    if i in old_confs_list:
                        old_confs_list.remove(i)
                    if i in all_confs_list:
                        all_confs_list.remove(i)
                except Exception as e:
                    print(
                        "ERROR: CONF{} could not be removed from user-input".format(i)
                    )
                    save_errors.append(
                        "ERROR: CONF{} could not be removed from user-input".format(i)
                    )
        # exit if there are no conformers at all
        if len(all_confs_list) == 0:
            print("ERROR: There are no conformers!")
            write_json("save_and_exit", json_dict, jsonfile)
        # which conformers are considered now and which not
        if len(old_confs_list) >= 1:
            print("The crude optimization was calculated before for:")
            print_block(["CONF" + str(x) for x in old_confs_list])
        if len(new_confs_list) >= 1:
            print("The crude optimization is calculated for:")
            print_block(["CONF" + str(x) for x in new_confs_list])
        else:
            print("No conformers are considered additionally.")

        # check if directories for conformers exist
        check_for_folder(
            ["CONF" + str(i) for i in old_confs_list], args.func, args.debug
        )

        if len(new_confs_list) > 0:
            directories = new_folders(new_confs_list, args.func, save_errors)
            # check if at least one conformer is left
            if len(new_confs_list) == 0 and len(old_confs_list) == 0:
                print("ERROR: No conformers left!")
                write_json("save_and_exit", json_dict, jsonfile)

            # get GFNx-xTB energies
            gfne = conformersxyz2coord(conformersxyz, nat, args.func, new_confs_list)

            # setup queues
            q = Queue()
            resultq = Queue()

            # preparations (necessary because of cefine)
            instructprep = {
                "jobtype": "prep",
                "chrg": args.chrg,
                "unpaired": args.unpaired,
                "func": args.func,
                "solv": args.solv,
                "sm": args.sm,
                "omp": args.omp,
                "progsettings": {
                    "tempprogpath": "",
                    "xtbpath": xtbpath,
                    "orca_old": orca_old,
                    "omp": args.omp,
                },
            }
            results = run_in_parallel(
                q, resultq, job, maxthreads, directories, instructprep, args.func
            )
            exit_log, fail_rate = check_tasks(results, args)
            for conf in list(results):
                if not conf.success:
                    print(
                        "\nERROR: A problem has occurred in the preparation of "
                        "the optimization of {}! The "
                        "conformer is removed.\n".format(conf.name)
                    )
                    json_dict[conf.name]["crude_opt"] = "failed"
                    json_dict[conf.name]["consider_for_part2"] = False
                    json_dict[conf.name]["backup_for_part2"] = False
                    save_errors.append(
                        "Conformer {} was removed, because preparation failed!".format(
                            conf.name
                        )
                    )
                    results.remove(conf)
            if exit_log:
                print(
                    "\nERROR: too many preparations failed ({:.2f} %)!".format(
                        fail_rate
                    )
                )
                write_json("save_and_exit", json_dict, jsonfile)
            print("Preparation completed.")

            # Crude optimization of part1:
            instructopt = {
                "jobtype": "opt",
                "chrg": args.chrg,
                "unpaired": args.unpaired,
                "func": args.func,
                "solv": args.solv,
                "sm": args.sm,
                "full": False,
                "progsettings": {"omp": args.omp, "tempprogpath": ""},
            }
            if args.ancopt == "on":
                instructopt["jobtype"] = "xtbopt"
                instructopt["progsettings"]["xtbpath"] = xtbpath
            else:
                instructopt["jobtype"] = "opt"
            if job == tm_job:
                instructopt["progsettings"]["tempprogpath"] = ""
            elif job == orca_job:
                instructopt["progsettings"]["tempprogpath"] = orcapath
                instructopt["progsettings"]["orca_old"] = orca_old

            results = run_in_parallel(
                q, resultq, job, maxthreads, results, instructopt, args.func
            )
            exit_log, fail_rate = check_tasks(results, args)
            # sort out conformers with failed optimizations
            for conf in list(results):
                if not conf.success:
                    print(
                        "\nERROR: A problem has occurred in the optimization "
                        "of {}! The conformer is removed.\n".format(conf.name)
                    )
                    json_dict[conf.name]["crude_opt"] = "failed"
                    json_dict[conf.name]["consider_for_part2"] = False
                    json_dict[conf.name]["backup_for_part2"] = False
                    save_errors.append(
                        "Conformer {} was removed, because the optimization "
                        "failed!".format(conf.name)
                    )
                    results.remove(conf)
            # update json_dict
            for conf in results:
                json_dict[conf.name]["crude_opt"] = "calculated"
                json_dict[conf.name]["energy_crude_opt"] = conf.energy

            if exit_log:
                print(
                    "\nERROR: too many optimizations failed ({:.2f} %)!".format(
                        fail_rate
                    )
                )
                write_json("save_and_exit", json_dict, jsonfile)
        # END of if args.part1 == 'on':
        else:
            # all conformers were already calculated
            results = []
        # adding conformers calculated before to results
        try:
            length = max([len(str(x)) for x in old_confs_list]) + 4
        except ValueError:
            length = 4 + int(args.nstruc)
        for item in old_confs_list:
            tmp = job()
            tmp.name = "".join("CONF" + str(item))
            tmp.success = True
            tmp.energy = json_dict["".join("CONF" + str(item))]["energy_crude_opt"]
            results.append(tmp)
            print(
                "Crude optimization energy of {:>{digits}}: {:.7f}".format(
                    tmp.name, float(tmp.energy), digits=length
                )
            )
        results.sort(key=lambda x: int(x.name[4:]))

        # check if at least one calculation was successful
        if len(results) == 0:
            print("ERROR: No conformer left!")
            write_json("save_and_exit", json_dict, jsonfile)
        # get GFN-xTB energies
        gfne = conformersxyz2coord(conformersxyz, nat, args.func, all_confs_list, True)
        # calculate relative energies for each conformer with GFN-xTB and DFT
        for i in results:
            for item in gfne:
                if item[0] == i.name:
                    i.xtb_energy = item[1]

        if len(results) == 1 or len(
            [i.xtb_energy for i in results if i.xtb_energy is not None]
        ) in [0, 1]:
            mingfn = 0.0
        else:
            mingfn = float(
                min([i.xtb_energy for i in results if i.xtb_energy is not None])
            )

        mindft = float(min([i.energy for i in results if i.energy is not None]))

        for i in results:
            try:
                i.rel_xtb_energy = float(
                    "{:.2f}".format((i.xtb_energy - mingfn) * 627.50947428)
                )
            except:
                print(
                    "Problem with GFN-xTB energy read in from {} for {},"
                    " relative energy is set to 0.00.".format(conformersxyz, i.name)
                )
                i.xtb_energy = 0.00
                i.rel_xtb_energy = 0.00
            try:
                i.rel_energy = float(
                    "{:.2f}".format((i.energy - mindft) * 627.50947428)
                )
            except:
                print(
                    "\nERROR: A problem has occurred in the optimization of {}!"
                    " The conformer is removed.\n".format(i.name)
                )
                json_dict[i.name]["crude_opt"] = "failed"
                json_dict[i.name]["energy_crude_opt"] = None
                save_errors.append(
                    "Conformer {} was removed, because the optimization "
                    "failed!".format(i.name)
                )
                results.remove(i)

        if not len([i.rel_energy for i in results if i.rel_energy is not None]) > 0:
            print("All DFT crude optimization energies are None! Going to exit!!!")
            write_json("save_and_exit", json_dict, jsonfile)
        else:
            if len([i.rel_energy for i in results if i.rel_energy is not None]) > 0:
                maxreldft = max(
                    [i.rel_energy for i in results if i.rel_energy is not None]
                )
            else:
                print(
                    "ERROR: The maximal relative DFT energy could not be "
                    "calculated! Going to exit!"
                )
                write_json("save_and_exit", json_dict, jsonfile)

        print("\n*********************************")
        print("* {:^29} *".format("energies of part 1"))
        print("*********************************")
        try:
            length = max([len(str(i.name)) for i in results])
        except ValueError:
            length = 4 + int(args.nstruc)
        print(
            "{:{digits}}     E(GFN-xTB)    Erel(GFN-xTB)   E({})     Erel({})".format(
                "CONF#", args.func, args.func, digits=length
            )
        )
        print(
            "{:{digits}}     [Eh]          [kcal/mol]      [Eh]           [kcal/mol]".format(
                "", digits=length
            )
        )
        lowestenergy = min([item.energy for item in results if item.energy is not None])
        for i in results:
            if i.energy != lowestenergy:
                print(
                    "{:{digits}}    {:>10.7f}     {:>5.2f}         {:>10.7f}     {:>5.2f}".format(
                        i.name,
                        i.xtb_energy,
                        i.rel_xtb_energy,
                        i.energy,
                        i.rel_energy,
                        digits=length,
                    )
                )
            else:
                print(
                    "{:{digits}}    {:>10.7f}     {:>5.2f}         {:>10.7f}     {:>5.2f}  <----- lowest".format(
                        i.name,
                        i.xtb_energy,
                        i.rel_xtb_energy,
                        i.energy,
                        i.rel_energy,
                        digits=length,
                    )
                )

        backuplist = []
        if maxreldft > args.thr1:
            print("\n*********************************")
            print("* conformers considered further *")
            print("*********************************")
            tmp_consider = []
            for i in list(results):
                json_dict[i.name]["crude_opt"] = "calculated"
                json_dict[i.name]["energy_crude_opt"] = i.energy
                if i.rel_energy <= args.thr1:
                    tmp_consider.append(i.name)
                    json_dict[i.name]["consider_for_part2"] = True
                    json_dict[i.name]["backup_for_part2"] = False
                elif args.thr1 < i.rel_energy <= (args.thr1 + 2):
                    backuplist.append(i)
                    json_dict[i.name]["consider_for_part2"] = False
                    json_dict[i.name]["backup_for_part2"] = True
                    results.remove(i)
                else:
                    json_dict[i.name]["consider_for_part2"] = False
                    json_dict[i.name]["backup_for_part2"] = False
                    results.remove(i)
            print_block(tmp_consider)
        else:
            print(
                "\nAll relative energies are below the initial threshold of {} "
                "kcal/mol.\nAll conformers are "
                "considered further.".format(str(args.thr1))
            )
            for i in list(results):
                json_dict[i.name]["crude_opt"] = "calculated"
                json_dict[i.name]["energy_crude_opt"] = i.energy
                json_dict[i.name]["consider_for_part2"] = True
                json_dict[i.name]["backup_for_part2"] = False

        if len(backuplist) >= 1:
            print("\n**********************")
            print("* backup conformers  *")
            print("**********************")
            print(
                "Conformers with a relative energy between {} and {} kcal/mol are\n"
                "sorted to the backup conformers.".format(
                    str(args.thr1), str(args.thr1 + 2)
                )
            )
            print("CONF#     E({})       Erel({})".format(args.func, args.func))
            print("          [Eh]              [kcal/mol]")
            for item in backuplist:
                print(
                    "{:8} {:5f}      {:5.2f}".format(
                        item.name, item.energy, item.rel_energy
                    )
                )
        else:
            print("There are no backup conformers.")

        write_json("save", json_dict, jsonfile)

        if len(save_errors) > 0:
            print("***---------------------------------------------------------***")
            print("Printing most relevant errors again, just for user convenience:")
            for error in list(save_errors):
                print(save_errors.pop())
            print("***---------------------------------------------------------***")
        print("\nEND of part1.\n")
    else:  # if part 1 is switched off
        print("PART1 has been skipped by user.")
        print("Reading from file {} if necessary.\n".format(jsonfile))
        try:
            results
        except NameError:
            results = []
    return results, json_dict


def part2(
    args, results, jsonfile, cwd, json_dict, maxthreads, xtbpath, environsettings
):
    """ """
    print("-----------------------------------------------------------")
    print(" PART 2 - optimizations and low level free energy")
    print("-----------------------------------------------------------\n")
    save_errors = []
    if args.part2 == "on":
        digilen = 60
        print(
            "program: {:{digits}} {}".format(
                "", args.prog, digits=digilen - len("program")
            )
        )
        print(
            "using ANCOPT algorithm for optimization: {:{digits}} {}".format(
                "",
                args.ancopt,
                digits=digilen - len("using ANCOPT algorithm for optimization"),
            )
        )
        print(
            "functional for part2: {:{digits}} {}".format(
                "", args.func, digits=digilen - len("functional for part2")
            )
        )
        if args.solv:
            print(
                "solvent model for optimizations: {:{digits}} {}".format(
                    "", args.sm, digits=digilen - len("solvent model for optimizations")
                )
            )
            print(
                "solvent model for Gsolv contribution: {:{digits}} {}".format(
                    "",
                    args.gsolv2,
                    digits=digilen - len("solvent model for Gsolv contribution"),
                )
            )
        print(
            "program for RRHO contribution: {:{digits}} {}".format(
                "", args.rrhoprog, digits=digilen - len("program for RRHO contribution")
            )
        )
        print(
            "current threshold: {:{digits}} {} kcal/mol".format(
                "", args.thr2, digits=digilen - len("current threshold")
            )
        )

        if args.prog == "tm":
            job = tm_job
        elif args.prog == "orca":
            job = orca_job

        # tmp_results = list of conformers calculated before,
        # results is the list with conformers that are calculated now
        tmp_results = []

        if args.part1 == "on":  # normal run
            # use all conformers passed on from part 1
            for item in list(results):
                if json_dict[item.name]["opt"] == "calculated":
                    # this conformer was already calculated successfully,
                    # add energy from json_dict and move to tmp_results
                    item.energy = json_dict[item.name]["energy_opt"]
                    tmp_results.append(results.pop(results.index(item)))
                elif json_dict[item.name]["opt"] == "failed":
                    # this conformer was already calculated, but the calculation failed,
                    # remove from results
                    print(
                        "The calculation of part 2 failed for {} in the previous"
                        " run. The conformer is sorted out.".format(item.name)
                    )
                    results.remove(item)
                elif json_dict[item.name]["opt"] == "not_calculated":
                    # this conformer has to be calculated now, stays in results
                    pass
        elif args.part1 == "off" and not args.backup:
            # for both, firstrun and second run; start this run with part 2 and no backup
            # take all conformers into account that have 'consider_for_part2' set to True
            # (default if not calculated)
            # create tmp_results or results item for each
            results = []
            for i in range(1, args.nstruc + 1):
                conf = "".join("CONF" + str(i))
                if json_dict[conf]["consider_for_part2"]:
                    if json_dict[conf]["opt"] == "calculated":
                        # this conformer was already calculated successfully,
                        # create tmp_results item
                        tmp = job()
                        tmp.name = conf
                        tmp.success = True
                        tmp.energy = json_dict[conf]["energy_opt"]
                        tmp_results.append(tmp)
                    elif json_dict[conf]["opt"] == "failed":
                        # this conformer was already calculated, but the calculation failed
                        print(
                            "The calculation of part 2 failed for {} in the previous run. "
                            "The conformer is sorted out.".format(conf)
                        )
                    elif json_dict[conf]["opt"] == "not_calculated":
                        # this conformer has to be calculated now, create results environment
                        tmp = job()
                        tmp.name = conf
                        results.append(tmp)
        elif args.part1 == "off" and args.backup:
            # not first run, start this run with part2 and calculate backup
            # calculate all conformers with 'backup_for_part2' = True
            results = []
            for item in json_dict:
                if "CONF" in item:
                    if json_dict[item]["backup_for_part2"]:
                        if json_dict[item]["opt"] == "not_calculated":
                            # this conformer has to be calculated now, create results item
                            tmp = job()
                            tmp.name = item
                            results.append(tmp)
                        elif json_dict[item]["opt"] == "failed":
                            # this conformer was calculated before and the calculation failed
                            print(
                                "Though {} is declared as backup conformer, the "
                                "calculation of part2 was performed before and "
                                "failed. The conformer is sorted out.".format(item)
                            )
                        else:
                            # this conformer was already calculated successfully, create tmp_results item
                            tmp = job()
                            tmp.name = item
                            tmp.success = True
                            tmp.energy = json_dict[item]["energy_opt"]
                            tmp_results.append(tmp)
                            print(
                                "Though {} is declared as backup conformer, the "
                                "calculation of part2 was performed before. "
                                "Using the energy from the previous "
                                "run.".format(item)
                            )
                    elif (
                        json_dict[item]["consider_for_part2"]
                        and json_dict[item]["opt"] == "not_calculated"
                    ):
                        # this conformer has to be calculated now, create results item
                        tmp = job()
                        tmp.name = item
                        results.append(tmp)
                    elif (
                        json_dict[item]["consider_for_part2"]
                        and json_dict[item]["opt"] == "calculated"
                    ):
                        # conformer was already calculated successfully in
                        # run without backup, create tmp_results item
                        tmp = job()
                        tmp.name = item
                        tmp.success = True
                        tmp.energy = json_dict[item]["energy_opt"]
                        tmp_results.append(tmp)

        for item in json_dict:
            if "CONF" in item:
                if json_dict[item]["removed_by_user"]:
                    for conf in tmp_results:
                        if item == conf.name:
                            tmp_results.remove(conf)
                            print(
                                "{} is removed as requested by the user!".format(item)
                            )
                    for conf in results:
                        if item == conf.name:
                            results.remove(conf)
                            print(
                                "{} is removed as requested by the user!".format(item)
                            )

        # check if there is at least one conformer for part 2
        if len(tmp_results) == 0 and len(results) == 0:
            print("ERROR: There are no conformers for part 2!")
            write_json("save_and_exit", json_dict, jsonfile)
        print("\nOptimization:")
        print(
            "number of further considered conformers: {:{digits}} {}".format(
                "",
                str(len(results) + len(tmp_results)),
                digits=digilen - len("number of further considered conformers"),
            )
        )
        if len(tmp_results) > 0:
            print(
                "number of conformers which were optimized before: {:{digits}} {}".format(
                    "",
                    str(len(tmp_results)),
                    digits=digilen
                    - len("number of conformers which were optimized before"),
                )
            )
            print_block([i.name for i in tmp_results])
            if len(results) > 0:
                print(
                    "number of conformers that are optimized now: {:{digits}} {}".format(
                        "",
                        str(len(results)),
                        digits=digilen
                        - len("number of conformers that are optimized now"),
                    )
                )
                print_block([i.name for i in results])
            else:
                print("No conformers are considered additionally.")
        else:
            print("Considered conformers:")
            print_block([i.name for i in results])

        check_for_folder([i.name for i in tmp_results], args.func, args.debug)

        # setup queues
        q = Queue()
        resultq = Queue()

        if len(results) > 0:
            # create new folders
            print("Setting up new directories.")
            for conf in list(results):
                tmp_dir = os.path.join(conf.name, args.func)
                if not os.path.isdir(tmp_dir):
                    try:
                        mkdir_p(tmp_dir)
                    except:
                        print(
                            "ERROR: Could not create folder for {}!".format(conf.name)
                        )
                        print("ERROR: Removing {}!".format(conf.name))
                        save_errors.append(
                            "Conformer {} was removed, because IO failed!".format(
                                conf.name
                            )
                        )
                        results.remove(conf)
            # create coord from crest_conformers.xyz if there is no coord file in the directory CONFx/args.func
            new_confs_list = [conf.name for conf in results]
            conformersxyz2coord(conformersxyz, nat, args.func, new_confs_list)
            print("Constructed all folders.")

            if len(results) == 0:
                print("ERROR: No conformers left!")
                write_json("save_and_exit", json_dict, jsonfile)

            if args.part1 == "off":
                instructprep = {
                    "jobtype": "prep",
                    "chrg": args.chrg,
                    "unpaired": args.unpaired,
                    "func": args.func,
                    "solv": args.solv,
                    "sm": args.sm,
                    "progsettings": {"omp": args.omp, "tempprogpath": ""},
                }

                results = run_in_parallel(
                    q, resultq, job, maxthreads, results, instructprep, args.func
                )
                exit_log, fail_rate = check_tasks(results, args)
                for i in list(results):
                    if not i.success:
                        print(
                            "\nERROR: A problem has occurred in the optimization "
                            "preparation of {}! The conformer "
                            "is removed.\n".format(i.name)
                        )
                        json_dict[i.name]["opt"] = "failed"
                        json_dict[i.name]["consider_for_part3"] = False
                        json_dict[i.name]["backup_for_part3"] = False
                        json_dict[i.name]["energy_opt"] = None
                        save_errors.append(
                            "Conformer {} was removed, because preparation failed!".format(
                                i.name
                            )
                        )
                        results.remove(i)

                if exit_log:
                    print(
                        "\nERROR: too many preparations failed ({:.2f} %)!".format(
                            fail_rate
                        )
                    )
                    write_json("save_and_exit", json_dict, jsonfile)
                print("Preparation completed.")

            instructopt = {
                "jobtype": "opt",
                "chrg": args.chrg,
                "unpaired": args.unpaired,
                "func": args.func,
                "solv": args.solv,
                "sm": args.sm,
                "full": True,
                "progsettings": {"omp": args.omp, "tempprogpath": ""},
            }
            if args.ancopt == "on":
                instructopt["jobtype"] = "xtbopt"
                instructopt["progsettings"]["xtbpath"] = xtbpath
            else:
                instructopt["jobtype"] = "opt"
            if job == tm_job:
                instructopt["progsettings"]["tempprogpath"] = ""
            elif job == orca_job:
                instructopt["progsettings"]["tempprogpath"] = orcapath
                instructopt["progsettings"]["orca_old"] = orca_old

            results = run_in_parallel(
                q, resultq, job, maxthreads, results, instructopt, args.func
            )
            exit_log, fail_rate = check_tasks(results, args)
            if args.prog == "tm":  # copy control --> control_opt
                for item in results:
                    try:
                        shutil.copy(
                            os.path.join(item.workdir, "control"),
                            os.path.join(item.workdir, "control_opt"),
                        )
                    except:
                        pass

            # sort out conformers with crashed optimizations
            for i in list(results):
                if not i.success:
                    print(
                        "\nERROR: A problem has occurred in the optimization of "
                        "{}! The conformer is removed.\n".format(i.name)
                    )
                    json_dict[i.name]["opt"] = "failed"
                    json_dict[i.name]["consider_for_part3"] = False
                    json_dict[i.name]["backup_for_part3"] = False
                    json_dict[i.name]["energy_opt"] = None
                    save_errors.append(
                        "Conformer {} was removed, because optimization failed!".format(
                            i.name
                        )
                    )
                    results.remove(i)

            # write energies to json_dict
            for i in results:
                json_dict[i.name]["opt"] = "calculated"
                json_dict[i.name]["energy_opt"] = i.energy

            if exit_log:
                print(
                    "\nERROR: too many optimizations failed ({:.2f} %)!".format(
                        fail_rate
                    )
                )
                write_json("save_and_exit", json_dict, jsonfile)
            # always save json after optimization!
            write_json("save", json_dict, jsonfile)
            # end of optimization for new conformers

        # adding conformers calculated before to results
        try:
            length = max([len(str(item.name)) for item in tmp_results])
        except ValueError:
            length = 4 + int(args.nstruc)
        for item in tmp_results:
            print(
                "Optimization energy of {:>{digits}}: {:.7f}".format(
                    item.name, item.energy, digits=length
                )
            )
            results.append(item)
        results.sort(key=lambda x: int(x.name[4:]))

        if len(results) == 0:
            print("ERROR: No conformers left!")
            write_json("save_and_exit", json_dict, jsonfile)

        if args.gsolv2 in ("cosmors", "gbsa_gsolv") and args.solv:
            # optimized in implicit solvent therefore gas phase single-point needed for free energy with COSMO-RS
            print("\nGSOLV")
            print(
                "Gas Phase Single-point at {} level,\n(to be used in combination"
                " with the solvation correction).".format(args.func)
            )
            tmp_results = []
            for i in list(results):
                if json_dict[i.name]["sp_part2"] == "calculated":
                    # this conformer was already calculated successfully, move it tmp_results
                    i.energy = json_dict[i.name]["energy_sp_part2"]
                    tmp_results.append(results.pop(results.index(i)))
                elif json_dict[i.name]["sp_part2"] == "failed":
                    # this conformer was already calculated but the calculation failed, remove this conformer
                    print(
                        "The calculation of part 1 failed for {} in the previous"
                        " run. The conformer is sorted "
                        "out.".format(i.name)
                    )
                    results.remove(i)
                elif json_dict[i.name]["sp_part2"] == "not_calculated":
                    # this conformer has to be calculated now
                    pass
            print(
                "number of further considered conformers: {:{digits}} {}".format(
                    "",
                    str(len(results) + len(tmp_results)),
                    digits=digilen - len("number of further considered conformers"),
                )
            )
            if len(tmp_results) > 0:
                print(
                    "The gas phase single-point was already calculated for {} "
                    "conformers:".format(str(len(tmp_results)))
                )
                print_block([i.name for i in tmp_results])
                if len(results) > 0:
                    print(
                        "The gas phase single-point is calculated now for {} conformers:".format(
                            str(len(results))
                        )
                    )
                    print_block([i.name for i in results])
                else:
                    print("No conformers are considered additionally.")
            else:
                print("Considered conformers:")
                print_block([i.name for i in results])
            # calculate SP in gas phase for all conformers that remain in list results
            if len(results) > 0:
                if args.prog == "tm":
                    for item in results:
                        try:
                            shutil.copy(
                                os.path.join(item.workdir, "mos"),
                                os.path.join(item.workdir, "mos-keep"),
                            )
                        except FileNotFoundError:
                            pass

                instructprep = {
                    "jobtype": "prep",
                    "chrg": args.chrg,
                    "unpaired": args.unpaired,
                    "func": args.func,
                    "solv": args.solv,
                    "sm": "gas",
                    "progsettings": {"omp": args.omp, "tempprogpath": ""},
                }
                results = run_in_parallel(
                    q, resultq, job, maxthreads, results, instructprep, args.func
                )
                exit_log, fail_rate = check_tasks(results, args)

                for i in list(results):
                    if not i.success:
                        print(
                            "\nERROR: A problem has occurred in the gas phase "
                            "single-point of {}! The conformer is "
                            "removed.\n".format(i.name)
                        )
                        json_dict[i.name]["opt"] = "failed"
                        json_dict[i.name]["consider_for_part3"] = False
                        json_dict[i.name]["backup_for_part3"] = False
                        json_dict[i.name]["energy_opt"] = None
                        save_errors.append(
                            "Conformer {} was removed, because preparation of gas"
                            " phase single-point failed (because of COSMO-RS"
                            ")!".format(i.name)
                        )
                        results.remove(i)

                if exit_log:
                    print(
                        "\nERROR: too many preparations failed ({:.2f} %)!".format(
                            fail_rate
                        )
                    )
                    write_json("save_and_exit", json_dict, jsonfile)
                print("Preparation completed.")

                if args.prog == "tm":
                    for item in results:
                        try:
                            shutil.copy(
                                os.path.join(item.workdir, "mos-keep"),
                                os.path.join(item.workdir, "mos"),
                            )
                        except FileNotFoundError:
                            pass

                # place single-points in queue:
                instructsp = {
                    "jobtype": "sp",
                    "chrg": args.chrg,
                    "unpaired": args.unpaired,
                    "func": args.func,
                    "solv": args.solv,
                    "sm": "gas",
                    "progsettings": {"omp": args.omp, "tempprogpath": ""},
                }
                if job == tm_job:
                    instructsp["progsettings"]["tempprogpath"] = ""
                elif job == orca_job:
                    instructsp["progsettings"]["tempprogpath"] = orcapath
                    instructsp["progsettings"]["orca_old"] = orca_old

                results = run_in_parallel(
                    q, resultq, job, maxthreads, results, instructsp, args.func
                )
                exit_log, fail_rate = check_tasks(results, args)

                if args.prog == "tm":
                    for item in results:
                        try:
                            shutil.copy(
                                os.path.join(item.workdir, "control"),
                                os.path.join(item.workdir, "control_sp_gas"),
                            )
                        except:
                            pass

                # sort out conformers with crashed single-points
                for i in list(results):
                    if not i.success:
                        print(
                            "\nERROR: A problem has occurred in the single-point "
                            "calculation of {}! The conformer is "
                            "removed.\n".format(i.name)
                        )
                        json_dict[i.name]["sp_part2"] = "failed"
                        json_dict[i.name]["consider_for_part3"] = False
                        json_dict[i.name]["backup_for_part3"] = False
                        json_dict[i.name]["energy_sp_part2"] = None
                        save_errors.append(
                            "Error in single-point calculation for conformer {}".format(
                                i.name
                            )
                        )
                        results.remove(i)
                for i in results:
                    json_dict[i.name]["sp_part2"] = "calculated"
                    json_dict[i.name]["energy_sp_part2"] = i.energy

                if exit_log:
                    print(
                        "\nERROR: too many single-point calculations failed ({:.2f} %)!".format(
                            fail_rate
                        )
                    )
                    write_json("save_and_exit", json_dict, jsonfile)

            else:  # all conformers were calculated before
                results = []
            # adding conformers calculated before to results
            try:
                length = max([len(str(i.name)) for i in tmp_results])
            except ValueError:
                length = 4 + int(args.nstruc)
            for item in tmp_results:
                print(
                    "single-point energy of {:{digits}}: {:.7f}".format(
                        item.name, item.energy, digits=length
                    )
                )
                results.append(item)
            results.sort(key=lambda x: int(x.name[4:]))

            # check if at least one conformer is left:
            if len(results) == 0:
                print("ERROR: No conformers left!")
                write_json("save_and_exit", json_dict, jsonfile)

            # calculate only GBSA-Gsolv or COSMO-RS
            results = additive_gsolv(
                args,
                q,
                resultq,
                job,
                True,
                False,
                results,
                json_dict,
                jsonfile,
                save_errors,
                digilen,
            )
            #
        elif args.gsolv2 not in ("cosmors", "gbsa_gsolv") and args.solv:
            print("Gsolv values are included in the optimizations.")
        elif not args.solv:
            print("Since calculation is in the gas phase, Gsolv is not required.")
        if not args.rrhoprog == "off":
            # run RRHO
            results, save_errors = rrho_part23(
                args,
                q,
                resultq,
                job,
                True,
                False,
                results,
                json_dict,
                jsonfile,
                save_errors,
                digilen,
                cwd,
            )
        else:
            # skipp RRHO calculation:
            for conf in results:
                conf.rrho = 0.0
            # enso.json will not be updated with rrho information!
        for i in results:  # calculate free energy
            if (
                i.energy is not None
                and i.rrho is not None
                and isinstance(i.gsolv, float)
            ):
                i.free_energy = i.energy + i.gsolv + i.rrho
            else:
                # something went wrong and we do not know what, so reset everything
                print(
                    "The conformer {} was removed because the free energy could "
                    "not be calculated! {}: energy: {}, rrho: {}, gsolv: {}.".format(
                        i.name, i.name, i.energy, i.rrho, i.gsolv
                    )
                )
                json_dict[i.name]["consider_for_part2"] = False
                json_dict[i.name]["energy_opt"] = None
                json_dict[i.name]["sp_part2"] = "failed"
                json_dict[i.name]["energy_sp_part2"] = None
                if args.gsolv2 == "cosmors" and args.solv:
                    json_dict[i.name]["cosmo-rs"] = "failed"
                    json_dict[i.name]["energy_cosmo-rs"] = None
                elif args.gsolv2 == "gbsa_gsolv" and args.solv:
                    json_dict[i.name]["gbsa_gsolv"] = "failed"
                    json_dict[i.name]["energy_gbsa_gsolv"] = None
                json_dict[i.name]["rrho"] = "failed"
                json_dict[i.name]["energy_rrho"] = None
                json_dict[i.name]["consider_for_part3"] = False
                json_dict[i.name]["backup_for_part3"] = False
                save_errors.append(
                    "Error in summing up free energy for conformer {}".format(i.name)
                )
                results.remove(i)

        # evaluate calculated energies
        try:
            min2 = min(
                [item.free_energy for item in results if item.free_energy is not None]
            )
        except:
            print("\nERROR: no value to calculate difference from PART2!")
            write_json("save_and_exit", json_dict, jsonfile)

        for item in results:
            if item.success:
                try:
                    item.rel_free_energy = float(
                        "{:.2f}".format((item.free_energy - min2) * 627.50947428)
                    )
                except ValueError:
                    print(
                        "Could not calculate relative free energy for conformer "
                        "{}. rel free energy is set to a "
                        "very high number!".format(item.name)
                    )
                    item.rel_free_energy = 1000.00

        # check if at least one calculation was successful
        if len(results) == 0:
            print("ERROR: No conformers left!")
            write_json("save_and_exit", json_dict, jsonfile)

        print("\n*********************************")
        print("* {:^29} *".format("Gibbs free energies of part 2"))
        print("*********************************")
        try:
            length = max([len(str(i.name)) for i in results])
        except ValueError:
            length = 4 + int(args.nstruc)
        print(
            "{:{digits}}  {:12}    Gsolv [Eh]   RRHO [Eh]    Gtot [Eh]    rel. Gtot [kcal/mol]]\n".format(
                "CONF#", "E [Eh]", digits=length
            )
        )
        lowestfree = min(
            [item.free_energy for item in results if item.free_energy is not None]
        )
        for i in results:
            if i.free_energy != lowestfree:
                print(
                    "{:{digits}}  {:>10.7f},  {:>10.7f},  {:>10.7f},  {:>10.7f} {:>10.2f}".format(
                        i.name,
                        i.energy,
                        i.gsolv,
                        i.rrho,
                        i.free_energy,
                        i.rel_free_energy,
                        digits=length,
                    )
                )
            else:
                print(
                    "{:{digits}}  {:>10.7f},  {:>10.7f},  {:>10.7f},  {:>10.7f} {:>10.2f} <----- lowest".format(
                        i.name,
                        i.energy,
                        i.gsolv,
                        i.rrho,
                        i.free_energy,
                        i.rel_free_energy,
                        digits=length,
                    )
                )

        ### check if conformers are rotamers or identical, possibly sort them out
        rotdict = crest_routine(results, args.func, args.crestcheck, json_dict)
        results.sort(key=lambda x: int(x.name[4:]))
        # evaluate which conformers to consider further
        # energy from part2 should be smaller than threshold
        backuplist2 = []
        print("\n*********************************")
        print("* conformers considered further *")
        print("*********************************")
        tmp_consider = []
        for item in list(results):
            if item.rel_free_energy is not None and item.rel_free_energy <= args.thr2:
                tmp_consider.append(item.name)
                json_dict[item.name]["consider_for_part3"] = True
                json_dict[item.name]["backup_for_part3"] = False
            elif (
                item.rel_free_energy is not None
                and args.thr2 < item.rel_free_energy <= (args.thr2 + 2)
            ):
                json_dict[item.name]["consider_for_part3"] = False
                json_dict[item.name]["backup_for_part3"] = True
                backuplist2.append(item)
                results.remove(item)
            else:
                json_dict[item.name]["consider_for_part3"] = False
                json_dict[item.name]["backup_for_part3"] = False
                results.remove(item)
        print_block(tmp_consider)

        if len(backuplist2) >= 1:
            print("\n**********************")
            print("* backup conformers  *")
            print("**********************")
            print(
                "Conformers with a relative free energy between {} and {} kcal/mol\n"
                "are sorted to the backup conformers.".format(
                    str(args.thr2), str(args.thr2 + 2)
                )
            )
            print("CONF#     G({})       Grel({})".format(args.func, args.func))
            print("          [Eh]              [kcal/mol]")
            for item in backuplist2:
                print(
                    "{:8} {:5f}      {:5.2f}".format(
                        str(item.name),
                        float(item.free_energy),
                        float(item.rel_free_energy),
                    )
                )

        write_json("save", json_dict, jsonfile)

        if len(save_errors) > 0:
            print("***---------------------------------------------------------***")
            print("Printing most relevant errors again, just for user convenience:")
            for error in list(save_errors):
                print(save_errors.pop())
            print("***---------------------------------------------------------***")

        print("\n END of part2.\n")
    else:  # if part 2 is switched off
        print("PART2 has been skipped by user!")
        print("Reading from file {} if necessary.\n".format(jsonfile))
        try:
            results
        except NameError:
            results = []
        try:
            rotdict
        except NameError:
            rotdict = []
    return results, json_dict, rotdict


def part3(
    args, results, jsonfile, cwd, json_dict, maxthreads, xtbpath, environsettings
):
    """ """
    save_errors = []
    print("-----------------------------------------------------------")
    print(" PART 3 - high level free energy calculation")
    print("-----------------------------------------------------------\n")
    if args.part3 == "on":
        digilen = 60
        print(
            "program for single-point: {:{digits}} {}".format(
                "", str(args.prog), digits=digilen - len("program for single-point")
            )
        )
        print(
            "functional for part 3: {:{digits}} {}".format(
                "", args.func3, digits=digilen - len("functional for part 3")
            )
        )
        if args.func3 == "dsd-blyp" and args.basis3 != "def2-TZVPP":
            print("WARNING: only def2-TZVPP is available with DSD-BLYP!")
            args.basis3 = "def2-TZVPP"
        print(
            "basis set for part 3: {:{digits}} {}".format(
                "", args.basis3, digits=digilen - len("basis set for part 3")
            )
        )
        print(
            "program for RRHO contribution: {:{digits}} {}".format(
                "", args.rrhoprog, digits=digilen - len("program for RRHO contribution")
            )
        )
        if args.rrhoprog == "xtb":
            # os.environ["OMP_NUM_THREADS"] = "{:d}".format(args.omp)
            environsettings["OMP_NUM_THREADS"] = "{:d}".format(args.omp)
            print(
                "number of cores used by xtb: {:{digits}} {}".format(
                    "", args.omp, digits=digilen - len("number of cores used by xtb")
                )
            )
        if args.solv:
            print(
                "solvent model: {:{digits}} {}".format(
                    "", args.sm3, digits=digilen - len("solvent model")
                )
            )

        if not args.boltzmann:
            print("\nSingle-points")
        if args.prog3 == "tm":
            job = tm_job
        elif args.prog3 == "orca":
            job = orca_job

        if args.part2 == "on":
            # normal run or backup run
            # use all conformers passed on from part 1
            tmp_results = []
            # list of conformers calculated before, results is the list with conformers that are calculated now
            for item in list(results):
                if json_dict[item.name]["sp_part3"] == "calculated":
                    # this conformer was already calculated successfully, add energy from json_dict
                    # and move to tmp_results
                    item.energy = json_dict[item.name]["energy_sp_part3"]
                    tmp_results.append(results.pop(results.index(item)))
                elif json_dict[item.name]["sp_part3"] == "failed":
                    # this conformer was already calculated, but the calculation failed, remove from results
                    print(
                        "The calculation of part 3 failed for {} in the previous run. The conformer is sorted "
                        "out.".format(item.name)
                    )
                    results.remove(item)
                elif json_dict[item.name]["sp_part3"] == "not_calculated":
                    # this conformer has to be calculated now, stays in results
                    pass
            if args.backup:
                # add backup conformers
                for item in json_dict:
                    if "CONF" in item and json_dict[item]["backup_for_part3"]:
                        if json_dict[item]["sp_part3"] == "not_calculated":
                            # this conformer has to be calculated now, create results item
                            tmp = job()
                            tmp.name = item
                            results.append(tmp)
                        elif json_dict[item]["backup_for_part3"] == "failed":
                            # this conformer was calculated before and the calculation failed
                            print(
                                "Though {} is declared as backup conformer, the calculation of part3 was performed "
                                "before and failed. The conformer is sorted out.".format(
                                    item
                                )
                            )
                        elif json_dict[item]["sp_part3"] == "calculated":
                            # this conformer was already calculated successfully, create tmp_results item
                            tmp = job()
                            tmp.name = item
                            tmp.success = True
                            tmp.energy = json_dict[item]["energy_sp_part3"]
                            tmp_results.append(tmp)
                            print(
                                "Though {} is declared as backup conformer, the calculation of part3 was performed "
                                "before. Using the energy from the previous run.".format(
                                    item
                                )
                            )
        elif args.part2 == "off" and not args.boltzmann:
            # part 2 is off
            if firstrun:
                # not possible
                print(
                    "ERROR: Part 2 is switched off but the file {} containing all information of the previous run "
                    "cannot be found!".format(jsonfile)
                )
                write_json("save_and_exit", json_dict, jsonfile)
            # normal restart
            # take all conformers into account that have 'consider_for_part3' set to True (default is False)
            # create tmp_results or results item for each
            results = []
            tmp_results = []
            # list of conformers calculated before, results is the list with conformers that are calculated now
            for item in json_dict:
                if "CONF" in item and json_dict[item]["consider_for_part3"]:
                    if json_dict[item]["sp_part3"] == "calculated":
                        # this conformer was already calculated successfully, create tmp_results item
                        tmp = job()
                        tmp.name = item
                        tmp.success = True
                        tmp.energy = json_dict[item]["energy_sp_part3"]
                        tmp_results.append(tmp)
                    elif json_dict[item]["sp_part3"] == "failed":
                        # this conformer was already calculated, but the calculation failed
                        print(
                            "The calculation of part 3 failed for {} in the previous run. The conformer is sorted "
                            "out.".format(item)
                        )
                    elif json_dict[item]["sp_part3"] == "not_calculated":
                        # this conformer has to be calculated now, create results environment
                        tmp = job()
                        tmp.name = item
                        results.append(tmp)
            if args.backup:
                # add backup conformers
                for item in json_dict:
                    if "CONF" in item and json_dict[item]["backup_for_part3"]:
                        if json_dict[item]["sp_part3"] == "not_calculated":
                            # this conformer has to be calculated now, create results item
                            tmp = job()
                            tmp.name = item
                            results.append(tmp)
                        elif json_dict[item]["backup_for_part3"] == "failed":
                            # this conformer was calculated before and the calculation failed
                            print(
                                "Though {} is declared as backup conformer, the calculation of part3 was performed "
                                "before and failed. The conformer is sorted out.".format(
                                    item
                                )
                            )
                        elif json_dict[item]["sp_part3"] == "calculated":
                            # this conformer was already calculated successfully, create tmp_results item
                            tmp = job()
                            tmp.name = item
                            tmp.success = True
                            tmp.energy = json_dict[item]["energy_sp_part3"]
                            tmp_results.append(tmp)
                            print(
                                "Though {} is declared as backup conformer, the calculation of part3 was performed "
                                "before. Using the energy from the previous run.".format(
                                    item
                                )
                            )
        if args.boltzmann:
            # take all conformers into account that have consider_for_part3 = True and all required energies
            results = []
            tmp_results = []
            if args.sm3 in ("cosmors", "gbsa_gsolv") and args.solv:
                if args.sm3 == "cosmors":
                    js_smodel = "cosmo-rs"
                    js_sm_energy = "energy_cosmo-rs"
                elif args.sm3 == "gbsa_gsolv":
                    js_smodel = "gbsa_gsolv"
                    js_sm_energy = "energy_gbsa_gsolv"
                for item in json_dict:
                    if "CONF" in item and json_dict[item]["consider_for_part3"]:
                        if (
                            json_dict[item]["sp_part3"] == "calculated"
                            and json_dict[item]["rrho"] == "calculated"
                            and json_dict[item][js_smodel] == "calculated"
                        ):
                            tmp = job()
                            tmp.name = item
                            tmp.success = True
                            tmp.energy = json_dict[item]["energy_sp_part3"]
                            tmp.rrho = json_dict[item]["energy_rrho"]
                            tmp.gsolv = json_dict[item][js_sm_energy]
                            results.append(tmp)
            else:
                for item in json_dict:
                    if "CONF" in item and json_dict[item]["consider_for_part3"]:
                        if (
                            json_dict[item]["sp_part3"] == "calculated"
                            and json_dict[item]["rrho"] == "calculated"
                        ):
                            tmp = job()
                            tmp.name = item
                            tmp.success = True
                            tmp.energy = json_dict[item]["energy_sp_part3"]
                            tmp.rrho = json_dict[item]["energy_rrho"]
                            tmp.gsolv = 0.0
                            results.append(tmp)
        for item in json_dict:
            if "CONF" in item:
                if json_dict[item]["removed_by_user"]:
                    for conf in tmp_results:
                        if item == conf.name:
                            tmp_results.remove(conf)
                            print(
                                "{} is removed as requested by the user!".format(item)
                            )
                    for conf in results:
                        if item == conf.name:
                            results.remove(conf)
                            print(
                                "{} is removed as requested by the user!".format(item)
                            )
        # check if there is at least one conformer for part 3
        if len(tmp_results) == 0 and len(results) == 0:
            print("ERROR: There are no conformers left for part 3!")
            write_json("save_and_exit", json_dict, jsonfile)
        elif (len(tmp_results) + len(results)) == 1:
            print(
                "There is only one conformer left. Part 3 is skipped and Part 4 is calculated directly."
            )
            for (
                i
            ) in (
                tmp_results
            ):  # the conformer could either be in results or in tmp_results
                results.append(i)
            for i in results:
                if not i.energy:
                    i.energy = 0.0
                if not i.gsolv:
                    i.gsolv = 0.0
                if not i.rrho:
                    i.rrho = 0.0
                i.bm_weight = 1.00
                i.new_bm_weight = 1.00
                json_dict[i.name]["consider_for_part4"] = True
        else:
            print(
                "number of further considered conformers: {:{digits}} {}".format(
                    "",
                    str(len(results) + len(tmp_results)),
                    digits=digilen - len("number of further considered conformers"),
                )
            )
            if len(tmp_results) > 0:
                print(
                    "The single-point was calculated for {} conformers before:".format(
                        str(len(tmp_results))
                    )
                )
                print_block([i.name for i in tmp_results])
                if len(results) > 0:
                    print(
                        "The single-point is calculated for {} conformers now:".format(
                            str(len(results))
                        )
                    )
                    print_block([i.name for i in results])
                else:
                    print("No conformers are considered additionally.")
            else:
                print("Considered conformers:")
                print_block([i.name for i in results])
            if not args.boltzmann:
                # check if directories of conformers calculated before exist,
                check_for_folder([i.name for i in tmp_results], args.func, args.debug)

                # setup queues
                q = Queue()
                resultq = Queue()

                if results:
                    # create new folder and copy coord from CONFx/args.func to CONFx/args.func3
                    print("Setting up new directories.")
                    if not args.func == args.func3:
                        # if args.func == args.func3 no new directory or copying of file coord is needed
                        for conf in list(results):
                            tmp_from = os.path.join(cwd, conf.name, args.func)
                            tmp_to = os.path.join(cwd, conf.name, args.func3)
                            if not os.path.isdir(tmp_to):
                                mkdir_p(tmp_to)
                            try:
                                shutil.copy(
                                    os.path.join(tmp_from, "coord"),
                                    os.path.join(tmp_to, "coord"),
                                )
                            except FileNotFoundError:
                                if not os.path.isfile(os.path.join(tmp_from, "coord")):
                                    print(
                                        "ERROR: while copying the coord file from {}! The corresponding "
                                        "file does not exist.".format(tmp_from)
                                    )
                                elif not os.path.isdir(tmp_to):
                                    print(
                                        "ERROR: Could not create folder {}!".format(
                                            tmp_to
                                        )
                                    )
                                print("ERROR: Removing conformer {}!".format(conf.name))
                                json_dict[conf.name]["sp_part3"] = "failed"
                                json_dict[conf.name]["consider_for_part4"] = False
                                json_dict[conf.name]["energy_sp_part3"] = None
                                save_errors.append(
                                    "Conformer {} was removed, because IO failed!".format(
                                        conf.name
                                    )
                                )
                                results.remove(conf)
                    else:  # args.func == args.func3
                        for conf in list(results):
                            tmp_from = os.path.join(cwd, conf.name, args.func)
                            if not os.path.isdir(tmp_from) or not os.path.isfile(
                                os.path.join(tmp_from, "coord")
                            ):
                                print(
                                    "ERROR: Either the folder {} or the coord file does not exist!".format(
                                        last_folders(tmp_from, 2)
                                    )
                                )
                                print("ERROR: Removing conformer {}!".format(conf.name))
                                json_dict[conf.name]["sp_part3"] = "failed"
                                json_dict[conf.name]["consider_for_part4"] = False
                                json_dict[conf.name]["energy_sp_part3"] = None
                                save_errors.append(
                                    "Conformer {} was removed, because IO failed!".format(
                                        conf.name
                                    )
                                )
                                results.remove(conf)
                    print("Constructed all new folders.")

                    if not results:
                        print("ERROR: No conformers left!")
                        write_json("save_and_exit", json_dict, jsonfile)

                    # place prep in queue:
                    instructprep = {
                        "jobtype": "prep",
                        "chrg": args.chrg,
                        "unpaired": args.unpaired,
                        "func": args.func3,
                        "basis": args.basis3,
                        "solv": args.solv,
                        "sm": "gas",
                        "boltzmann": args.boltzmann,
                        "progsettings": {"omp": args.omp, "tempprogpath": ""},
                    }
                    if args.sm3 == "smd" and args.solv:
                        instructprep["sm"] = "smd"  # solvation consider during SP
                    elif args.sm3 == "dcosmors" and args.solv:
                        instructprep["sm"] = "dcosmors"  # solvation consider during SP
                    else:
                        instructprep["sm"] = "gas"
                        # solvation calculated separately with COMSO-RS

                    results = run_in_parallel(
                        q, resultq, job, maxthreads, results, instructprep, args.func3
                    )
                    exit_log, fail_rate = check_tasks(results, args)

                    # sort out conformers with crashed preparations
                    for i in list(results):
                        if not i.success:
                            print(
                                "\nERROR: A problem has occurred in the single-point preperation of {}! The conformer"
                                " is removed.\n".format(i.name)
                            )
                            json_dict[i.name]["sp_part3"] = "failed"
                            json_dict[i.name]["consider_for_part4"] = False
                            json_dict[i.name]["energy_sp_part3"] = None
                            save_errors.append(
                                "Conformer {} was removed, because preparation failed!".format(
                                    i.name
                                )
                            )
                            results.remove(i)

                    if exit_log:
                        print(
                            "\nERROR: too many preparations failed ({:.2f} %)!".format(
                                fail_rate
                            )
                        )
                        write_json("save_and_exit", json_dict, jsonfile)
                    print("Preparation completed.")

                    # place SP work in queue:
                    instructsp = {
                        "jobtype": "sp",
                        "chrg": args.chrg,
                        "unpaired": args.unpaired,
                        "func": args.func3,
                        "solv": args.solv,
                        "basis": args.basis3,
                        "boltzmann": args.boltzmann,
                        "progsettings": {"omp": args.omp, "tempprogpath": ""},
                    }
                    if job == tm_job:
                        instructsp["progsettings"]["tempprogpath"] = ""
                    elif job == orca_job:
                        instructsp["progsettings"]["tempprogpath"] = orcapath
                        instructsp["progsettings"]["orca_old"] = orca_old

                    results = run_in_parallel(
                        q, resultq, job, maxthreads, results, instructsp, args.func3
                    )
                    exit_log, fail_rate = check_tasks(results, args)

                    # sort out conformers with crashed single-points
                    for i in list(results):
                        if not i.success:
                            print(
                                "\nERROR: A problem has occurred in the single-point calculation of {}! The conformer"
                                " is removed.\n".format(i.name)
                            )
                            json_dict[i.name]["sp_part3"] = "failed"
                            json_dict[i.name]["consider_for_part4"] = False
                            json_dict[i.name]["energy_sp_part3"] = None
                            save_errors.append(
                                "Conformer {} was removed, because the single-point calculation "
                                "failed!".format(i.name)
                            )
                            results.remove(i)
                    # write energies to json_dict
                    for i in results:
                        json_dict[i.name]["sp_part3"] = "calculated"
                        json_dict[i.name]["energy_sp_part3"] = i.energy
                        # write energy to sp3_energy
                        i.sp3_energy = i.energy

                    if exit_log:
                        print(
                            "\nERROR: too many single-points failed ({:.2f} %)!".format(
                                fail_rate
                            )
                        )
                        write_json("save_and_exit", json_dict, jsonfile)
                    # end of sp for new conformers

                # adding conformers calculated before to results
                try:
                    length = max([len(str(i.name)) for i in tmp_results])
                except ValueError:
                    length = 4 + int(args.nstruc)
                for item in tmp_results:
                    print(
                        "single-point energy of {:{digits}}: {:.7f}".format(
                            item.name, item.energy, digits=length
                        )
                    )
                    results.append(item)
                results.sort(key=lambda x: int(x.name[4:]))

                if len(results) == 0:
                    print("ERROR: No conformers left!")
                    write_json("save_and_exit", json_dict, jsonfile)

                print("\nGSOLV")
                if args.sm3 in ("cosmors", "gbsa_gsolv") and args.solv:
                    # calculate only GBSA-Gsolv or COSMO-RS
                    results = additive_gsolv(
                        args,
                        q,
                        resultq,
                        job,
                        True,
                        False,
                        results,
                        json_dict,
                        jsonfile,
                        save_errors,
                        digilen,
                    )
                else:
                    for conf in results:
                        # if COSMO-RS was used in part 2 but not in part 3, gsolv has to be reset here
                        conf.gsolv = 0.00
                    if args.sm3 == "smd" and args.solv:
                        print(
                            "Gsolv values are included in the single-point calculations with SMD.\n"
                        )
                    elif args.sm3 == "dcosmors" and args.solv:
                        print(
                            "Gsolv values are included in the single-point calculations with DCOSMO-RS.\n"
                        )
                    elif not args.solv:
                        print(
                            "Since the calculation is performed in gas phase, Gsolv is not required.\n"
                        )

                if not args.rrhoprog == "off":
                    # run RRHO
                    results, save_errors = rrho_part23(
                        args,
                        q,
                        resultq,
                        job,
                        False,
                        True,
                        results,
                        json_dict,
                        jsonfile,
                        save_errors,
                        digilen,
                        cwd,
                    )
                else:
                    # skipp RRHO calculation:
                    for conf in results:
                        conf.rrho = 0.0
                    # enso.json will not be updated with rrho information!
            for i in results:  # calculate free energy
                if (
                    i.energy is not None
                    and i.rrho is not None
                    and isinstance(i.gsolv, float)
                ):
                    i.free_energy = i.energy + i.gsolv + i.rrho
                else:
                    # something went wrong and we do not know what, so reset everything
                    print(
                        "The conformer {} was removed because the free energy could not be calculated!"
                        "{}: energy: {}, rrho: {}, gsolv: {}.".format(
                            i.name, i.name, i.energy, i.rrho, i.gsolv
                        )
                    )
                    json_dict[i.name]["sp_part3"] = "failed"
                    json_dict[i.name]["energy_sp_part3"] = None
                    if args.gsolv2 == "cosmors" and args.solv:
                        json_dict[i.name]["cosmo-rs"] = "failed"
                        json_dict[i.name]["energy_cosmo-rs"] = None
                    if args.gsolv2 == "gbsa_gsolv" and args.solv:
                        json_dict[i.name]["gbsa_gsolv"] = "failed"
                        json_dict[i.name]["energy_gbsa_gsolv"] = None
                    json_dict[i.name]["rrho"] = "failed"
                    json_dict[i.name]["energy_rrho"] = None
                    json_dict[i.name]["consider_for_part4"] = False
                    save_errors.append(
                        "Error in summing up free energy for conformer {}".format(
                            i.name
                        )
                    )
                    results.remove(i)

            # check for None in free_energy:
            counter = 0
            allitems = len(results)
            for conf in list(results):
                if conf.free_energy is None:
                    counter += 1
                    save_errors.append(
                        "Error in free energy for conformer {}".format(conf.name)
                    )
                    results.remove(conf)
                    print("removed {}, because free energy is None.".format(conf.name))
            if len(results) == 0:
                print("\nERROR: all free energies are None!")
                write_json("save_and_exit", json_dict, jsonfile)

            # calculate delta G
            if len([item for item in results if item.free_energy is not None]) == 1:
                minfree = [
                    item.free_energy for item in results if item.free_energy is not None
                ][0]
            else:
                try:
                    minfree = min(
                        [
                            item.free_energy
                            for item in results
                            if item.free_energy is not None
                        ]
                    )
                except:
                    print("No value to calculate difference for PART3!")
                    write_json("save_and_exit", json_dict, jsonfile)
            for item in results:
                item.rel_free_energy = (item.free_energy - minfree) * 627.50947428

            print("\n*********************************")
            print("* {:^29} *".format("Gibbs free energies of part 3"))
            print("*********************************")
            print(
                "CONF#          E [Eh]       Gsolv [Eh]   RRHO [Eh]    Gtot [Eh]    rel. Gtot [kcal/mol]]\n"
            )
            lowestfree = min(
                [item.free_energy for item in results if item.free_energy is not None]
            )
            for i in results:
                if i.free_energy != lowestfree:
                    print(
                        "{:10}  {:>10.7f},  {:>10.7f},  {:>10.7f},  {:>10.7f} {:>10.2f}".format(
                            i.name,
                            i.energy,
                            i.gsolv,
                            i.rrho,
                            i.free_energy,
                            i.rel_free_energy,
                        )
                    )
                else:
                    print(
                        "{:10}  {:>10.7f},  {:>10.7f},  {:>10.7f},  {:>10.7f} {:>10.2f} <----- lowest".format(
                            i.name,
                            i.energy,
                            i.gsolv,
                            i.rrho,
                            i.free_energy,
                            i.rel_free_energy,
                        )
                    )

            if float(counter) / allitems >= 0.25 and args.check:
                print(
                    "\nERROR: too many free energies are None ({:.2f} %)!".format(
                        float(counter) / float(allitems) * 100
                    )
                )
                write_json("save_and_exit", json_dict, jsonfile)

            print("Calculate Boltzmann populations")
            au2J = 4.3597482e-18  # a.u.(hartree/mol) to J
            # kcalmol2J = 7.993840458653155e-21  # kcal/mol to J  =   *4.184 / N_{A}
            kb = 1.3806485279e-23  # J/K
            try:
                T = float(args.temperature)
            except:
                print("WARNING:Temperature could not be converted using T=298.15 K")
                T = 298.15  # K
            ## avoid division by zero
            if T == 0:
                T += 0.00001
            args.temperature
            sum = 0.0
            for item in results:
                sum += math.exp(-((item.free_energy - minfree) * au2J) / (kb * T))
            for item in results:
                item.bm_weight = (
                    math.exp(-((item.free_energy - minfree) * au2J) / (kb * T)) / sum
                )

            # sort out all conformers  with a Boltzmann population smaller
            # than 1% and recalculate BW
            print(
                "Recalculate Boltzmann populations (boltzmann new) by neglecting "
                "all conformers with an\noriginal Boltzmann population "
                "(boltzmann old) below 1%"
            )
            bw_list = []
            for item in list(results):
                if item.bm_weight >= 0.01:
                    bw_list.append(item.name)
            sum = 0.0
            for item in results:
                if item.name in bw_list:
                    sum += math.exp(-((item.free_energy - minfree) * au2J) / (kb * T))
            for item in results:
                if item.name in bw_list:
                    item.new_bm_weight = (
                        math.exp(-((item.free_energy - minfree) * au2J) / (kb * T))
                        / sum
                    )
                else:
                    item.new_bm_weight = 0.0

            print("\n*********************************")
            print("* {:^29} *".format("boltzmann weight of part 3"))
            print("*********************************")
            print(
                "#CONF         Gtot [Eh]   rel Gtot[kcal/mol]   boltzmann old [%]   boltzmann new [%] \n"
            )
            for item in results:
                print(
                    "{:10} {:>12.7f}  {:>12.2f}  {:>18.2f}  {:>17.2f}".format(
                        str(item.name),
                        item.free_energy,
                        item.rel_free_energy,
                        100 * item.bm_weight,
                        100 * item.new_bm_weight,
                    )
                )

            # evaluate which conformers to consider further
            print("\n*********************************")
            print("* conformers considered further *")
            print("*********************************")
            thr3 = 0.98
            results.sort(key=lambda x: float(x.new_bm_weight), reverse=True)
            sum_bw = 0.0
            for conf in results:
                if conf.new_bm_weight == 0.0:
                    continue
                print(conf.name)
                json_dict[conf.name]["consider_for_part4"] = True
                sum_bw += conf.new_bm_weight
                if sum_bw >= thr3:
                    break
                else:
                    continue
            # remove conformers from list results if they are not considered further
            for conf in list(results):
                if not json_dict[conf.name]["consider_for_part4"]:
                    results.remove(conf)
            results.sort(key=lambda x: int(x.name[4:]))

            # check if at least one calculation was successful
            if not results:
                print("\nERROR: No conformers left!")
                write_json("save_and_exit", json_dict, jsonfile)

        # the following is executed even if there is only one conformer left
        #  write anmr files, required by anmr program
        # move old anmr_enso files to anmr_enso.x+1
        move_recursively(cwd, "anmr_enso")
        write_anmr_enso(cwd, results)

        # write populated confs from part3
        move_recursively(cwd, "populated-conf-part3.xyz")
        write_trj(
            sorted(results, key=lambda x: float(x.free_energy)),
            "populated-conf-part3.xyz",
            args.func,
        )

        write_json("save", json_dict, jsonfile)

        if save_errors:
            print("***---------------------------------------------------------***")
            print("Printing most relevant errors again, just for user convenience:")
            for error in list(save_errors):
                print(save_errors.pop())
            print("***---------------------------------------------------------***")

        print("\n END of part3.\n")
    else:  # if part 3 is switched off
        print("PART3 has been skipped by user.")
        print("Reading from file {} if necessary.\n".format(jsonfile))
        try:
            results
        except NameError:
            results = []
    return results, json_dict


def part4(
    args,
    results,
    jsonfile,
    cwd,
    json_dict,
    maxthreads,
    xtbpath,
    environsettings,
    rotdict,
    spectrumlist,
    escfpath,
    mpshiftpath,
):
    """ """
    save_errors = []
    removelist = []
    print("\n-----------------------------------------------------------")
    print(" PART 4 - couplings and shieldings")
    print("-----------------------------------------------------------\n")
    digilen = 60
    if args.part4 == "on":
        if args.unpaired > 0:
            print(
                "ERROR: Coupling and shift calculation is only available for "
                "closed-shell systems!"
            )
            write_json("save_and_exit", json_dict, jsonfile)
        print(
            "main program: {:{digits}} {}".format(
                "", args.prog, digits=digilen - len("main program")
            )
        )
        if args.prog4:
            print(
                "program for part 4: {:{digits}} {}".format(
                    "", args.prog4, digits=digilen - len("program for part 4")
                )
            )
            if args.prog4 == "tm":
                job = tm_job
            elif args.prog4 == "orca":
                job = orca_job
        else:
            args.prog4 = args.prog
            print(
                "program for part 4 was automatically set to: "
                "{:{digits}} {}".format(
                    "",
                    args.prog4,
                    digits=digilen - len("program for part 4 was automatically set to"),
                )
            )
            if args.prog4 == "tm":
                job = tm_job
            elif args.prog4 == "orca":
                job = orca_job
        print(
            "calculate couplings: {:{digits}} {}".format(
                "", args.calcJ, digits=digilen - len("calculate couplings")
            )
        )
        print(
            "functional for coupling calculation: {:{digits}} {}".format(
                "",
                args.funcJ,
                digits=digilen - len("functional for coupling calculation"),
            )
        )
        if args.basisJ == "default":
            if job == tm_job:
                args.basisJ = "def2-TZVP"
            elif job == orca_job:
                args.basisJ = "pcJ-0"
        print(
            "basis set for coupling calculation: {:{digits}} {}".format(
                "",
                args.basisJ,
                digits=digilen - len("basis set for coupling calculation"),
            )
        )
        print(
            "calculate shieldings: {:{digits}} {}".format(
                "", args.calcS, digits=digilen - len("calculate shieldings")
            )
        )
        print(
            "functional for shielding calculation: {:{digits}} {}".format(
                "",
                args.funcS,
                digits=digilen - len("functional for shielding calculation"),
            )
        )
        if args.basisS == "default":
            if job == tm_job:
                args.basisS = "def2-TZVP"
            elif job == orca_job:
                args.basisS = "pcSseg-2"
        print(
            "basis set for shielding calculation: {:{digits}} {}".format(
                "",
                args.basisS,
                digits=digilen - len("basis set for shielding calculation"),
            )
        )
        if args.prog4 == "tm" and args.solv:
            if args.sm4 != "cosmo":
                print(
                    "WARNING: With TM, only COSMO is available as solvent model for"
                    " the shielding and coupling constants."
                )
                args.sm4 = "cosmo"
        elif args.prog4 == "orca" and args.solv:
            if args.sm4 != "cpcm" and args.sm4 != "smd":
                print(
                    "WARNING: With ORCA, only CPCM or SMD are avaible as solvent"
                    " model for the shielding and coupling constants.\n"
                    " CPCM is set as default."
                )
                args.sm4 = "cpcm"
        if args.solv:
            print(
                "solvent model: {:{digits}} {}".format(
                    "", args.sm4, digits=digilen - len("solvent model")
                )
            )
        print(
            "calculating spectra for: {:{digits}} {}".format(
                "",
                ", ".join(spectrumlist),
                digits=digilen - len("calculating spectra for"),
            )
        )
        if "1H" in spectrumlist:
            print(
                "reference for 1H: {:{digits}} {}".format(
                    "", str(args.href), digits=digilen - len("reference for 1H")
                )
            )
        if "13C" in spectrumlist:
            print(
                "reference for 13C: {:{digits}} {}".format(
                    "", str(args.cref), digits=digilen - len("reference for 13C")
                )
            )
        if "19F" in spectrumlist:
            print(
                "reference for 19F: {:{digits}} {}".format(
                    "", str(args.fref), digits=digilen - len("reference for 19F")
                )
            )
        if "31P" in spectrumlist:
            print(
                "reference for 31P: {:{digits}} {}".format(
                    "", str(args.pref), digits=digilen - len("reference for 31P")
                )
            )
        if "29Si" in spectrumlist:
            print(
                "reference for 29Si: {:{digits}} {}".format(
                    "", str(args.siref), digits=digilen - len("reference for 29Si")
                )
            )

        print(
            "resonance frequency: {:{digits}} {}".format(
                "", str(args.mf), digits=digilen - len("resonance frequency")
            )
        )
        print("")

        # setup queues
        q = Queue()
        resultq = Queue()

        # COUPLINGS
        if args.calcJ == "on":
            tmp_results = []
            # list with conformers calculated before
            if args.part3 == "on":  # use all conformers of list results
                for item in list(results):
                    # it is possible to calculate multiple nuclei at once,
                    # so check for each nuc if it was
                    # calculated before or if it failed
                    tmp_failed = False
                    tmp_calculate = False
                    for nuc in spectrumlist:
                        nuc_string = "".join(nuc + "_J")
                        if json_dict[item.name][nuc_string] == "not_calculated":
                            tmp_calculate = True
                        elif json_dict[item.name][nuc_string] == "failed":
                            tmp_failed = True
                            print(
                                "ERROR: The {} coupling calculation for {} failed"
                                " in the previous run!\n"
                                "The conformer is sorted out.".format(nuc, item.name)
                            )
                            break
                    if tmp_failed:
                        results.remove(item)
                        removelist.append(item.name)
                    elif not tmp_calculate and not tmp_failed:
                        # for this conformer, all activated coupling calculation
                        # were performed in the previous run
                        tmp_results.append(results.pop(results.index(item)))
            if args.part3 == "off":
                # use all conformers with 'consider_for_part4' = True
                if firstrun:  # not possible
                    print(
                        "ERROR: Part 3 is switched off but the file {} containing"
                        " all information of the previous "
                        "run cannot be found!".format(jsonfile)
                    )
                    write_json("save_and_exit", json_dict, jsonfile)
                results = []  # list of conformers that have to be calculated now
                for item in json_dict:
                    if "CONF" in item and json_dict[item]["consider_for_part4"]:
                        # it is possible to calculate multiple nuclei at once,
                        # so check for each nuc if it was
                        # calculated before or if it failed
                        tmp_failed = False
                        tmp_calculate = False
                        for nuc in spectrumlist:
                            nuc_string = "".join(nuc + "_J")
                            if json_dict[item][nuc_string] == "not_calculated":
                                tmp_calculate = True
                            elif json_dict[item][nuc_string] == "failed":
                                tmp_failed = True
                                print(
                                    "ERROR: The {} coupling calculation for {} "
                                    "failed in the previous run! The "
                                    "conformer is sorted out.".format(nuc, item.name)
                                )
                                removelist.append(item.name)
                                break
                        if tmp_calculate and not tmp_failed:
                            # this conformer has be calculated now
                            tmp = job()
                            tmp.name = item
                            tmp.success = True
                            results.append(tmp)
                        elif not tmp_calculate and not tmp_failed:
                            # for this conformer, all activated coupling
                            # calculation were performed in the previous run
                            tmp = job()
                            tmp.name = item
                            tmp.success = True
                            tmp_results.append(tmp)

            for item in json_dict:
                if "CONF" in item:
                    if json_dict[item]["removed_by_user"]:
                        for conf in tmp_results:
                            if item == conf.name:
                                tmp_results.remove(conf)
                                removelist.append(item.name)
                                print(
                                    "{} is removed as requested by the user!".format(
                                        item
                                    )
                                )
                        for conf in results:
                            if item == conf.name:
                                results.remove(conf)
                                removelist.append(conf.name)
                                print(
                                    "{} is removed as requested by the "
                                    "user!".format(item)
                                )
            if not tmp_results and not results:
                print("ERROR: There are no conformers for part 4!\n" " Going to exit!")
                sys.exit()
            print(
                "number of further considered conformers: {:{digits}}"
                " {}".format(
                    "",
                    str(len(results) + len(tmp_results)),
                    digits=digilen - len("number of further considered conformers"),
                )
            )
            if tmp_results:
                print(
                    "The coupling constants were calculated before for {} "
                    "conformers:".format(str(len(tmp_results)))
                )
                print_block([i.name for i in tmp_results])
                if results:
                    print(
                        "The coupling constants are calculated now for {} "
                        "conformers:".format(str(len(results)))
                    )
                    print_block([i.name for i in results])
                else:
                    print("No conformers are considered additionally.")
            else:
                print("Considered conformers:")
                print_block([i.name for i in results])

            if results:
                # create NMR folders for shieldings and
                # couplings within the CONFX folders!
                print("Setting up new directories.")
                cwd = os.getcwd()
                for conf in list(results):
                    tmp_from = os.path.join(cwd, conf.name, args.func)
                    tmp_to = os.path.join(cwd, conf.name, "NMR")
                    if not os.path.isdir(tmp_to):
                        mkdir_p(tmp_to)
                    try:
                        shutil.copy(
                            os.path.join(tmp_from, "coord"),
                            os.path.join(tmp_to, "coord"),
                        )
                    except FileNotFoundError:
                        if not os.path.isfile(os.path.join(tmp_from, "coord")):
                            print(
                                "ERROR: while copying the coord file from {}! "
                                "The corresponding file does not exist.".format(
                                    tmp_from
                                )
                            )
                        elif not os.path.isdir(tmp_to):
                            print("ERROR: Could not create folder {}!".format(tmp_to))
                        print("ERROR: Removing conformer {}!".format(conf.name))
                        for nuc in spectrumlist:
                            json_dict[conf.name]["".join(nuc + "_J")] = "failed"
                        save_errors.append(
                            "Conformer {} was removed, because IO failed!".format(
                                conf.name
                            )
                        )
                        results.remove(conf)
                        removelist.append(conf.name)
                print("Constructed all new folders.")

                if not results:
                    print("ERROR: No conformers left!")
                    write_json("save_and_exit", json_dict, jsonfile)

                if job == tm_job:
                    # need converged mos, therefore,
                    # cefine and single-point have to be performed
                    # before NMR calculations
                    # place work in queue:

                    instructprep = {
                        "jobtype": "prep",
                        "chrg": args.chrg,
                        "unpaired": args.unpaired,
                        "func": args.funcJ,
                        "basis": args.basisJ,
                        "solv": args.solv,
                        "sm": args.sm4,
                        "NMR": True,
                        "progsettings": {"omp": args.omp, "tempprogpath": ""},
                    }
                    if args.hactive == "on":
                        instructprep["hactive"] = True
                    if args.cactive == "on":
                        instructprep["cactive"] = True
                    if args.factive == "on":
                        instructprep["factive"] = True
                    if args.pactive == "on":
                        instructprep["pactive"] = True
                    if args.siactive == "on":
                        instructprep["siactive"] = True
                    results = run_in_parallel(
                        q, resultq, job, maxthreads, results, instructprep, "NMR"
                    )
                    exit_log, fail_rate = check_tasks(results, args)

                    # sort out conformers for which the single-point calculation failed
                    for item in list(results):
                        if not item.success:
                            print(
                                "\nERROR: A problem has occurred in the "
                                "single-point preparation of {} required for "
                                "the coupling calculation with TM! The conformer "
                                "is removed.\n".format(item.name)
                            )
                            for nuc in spectrumlist:
                                # update json_dict
                                nuc_string = "".join(nuc + "_J")
                                json_dict[item.name][nuc_string] = "failed"
                            save_errors.append(
                                "Conformer {} was removed, because single-point "
                                "preparation failed!".format(item.name)
                            )
                            results.remove(item)
                            removelist.append(item.name)
                    if exit_log:
                        print(
                            "\nERROR: too many preparations failed ({:.2f} %)!".format(
                                fail_rate
                            )
                        )
                        write_json("save_and_exit", json_dict, jsonfile)
                    print("Preparation completed.")

                    # place SP work in queue:
                    instructsp = {
                        "jobtype": "sp",
                        "chrg": args.chrg,
                        "unpaired": args.unpaired,
                        "func": args.funcJ,
                        "solv": args.solv,
                        "basis": args.basisJ,
                        "progsettings": {"omp": args.omp, "tempprogpath": ""},
                    }
                    if job == tm_job:
                        instructsp["progsettings"]["tempprogpath"] = ""
                    elif job == orca_job:
                        instructsp["progsettings"]["tempprogpath"] = orcapath
                        instructsp["progsettings"]["orca_old"] = orca_old

                    results = run_in_parallel(
                        q, resultq, job, maxthreads, results, instructsp, "NMR"
                    )
                    exit_log, fail_rate = check_tasks(results, args)

                    # move control to control_J
                    for item in results:
                        try:
                            shutil.copy(
                                os.path.join(item.workdir, "control"),
                                os.path.join(item.workdir, "control_J"),
                            )
                        except FileNotFoundError:
                            pass

                    # sort out conformers for which the single-point calculation failed
                    for item in list(results):
                        if not item.success:
                            print(
                                "\nERROR: A problem has occurred in the "
                                "single-point calculation of {} required for "
                                "the coupling calculation with TM! The conformer "
                                "is removed.\n".format(item.name)
                            )
                            for nuc in spectrumlist:
                                # update json_dict
                                nuc_string = "".join(nuc + "_J")
                                json_dict[item.name][nuc_string] = "failed"
                            save_errors.append(
                                "Conformer {} was removed, because the "
                                "single-point calculation failed!".format(item.name)
                            )
                            results.remove(item)
                            removelist.append(item.name)
                    if exit_log:
                        print(
                            "\nERROR: too many single-points failed ({:.2f} %)!".format(
                                fail_rate
                            )
                        )
                        write_json("save_and_exit", json_dict, jsonfile)

                    # move mos to mos-keep_J
                    for item in results:
                        try:
                            shutil.copy(
                                os.path.join(item.workdir, "mos"),
                                os.path.join(item.workdir, "mos-keep_J"),
                            )
                        except FileNotFoundError:
                            pass
                    print("single-points completed.")

                # calculate J
                # place work in queue
                instructJ = {
                    "jobtype": "nmrJ",
                    "chrg": args.chrg,
                    "unpaired": args.unpaired,
                    "func": args.funcJ,
                    "basis": args.basisJ,
                    "solv": args.solv,
                    "sm": args.sm4,
                    "progsettings": {"omp": args.omp},
                }
                if args.hactive == "on":
                    instructJ["hactive"] = True
                if args.cactive == "on":
                    instructJ["cactive"] = True
                if args.factive == "on":
                    instructJ["factive"] = True
                if args.pactive == "on":
                    instructJ["pactive"] = True
                if args.siactive == "on":
                    instructJ["siactive"] = True
                if job == tm_job:
                    instructJ["progsettings"]["tempprogpath"] = escfpath
                elif job == orca_job:
                    instructJ["progsettings"]["tempprogpath"] = orcapath
                    instructJ["progsettings"]["orca_old"] = orca_old

                results = run_in_parallel(
                    q, resultq, job, maxthreads, results, instructJ, "NMR"
                )
                exit_log, fail_rate = check_tasks(results, args)

                # sort out conformers for which the coupling calculation failed
                for item in list(results):
                    if not item.success:
                        print(
                            "\nERROR: A problem has occurred in the coupling "
                            "calculation of {}! The conformer is removed.\n".format(
                                item.name
                            )
                        )
                        for nuc in spectrumlist:
                            json_dict[item.name]["".join(nuc + "_J")] = "failed"
                        save_errors.append(
                            "Conformer {} was removed, because the couplings calculation "
                            "failed!".format(item.name)
                        )
                        results.remove(item)
                        removelist.append(item.name)
                # update json_dict
                for item in results:
                    for nuc in spectrumlist:
                        nuc_string = "".join(nuc + "_J")
                        json_dict[item.name][nuc_string] = "calculated"

                if exit_log:
                    print(
                        "\nERROR: too many coupling calculations failed ({:.2f} %)!".format(
                            fail_rate
                        )
                    )
                    write_json("save_and_exit", json_dict, jsonfile)

            # adding conformers calculated before to results
            for item in tmp_results:
                results.append(item)
            results.sort(key=lambda x: int(x.name[4:]))

            # check if at least one conformer is left:
            if len(results) == 0:
                print("ERROR: No conformers left!")
                write_json("save_and_exit", json_dict, jsonfile)

            print("Coupling calculations completed.\n")

        # SHIELDINGS
        if args.calcS == "on":
            tmp_results = []
            # list with conformers calculated before
            if args.calcJ == "on" or args.part3 == "on":
                for item in list(results):
                    # it is possible to calculate multiple nuclei at once,
                    # so check for each nuc if it was calculated before or if it failed
                    tmp_failed = False
                    tmp_calculate = False
                    for nuc in spectrumlist:
                        nuc_string = "".join(nuc + "_S")
                        if json_dict[item.name][nuc_string] == "not_calculated":
                            tmp_calculate = True
                        elif json_dict[item.name][nuc_string] == "failed":
                            tmp_failed = True
                            print(
                                "ERROR: The {} shielding calculation for {} failed in the previous run! The "
                                "conformer is sorted out.".format(nuc, item.name)
                            )
                            break
                    if tmp_failed:
                        results.remove(item)
                        removelist.append(item.name)
                    elif not tmp_calculate and not tmp_failed:
                        # for this conformer, all activated shielding calculation
                        # were performed in the previous run
                        tmp_results.append(results.pop(results.index(item)))
                    elif tmp_calculate and not tmp_failed:
                        # this conformer has to be calculated now
                        pass
            elif args.part3 == "off" and args.calcJ == "off":
                # use all conformers with 'consider_for_part4' = True
                if firstrun:  # not possible
                    print(
                        "ERROR: Part 3 is switched off but the file {} containing"
                        "\n       all information of the previous "
                        "run cannot be found!".format(jsonfile)
                    )
                    write_json("save_and_exit", json_dict, jsonfile)
                results = []
                # list of conformers that have to be calculated now
                for item in json_dict:
                    if "CONF" in item and json_dict[item]["consider_for_part4"]:
                        # it is possible to calculate multiple nuclei at once,
                        # so check for each nuc if it was calculated before or if it failed
                        tmp_failed = False
                        tmp_calculate = False
                        for nuc in spectrumlist:
                            nuc_string = "".join(nuc + "_S")
                            if json_dict[item][nuc_string] == "not_calculated":
                                tmp_calculate = True
                            elif json_dict[item][nuc_string] == "failed":
                                tmp_failed = True
                                print(
                                    "ERROR: The {} shielding calculation for {} "
                                    "failed in the previous run!\n The "
                                    "conformer is sorted out.".format(nuc, item.name)
                                )
                                removelist.append(item.name)
                                break
                        if tmp_calculate and not tmp_failed:
                            # this conformer has be calculated now
                            tmp = job()
                            tmp.name = item
                            tmp.success = True
                            results.append(tmp)
                        elif not tmp_calculate and not tmp_failed:
                            # for this conformer, all activated coupling
                            # calculation were performed in the previous run
                            tmp = job()
                            tmp.name = item
                            tmp.success = True
                            tmp_results.append(tmp)

            if not tmp_results and not results:
                print("ERROR: There are no conformers for part 4!")
                sys.exit()
            print(
                "number of further considered conformers: {:{digits}} "
                "{}".format(
                    "",
                    str(len(results) + len(tmp_results)),
                    digits=digilen - len("number of further considered conformers"),
                )
            )
            if tmp_results:
                print(
                    "The shielding constants were calculated before for {} "
                    "conformers:".format(str(len(tmp_results)))
                )
                print_block([i.name for i in tmp_results])
                if results:
                    print(
                        "The shielding constants are calculated now for {} "
                        "conformers::".format(str(len(results)))
                    )
                    for i in results:
                        print(i.name)
                else:
                    print("No conformers are considered additionally.")
            else:
                print("Considered conformers:")
                print_block([i.name for i in results])

            if results:
                print("Setting up new directories.")
                cwd = os.getcwd()
                for conf in list(results):
                    tmp_from = os.path.join(cwd, conf.name, args.func)
                    tmp_to = os.path.join(cwd, conf.name, "NMR")
                    if not os.path.isdir(tmp_to):
                        mkdir_p(tmp_to)
                    try:
                        shutil.copy(
                            os.path.join(tmp_from, "coord"),
                            os.path.join(tmp_to, "coord"),
                        )
                    except FileNotFoundError:
                        if not os.path.isfile(os.path.join(tmp_from, "coord")):
                            print(
                                "ERROR: while copying the coord file from {}! The corresponding "
                                "file does not exist.".format(tmp_from)
                            )
                        elif not os.path.isdir(tmp_to):
                            print("ERROR: Could not create folder {}!".format(tmp_to))
                        print("ERROR: Removing conformer {}!".format(conf.name))
                        for nuc in spectrumlist:
                            json_dict[conf.name]["".join(nuc + "_S")] = "failed"
                        save_errors.append(
                            "Conformer {} was removed, because IO failed!".format(
                                conf.name
                            )
                        )
                        results.remove(conf)
                        removelist.append(conf.name)
                print("Constructed all new folders.")

                if not results:
                    print("ERROR: No conformers left!")
                    write_json("save_and_exit", json_dict, jsonfile)

                if job == tm_job:
                    do_cefine = False
                    if args.calcJ == "off":
                        do_cefine = True
                    else:
                        if args.funcJ != args.funcS and args.basisS != args.basisJ:
                            do_cefine = True
                    if do_cefine:
                        # need converged mos, therefore,
                        # cefine and single-point have to be performed before NMR calculations
                        # place work in queue:
                        instructprep = {
                            "jobtype": "prep",
                            "chrg": args.chrg,
                            "unpaired": args.unpaired,
                            "func": args.funcS,
                            "basis": args.basisS,
                            "solv": args.solv,
                            "sm": args.sm4,
                            "NMR": True,
                            "progsettings": {"omp": args.omp, "tempprogpath": ""},
                        }
                        if args.hactive == "on":
                            instructprep["hactive"] = True
                        if args.cactive == "on":
                            instructprep["cactive"] = True
                        if args.factive == "on":
                            instructprep["factive"] = True
                        if args.pactive == "on":
                            instructprep["pactive"] = True
                        if args.siactive == "on":
                            instructprep["siactive"] = True
                        results = run_in_parallel(
                            q, resultq, job, maxthreads, results, instructprep, "NMR"
                        )
                        exit_log, fail_rate = check_tasks(results, args)

                        # sort out conformers for which the single-point calculation failed
                        for item in list(results):
                            if not item.success:
                                print(
                                    "\nERROR: A problem has occurred in the single-point "
                                    "preparation of {} required for the shielding "
                                    "calculation with TM! The conformer is removed.\n".format(
                                        item.name
                                    )
                                )
                                for nuc in spectrumlist:
                                    # update json_dict
                                    nuc_string = "".join(nuc + "_S")
                                    json_dict[item.name][nuc_string] = "failed"
                                save_errors.append(
                                    "Conformer {} was removed, because single-point"
                                    " preparation failed!".format(item.name)
                                )
                                results.remove(item)
                                removelist.append(item.name)

                        if exit_log:
                            print(
                                "\nERROR: too many preparations failed ({:.2f} %)!".format(
                                    fail_rate
                                )
                            )
                            write_json("save_and_exit", json_dict, jsonfile)

                        print("Preparation completed.")

                        # place SP work in queue:
                        instructsp = {
                            "jobtype": "sp",
                            "chrg": args.chrg,
                            "unpaired": args.unpaired,
                            "func": args.funcS,
                            "solv": args.solv,
                            "basis": args.basisS,
                            "progsettings": {"omp": args.omp, "tempprogpath": ""},
                        }
                        if job == tm_job:
                            instructsp["progsettings"]["tempprogpath"] = ""
                        elif job == orca_job:
                            instructsp["progsettings"]["tempprogpath"] = orcapath
                            instructsp["progsettings"]["orca_old"] = orca_old
                        results = run_in_parallel(
                            q, resultq, job, maxthreads, results, instructsp, "NMR"
                        )
                        exit_log, fail_rate = check_tasks(results, args)

                        # move control to control_S
                        for item in results:
                            try:
                                shutil.copy(
                                    os.path.join(item.workdir, "control"),
                                    os.path.join(item.workdir, "control_S"),
                                )
                            except FileNotFoundError:
                                pass

                        # sort out conformers for which the single-point calculation failed
                        for item in list(results):
                            if not item.success:
                                print(
                                    "\nERROR: A problem has occurred in the single-point preparation of {} required "
                                    "for the shielding calculation with TM! The conformer is removed.\n".format(
                                        item.name
                                    )
                                )
                                for nuc in spectrumlist:
                                    # update json_dict
                                    nuc_string = "".join(nuc + "_S")
                                    json_dict[item.name][nuc_string] = "failed"
                                save_errors.append(
                                    "Conformer {} was removed, because single-point calculation "
                                    "failed!".format(item.name)
                                )
                                results.remove(item)
                                removelist.append(item.name)

                        # move mos to mos-keep_S
                        for item in results:
                            try:
                                shutil.copy(
                                    os.path.join(item.workdir, "mos"),
                                    os.path.join(item.workdir, "mos-keep_S"),
                                )
                            except FileNotFoundError:
                                pass

                        if exit_log:
                            print(
                                "\nERROR: too many single-points failed ({:.2f} %)!".format(
                                    fail_rate
                                )
                            )
                            write_json("save_and_exit", json_dict, jsonfile)

                        print("single-points completed.")

                instructS = {
                    "jobtype": "nmrS",
                    "chrg": args.chrg,
                    "unpaired": args.unpaired,
                    "func": args.funcS,
                    "basis": args.basisS,
                    "solv": args.solv,
                    "sm": args.sm4,
                    "progsettings": {"omp": args.omp, "tempprogpath": ""},
                }
                if args.hactive == "on":
                    instructS["hactive"] = True
                if args.cactive == "on":
                    instructS["cactive"] = True
                if args.factive == "on":
                    instructS["factive"] = True
                if args.pactive == "on":
                    instructS["pactive"] = True
                if args.siactive == "on":
                    instructS["siactive"] = True
                if job == tm_job:
                    instructS["progsettings"]["tempprogpath"] = mpshiftpath
                elif job == orca_job:
                    instructS["progsettings"]["tempprogpath"] = orcapath
                    instructS["progsettings"]["orca_old"] = orca_old

                results = run_in_parallel(
                    q, resultq, job, maxthreads, results, instructS, "NMR"
                )
                exit_log, fail_rate = check_tasks(results, args)

                # sort out conformers for which the coupling calculation failed
                for conf in list(results):
                    if not conf.success:
                        print(
                            "\nERROR: A problem has occurred in the shielding calculation of {}!\n".format(
                                conf.name
                            )
                        )
                        for nuc in spectrumlist:
                            json_dict[conf.name]["".join(nuc + "_S")] = "failed"
                        save_errors.append(
                            "Conformer {} was removed, because the shielding calculation "
                            "failed!".format(conf.name)
                        )
                        results.remove(conf)
                        removelist.append(conf.name)
                    else:
                        for nuc in spectrumlist:
                            json_dict[conf.name]["".join(nuc + "_S")] = "calculated"

                if exit_log:
                    print(
                        "ERROR: {:.2f} % of the shielding calculations failed!".format(
                            fail_rate
                        )
                    )
            print("Shielding calculations completed.")
        write_json("save", json_dict, jsonfile)
        # write generic:
        instructgeneric = {
            "jobtype": "genericout",
            "nat": nat,
            "progsettings": {"omp": args.omp, "tempprogpath": ""},
        }
        results = run_in_parallel(
            q, resultq, job, maxthreads, results + tmp_results, instructgeneric, "NMR"
        )
        # exit_log, fail_rate = check_tasks(results, args)

        # write .anmrrc
        print("\nwriting .anmrrc")
        writing_anmrrc(
            args.prog,
            args.func,
            args.sm,
            args.sm4,
            args.funcS,
            args.basisS,
            args.solv,
            args.mf,
            args.href,
            args.cref,
            args.fref,
            args.pref,
            args.siref,
            args.hactive,
            args.cactive,
            args.factive,
            args.pactive,
            args.siactive,
            args.calcJ,
            args.calcS,
            save_errors,
        )
        # check if rotamers of each other are still remaining in the final ensemble
        if rotdict:
            try:
                length = max([len(j) for i in rotdict.items() for j in i])
                for item in results + tmp_results:
                    if item.name in rotdict.keys():
                        if rotdict[item.name] in [
                            x.name for x in results + tmp_results
                        ]:
                            save_errors.append(
                                "Possible rotamers of each other still in final "
                                "ensemble: {:{digits}} <--> {:{digits}}. Please "
                                "check by hand!".format(
                                    item.name, rotdict[item.name], digits=length
                                )
                            )
            except:
                pass

        if save_errors:
            print("***---------------------------------------------------------***")
            print("Printing most relevant errors again, just for user convenience:")
            for error in list(save_errors):
                print(save_errors.pop())
            print("***---------------------------------------------------------***")

        ### if couplings or shieldings failed:
        if removelist:
            correct_anmr_enso(cwd, removelist)

        print("\n END of part4.\n")
    else:  # if part 4 is switched off
        print("PART4 has been skipped by user\n")
    print("-----------------------------------------------------------")
    if (
        args.part1 == "off"
        and args.part2 == "off"
        and args.part3 == "off"
        and args.part4 == "off"
    ):
        write_json("save", json_dict, jsonfile)
    return


# start of ENSO:
if __name__ == "__main__":

    cwd = os.getcwd()
    # RUNNING ENSO STARTUP
    (
        args,
        jsonfile,
        json_dict,
        firstrun,
        conformersxyz,
        nat,
        maxthreads,
        xtbpath,
        environsettings,
        dbpath,
        cosmothermversion,
        cosmorssetup,
        orcapath,
        orca_old,
        crestpath,
        mpshiftpath,
        escfpath,
        spectrumlist,
    ) = enso_startup(cwd)

    # RUNNING PART1
    try:
        results, json_dict = part1(
            args,
            jsonfile,
            json_dict,
            conformersxyz,
            nat,
            maxthreads,
            xtbpath,
            environsettings,
        )
    except Exception as e:
        print("ERROR in part1!")
        print("The error-message is {}\n\n".format(e))
        traceback.print_exc()
        print("Going to exit!")
        sys.exit(1)

    # RUNNING PART2
    try:
        results, json_dict, rotdict = part2(
            args,
            results,
            jsonfile,
            cwd,
            json_dict,
            maxthreads,
            xtbpath,
            environsettings,
        )
    except Exception as e:
        print("ERROR in part2!")
        print("The error-message is {}\n\n".format(e))
        traceback.print_exc()
        print("Going to exit!")
        sys.exit(1)

    # RUNNING PART3
    try:
        results, json_dict = part3(
            args,
            results,
            jsonfile,
            cwd,
            json_dict,
            maxthreads,
            xtbpath,
            environsettings,
        )
    except Exception as e:
        print("ERROR in part3!")
        print("The error-message is {}\n\n".format(e))
        traceback.print_exc()
        print("Going to exit!")
        sys.exit(1)

    # RUNNING PART4
    try:
        part4(
            args,
            results,
            jsonfile,
            cwd,
            json_dict,
            maxthreads,
            xtbpath,
            environsettings,
            rotdict,
            spectrumlist,
            escfpath,
            mpshiftpath,
        )
    except Exception as e:
        print("ERROR in part4!")
        print("The error-message is {}\n\n".format(e))
        traceback.print_exc()
        print("Going to exit!")
        sys.exit(1)

    # ALL DONE
    print("\nENSO all done!")
