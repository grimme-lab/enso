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


def cml(descr, solvents, func, func3, funcJ, funcS, gfnv, href, cref, fref, pref):
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
        J_list = ["1H_J", "13C_J", "19F_J", "31P_J"]
        S_list = ["1H_S", "13C_S", "19F_S", "31P_S"]
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
            "1H_S",
            "13C_S",
            "19F_S",
            "31P_S",
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
                    "with respect to the previous run.\n GRRHO is reset to "
                    "not_calculated for all conformers."
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
                if json_dict["flags"]["sm3"] in ("cosmors", "gbsa_gsolv") and flags_dict["sm3"] in ("cosmors", "gbsa_gsolv"):
                    # dont recalculate SP
                    pass
                else:
                    print(
                        "WARNING: The density functional / the basis set and/or"
                        " the solvent model for part 3 was changed \n"
                        "         with respect to the previous run. The "
                        "single-point of part 3 is reset to not_calculated"
                        " for all conformers."
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
                        "to the previous run. The coupling calculation is reset "
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
                        "to the previous run. The shielding calculation is "
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


def conformersxyz2coord(conformersxyz, nat, directory, conflist):
    """read conformersxyz or xtb_confg.xyz and write coord into 
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
        atom = []
        x = []
        y = []
        z = []
        gfne[counter][0] = "".join(("CONF", str(i)))
        gfne[counter][1] = check_for_float(stringfile_lines[(i - 1) * (nat + 2) + 1])
        if gfne[counter][1] is None:
            print(
                "Error in float conversion while reading file"
                " {}!".format(conformersxyz)
            )
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
                "{: 09.7f} {: 09.7f}  {: 09.7f}  {}".format(x[j], y[j], z[j], atom[j])
            )
        if not os.path.isfile(os.path.join("CONF{}".format(i), directory, "coord")):
            print(
                "Write new coord file in {}".format(
                    os.path.join("CONF{}".format(i), directory)
                )
            )
            with open(
                os.path.join("CONF{}".format(i), directory, "coord"), "w", newline=None
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
            inp.write("  {}\n".format(nat))  ### number of atoms
            inp.write("{:20.8f}        !{}\n".format(i.energy, i.name))  ### energy
            conf_xyz = coord2xyz(
                os.path.join(cwd, i.name, func)
            )  ### coordinates in xyz
            for j in conf_xyz:
                inp.write(
                    " {:3} {:19.10f} {:19.10f} {:19.10f}\n".format(
                        j[0].upper(), float(j[1]), float(j[2]), float(j[3])
                    )
                )

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
                out.write("  {}\n".format(nat))  ### number of atoms
                out.write("E= {:20.8f}  G= {:20.8f}      !{}\n".format(i.energy, i.free_energy, i.name))  ### energy
                conf_xyz = coord2xyz(
                    os.path.join(cwd, i.name, optfolder)
                )  ### coordinates in xyz
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
        coordxyz = []
        for i in range(len(x)):
            coordxyz.append([atom[i], x[i], y[i], z[i]])
    return coordxyz


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
        fmtlenght = 7
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
    hactive,
    cactive,
    factive,
    pactive,
    calcJ,
    calcS,
):

    if func == "pbe0":
        func = "pbeh-3c"

    # 1H
    h_tm_shieldings = {
        "pbeh-3c": {
            "pbe0": {
                "TMS": {
                    "acetone": ["31.778"],
                    "chcl3": ["31.786"],
                    "ch2cl2": ["31.785"],
                    "dmso": ["31.776"],
                    "h2o": ["31.775"],
                    "methanol": ["31.780"],
                    "thf": ["31.783"],
                    "toluene": ["31.790"],
                    "gas": ["32.022"],
                },
                "DSS": {
                    "acetone": ["31.763"],
                    "chcl3": ["31.767"],
                    "ch2cl2": ["31.767"],
                    "dmso": ["31.763"],
                    "h2o": ["31.762"],
                    "methanol": ["31.766"],
                    "thf": ["31.765"],
                    "toluene": ["31.768"],
                    "gas": ["31.793"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["32.015"],
                    "chcl3": ["32.022"],
                    "ch2cl2": ["32.022"],
                    "dmso": ["32.014"],
                    "h2o": ["32.013"],
                    "methanol": ["32.018"],
                    "thf": ["32.019"],
                    "toluene": ["32.023"],
                    "gas": ["32.051"],
                },
                "DSS": {
                    "acetone": ["31.998"],
                    "chcl3": ["32.000"],
                    "ch2cl2": ["32.002"],
                    "dmso": ["31.998"],
                    "h2o": ["31.997"],
                    "methanol": ["32.002"],
                    "thf": ["31.999"],
                    "toluene": ["32.000"],
                    "gas": ["32.021"],
                },
            },
        },
        "b97-3c": {
            "pbe0": {
                "TMS": {
                    "acetone": ["31.825"],
                    "chcl3": ["31.832"],
                    "ch2cl2": ["31.830"],
                    "dmso": ["31.824"],
                    "h2o": ["31.821"],
                    "methanol": ["31.826"],
                    "thf": ["31.829"],
                    "toluene": ["31.836"],
                    "gas": ["31.864"],
                },
                "DSS": {
                    "acetone": ["31.811"],
                    "chcl3": ["31.811"],
                    "ch2cl2": ["31.812"],
                    "dmso": ["31.810"],
                    "h2o": ["31.809"],
                    "methanol": ["31.813"],
                    "thf": ["31.814"],
                    "toluene": ["31.815"],
                    "gas": ["31.793"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["32.061"],
                    "chcl3": ["32.067"],
                    "ch2cl2": ["32.066"],
                    "dmso": ["32.061"],
                    "h2o": ["32.058"],
                    "methanol": ["32.063"],
                    "thf": ["32.065"],
                    "toluene": ["32.069"],
                    "gas": ["32.093"],
                },
                "DSS": {
                    "acetone": ["32.045"],
                    "chcl3": ["32.044"],
                    "ch2cl2": ["32.046"],
                    "dmso": ["32.045"],
                    "h2o": ["32.043"],
                    "methanol": ["32.048"],
                    "thf": ["32.048"],
                    "toluene": ["32.046"],
                    "gas": ["32.021"],
                },
            },
        },
        "tpss": {
            "pbe0": {
                "TMS": {
                    "acetone": ["31.598"],
                    "chcl3": ["31.605"],
                    "ch2cl2": ["31.604"],
                    "dmso": ["31.629"],
                    "h2o": ["31.596"],
                    "methanol": ["31.599"],
                    "thf": ["31.603"],
                    "toluene": ["31.609"],
                    "gas": ["31.638"],
                },
                "DSS": {
                    "acetone": ["31.584"],
                    "chcl3": ["31.587"],
                    "ch2cl2": ["31.589"],
                    "dmso": ["31.585"],
                    "h2o": ["31.582"],
                    "methanol": ["31.588"],
                    "thf": ["31.588"],
                    "toluene": ["31.590"],
                    "gas": ["31.793"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["31.836"],
                    "chcl3": ["31.841"],
                    "ch2cl2": ["31.841"],
                    "dmso": ["31.835"],
                    "h2o": ["31.834"],
                    "methanol": ["31.838"],
                    "thf": ["31.840"],
                    "toluene": ["31.843"],
                    "gas": ["31.869"],
                },
                "DSS": {
                    "acetone": ["31.820"],
                    "chcl3": ["31.821"],
                    "ch2cl2": ["31.824"],
                    "dmso": ["31.821"],
                    "h2o": ["31.818"],
                    "methanol": ["31.824"],
                    "thf": ["31.823"],
                    "toluene": ["31.823"],
                    "gas": ["32.021"],
                },
            },
        },
    }
    h_orca_shieldings = {
        "pbeh-3c": {
            "pbe0": {
                "TMS": {
                    "acetone": ["31.506"],
                    "chcl3": ["31.516"],
                    "ch2cl2": ["31.510"],
                    "dmso": ["31.506"],
                    "h2o": ["31.517"],
                    "methanol": ["31.504"],
                    "thf": ["31.512"],
                    "toluene": ["31.530"],
                    "gas": ["31.576"],
                },
                "DSS": {
                    "acetone": ["31.491"],
                    "chcl3": ["31.500"],
                    "ch2cl2": ["31.493"],
                    "dmso": ["31.492"],
                    "h2o": ["31.499"],
                    "methanol": ["31.483"],
                    "thf": ["31.497"],
                    "toluene": ["31.509"],
                    "gas": ["31.541"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["31.837"],
                    "chcl3": ["31.845"],
                    "ch2cl2": ["31.840"],
                    "dmso": ["31.838"],
                    "h2o": ["31.848"],
                    "methanol": ["31.835"],
                    "thf": ["31.842"],
                    "toluene": ["31.857"],
                    "gas": ["31.896"],
                },
                "DSS": {
                    "acetone": ["31.818"],
                    "chcl3": ["31.823"],
                    "ch2cl2": ["31.819"],
                    "dmso": ["31.819"],
                    "h2o": ["31.828"],
                    "methanol": ["31.809"],
                    "thf": ["31.822"],
                    "toluene": ["31.831"],
                    "gas": ["31.542"],
                },
            },
        },
        "b97-3c": {
            "pbe0": {
                "TMS": {
                    "acetone": ["31.556"],
                    "chcl3": ["31.566"],
                    "ch2cl2": ["31.560"],
                    "dmso": ["31.558"],
                    "h2o": ["31.571"],
                    "methanol": ["31.554"],
                    "thf": ["31.563"],
                    "toluene": ["31.581"],
                    "gas": ["31.623"],
                },
                "DSS": {
                    "acetone": ["31.539"],
                    "chcl3": ["31.545"],
                    "ch2cl2": ["31.542"],
                    "dmso": ["31.540"],
                    "h2o": ["31.554"],
                    "methanol": ["31.540"],
                    "thf": ["31.542"],
                    "toluene": ["31.553"],
                    "gas": ["31.586"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["31.887"],
                    "chcl3": ["31.895"],
                    "ch2cl2": ["31.890"],
                    "dmso": ["31.889"],
                    "h2o": ["31.903"],
                    "methanol": ["31.885"],
                    "thf": ["31.893"],
                    "toluene": ["31.908"],
                    "gas": ["31.944"],
                },
                "DSS": {
                    "acetone": ["31.862"],
                    "chcl3": ["31.865"],
                    "ch2cl2": ["31.863"],
                    "dmso": ["31.863"],
                    "h2o": ["31.876"],
                    "methanol": ["31.860"],
                    "thf": ["31.863"],
                    "toluene": ["31.871"],
                    "gas": ["31.898"],
                },
            },
        },
        "tpss": {
            "pbe0": {
                "TMS": {
                    "acetone": ["31.321"],
                    "chcl3": ["31.330"],
                    "ch2cl2": ["31.324"],
                    "dmso": ["31.322"],
                    "h2o": ["31.333"],
                    "methanol": ["31.318"],
                    "thf": ["31.327"],
                    "toluene": ["31.343"],
                    "gas": ["31.385"],
                },
                "DSS": {
                    "acetone": ["31.305"],
                    "chcl3": ["31.313"],
                    "ch2cl2": ["31.307"],
                    "dmso": ["31.304"],
                    "h2o": ["31.316"],
                    "methanol": ["31.302"],
                    "thf": ["31.310"],
                    "toluene": ["31.305"],
                    "gas": ["31.354"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["31.652"],
                    "chcl3": ["31.659"],
                    "ch2cl2": ["31.655"],
                    "dmso": ["31.653"],
                    "h2o": ["31.665"],
                    "methanol": ["31.649"],
                    "thf": ["31.657"],
                    "toluene": ["31.669"],
                    "gas": ["31.705"],
                },
                "DSS": {
                    "acetone": ["31.631"],
                    "chcl3": ["31.633"],
                    "ch2cl2": ["31.632"],
                    "dmso": ["31.631"],
                    "h2o": ["31.644"],
                    "methanol": ["31.629"],
                    "thf": ["31.633"],
                    "toluene": ["31.619"],
                    "gas": ["31.667"],
                },
            },
        },
    }
    # 13C
    c_tm_shieldings = {
        "pbeh-3c": {
            "pbe0": {
                "TMS": {
                    "acetone": ["189.942"],
                    "chcl3": ["189.674"],
                    "ch2cl2": ["189.860"],
                    "dmso": ["189.986"],
                    "h2o": ["189.999"],
                    "methanol": ["189.997"],
                    "thf": ["189.798"],
                    "toluene": ["189.348"],
                    "gas": ["188.900"],
                },
                "DSS": {
                    "acetone": ["191.817"],
                    "chcl3": ["191.620"],
                    "ch2cl2": ["191.802"],
                    "dmso": ["191.822"],
                    "h2o": ["191.906"],
                    "methanol": ["191.882"],
                    "thf": ["191.664"],
                    "toluene": ["191.290"],
                    "gas": ["190.866"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["187.647"],
                    "chcl3": ["187.407"],
                    "ch2cl2": ["187.573"],
                    "dmso": ["187.689"],
                    "h2o": ["187.700"],
                    "methanol": ["187.701"],
                    "thf": ["187.513"],
                    "toluene": ["187.103"],
                    "gas": ["186.697"],
                },
                "DSS": {
                    "acetone": ["189.361"],
                    "chcl3": ["189.180"],
                    "ch2cl2": ["189.351"],
                    "dmso": ["189.364"],
                    "h2o": ["189.443"],
                    "methanol": ["189.423"],
                    "thf": ["189.219"],
                    "toluene": ["188.873"],
                    "gas": ["188.488"],
                },
            },
        },
        "b97-3c": {
            "pbe0": {
                "TMS": {
                    "acetone": ["190.235"],
                    "chcl3": ["189.958"],
                    "ch2cl2": ["190.144"],
                    "dmso": ["190.290"],
                    "h2o": ["190.295"],
                    "methanol": ["190.288"],
                    "thf": ["190.092"],
                    "toluene": ["189.633"],
                    "gas": ["189.177"],
                },
                "DSS": {
                    "acetone": ["192.102"],
                    "chcl3": ["191.907"],
                    "ch2cl2": ["192.087"],
                    "dmso": ["192.108"],
                    "h2o": ["192.206"],
                    "methanol": ["192.062"],
                    "thf": ["191.970"],
                    "toluene": ["191.588"],
                    "gas": ["190.866"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["187.957"],
                    "chcl3": ["187.708"],
                    "ch2cl2": ["187.881"],
                    "dmso": ["188.008"],
                    "h2o": ["188.012"],
                    "methanol": ["188.007"],
                    "thf": ["187.832"],
                    "toluene": ["187.404"],
                    "gas": ["186.989"],
                },
                "DSS": {
                    "acetone": ["189.665"],
                    "chcl3": ["189.486"],
                    "ch2cl2": ["189.655"],
                    "dmso": ["189.669"],
                    "h2o": ["189.761"],
                    "methanol": ["189.625"],
                    "thf": ["189.545"],
                    "toluene": ["189.189"],
                    "gas": ["188.488"],
                },
            },
        },
        "tpss": {
            "pbe0": {
                "TMS": {
                    "acetone": ["188.783"],
                    "chcl3": ["188.500"],
                    "ch2cl2": ["188.692"],
                    "dmso": ["187.554"],
                    "h2o": ["188.851"],
                    "methanol": ["188.839"],
                    "thf": ["188.636"],
                    "toluene": ["188.158"],
                    "gas": ["187.671"],
                },
                "DSS": {
                    "acetone": ["190.655"],
                    "chcl3": ["190.446"],
                    "ch2cl2": ["190.643"],
                    "dmso": ["190.670"],
                    "h2o": ["190.757"],
                    "methanol": ["190.714"],
                    "thf": ["190.516"],
                    "toluene": ["190.115"],
                    "gas": ["190.866"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["186.481"],
                    "chcl3": ["186.216"],
                    "ch2cl2": ["186.397"],
                    "dmso": ["186.530"],
                    "h2o": ["186.544"],
                    "methanol": ["186.534"],
                    "thf": ["186.343"],
                    "toluene": ["185.907"],
                    "gas": ["185.463"],
                },
                "DSS": {
                    "acetone": ["188.192"],
                    "chcl3": ["188.001"],
                    "ch2cl2": ["188.186"],
                    "dmso": ["188.204"],
                    "h2o": ["188.287"],
                    "methanol": ["188.248"],
                    "thf": ["188.066"],
                    "toluene": ["187.693"],
                    "gas": ["188.488"],
                },
            },
        },
    }
    c_orca_shieldings = {
        "pbeh-3c": {
            "pbe0": {
                "TMS": {
                    "acetone": ["187.333"],
                    "chcl3": ["187.153"],
                    "ch2cl2": ["187.261"],
                    "dmso": ["187.383"],
                    "h2o": ["187.459"],
                    "methanol": ["187.350"],
                    "thf": ["187.234"],
                    "toluene": ["186.919"],
                    "gas": ["186.373"],
                },
                "DSS": {
                    "acetone": ["189.411"],
                    "chcl3": ["189.270"],
                    "ch2cl2": ["189.348"],
                    "dmso": ["189.512"],
                    "h2o": ["189.583"],
                    "methanol": ["189.249"],
                    "thf": ["189.365"],
                    "toluene": ["189.075"],
                    "gas": ["188.510"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["188.595"],
                    "chcl3": ["188.442"],
                    "ch2cl2": ["188.535"],
                    "dmso": ["188.639"],
                    "h2o": ["188.714"],
                    "methanol": ["188.608"],
                    "thf": ["188.513"],
                    "toluene": ["188.246"],
                    "gas": ["187.761"],
                },
                "DSS": {
                    "acetone": ["190.490"],
                    "chcl3": ["190.362"],
                    "ch2cl2": ["190.440"],
                    "dmso": ["190.582"],
                    "h2o": ["190.654"],
                    "methanol": ["190.323"],
                    "thf": ["190.451"],
                    "toluene": ["190.212"],
                    "gas": ["188.516"],
                },
            },
        },
        "b97-3c": {
            "pbe0": {
                "TMS": {
                    "acetone": ["187.699"],
                    "chcl3": ["187.527"],
                    "ch2cl2": ["187.626"],
                    "dmso": ["187.764"],
                    "h2o": ["187.847"],
                    "methanol": ["187.724"],
                    "thf": ["187.609"],
                    "toluene": ["187.316"],
                    "gas": ["186.704"],
                },
                "DSS": {
                    "acetone": ["189.731"],
                    "chcl3": ["189.573"],
                    "ch2cl2": ["189.667"],
                    "dmso": ["189.812"],
                    "h2o": ["189.912"],
                    "methanol": ["189.761"],
                    "thf": ["189.671"],
                    "toluene": ["189.372"],
                    "gas": ["188.815"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["188.961"],
                    "chcl3": ["188.813"],
                    "ch2cl2": ["188.897"],
                    "dmso": ["189.019"],
                    "h2o": ["189.100"],
                    "methanol": ["188.981"],
                    "thf": ["188.884"],
                    "toluene": ["188.638"],
                    "gas": ["188.089"],
                },
                "DSS": {
                    "acetone": ["190.792"],
                    "chcl3": ["190.653"],
                    "ch2cl2": ["190.737"],
                    "dmso": ["190.861"],
                    "h2o": ["190.957"],
                    "methanol": ["190.800"],
                    "thf": ["190.745"],
                    "toluene": ["190.490"],
                    "gas": ["190.037"],
                },
            },
        },
        "tpss": {
            "pbe0": {
                "TMS": {
                    "acetone": ["186.129"],
                    "chcl3": ["185.938"],
                    "ch2cl2": ["186.051"],
                    "dmso": ["186.190"],
                    "h2o": ["186.267"],
                    "methanol": ["186.151"],
                    "thf": ["186.027"],
                    "toluene": ["185.693"],
                    "gas": ["185.081"],
                },
                "DSS": {
                    "acetone": ["187.985"],
                    "chcl3": ["187.816"],
                    "ch2cl2": ["187.897"],
                    "dmso": ["188.044"],
                    "h2o": ["188.171"],
                    "methanol": ["188.012"],
                    "thf": ["187.876"],
                    "toluene": ["186.552"],
                    "gas": ["187.040"],
                },
            },
            "tpss": {
                "TMS": {
                    "acetone": ["187.408"],
                    "chcl3": ["187.246"],
                    "ch2cl2": ["187.342"],
                    "dmso": ["187.464"],
                    "h2o": ["187.541"],
                    "methanol": ["187.428"],
                    "thf": ["187.324"],
                    "toluene": ["187.039"],
                    "gas": ["186.493"],
                },
                "DSS": {
                    "acetone": ["189.091"],
                    "chcl3": ["188.931"],
                    "ch2cl2": ["189.014"],
                    "dmso": ["189.143"],
                    "h2o": ["189.270"],
                    "methanol": ["189.107"],
                    "thf": ["188.991"],
                    "toluene": ["187.788"],
                    "gas": ["188.298"],
                },
            },
        },
    }
    # 19F
    f_tm_shieldings = {
        "pbeh-3c": {
            "pbe0": {
                "CFCl3": {
                    "acetone": ["181.52"],
                    "chcl3": ["182.57"],
                    "ch2cl2": ["182.38"],
                    "dmso": ["182.48"],
                    "h2o": ["181.76"],
                    "methanol": ["181.04"],
                    "thf": ["181.25"],
                    "toluene": ["180.84"],
                    "gas": ["168.09"],
                }
            },
            "tpss": {
                "CFCl3": {
                    "acetone": ["165.40"],
                    "chcl3": ["166.49"],
                    "ch2cl2": ["166.29"],
                    "dmso": ["166.35"],
                    "h2o": ["165.63"],
                    "methanol": ["164.92"],
                    "thf": ["165.15"],
                    "toluene": ["164.81"],
                    "gas": ["152.16"],
                }
            },
        },
        "b97-3c": {
            "pbe0": {
                "CFCl3": {
                    "acetone": ["169.13"],
                    "chcl3": ["170.11"],
                    "ch2cl2": ["169.97"],
                    "dmso": ["170.09"],
                    "h2o": ["168.87"],
                    "methanol": ["168.24"],
                    "thf": ["168.84"],
                    "toluene": ["168.21"],
                    "gas": ["155.60"],
                }
            },
            "tpss": {
                "CFCl3": {
                    "acetone": ["151.51"],
                    "chcl3": ["152.57"],
                    "ch2cl2": ["152.38"],
                    "dmso": ["152.45"],
                    "h2o": ["151.25"],
                    "methanol": ["150.62"],
                    "thf": ["151.26"],
                    "toluene": ["150.74"],
                    "gas": ["138.28"],
                }
            },
        },
        "tpss": {
            "pbe0": {
                "CFCl3": {
                    "acetone": ["165.42"],
                    "chcl3": ["166.38"],
                    "ch2cl2": ["166.26"],
                    "dmso": ["166.36"],
                    "h2o": ["165.55"],
                    "methanol": ["164.95"],
                    "thf": ["165.10"],
                    "toluene": ["164.51"],
                    "gas": ["151.77"],
                }
            },
            "tpss": {
                "CFCl3": {
                    "acetone": ["148.14"],
                    "chcl3": ["149.16"],
                    "ch2cl2": ["149.00"],
                    "dmso": ["149.08"],
                    "h2o": ["148.27"],
                    "methanol": ["147.67"],
                    "thf": ["147.85"],
                    "toluene": ["147.35"],
                    "gas": ["151.77"],
                }
            },
        },
    }
    f_orca_shieldings_old = {
        "pbeh-3c": {
            "pbe0": {
                "CFCl3": {
                    "acetone": ["171.51"],
                    "chcl3": ["172.63"],
                    "ch2cl2": ["172.37"],
                    "dmso": ["172.47"],
                    "h2o": ["171.78"],
                    "methanol": ["171.00"],
                    "thf": ["171.31"],
                    "toluene": ["170.95"],
                    "gas": ["164.70"],
                }
            },
            "tpss": {
                "CFCl3": {
                    "acetone": ["159.20"],
                    "chcl3": ["160.37"],
                    "ch2cl2": ["160.08"],
                    "dmso": ["160.15"],
                    "h2o": ["159.46"],
                    "methanol": ["158.69"],
                    "thf": ["159.03"],
                    "toluene": ["158.77"],
                    "gas": ["153.18"],
                }
            },
        },
        "b97-3c": {
            "pbe0": {
                "CFCl3": {
                    "acetone": ["171.54"],
                    "chcl3": ["172.65"],
                    "ch2cl2": ["172.40"],
                    "dmso": ["172.50"],
                    "h2o": ["171.81"],
                    "methanol": ["171.03"],
                    "thf": ["171.33"],
                    "toluene": ["172.42"],
                    "gas": ["151.38"],
                }
            },
            "tpss": {
                "CFCl3": {
                    "acetone": ["159.15"],
                    "chcl3": ["160.32"],
                    "ch2cl2": ["160.03"],
                    "dmso": ["160.10"],
                    "h2o": ["159.40"],
                    "methanol": ["158.63"],
                    "thf": ["158.98"],
                    "toluene": ["160.03"],
                    "gas": ["138.67"],
                }
            },
        },
        "tpss": {
            "pbe0": {
                "CFCl3": {
                    "acetone": ["171.54"],
                    "chcl3": ["172.65"],
                    "ch2cl2": ["172.40"],
                    "dmso": ["172.50"],
                    "h2o": ["171.81"],
                    "methanol": ["171.03"],
                    "thf": ["171.33"],
                    "toluene": ["170.97"],
                    "gas": ["138.67"],
                }
            },
            "tpss": {
                "CFCl3": {
                    "acetone": ["159.15"],
                    "chcl3": ["160.32"],
                    "ch2cl2": ["160.03"],
                    "dmso": ["160.10"],
                    "h2o": ["159.40"],
                    "methanol": ["158.63"],
                    "thf": ["158.98"],
                    "toluene": ["158.72"],
                    "gas": ["134.81"],
                }
            },
        },
    }
    f_orca_shieldings_new = {
        "pbeh-3c": {
            "pbe0": {
                "CFCl3": {
                    "acetone": ["180.16"],
                    "chcl3": ["177.87"],
                    "ch2cl2": ["178.81"],
                    "dmso": ["179.48"],
                    "h2o": ["180.45"],
                    "methanol": ["180.80"],
                    "thf": ["179.70"],
                    "toluene": ["178.06"],
                    "gas": ["188.02"],
                }
            },
            "tpss": {
                "CFCl3": {
                    "acetone": ["168.35"],
                    "chcl3": ["166.10"],
                    "ch2cl2": ["167.02"],
                    "dmso": ["167.67"],
                    "h2o": ["168.65"],
                    "methanol": ["168.99"],
                    "thf": ["167.92"],
                    "toluene": ["166.35"],
                    "gas": ["176.45"],
                }
            },
        },
        "b97-3c": {
            "pbe0": {
                "CFCl3": {
                    "acetone": ["166.76"],
                    "chcl3": ["164.46"],
                    "ch2cl2": ["165.41"],
                    "dmso": ["166.12"],
                    "h2o": ["167.15"],
                    "methanol": ["167.40"],
                    "thf": ["166.31"],
                    "toluene": ["164.67"],
                    "gas": ["174.68"],
                }
            },
            "tpss": {
                "CFCl3": {
                    "acetone": ["153.68"],
                    "chcl3": ["151.45"],
                    "ch2cl2": ["152.36"],
                    "dmso": ["153.04"],
                    "h2o": ["154.09"],
                    "methanol": ["154.32"],
                    "thf": ["153.28"],
                    "toluene": ["151.74"],
                    "gas": ["161.97"],
                }
            },
        },
        "tpss": {
            "pbe0": {
                "CFCl3": {
                    "acetone": ["163.16"],
                    "chcl3": ["160.80"],
                    "ch2cl2": ["161.79"],
                    "dmso": ["162.53"],
                    "h2o": ["163.57"],
                    "methanol": ["163.81"],
                    "thf": ["162.68"],
                    "toluene": ["160.91"],
                    "gas": ["170.71"],
                }
            },
            "tpss": {
                "CFCl3": {
                    "acetone": ["150.23"],
                    "chcl3": ["147.93"],
                    "ch2cl2": ["148.88"],
                    "dmso": ["149.60"],
                    "h2o": ["150.66"],
                    "methanol": ["150.88"],
                    "thf": ["149.79"],
                    "toluene": ["148.12"],
                    "gas": ["158.11"],
                }
            },
        },
    }

    # 31P
    p_tm_shieldings = {
        "pbeh-3c": {
            "pbe0": {
                "TMP": {
                    "acetone": ["292.9"],
                    "chcl3": ["291.9"],
                    "ch2cl2": ["291.7"],
                    "dmso": ["291.7"],
                    "h2o": ["292.4"],
                    "methanol": ["291.3"],
                    "thf": ["292.9"],
                },
                "PH3": {"toluene": ["328.1"], "gas": ["308.8"]},
            },
            "tpss": {
                "TMP": {
                    "acetone": ["297.9"],
                    "chcl3": ["296.9"],
                    "ch2cl2": ["296.7"],
                    "dmso": ["296.6"],
                    "h2o": ["297.6"],
                    "methanol": ["296.3"],
                    "thf": ["297.9"],
                },
                "PH3": {"toluene": ["315.2"], "gas": ["295.8"]},
            },
        },
        "b97-3c": {
            "pbe0": {
                "TMP": {
                    "acetone": ["285.9"],
                    "chcl3": ["284.7"],
                    "ch2cl2": ["284.5"],
                    "dmso": ["284.1"],
                    "h2o": ["285.4"],
                    "methanol": ["284.5"],
                    "thf": ["285.7"],
                },
                "PH3": {"toluene": ["326.2"], "gas": ["306.6"]},
            },
            "tpss": {
                "TMP": {
                    "acetone": ["290.8"],
                    "chcl3": ["289.6"],
                    "ch2cl2": ["289.4"],
                    "dmso": ["288.9"],
                    "h2o": ["290.5"],
                    "methanol": ["289.5"],
                    "thf": ["290.6"],
                },
                "PH3": {"toluene": ["292.6"], "gas": ["313.2"]},
            },
        },
        "tpss": {
            "pbe0": {
                "TMP": {
                    "acetone": ["289.0"],
                    "chcl3": ["287.9"],
                    "ch2cl2": ["288.4"],
                    "dmso": ["287.7"],
                    "h2o": ["289.0"],
                    "methanol": ["288.2"],
                    "thf": ["288.6"],
                },
                "PH3": {"toluene": ["325.2"], "gas": ["305.7"]},
            },
            "tpss": {
                "TMP": {
                    "acetone": ["293.7"],
                    "chcl3": ["292.8"],
                    "ch2cl2": ["293.3"],
                    "dmso": ["292.5"],
                    "h2o": ["294.1"],
                    "methanol": ["293.1"],
                    "thf": ["294.4"],
                },
                "PH3": {"toluene": ["312.2"], "gas": ["292.6"]},
            },
        },
    }
    p_orca_shieldings = {
        "pbeh-3c": {
            "pbe0": {
                "TMP": {
                    "acetone": ["307.9"],
                    "chcl3": ["306.6"],
                    "ch2cl2": ["306.4"],
                    "dmso": ["306.7"],
                    "h2o": ["307.0"],
                    "methanol": ["306.3"],
                    "thf": ["307.7"],
                },
                "PH3": {"toluene": ["338.4"], "gas": ["316.6"]},
            },
            "tpss": {
                "TMP": {
                    "acetone": ["318.8"],
                    "chcl3": ["317.4"],
                    "ch2cl2": ["317.2"],
                    "dmso": ["317.4"],
                    "h2o": ["317.8"],
                    "methanol": ["317.1"],
                    "thf": ["318.5"],
                },
                "PH3": {"toluene": ["343.5"], "gas": ["322.0"]},
            },
        },
        "b97-3c": {
            "pbe0": {
                "TMP": {
                    "acetone": ["300.8"],
                    "chcl3": ["300.1"],
                    "ch2cl2": ["298.6"],
                    "dmso": ["299.9"],
                    "h2o": ["300.2"],
                    "methanol": ["299.5"],
                    "thf": ["300.6"],
                },
                "PH3": {"toluene": ["336.0"], "gas": ["313.9"]},
            },
            "tpss": {
                "TMP": {
                    "acetone": ["311.6"],
                    "chcl3": ["310.9"],
                    "ch2cl2": ["309.4"],
                    "dmso": ["310.7"],
                    "h2o": ["311.1"],
                    "methanol": ["310.3"],
                    "thf": ["311.5"],
                },
                "PH3": {"toluene": ["341.1"], "gas": ["319.3"]},
            },
        },
        "tpss": {
            "pbe0": {
                "TMP": {
                    "acetone": ["301.8"],
                    "chcl3": ["301.3"],
                    "ch2cl2": ["300.5"],
                    "dmso": ["300.6"],
                    "h2o": ["301.3"],
                    "methanol": ["301.1"],
                    "thf": ["301.7"],
                },
                "PH3": {"toluene": ["334.8"], "gas": ["312.7"]},
            },
            "tpss": {
                "TMP": {
                    "acetone": ["312.6"],
                    "chcl3": ["312.0"],
                    "ch2cl2": ["311.2"],
                    "dmso": ["311.3"],
                    "h2o": ["312.1"],
                    "methanol": ["311.8"],
                    "thf": ["312.4"],
                },
                "PH3": {"toluene": ["339.9"], "gas": ["318.1"]},
            },
        },
    }

    # shieldings and solvent models
    if solv is None:
        solv = "gas"
    if prog == "tm":
        # print('NMR data: func {}, funcS {}, href {}, solv {}'.format(str(func), str(funcS), str(href), str(solv)))
        hshielding = h_tm_shieldings[func][funcS][href][solv][0]
        cshielding = c_tm_shieldings[func][funcS][cref][solv][0]
        fshielding = f_tm_shieldings[func][funcS][fref][solv][0]
        pshielding = p_tm_shieldings[func][funcS][pref][solv][0]
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
        sm4 = "COSMO"
        basisS = "def2-TZVP"
    elif prog == "orca":
        hshielding = h_orca_shieldings[func][funcS][href][solv][0]
        cshielding = c_orca_shieldings[func][funcS][cref][solv][0]
        if orca_old:
            fshielding = f_orca_shieldings_old[func][funcS][fref][solv][0]
        else:
            fshielding = f_orca_shieldings_new[func][funcS][fref][solv][0]
        pshielding = p_orca_shieldings[func][funcS][pref][solv][0]
        if sm == "cpcm":
            print(
                "WARNING: The geometry optimization of the reference molecule "
                "was calculated with SMD instead of CPCM as solvent model (sm)!"
            )
        if sm4 == "smd":
            print(
                "WARNING: The reference shielding was calculated with CPCM "
                "instead of SMD as solvent model (sm)!"
            )
        if basisS != "pcSseg-2":
            print(
                "WARNING: The reference shielding was calculated with the "
                "basis pcSseg-2 (basisS)!"
            )
        sm = "SMD"
        sm4 = "CPCM"
        basisS = "pcSseg-2"
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
    if hactive == "on":
        ha = 1
    else:
        ha = 0
    if cactive == "on":
        ca = 1
    else:
        ca = 0
    if factive == "on":
        fa = 1
    else:
        fa = 0
    if pactive == "on":
        pa = 1
    else:
        pa = 0

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
                [len(i) for i in [hshielding, cshielding, fshielding, pshielding]]
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
    factive = False
    symmetry = "C1"
    temperature = 298.15
    nat = 0
    progsettings = {
                    "tempprogpath": "", 
                    "xtbpath": "",
                    "orca_old" : "",
                    "omp": 1,
                    "cosmorssetup": None,
                    }

    def execute(self):
        pass

    def _sp(self):
        pass

    def _opt(self):
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
                "dmso": "f = dimethylsulfoxide.cosmo ",
                "methanol": "f = methanol.cosmo ",
                "thf": "f = thf.cosmo ",
                "toluene": "f = toluene_c0.cosmo ",
            }
            if fine:
                solv_data = os.path.join(
                    os.path.split(self.progsettings["cosmorssetup"].split()[5].strip('"'))[0],
                    "DATABASE-COSMO/BP-TZVPD-FINE",
                )
            else:
                solv_data = os.path.join(
                    os.path.split(self.progsettings["cosmorssetup"].split()[5].strip('"'))[0],
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
        removegf=False
        if self.basis == 'def2-QZVP(-gf)':
            self.basis = 'def2-QZVP'
            removegf=True
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
        }
        # remove -fg functions from def2-QZVP basis set
        if removegf:
            cef_calls[self.func] =  cef_calls[self.func] + ["-gf"]
        for k in range(2):
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
                    if self.func == 'b97-3c':
                        testfunc = 'b973c'
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
                    newcontrol.write('$nucsel "h" "c" "f" "p"\n')
                    newcontrol.write('$nucsel2 "h" "c" "f" "p"\n')
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
                    callargs = [self.progsettings["xtbpath"], "coord", "-opt", "lax", "-tm"]
                else:  # gas phase
                    callargs = [self.progsettings["xtbpath"], "coord", "-opt", "-tm"]
            else:  # only crude
                callargs = [self.progsettings["xtbpath"], "coord", "-opt", "crude", "-tm"]

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
            # if optimization was carried out in solution, need new optimization in gas phase
            if opt:
                self._opt(environsettings)
            else:  # if optimization was carried out in gas phase, coord is sufficient without reoptimization
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
                [self.progsettings["tempprogpath"], "-smpcpus", str(self.progsettings["omp"])],
                shell=False,
                stdin=None,
                stderr=subprocess.STDOUT,
                universal_newlines=False,
                cwd=self.workdir,
                stdout=outputfile,
                #env=environsettings,
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
                [self.progsettings["tempprogpath"], "-smpcpus", str(self.progsettings["omp"])],
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
                        "ERROR: shiedling calculation failed in {:18}".format(
                            last_folders(self.workdir, 1)
                        ),
                        file=sys.stderr,
                    )
                    self.success = False
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


class orca_job(qm_job):
    coord = []
    solvent_smd_new = {
        "acetone": ["%cpcm", " smd     true", ' smdsolvent "acetone"', "end"],
        "chcl3": ["%cpcm", " smd     true", ' smdsolvent "chloroform"', "end"],
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
        "ch2cl2": ["!CPCM(CH2Cl2)"],
        "dmso": ["!CPCM(DMSO)"],
        "h2o": ["!CPCM(Water)"],
        "methanol": ["!CPCM(Methanol)"],
        "thf": ["!CPCM(THF)"],
        "toluene": ["!CPCM(Toluene)"],
    }

    def coord2xyz(self):
        """convert TURBOMOLE coord file to xyz"""
        time.sleep(0.1)
        bohr2ang = 0.52917721067
        with open(
            os.path.join(self.workdir, "coord"), "r", encoding=coding, newline=None
        ) as f:
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
                "{} {: 09.7f}  {: 09.7f}  {: 09.7f}".format(atom[i], x[i], y[i], z[i])
            )
        self.coord = coordxyz
        self.nat = len(coordxyz)
        return

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
            self.coord2xyz()
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
            self.coord2xyz()
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
                newcoord.write("   orca bin= {} \n".format(os.path.join(self.progsettings["tempprogpath"], 'orca')))
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
                    callargs = [self.progsettings["xtbpath"], "coord", "--opt", "lax", "--orca"]
                else:  # gas phase
                    callargs = [self.progsettings["xtbpath"], "coord", "--opt", "--orca"]
            else:  # only crude
                callargs = [self.progsettings["xtbpath"], "coord", "--opt", "crude", "--orca"]
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
            self.coord2xyz()
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
            self.coord2xyz()
            # input generation
            ofreq_calls = {
                "b97-3c": [
                    "%MaxCore 8000",
                    "! def2-mTZVP b97-3c grid4 NumFreq",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "pbeh-3c": [
                    "%MaxCore 8000",
                    "! def2-mSVP pbeh-3c grid3 NumFreq",
                    "!     smallprint printgap noloewdin",
                    "%output",
                    "       print[P_BondOrder_M] 1",
                    "       print[P_Mayer] 1",
                    "       print[P_basis] 2",
                    "end",
                ],
                "tpss": [
                    "%MaxCore 8000",
                    "! def2-TZVP(-f) tpss grid4 NumFreq",
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
                    "! def2-mSVP pbeh-3c Opt grid3 NumFreq",
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
        self.coord2xyz()
        # generate input   # double hybrids not implemented  # tpss needs for some reason additionally def2/j
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
        }

        # 'dsd-blyp': ['%MaxCore 8000', '! dsd-blyp pcJ-0 grid5 rijk def2/jk def2-TZVPP/C nofinalgrid nososcf',
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
        self.coord2xyz()
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
        }
        # 'dsd-blyp': ['%MaxCore 8000',
        #          '! dsd-blyp pcSseg-2 grid5 rijcosx def2/j def2-TZVPP/C nofinalgrid nososcf gridx4',
        #          '!     smallprint printgap noloewdin', '%output', '       print[P_BondOrder_M] 1',
        #          '       print[P_Mayer] 1', '       print[P_basis] 2', 'end', '%mp2', '    RI true',
        #          '    density relaxed', 'end'] }

        with open(os.path.join(self.workdir, "inpS"), "w", newline=None) as inp:
            # functional
            for line in nmrs_calls[self.func]:
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
    try:
        tmp = instructdict["jobtype"]
    except KeyError:
        raise ("jobtype is missing in instructdict!")

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

    known_keys = (
        "reference for 1H",
        "reference for 13C",
        "reference for 19F",
        "reference for 31P",
        "1H active",
        "13C active",
        "19F active",
        "31P active",
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
        "resonance frequency": "mf",
        "reference for 1H": "href",
        "reference for 13C": "cref",
        "reference for 31P": "pref",
        "reference for 19F": "fref",
    }
    # später für unittest!
    # print(list(set(known_keys) - set(key_args_dict.keys())))

    enso_internal_defaults = {
        "reference for 1H": "TMS",
        "reference for 13C": "TMS",
        "reference for 19F": "default",
        "reference for 31P": "default",
        "1H active": "on",
        "13C active": "off",
        "19F active": "off",
        "31P active": "off",
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
        "basisJ": "default",
        "funcS": "pbe0",
        "basisS": "default",
        "couplings": "on",
        "shieldings": "on",
        "part1_threshold": "4.0",
        "part2_threshold": "2.0",
        "sm": "default",
        "sm3": "default",
        "sm4": "default",
        "check": "on",
        "crestcheck": "off",
        "maxthreads": "1",
        "omp": "1",
        "smgsolv2": "cosmors",
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

    def __init__(self, solvents, gfnv, func, func3, funcJ, funcS, href, cref):
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
            "prog_rrho": ("xtb", "prog"),
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
            "reference for 19F": ("default"),
            "reference for 31P": ("default"),
            "1H active": ("on", "off"),
            "13C active": ("on", "off"),
            "19F active": ("on", "off"),
            "31P active": ("on", "off"),
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
        patherror_logic = False
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
                outdata.write("ORCA version: 4.1\n")
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
        self, args, solvents, func, func3, funcJ, funcS, gfnv, href, cref, fref, pref
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
                            "WARNING: DCOSMO-RS is used as default solvent model with TM."
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
            # inp.write('reference for 19F: {:{digits}} #{}\n'.format(str(args.fref),', '.join(fref),digits=40-len('reference for 19F')))
            # inp.write('reference for 31P: {:{digits}} #{}\n'.format(str(args.pref),', '.join(pref),digits=40-len('reference for 31P')))
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
                            "\nERROR: Can not convert number of conformers (nconf)! Number of conformers has to be "
                            "an integer or set to all."
                        )
                        error_logical = True
        if "charge" in self.configdata:
            try:
                args.chrg = int(self.configdata["charge"])
            except ValueError:
                print(
                    "\nERROR: Can not read the charge information! Charge has to be an integer."
                )
                error_logical = True
        if "unpaired" in self.configdata:
            try:
                args.unpaired = int(self.configdata["unpaired"])
            except ValueError:
                print(
                    "\nERROR: Can not convert number of unpaired electrons! Number of unpaired electrons has to be "
                    "an integer."
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
                    "\nERROR: Prog is not recognized! Options for the main program are ORCA (orca) or TURBOMOLE (tm)."
                )
                error_logical = True
        if "prog3" in self.configdata:
            if self.configdata["prog3"] in ["tm", "orca"]:
                args.prog3 = self.configdata["prog3"]
            elif self.configdata["prog3"] == "prog":
                args.prog3 = args.prog
            else:
                print(
                    "\nERROR: Prog3 is not recognized! Options for the main program are ORCA (orca) or TURBOMOLE (tm)."
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
                        "\nERROR: Program for part 4 is not recognized! Options are ORCA (orca), TURBOMOLE (tm), "
                        "or the main program (prog)."
                    )
                    error_logical = True
                else:
                    print(
                        "\nWARNING: Program for part 4 is not recognized! Options are ORCA (orca), TURBOMOLE (tm), "
                        "or the main program (prog)."
                    )
        if "ancopt" in self.configdata:
            if self.configdata["ancopt"] in ["on", "off"]:
                args.ancopt = self.configdata["ancopt"]
            else:
                print(
                    "\nWARNING: The use of the ANCOPT-optimizer implemented in GFN-xTB is not recognized! "
                    "Options are on or off."
                )
                error_logical = True
        if "prog_rrho" in self.configdata:
            if self.configdata["prog_rrho"] == "xtb":
                args.rrhoprog = self.configdata["prog_rrho"]
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
                        "WARNING: Currently are only GFN-xTB hessians possible and no TM hessians"
                    )
            else:
                print(
                    "\nERROR: Program for RRHO contribution in part 2 and 3 (rrhoprog) is not recognized! Options are "
                    "GFN-xTB (xtb) or the main program (prog)."
                )
                error_logical = True
        if "gfn_version" in self.configdata:
            if self.configdata["gfn_version"] in gfnv:
                args.gfnv = self.configdata["gfn_version"]
            else:
                print(
                    "\nWARNING: GFN-xTB version is not recognized! Options are {}.\nGFN2-xTB is used as "
                    "default.".format(", ".join(gfnv))
                )
                args.gfnv = "gfn2"
        if "temperature" in self.configdata:
            try:
                tmp = float(self.configdata["temperature"])
                args.temperature = self.configdata["temperature"]
            except ValueError:
                print(
                    '\nWARNING: Temperature "{}" could not be converted to a float! Using 298.15 K instead as default!'.format(
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
                    "\nERROR: Boltzmann is not recognized! Options are off for calculating everything or on for "
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
                    "\nERROR: Backup is not recognized! Options are off for the normal procedure or on for "
                    "reevaluating the conformers sorted out."
                    ""
                )
                error_logical = True
        if "func" in self.configdata:
            if self.configdata["func"] in func:
                args.func = self.configdata["func"]
            else:
                print(
                    "\nERROR: Chosen functional for part 1 and 2 (func) is not implemented! Options are {}.".format(
                        ", ".join(func)
                    )
                )
                error_logical = True
        if "func3" in self.configdata:
            if self.configdata["func3"] in func3:
                args.func3 = self.configdata["func3"]
                if args.prog3 == "tm" and args.func3 != "pw6b95":
                    args.func3 = "pw6b95"
                    print(
                        "\nWARNING: Only PW6B95 is implemented as functional for part 3 for TM!"
                    )
            else:
                print(
                    "\nERROR: Chosen functional for part 3 (func3) is not implemented! Options are {}.".format(
                        ", ".join(func3)
                    )
                )
                error_logical = True
        if "basis3" in self.configdata:
            try:
                args.basis3 = str(self.configdata["basis3"])
                if args.basis3 not in self.knownbasissets3:
                    print(
                        "WARNING! Basis for part3: {} is used but could not be checked for correct input!".format(
                            args.basis3
                        )
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
        if "funcJ" in self.configdata:
            if self.configdata["funcJ"] in funcJ:
                args.funcJ = self.configdata["funcJ"]
            else:
                print(
                    "\nERROR: Chosen functional for coupling calculation (funcJ) is not implemented! Options are {}.".format(
                        ", ".join(funcJ)
                    )
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
                if args.prog4 == "tm":
                    args.basisJ = "def2-TZVP"
                elif args.prog4 == "orca":
                    args.basisJ = "pcJ-0"
            if args.basisJ not in self.knownbasissetsJ:
                print(
                    "WARNING! Basis for coupling calculation: {} is used but could not be checked for correct "
                    "input!".format(args.basisJ)
                )
        if "shieldings" in self.configdata:
            if self.configdata["shieldings"] in ["on", "off"]:
                args.calcS = self.configdata["shieldings"]
            else:
                print(
                    '\nERROR: Shieldings is not recognized! Options are "on" or "off".'
                )
                error_logical = True
        if "funcS" in self.configdata:
            if self.configdata["funcS"] in funcS:
                args.funcS = self.configdata["funcS"]
            else:
                print(
                    "\nERROR: Chosen functional for shielding calculation (funcS) is not implemented! Options are {}.".format(
                        ", ".join(funcS)
                    )
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
                if args.prog4 == "tm":
                    args.basisS = "def2-TZVP"
                elif args.prog4 == "orca":
                    args.basisS = "pcSseg-2"
            if args.basisS not in self.knownbasissetsS:
                print(
                    "WARNING! Basis for shielding calculation: {} is used but could not be checked for correct "
                    "input!".format(args.basisS)
                )
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
            elif args.prog == "orca":
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
            elif args.prog == "tm":
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
                args.sm4 = "cosmo"
            elif args.prog4 == "orca":
                if self.configdata["sm4"] == "smd":
                    args.sm4 = "smd"
                else:
                    args.sm4 = "cpcm"
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
        # here 19 F and 31 P
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
        if "resonance frequency" in self.configdata:
            try:
                args.mf = int(self.configdata["resonance frequency"])
            except ValueError:
                args.mf = None
                print(
                    "\nWARNING: Can not convert resonance frequency! Flag -mf has to be set while executing the "
                    "anmr program."
                )
        if error_logical:
            print("Going to exit due to input errors!")
            sys.exit(1)
        return spectrumlist


def check_for_folder(conflist, functional):
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
    if error_logical:
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
            "implemented! For a RRHO contribution in solution, please use GFN-xTB."
            "\nRRHO is calculated in the gas phase.\n"
        )
    elif args.solv and args.rrhoprog == "tm":
        print(
            "WARNING: TM RRHO calculation with solvation correction is currently not "
            "implemented! For a RRHO contribution in solution, please use GFN-xTB."
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
            #os.environ["OMP_NUM_THREADS"] = "{:d}".format(args.omp)
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
            "temperature" : args.temperature,
            "progsettings": {"omp": args.omp,
                             "tempprogpath": "",
                            },
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
            instructrrho["progsettings"]["tempprogpath"] = ''

        results = run_in_parallel(
            q, resultq, job, maxthreads, results, instructrrho,"rrho"
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
                "progsettings": {"omp": args.omp,
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
                "progsettings": {"omp": args.omp,
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
        "ch2cl2",
        "dmso",
        "h2o",
        "methanol",
        "thf",
        "toluene",
        "gas",
    )
    impfunc = ("pbeh-3c", "b97-3c", "tpss")  # pbe0
    impfunc3 = ("pw6b95", "wb97x", "dsd-blyp")  # pbe0
    impfuncJ = ("tpss", "pbe0")
    impfuncS = ("tpss", "pbe0")
    impgfnv = ("gfn1", "gfn2")
    imphref = ("TMS", "DSS")
    impcref = ("TMS", "DSS")
    impfref = ("CFCl3",)
    imppref = ("TMP/PH3",)
    orca_old = False

    descr = """
     __________________________________________________
    |                                                  |
    |                                                  |
    |                       ENSO -                     |
    |         energetic sorting of CREST CRE           |
    |          for automated NMR calculations          |
    |             University of Bonn, MCTC             |
    |                    July 2018                     |
    |                   version 1.272                  |
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
    hactive = "off"
    cactive = "off"
    factive = "off"
    pactive = "off"
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
            if args.solv is None or args.solv == "toluene":
                args.pref = "PH3"
            else:
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


def part1(args, jsonfile, json_dict, conformersxyz, nat, maxthreads, xtbpath, environsettings):
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
        check_for_folder(["CONF" + str(i) for i in old_confs_list], args.func)

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
                "progsettings" :{
                                "tempprogpath": "", 
                                "xtbpath": xtbpath,
                                "orca_old" : orca_old,
                                "omp": args.omp,
                                }
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
                "progsettings": {"omp": args.omp,
                                 "tempprogpath": "",
                                },
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
        # get GFN-xTB energies
        gfne = conformersxyz2coord(conformersxyz, nat, args.func, all_confs_list)
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


def part2(args, results, jsonfile, cwd, json_dict, maxthreads, xtbpath, environsettings):
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

        check_for_folder([i.name for i in tmp_results], args.func)

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
                    "progsettings": {"omp": args.omp,
                                     "tempprogpath": "",
                                    },
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
                "progsettings": {"omp": args.omp,
                                 "tempprogpath": "",
                                },
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
                    "progsettings": {"omp": args.omp,
                                     "tempprogpath": "",
                                    },
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
                    "progsettings": {"omp": args.omp,
                                     "tempprogpath": "",
                                    },
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
        lowestfree = min([item.free_energy for item in results if item.free_energy is not None])
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


def part3(args, results, jsonfile, cwd, json_dict, maxthreads, xtbpath, environsettings):
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
            #os.environ["OMP_NUM_THREADS"] = "{:d}".format(args.omp)
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
                check_for_folder([i.name for i in tmp_results], args.func)

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
                        "progsettings": {"omp": args.omp,
                                         "tempprogpath": "",
                                        },
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
                        "progsettings": {"omp": args.omp,
                                         "tempprogpath": "",
                                        },
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

                #
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
            lowestfree = min([item.free_energy for item in results if item.free_energy is not None])
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
            kcalmol2J = 7.993840458653155e-21  # kcal/mol to J  =   *4.184 / N_{A}
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
    mpshiftpath):
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
                        "progsettings": {"omp": args.omp,
                                         "tempprogpath": "",
                                        },
                    }
                    if args.hactive == "on":
                        instructprep["hactive"] = True
                    if args.cactive == "on":
                        instructprep["cactive"] = True
                    if args.factive == "on":
                        instructprep["factive"] = True
                    if args.pactive == "on":
                        instructprep["pactive"] = True

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
                        "progsettings": {"omp": args.omp,
                                         "tempprogpath": "",
                                        },
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
                    "progsettings": {"omp": args.omp,
                                    },
                }
                if args.hactive == "on":
                    instructJ["hactive"] = True
                if args.cactive == "on":
                    instructJ["cactive"] = True
                if args.factive == "on":
                    instructJ["factive"] = True
                if args.pactive == "on":
                    instructJ["pactive"] = True
                if job ==tm_job:
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
                            "progsettings": {"omp": args.omp,
                                             "tempprogpath": "",
                                            },
                        }
                        if args.hactive == "on":
                            instructprep["hactive"] = True
                        if args.cactive == "on":
                            instructprep["cactive"] = True
                        if args.factive == "on":
                            instructprep["factive"] = True
                        if args.pactive == "on":
                            instructprep["pactive"] = True

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
                            "progsettings": {"omp": args.omp,
                                             "tempprogpath": "",
                                            },

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
                    "progsettings": {"omp": args.omp,
                                     "tempprogpath": "",
                                    },
                }
                if args.hactive == "on":
                    instructS["hactive"] = True
                if args.cactive == "on":
                    instructS["cactive"] = True
                if args.factive == "on":
                    instructS["factive"] = True
                if args.pactive == "on":
                    instructS["pactive"] = True
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
            args.hactive,
            args.cactive,
            args.factive,
            args.pactive,
            args.calcJ,
            args.calcS,
        )
        write_json("save", json_dict, jsonfile)

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
            args, jsonfile, json_dict, conformersxyz, nat, maxthreads, xtbpath, environsettings
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
            args, results, jsonfile, cwd, json_dict, maxthreads, xtbpath, environsettings
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
            args, results, jsonfile, cwd, json_dict, maxthreads, xtbpath, environsettings
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
            mpshiftpath
        )
    except Exception as e:
        print("ERROR in part4!")
        print("The error-message is {}\n\n".format(e))
        traceback.print_exc()
        print("Going to exit!")
        sys.exit(1)

    # ALL DONE
    print("\nENSO all done!")
