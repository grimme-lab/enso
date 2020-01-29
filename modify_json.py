#!/usr/bin/env python3

coding = 'ISO-8859-1'

descr='''
     __________________________________________________
    |                                                  |
    |                      EJMT                        |
    |            enso.json modification tool           |
    |             University of Bonn, MCTC             |
    |                   July 2019                      |
    |                     v 1.01                       |
    |                  K. Schmitz,F.Bohle              |
    |__________________________________________________|
    '''

structure_input='''
structure:                     
    CONF1                      
    1H_S: calculated           
    energy_cosmors: -1.23      
    CONF3                      
    energy_sp_part3: -123.456  
'''

def cml(descr,structure_input):
    """ Get args object from commandline interface.
        Needs argparse module."""

    #parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter,usage=argparse.SUPPRESS,description=useit) #argparse.RawDescriptionHelpFormatter) #,
    parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter,usage=argparse.SUPPRESS)
    parser.add_argument('-i', "--input_file",
                        dest='inputfile',
                        action='store',
                        required=False,
                        type=str,
                        help='''Provide input file containing the modifications. If only conformers should be added, no input file is required. Defaulte: None.
                        {} '''.format(structure_input))
    parser.add_argument('-j', "--json_file",
                        dest='jsonfile',
                        action='store',
                        required=False,
                        default='enso.json',
                        type=str,
                        help='Provide json_file that should be modified. Default: enso.json')
    parser.add_argument('-o', "--output_file",
                        dest='outputfile',
                        action='store',
                        required=False,
                        default='enso.json',
                        type=str,
                        help='Provide output file (new json_file). If the output file and the json_file are the same, the json_file is '
                             'moved recursively. Default enso.json.')
    parser.add_argument('-a', "--add",
                        dest='add',
                        action='store',
                        required=False,
                        nargs='+',
                        help='Add new conformer(s). Requires only the number of each conformer that should be added. Default: None.')

    args = parser.parse_args()
    return args

def decomment(csvfile):
    ''' remove any comments from file before parsing with csv.DictReader'''
    for row in csvfile:
        raw = row.split('#')[0].strip()
        if raw: yield raw
    return

def read_json(cwd, jsonfile): # jsonfile
    '''Reading in jsonfile'''
    # check if jsonfile exists
    if not os.path.isfile(os.path.join(cwd, jsonfile)):
        print('ERROR: file {} does not exist!\nGoing to exit.'.format(jsonfile))
        sys.exit(1)
    json_dict = {}
    print('Reading in {}.'.format(jsonfile))
    try:
        with open(os.path.join(cwd, jsonfile), "r", encoding=coding, newline=None) as inp:
            json_dict = json.load(inp)
    except ValueError as e:
        print('Your Jsonfile is corrupted!\n{}'.format(e))
        sys.exit(1)
    return json_dict

def read_input(path, filename):
    changes = {}
    tmpname = None
    tmpkeys = None
    with open(os.path.join(path, filename), 'r') as inp:
        reader = csv.DictReader(decomment(inp), fieldnames=('key', 'value'), skipinitialspace=True,
                                    delimiter=':')
        for row in reader:
            if 'CONF' in row['key']:
                if tmpkeys is not None:
                    # dump conformer to changes
                    if tmpname in changes.keys():
                        print('ERROR: {} is listed more than one time in {}! Please write all changes concerning one conformer together. \nGoing to exit.'.format(tmpname, filename))
                        sys.exit(1)
                    changes[tmpname] = tmpkeys
                tmpname = row['key']
                tmpkeys = {}
            else:
                tmpkeys[row['key']] = row['value']
        if tmpkeys is None:
            print('ERROR: {} is empty! \nGoing to exit.'.format(filename))
            sys.exit()
        changes[tmpname] = tmpkeys
        # groß und kleinschreibung nervt sonst
        for conformer in changes:
            for key, value in changes[conformer].items():
                if value == 'false':
                    changes[conformer][key] = 'False'
                elif value == 'true':
                    changes[conformer][key] = 'True'
    return changes

def check_changes(changes_dict, inpfile):
    '''Check input'''
    print('Checking input.')
    error_logical = False
    partkey_list = ['crude_opt', 'opt', 'sp_part2', 'cosmo-rs', 'rrho', 'sp_part3', '1H_J', '13C_J', '19F_J',
                        '31P_J', '1H_S', '13C_S', '19F_S', '31P_S']
    partkeys = ['calculated', 'failed', 'not_calculated']
    numberkey_list = ['energy_crude_opt', 'energy_opt', 'energy_sp_part2', 'energy_sp_part3', 'energy_cosmo-rs',
                          'energy_rrho']
    true_false_keylist = ['consider_for_part2', 'consider_for_part3', 'consider_for_part4', 'backup_for_part2',
                              'backup_for_part3', 'removed_by_user']
    other_keylist = ['symmetry']
    for item in changes_dict.keys():
        if not 'CONF' in item:
            print('ERROR: {} is not recognized!\nGoing to exit.'.format(item))
            error_logical = True
        for key in changes_dict[item]:
            if key in partkey_list:
                if changes_dict[item][key] not in partkeys:
                    print('ERROR: information about {} for {} in the file {} are not recognized!'.format(key, item, inpfile))
                    error_logical = True
            elif key in numberkey_list:
                if changes_dict[item][key] == 'None':
                    changes_dict[item][key] = None
                else:
                    try:
                        changes_dict[item][key] = float(changes_dict[item][key] )
                    except:
                        print('ERROR: {} of {} has to be None or a number!'.format(key, item))
                        error_logical = True
            elif key in true_false_keylist:
                if changes_dict[item][key] == 'True':
                    changes_dict[item][key] =  True
                elif changes_dict[item][key] == 'False':
                    changes_dict[item][key] =  False
                else:
                    print('ERROR: {} of {} has to be True or False!'.format(key, item))
                    error_logical = True
            elif key in other_keylist:
                pass
            else:
                print('ERROR: {} of {} is not recognized!'.format(key, item))
                error_logical = True
    if error_logical == True:
        print('One or multiple errors were found in the file {}. The following values and keys are possible:'.format(inpfile))
        print('"calculated", "not_calculated" or "failed" are options for:')
        for i in partkey_list:
            print("    {}".format(i))
        print('None or a number are options for:')
        for i in numberkey_list:
            print("    {}".format(i))
        print('True or False are options for:')
        for i in true_false_keylist:
            print("    {}".format(i))
        print('There are no specifications for:')
        for i in other_keylist:
            print("    {}".format(i))
        print('Please check the flag "-h" for further information regarding the input.\nGoing to exit.')
        sys.exit(1)
    return

def modify_json(changes_dict, json_dict, new_conf):
    '''Apply changes to json_dict'''
    print('Applying changes.')

    energy_dict = {'energy_opt': 'opt',
                  'energy_crude_opt': 'crude_opt',
                  'energy_rrho': 'rrho',
                  'energy_cosmo-rs': 'cosmo-rs',
                  'energy_sp_part2': 'sp_part2',
                  'energy_sp_part3': 'sp_part3'}

    match_dict = {'consider_for_part2': ['crude_opt'],
                  'backup_for_part2': ['crude_opt'],
                  'consider_for_part3': ['opt', 'rrho'],
                  'backup_for_part3': ['opt', 'rrho'],
                  'consider_for_part4': ['sp_part3', 'rrho']}

    for item in changes_dict.keys():
        if item not in json_dict.keys():
            # create new conformer
            print('Create new conformer {}.'.format(item))
            json_dict[item] = new_conf.copy()
        for key in changes_dict[item]:
            try:
                json_dict[item][key] = changes_dict[item][key]
            except:
                print('ERROR: Could not change {} of {}!'.format(key, item))
                break
            # match keys if necessary
            if 'energy' in key:
                tmpkey = key.replace('energy_', '')
                if key == None and tmpkey not in changes_dict[item]:
                    json_dict[item][tmpkey] = 'not_calculated'
                elif key != None and tmpkey not in changes_dict[item]:
                    json_dict[item][tmpkey] = 'calculated'
        # check if status of calculation have to be changed
        for energy in energy_dict.keys():
            if energy in changes_dict[item] and energy_dict[energy] not in changes_dict[item]:
                if changes_dict[item][energy] != None:
                    json_dict[item][energy_dict[energy]] = "calculated"
                else:
                    json_dict[item][energy_dict[energy]] = "not_calculated"
        # check if considered and/or backup have to be reset to defaults
        for match_key in match_dict:
            # only change keys automatically if they are not in change_dict
            if match_key not in changes_dict[item]:
                for match_item in match_dict[match_key]:
                    if json_dict[item][match_item] in ['failed' , 'not_calculated']:
                        json_dict[item][match_key] = False
                    elif match_key == 'consider_for_part2' and json_dict[item][match_item] == 'not_calculated':
                        json_dict[item][match_key] = True

def add_conformers(changes_dict, json_dict, new_conf):
    '''Add conformers given as flag'''
    print('Add conformers passed as flag.')
    for item in args.add:
        tmpitem = 'CONF' + str(item)
        if tmpitem in json_dict.keys():
            if tmpitem in changes_dict.keys():
                pass
            else:
                print('WARNING: {} is already in the original json file. It is not added.'.format(tmpitem))
        else:
            print('Create new conformer {}.'.format(item))
            json_dict[tmpitem] = new_conf.copy()
    return json_dict

def splitting(item):
    '''Used in move recursively'''
    try:
        return int(item.rsplit('.', 1)[1])
    except ValueError:
        return 0

def move_recursively(path, filename):
    '''Check if file or file.x exists and move them to file.x+1
       ignores e.g. file.save'''
    files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path,f))] # list of all files in directory
    newfiles = [] # list of all files in directory that contain filename and '.'
    for item in files:
        if filename+'.' in item:
            newfiles.append(item)
    newfiles.sort(key=splitting, reverse=True)
    for item in newfiles:
        try:
            data = item.rsplit('.', 1) # splits only at last '.'
            int(data[1])
        except ValueError:
            continue
        tmp_from = os.path.join(path, item)
        newfilename = str(data[0]) + "." + str(int(data[1]) + 1)
        tmp_to   = os.path.join(path, newfilename)

        print('Backing up {} to {}.'.format(item, newfilename))
        shutil.move(tmp_from, tmp_to)

    if filename in files:
        print('Backing up {} to {}.'.format(filename, filename+'.1'))
        shutil.move(os.path.join(path, filename), os.path.join(path, filename + '.1'))
    return

def write_json(json_dict, path, jsonfile):
    '''Writing jsonfile'''
    outfile = jsonfile
    print('Modifications are written to {}.'.format(outfile))
    with open(os.path.join(path, outfile), 'w') as out:
        json.dump(json_dict, out, indent=4, sort_keys=False)
    return

if __name__ == "__main__":
    try:
        import argparse
    except ImportError:
        raise ImportError('EJMT requires the module argparse. Please install the module argparse.')
    try:
        import sys
    except ImportError:
        raise ImportError('EJMT requires the module sys. Please install the module sys.')
    try:
        import os
    except ImportError:
        raise ImportError('EJMT requires the module os. Please install the module sys.')
    from os.path import expanduser
    try:
        import json
    except ImportError:
        raise ImportError('EJMT requires the module json. Please install the module json.')
    try:
        import shutil
    except ImportError:
        raise ImportError('EJMT requires the module shutil. Please install the module shutil.')
    try:
        import csv
    except ImportError:
        raise ImportError('EJMT requires the module csv. Please install the module csv.')

    print(descr)  ### Program description
    args = cml(descr, structure_input)

    cwd = os.getcwd()

    #template for new conformer
    new_conf = {'crude_opt': 'not_calculated', 'energy_crude_opt': None, 'backup_for_part2': False,
            'consider_for_part2': True,
            'opt': 'not_calculated', 'energy_opt': None, 'backup_for_part3': False,
            'sp_part2': 'not_calculated', 'energy_sp_part2': None, 'consider_for_part3': False,
            'sp_part3': 'not_calculated', 'energy_sp_part3': None,
            'cosmo-rs': 'not_calculated', 'energy_cosmo-rs': None,
            'rrho': 'not_calculated', 'energy_rrho': None, 'symmetry': 'C1',
            'consider_for_part4': False,
            '1H_J': 'not_calculated', '1H_S': 'not_calculated', '13C_J': 'not_calculated',
            '13C_S': 'not_calculated', '19F_J': 'not_calculated', '19F_S': 'not_calculated',
            '31P_J': 'not_calculated', '31P_S': 'not_calculated', 'removed_by_user': False}

    # read json file
    json_dict = read_json(cwd, args.jsonfile)

    if args.inputfile:
        # read input
        changes_dict = read_input(cwd, args.inputfile)

        # check input
        check_changes(changes_dict, args.inputfile)

        # modify json file
        modify_json(changes_dict, json_dict, new_conf)

    # add conformers
    if args.add and args.inputfile:
        add_conformers(changes_dict, json_dict, new_conf)
    elif args.add and not args.inputfile:
        changes_dict = {}
        add_conformers(changes_dict, json_dict, new_conf)

    # move enso.json recursively
    if args.outputfile == args.jsonfile:
        move_recursively(cwd, args.jsonfile)

    # dump json_file
    write_json(json_dict, cwd, args.outputfile)

print('\nAll done!')

#CONF1
#energy_crudeopt: -123.231
#energy_cosmors: -0.2341
#CONF2