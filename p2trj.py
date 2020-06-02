#!/usr/bin/env python3

""" generate structure-ensemble from  ouputfiles and conformer folders"""



from sys import argv,exit
import os

print('Usage: ./p2trj.py part1_energies.dat pbeh-3c 0.8')

def coord2xyz(pathin, commentary, pathout):
    """convert TURBOMOLE coord file to xyz"""
    bohr2ang = 0.52917721067
    with open(os.path.join(pathin, "coord"), "r", newline=None) as f:
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
    nat = int(len(x))
    coordxyz = []
    for i in range(len(x)):
        coordxyz.append(
            "{:3} {: 19.10f}  {: 19.10f}  {: 19.10f}".format(
                atom[i][0].upper() + atom[i][1:], x[i], y[i], z[i]
            )
        )
    with open(pathout, "a", newline=None) as out:
        out.write("  {}\n".format(nat))  ### number of atoms
        out.write("{}\n".format(str(commentary)))
        for line in coordxyz:
            out.write(line + "\n")

if argv[1] in ('part1_energies.dat', 'part2_free_energies.dat') and os.path.isfile(argv[1]):
    with open(argv[1], 'r') as inp:
        data = inp.readlines()
else:
    print('Either need filename or file not found!')
    exit()
if argv[1] == 'part1_energies.dat':
    column = 4
elif argv[1] == 'part2_free_energies.dat':
    column = 5
for line in data[2:]:
    tmp = line.split()
    if float(tmp[column]) < float(argv[3]):
        print(tmp[0])
        coord2xyz(os.path.join(tmp[0], argv[2]), str(tmp[column])+ '  '+ tmp[0], 'all.xyz')

