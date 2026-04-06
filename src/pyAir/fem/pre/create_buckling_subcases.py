# -*- coding: utf-8 -*-
"""
Created on Mon May 26 11:26:23 2025

@author: U69432
"""

from Common.UI.filedialog import askOpenBDF
__version__ = '1.1'

import os

# NOTES:
# -----------------------------------------------
# Nastran has some rules related to the SUBASES
# - Must be an integer between 1 and 9999999
# - The must be ordered


def create_sol105_subcases(
        method_sid: int = 10,
        exlude_sid: int = None,
        offset: int = 1000000,
        verbose: bool = True
):

    # Asking for SUBCASE BDF file
    bdf_path = askOpenBDF(verbose)
    dirname = os.path.dirname(bdf_path)
    bdf_fn = os.path.basename(bdf_path)

    # Creating the BUCKLING SUBCASES BDF FILE
    out_fn = bdf_fn.split('.')[0] + '_SOL105.' + bdf_fn.split('.')[1]
    out_path = os.path.join(dirname, out_fn)

    # Reading the content of the BDF file
    with open(bdf_path, 'r') as input_file:
        input_lines = input_file.readlines()

    # Retrieving the static SID
    static_sids = []
    for line in input_lines:
        if line.upper().startswith('SUBCASE'):
            static_sid = line.upper().split('SUBCASE')[1].strip()
            static_sids.append(static_sid)

    # Autoselecting the Labeling rule for the buckling subcases
    max_sid = static_sids[-1]
    digits_max_sid = len(max_sid)
    if digits_max_sid <= 6:
        label_rule = 'Adding digit 1 at the leftside'
        buck_sids = ['1' + sid for sid in static_sids]
    else:
        if int(max_sid) + offset <= 9999999:
            label_rule = 'Offsting the SID by 1000000'
            buck_sids = [int(sid)+1000000 for sid in static_sids]
            for buck_sid in buck_sids:
                if buck_sid in static_sids:
                    buck_sids = None
                else:
                    pass
        else:
            print('Impossible to generate a labeling rule for the buckling subacases')
            print('The static SID have either more than 6 digitis or if the its is'
                  'increased the first digit of the maximum SID, it exceeds the '
                  ' 9999999 limit')
            buck_sids = None

    # Generating the file with the SOL105 subcases
    if buck_sids is None:
        pass
    else:
        static_buckling_sid_map = dict(zip(static_sids, buck_sids))
        output_lines = []
        with open(out_path, 'w') as outputfile:
            for line in input_lines:
                if line.startswith('$'):
                    output_lines += [line]
                else:
                    if line.upper().startswith('SUBCASE'):
                        static_sid = line.upper().split('SUBCASE')[1].strip()
                        buck_sid = static_buckling_sid_map[static_sid]
                        output_lines += [
                            'SUBCASE{0}{1}\n'.format(5*' ', buck_sid)]
                    if 'TITLE' in line.upper():
                        whitespaces = line.upper().split('TITLE')[0]
                        output_lines += [line]
                        output_lines += [
                            '{3}STATSUB{0}={1}{2}\n'.format(
                                ' ', ' ', static_sid, whitespaces),
                            '{3}METHOD{0}={1}{2}\n'.format(
                                ' ', ' ', method_sid, whitespaces),
                            '{2}VECTOR(PLOT){0}={1}ALL\n'.format(
                                ' ', ' ', whitespaces)
                        ]
                        if not (exlude_sid is None):
                            output_lines += [
                                '{1}PARAM,EXCLUDE = {0}\n'.format(
                                    exlude_sid, whitespaces)
                            ]
                    if 'SPC' in line.upper():
                        output_lines += [line]
                    if 'MPC' in line.upper():
                        output_lines += [line]

            outputfile.writelines(output_lines)
        if verbose:
            print('SOL105 SUBCASES file generated at:')
            print('---------------------------------')
            print(out_path)
            print('---------------------------------')
            print('Labelling rule : ', label_rule)
            print('Minimum Buckling Subcase ID : ', buck_sids[0])
            print('Maximum Buckling Subcase ID : ', buck_sids[-1])
            print('---------------------------------')
            print('Please include in the launcher file the following card:')
            print('EIGRL   {0:>8}{1:>8}{2:>8}{3:>8}'.format(
                method_sid, 1, '', 15))


create_sol105_subcases(
    method_sid=2,
    # exlude_sid=-1050
)
