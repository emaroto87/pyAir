# -*- coding: utf-8 -*-
"""
Created on Mon May 26 10:02:21 2025

@author: U69432
"""

from typing import Optional
from typing import Dict
from typing import List
from typing import Literal

from tqdm import tqdm

import numpy as np
import pandas as pd

from datetime import datetime

from common.ui.filedialogs import askOpenOP2
from common.ui.filedialogs import askSaveAsExcel

from common.ui.tui import TUI as ui_menu
from common.ui.tui import ask_entities_ids_list
from common.io.writers.excel_writer import write_dataframes_to_excel as to_excel
from common.general import load_op2_file
from common.general import get_result_from_op2
from common.general import filter_df_by_entity
from common.general import split_df_by_subcases


from pyNastran.op2.op2_geom import OP2Geom

from collections import namedtuple
from collections import defaultdict

from pathlib import Path

__version__ = '1.0.3'
__2Delem_labels = ['CQUAD4', 'CTRIA3']


def get_joint_info(op2_filename, verbose=False):
    """
    Return from a Nastran OP2 file a pandas DataFrame with the following info
    about the CBUSHES of the file. It only works if the CBUSH connections is 
    node by coincident nodes. In other words, no RBE3, RBE2 o intermediate
    elements between the CBUSH and the 2D element representing the plate.

    Data:
        - CBUSH ID
        - CBUSH NODE 
        - PLATE ELEMENT ID 
        - PLATE TYPE
        - PLATE PID
        - PLATE MID
        - PLATE MAT E MODULUS
        - PLATE THICKNESS
    """
    cbush_joint_labels = [
        'cbush_id',
        'cbush_node',
        'plate_eid',
        'plate_type',
        'plate_pid',
        'plate_mid',
        'plate_mat_E',
        'thickness'
    ]

    cbush_joint = namedtuple('cbush_joint', cbush_joint_labels)
    joint_dict = defaultdict(int)

    # Reading the OP2 file
    op2 = OP2Geom(make_geom=True, debug=False)
    op2.read_op2(op2_filename)

    # Retrieving all the CBUSH from OP2
    elems = op2.elements
    cbushs = [(eid, elem)
              for (eid, elem) in op2.elements.items() if elem.type == 'CBUSH']

    i = 0
    for eid, cbush in tqdm(cbushs, desc='Mapping CBUSHes...'):

        nid1, nid2 = cbush.node_ids
        if verbose:
            print(f"\nCBUSH {eid} conecta nodos {nid1} y {nid2}", end='r')

        for nid in [nid1, nid2]:

            # Getting the elements connected to the CBUSH node
            elems_c2_nid = [elem for eid,
                            elem in elems.items() if nid in elem.node_ids]

            j = None
            # Getting the 2D elements connected to the CBUSH node
            for elem in elems_c2_nid:

                if elem.type in __2Delem_labels:
                    prop = op2.properties[elem.pid]
                    mat = op2.materials[prop.mid1]
                    if verbose:
                        print(
                            f"{elem.type} {elem.eid}: PID={elem.pid},"
                            f" espesor={prop.t},"
                            f" material MID={prop.mid1}"
                        )
                    j = cbush_joint(
                        eid,
                        nid,
                        elem.eid,
                        elem.type,
                        elem.pid,
                        prop.mid1,
                        mat.E(),
                        prop.t
                    )
                else:
                    pass

            if not (j is None):
                joint_dict[i] = j
                i += 1

    # Storing the info into pandas DataFrame
    df = pd.DataFrame(joint_dict).transpose()
    df.columns = cbush_joint_labels
    df.reset_index(inplace=True, drop=True)
    df.set_index('cbush_id', inplace=True)

    return df


def get_cbush_forces(op2_file_path: 'str', shear_comp='YZ', verbose=False):
    '''
    Retrieves all the loads for all the cbushes within a OP2 file

    Parameters
    ----------
    op2_file_path : 'path-like' string
        Path of the NASTRAN OP2 file
    shear_comp : 'str', optional
        Directions in order to compute the shear and axial forces of the cbush 
        from the force components. One of three options are availabe:
            - 'XY', 'XZ' or 'YZ'

        The default is 'YZ'.

    Returns
    -------
    cbush_force : pandas DataFrame
        DataFrame of all the information of the cbushes in the OP2 file

    '''

    shear_opts = {'XY': ('fx', 'fy', 'fz'),
                  'XZ': ('fx', 'fz', 'fy'),
                  'YZ': ('fy', 'fz', 'fx'), }

    try:
        op2 = load_op2_file(op2_file_path)
    except:
        raise (IOError)

    if not (shear_comp in shear_opts.keys()):
        raise ('[Error] Invalid shear component string')
        return None

    else:
        # try:
        # print('Reading CBUSH forces...',end='')
        cbush_force = get_result_from_op2(op2, 'cbush_force')
        # print('Done')
        # except:
        #     print('Unable to extract CBUSH forces. Please check OP2 file')

        def shear_lambda(row, comp1, comp2):
            return ((row[comp1])**2+(row[comp2])**2)**(0.5)

        cbush_force['shear'] = shear_lambda(
            cbush_force,
            shear_opts['YZ'][0],
            shear_opts['YZ'][1])
        cbush_force['axial'] = cbush_force[shear_opts['YZ'][2]]

        max_shear_cbush = cbush_force[cbush_force.shear ==
                                      cbush_force.shear.max()]
        max_axial_cbush = cbush_force[cbush_force.axial ==
                                      cbush_force.axial.max()]
        if verbose:
            print('Max Shear : {0:.2f} at element {1} and subcase {2}'.format(
                max_shear_cbush.shear.values[0],
                max_shear_cbush.index[0],
                max_shear_cbush.SubcaseID.values[0]
            ))
        if verbose:
            print('Max Axial : {0:.2f} at element {1} and subcase {2}'.format(
                max_axial_cbush.axial.values[0],
                max_axial_cbush.index[0],
                max_axial_cbush.SubcaseID.values[0]
            ))

        return cbush_force


if __name__ == '__main__':
    # CREATING INTERACTIVE MENU
    # =========================================================================
    shear_comp_menu = ui_menu(
        options=['XY', 'XZ', 'YZ'],
        label='shear_comp_sel',
        is_main=False
    )

    subcases_selection = ui_menu(
        options=['Yes', 'No'],
        label='subcase_filter',
        is_main=False
    )

    add_joint_map = ui_menu(
        options=['Yes', 'No'],
        label='joint_map_sel',
        is_main=False
    )

    # REQUESTING INPUTS
    # =========================================================================
    # Asking the user the components of the shear force
    print('Select the directions that composes the CBUSH shear.')
    shear_comp_menu.loop()

    # Asking the user if wants filter the info by Subcase
    print(' Filter by Subcase ID list (Yes/No):')
    subcases_selection.loop()

    # Asking the user if wants to mapp all the joints info
    print(' Include a Mapping of the joints (Only for Node-Dependent Meshes) (Yes/No):')
    add_joint_map.loop()

    op2_path = askOpenOP2(multiple=False)
    excel_path = askSaveAsExcel()[0]

    # RUNNING SCRIPT
    # =========================================================================
    exit_cond1 = shear_comp_menu.selection.upper() == 'BACK'
    exit_cond2 = add_joint_map.selection.upper() == 'BACK'

    if not ((op2_path is None) or exit_cond1 or exit_cond2):

        # Joint Mapping
        if add_joint_map.selection.upper() == 'YES':
            cbushes_joint_table = get_joint_info(op2_path)

        # Joint Max Shear and Axial Force
        cbushes_forces_table = get_cbush_forces(
            op2_file_path=op2_path,
            shear_comp=shear_comp_menu.selection
        )

        if subcases_selection.selection.upper() == 'YES':
            subcases = ask_entities_ids_list()
            cbushes_forces_table = filter_df_by_entity(
                cbushes_forces_table,
                subcases,
                'SubcaseID',
                verbose=True)

        # Computing the maximun shear and axial force per cbush
        # Note : In terms of axial load, it gets the maximum absolute
        cbushes_max_shear = cbushes_forces_table.groupby(
            ['ElementID'])[['shear']].max()
        cbushes_max_axial = cbushes_forces_table.groupby(
            ['ElementID'])[['axial']].agg(
                lambda x: x.max() if abs(x.max()) > abs(x.min()) else x.min()
        )

        # Generating the table with MaxShear and MaxAxial
        # ( Shear and Axial LC correlated)
        table_cols = cbushes_forces_table.columns

        max_shear_cbushes_table = cbushes_max_shear.merge(
            cbushes_forces_table,
            on=['ElementID', 'shear']
        )
        max_axial_cbushes_table = cbushes_max_axial.merge(
            cbushes_forces_table,
            on=['ElementID', 'axial']
        )

        max_shear_cbushes_table.drop_duplicates(subset=['shear'], inplace=True)
        max_axial_cbushes_table.drop_duplicates(subset=['axial'], inplace=True)

        reorder_max_shear_cbushes_table = max_shear_cbushes_table[table_cols]
        reorder_max_axial_cbushes_table = max_axial_cbushes_table[table_cols]

        # Generating the table with the envelope (Max shear and Max)
        # ( Shear and Axial LC NOT correlated)

        # --> Renaming Subcases Columns
        max_shear_cbushes_table.rename(
            columns={'SubcaseID': 'SubcaseID shear'},
            inplace=True)
        max_axial_cbushes_table.rename(
            columns={'SubcaseID': 'SubcaseID axial'},
            inplace=True)

        max_shear_table_columns = max_shear_cbushes_table[[
            'SubcaseID shear', 'shear']]
        max_axial_table_columns = max_axial_cbushes_table[[
            'SubcaseID axial', 'axial']]

        envlp_table = pd.merge(
            left=max_shear_table_columns,
            right=max_axial_table_columns,
            on=['ElementID']
        )

        # Generating Info data:
        info_dict = {
            'Op2 file': op2_path,
            'Date': datetime.today().strftime('%d-%m-%Y'),
            'PyScript': 'get_cbush_forces',
            'Version': __version__
        }

        info_table = pd.DataFrame(
            data=info_dict.values(),
            index=info_dict.keys(),
            columns=[''],
        )

        # Exporting into excel file
        if add_joint_map.selection.upper() == 'YES':
            excel_sheets = {
                'Info': info_table,
                'Joint_Mapping': cbushes_joint_table,
                'Joint_Max_Shear': reorder_max_shear_cbushes_table,
                'Joint_Max_Axial': reorder_max_axial_cbushes_table,
                'Joint_Envelope': envlp_table,
                'Joint_ALL_Data ': cbushes_forces_table,
            }
        else:
            excel_sheets = {
                'Info': info_table,
                'Joint_Max_Shear': reorder_max_shear_cbushes_table,
                'Joint_Max_Axial': reorder_max_axial_cbushes_table,
                'Joint_Envelope': envlp_table,
                'Joint_ALL_Data ': cbushes_forces_table,
            }

        to_excel(
            filepath=excel_path,
            sheets=excel_sheets,
            include_index=True)
    else:
        print('[Error] : Not enough inputs to compute the CBUSH forces.')
