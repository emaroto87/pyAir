# -*- coding: utf-8 -*-

from typing import Optional, Dict, List, Literal
import re
import os
import numpy as np
import pandas as pd
import tempfile as tmpf

from common.ui.filedialogs import askOpenOP2
from common.ui.filedialogs import askOpenGroupFile
from common.ui.filedialogs import askSaveAsExcel
from common.ui.tui import ask_entities_fields_list
from common.ui.tui import TUI as ui_menu

from common.errors.errors import check_type
from common.errors.errors import check_option
from common.general import __ENTITY_TYPES
from common.general import load_op2_file
from common.general import get_result_from_op2
from common.general import split_df_by_subcases
from common.general import filter_df_by_entity


from common.io.writers.excel_writer import write_dataframes_to_excel as to_excel


# from general import get_envelope_from_df

from pyNastran.bdf.bdf import BDF
from pyNastran.op2.op2_geom import OP2Geom

from scipy.spatial import ConvexHull

from tqdm import tqdm

from collections import defaultdict
from collections import namedtuple

from pathlib import Path

__version__ = '1.0.8'
__releaseDateShort__ = '2025/12/03'
__releaseDateLong__ = 'Feb 05, 2026'

__author__ = 'Enrique Maroto'
__email__ = 'enrique.maroto-herreros@capgemini.com'
__desc__ = 'Nastran OP2 Post-Processing tool'


'''
Control de versiones:
    v1.0.0
    ------
    Primera release

    v1.0.1
    ------
        1. Errores solucionados:
            * read_sets_from_file:
                En la version 1.0.0 no era capaz de leer grupos de un solo
                miembro. Ahora si. Por otro lado, la funcion tenia problemas
                para detectar las etiquetas de los grupos si éstos incluian
                el caracter "-". Esta solucionado.

        2. Se ha incluido...
            * CROD_Forces

        3. Se ha eliminado...
            * N/A

    v1.0.2
    ------
        1. Errores solucionados
            * read_sets_from_file:
                Al meter la opción de leer grupos de un solo miembro se jodió
                la lectura de groupos con multiples miembros.
        2. Se ha incluido
            2.1 Interfaz grafica con Tkinter

    v1.0.3
    ------
        1. Errores solucionados: N/A

        2. Se ha incluido...
        2.1
        
    v1.0.7
    ------
        1. Errores solucionados: 
            * Había un error a en la función group_envelope que hacía que no
            no reconociera los ids en la lista de nodos o de elementos
        2. Errores    
'''


# Creating and definition of Links
link = namedtuple('link', ['entity_type', 'op2_results_labels', 'op2_results'])

__RealPlateForce = link(
    __ENTITY_TYPES[1],
    ['force.cquad4_force', 'force.ctria3_force'],
    ('mx', 'my', 'mxy', 'bmx', 'bmy', 'bmxy', 'tx', 'ty')
)
__RealPlateStress = link(
    __ENTITY_TYPES[1],
    ['stress.cquad4_stress', 'stress.ctria3_stress'],
    ('fiber_distance', 'oxx', 'oyy', 'txy',
     'angle', 'omax', 'omin', 'von_mises')
)
__RealPlateStrain = link(
    __ENTITY_TYPES[1],
    ['strain.cquad4_strain', 'strain.ctria3_strain'],
    ('fiber_curvature', 'exx', 'eyy', 'exy',
     'angle', 'emax', 'emin', 'von_mises')
)
__RealMPCForces = link(
    __ENTITY_TYPES[0],
    ['mpc_forces'],
    ('t1', 't2', 't3', 'r1', 'r2', 'r3')
)
__RealSPCForces = link(
    __ENTITY_TYPES[0],
    ['spc_forces'],
    ('t1', 't2', 't3', 'r1', 'r2', 'r3')
)
__RealDisplacement = link(
    __ENTITY_TYPES[0],
    ['displacements'],
    ('t1', 't2', 't3', 'r1', 'r2', 'r3')
)
__RealGriPointForces = link(
    __ENTITY_TYPES[0],
    ['grid_point_forces'],
    ('f1', 'f2', 'f3', 'm1', 'm2', 'm3')
)
__RealCBushForce = link(
    __ENTITY_TYPES[1],
    ['force.cbush_force'],
    ('fx', 'fy', 'fz', 'mx', 'my', 'mz')
)
__RealCBarForce = link(
    __ENTITY_TYPES[1],
    ['force.cbar_force'],
    ['bending_moment_a1',
     'bending_moment_a2',
     'bending_moment_b1',
     'bending_moment_b2',
     'shear1',
     'shear2',
     'axial',
     'torque']
)
__RealCRodForce = link(
    __ENTITY_TYPES[1],
    ['force.crod_force'],
    [
        'axial',
        'torsion'
    ]
)
__RealCBeamForce = link(
    __ENTITY_TYPES[1],
    ['force.cbeam_force'],
    ['sd',
     'bending_moment_1',
     'bending_moment_2',
     'shear1',
     'shear2',
     'axial_force',
     'total_torque',
     'warping_torque']
)

__Result_Labels_Links = {
    'Displacement': __RealDisplacement,
    'Plate_Force': __RealPlateForce,
    'Plate_Stress': __RealPlateStress,
    'Plate_Strain': __RealPlateStrain,
    'CBush_Force': __RealCBushForce,
    'CRod_Force': __RealCRodForce,
    'CBar_Force': __RealCBarForce,
    'CBeam_Force': __RealCBeamForce,
    'SPC_Forces': __RealSPCForces,
    'MPC_Forces': __RealMPCForces,
    'Grid_Point_Forces': __RealGriPointForces,
}

ENVELOPE_OPTIONS = (
    'MAX',
    'MIN',
    'EXTREME',
    'MAX_ABS',
    'QHULL',
    'AVERAGE',
    'NONE'
)


class ResultLink:
    def __init__(self, links_dict: Optional[Dict[str, str]] = None):
        if links_dict:
            for label, link in links_dict.items():
                self.add_link(label, link)

    def add_link(self, link_label: str, link):
        setattr(self, link_label, link)


def qhull(df: pd.DataFrame,
          columns: List[str],
          debug: Optional[bool] = False) -> pd.DataFrame:
    '''

    Parameters
    ----------
    df : pandas DataFrame
        Table of data.

    columns : List[str]
        List containing the labels of the columns from which extract and
        generate the envelope

    debug : Optional[bool], optional
        If True, it prints more information for debugging purposes. 
        The default is False.

    Returns
    -------
    envelope : pd.Dataframe
        pandas DataFrame with the Qhull algorythm applied.

    '''

    # Checkings
    check_type(df, pd.DataFrame)
    check_type(columns, list)
    check_type(debug, bool)

    # Filtering the data by the selected columns
    df.reset_index(inplace=True)
    data = df[columns].to_numpy()
    if debug:
        print(f'Points :\n{data}')

    hull = ConvexHull(data)
    if debug:
        print('QHull envelope Vertices:')
        print(hull.vertices)

    envelope_index = df.index[hull.vertices]
    envelope = df[df.index.isin(envelope_index)]
    if debug:
        print(f'Output Table:\n{envelope}')

    return envelope


def get_all_ids(op2_filename: Path, debug_geom: Optional[bool] = False):
    op2 = OP2Geom(make_geom=True, debug=debug_geom)
    geom = op2.read_op2(op2_filename)
    ids = {
        __ENTITY_TYPES[0]: geom.nodes,
        __ENTITY_TYPES[1]: geom.elements,
    }
    return ids


def gen_all_group(
        op2_path: Path,
        entity_type: str,
        debug_geom: Optional[bool] = False
) -> Dict[str, list]:
    '''
    Given a NASTRAN OP2 file, it generates a group with the label "ALL" 
    containing either all the ids of the nodes or elements present in the file.

    Parameters
    ----------
    op2_filename : Path
        Path of the NASTRAN OP2 file.
    entity_type : str
        Either one of the following options : NodeID or ElementID
    debug_geom : Optional[bool], optional
        Parameter inherited from OP2Geom class. If False, it skips some errors
        while reading the OP2 file. The default is False.

    Raises
    ------
    ValueError
        When the entity option is not one of the allowable ones NodeID or 
        Element ID.

    Returns
    -------
    Dict[str,list]
        Returns a dictionary with one unique key "ALL" whose value is a list
        containg all the ids.

    '''

    # Opening the OP2 to read geometry
    op2 = OP2Geom(make_geom=True, debug=debug_geom)
    op2.read_op2(op2_path)

    # Extracting all the ids
    groups = {}

    if entity_type == __ENTITY_TYPES[0]:
        groups['ALL'] = list(op2.node_ids)
    elif entity_type == __ENTITY_TYPES[1]:
        groups['ALL'] = list(op2.element_ids)
    else:
        raise ValueError(f'{entity_type} is not a valid ENTITY_TYPE')

    return groups


def get_result_ids_group(op2_path, result_link):
    o = load_op2_file(op2_path)

    # Extracting all the ids
    groups = {}
    ids = []
    for label in result_link.op2_results_labels:
        r = getattr(o, label)
        t = r[list(r.keys())[0]]
        if result_link.entity_type == __ENTITY_TYPES[0]:
            ids += list(t.node_gridtype[:, 0])

        elif result_link.entity_type == __ENTITY_TYPES[1]:
            ids += list(t.element)
        else:
            raise ValueError(
                f'{result_link.entity_type} is not a valid ENTITY_TYPE')

    groups['ALL'] = ids

    return groups


def envelope(df: pd.DataFrame,
             envelope_option: str,
             envelope_columns: List[str]
             ) -> pd.DataFrame:
    '''


    Parameters
    ----------
    df : pandas DataFrame
        Table of data.
    envelope_option : str
        One of the available options ('MAX','MIN','EXTREME','MAX_ABS','QHULL','AVERAGE','NONE').
    envelope_columns : List[str]
        List containing the labels of the columns from which extract and
        generate the envelope

    Raises
    ------
    ValueError
        If the any of the columns are not included in the table
    AssertionError
        For the En

    Returns
    -------
    pandas Dataframe
        Table with the points that generate the envelope

    '''

    # Checking arguments
    check_type(df, pd.DataFrame)
    check_type(envelope_option, str)
    check_type(envelope_columns, list)

    for envelope_column in envelope_columns:
        if not (envelope_column in df.columns):
            raise ValueError(
                f'Column "{envelope_column}" is not present in Pandas DataFrame'
                f'.Please review envelope_columns list'
            )

    if envelope_option in ENVELOPE_OPTIONS:
        # MAX OPTION
        if envelope_option.upper() == ENVELOPE_OPTIONS[0]:
            if len(envelope_columns) == 1:
                df_env_rows = []
                column = envelope_columns[0]
                max_val = df[column].max()
                max_row = df[df[column] == max_val]
                df_env_rows += [max_row]
                return pd.concat(df_env_rows)
            else:
                raise AssertionError(
                    f'When using envelope options {ENVELOPE_OPTIONS[0]}, '
                    f'{ENVELOPE_OPTIONS[1]},{ENVELOPE_OPTIONS[2]} and '
                    f'{ENVELOPE_OPTIONS[3]} the number of elements within '
                    f'envelope_columns list must be only one. Please check'
                )
        # MIN OPTION
        if envelope_option.upper() == ENVELOPE_OPTIONS[1]:
            if len(envelope_columns) == 1:
                df_env_rows = []
                column = envelope_columns[0]
                min_val = df[column].min()
                min_row = df[df[column] == min_val]
                df_env_rows += [min_row]
                return pd.concat(df_env_rows)
            else:
                raise AssertionError(
                    f'When using envelope options {ENVELOPE_OPTIONS[0]}, '
                    f'{ENVELOPE_OPTIONS[1]},{ENVELOPE_OPTIONS[2]} and '
                    f'{ENVELOPE_OPTIONS[3]} the number of elements within '
                    f'envelope_columns list must be only one. Please check'
                )

        # EXTREME OPTION
        if envelope_option.upper() == ENVELOPE_OPTIONS[2]:
            if len(envelope_columns) == 1:
                df_env_rows = []
                column = envelope_columns[0]
                max_val = df[column].max()
                min_val = df[column].min()
                max_row = df[df[column] == max_val]
                min_row = df[df[column] == min_val]
                df_env_rows += [max_row]
                df_env_rows += [min_row]
                return pd.concat(df_env_rows)
            else:
                raise AssertionError(
                    f'When using envelope option {ENVELOPE_OPTIONS[0]} or '
                    f'{ENVELOPE_OPTIONS[1]} the number of elements within '
                    f'envelope_columns list must be only one. Please check'
                )
        # MAX_ABS OPTION
        if envelope_option.upper() == ENVELOPE_OPTIONS[3]:
            if len(envelope_columns) == 1:
                column = envelope_columns[0]
                max_val = df[column].max()
                min_val = df[column].min()
                max_abs = max_val if abs(max_val) >= abs(min_val) else min_val
                max_abs_row = df[df[column] == max_abs]
                # print(max_abs_row)
                return max_abs_row
            else:
                raise AssertionError(
                    f'When using envelope option {ENVELOPE_OPTIONS[0]} or '
                    f'{ENVELOPE_OPTIONS[0]} the number of elements within '
                    f'envelope_columns list must be only one. Please check'
                )

        # QHULL OPTION
        if envelope_option.upper() == ENVELOPE_OPTIONS[4]:
            if len(envelope_columns) >= 2:
                return qhull(df, envelope_columns)

        # AVERAGE
        # if envelope_option.upper() == ENVELOPE_OPTIONS[5]):

    else:
        raise ValueError(f'{envelope_option} is not a valid option. Please sel'
                         'ect one of the valid options {ENVELOPE_OPTIONS}')


def group_envelope(dataframe,
                   groups: dict,
                   result_link,
                   envelope_type: str,
                   envelope_columns: list,
                   groups_summary: bool = False,
                   output_path: str = None,
                   debug: bool = False
                   ):

    tables = {}
    entity_type = result_link.entity_type

    # Flag para saber si la entidad esta en el indice o en la columna y
    # así exportarlo en el excel.
    index_flag = False
    try:
        dataframe_members = list(dataframe[entity_type])
    except:
        pass
    try:
        dataframe_members = list(dataframe.index)
        index_flag = True
    except:
        print('Quique dale una vuelta que falla algo...MELON!!!')
    print(index_flag)

    # --> BEGIN LOOP: for each group if list of groups is provided
    for label, group in groups.items():

        print(f'Analysing group {label}...')
        # Step 1 : the raw table is filtered per each group of the list of groups
        # The table contains all the results for each load case and for all
        # the entities of the group. In other words, is the raw table per
        # each group (stored in raw_group_table)
        raw_group_table = filter_df_by_entity(
            dataframe,
            group,
            entity_type=entity_type
        )
        print(raw_group_table)
        # Step 2 : Performing the envelope operations
        # In this part, the raw_group_table is taken and filtered per each
        # item of the contained in the group as preliminary step before
        # performing any "envelope" operation (max,min,qhull..). Once the
        # filtering is done, the envelope operation is done per each item
        # table.

        if not (envelope_type.upper() in ['NONE', 'AVERAGE']):
            items_env_table = []
            for item in tqdm(group):
                # print(item)
                if item in dataframe_members:
                    if debug:
                        print(f'\t{entity_type} : {item}')

                    # Look for the item in the column or index
                    if entity_type in raw_group_table.columns:
                        item_table = raw_group_table[
                            raw_group_table[entity_type].isin([item])]
                    else:
                        item_table = raw_group_table[
                            raw_group_table.index.isin([item])]

                    # Apply the envelope if necessary
                    try:
                        df_eid_env = envelope(
                            item_table,
                            envelope_type,
                            envelope_columns
                        )
                        items_env_table += [df_eid_env]
                    except:
                        print(f'Unable to get envelope type {envelope_type}'
                              f'for {entity_type} : {item}'
                              )
                else:
                    print(
                        f'{entity_type} {item} is not in the results. Please '
                        f'check if the ID or Type is correct'
                    )
                    continue

            # Concatenating all the envelopes per ID of the group
            if len(items_env_table) == 0:
                env_group_table = pd.DataFrame()
            else:
                env_group_table = pd.concat(items_env_table)

        # ENVELOPE : AVERAGE
        elif envelope_type.upper() in ['AVERAGE']:
            group_average_table = []
            subcases = getattr(dataframe, 'SubcaseID').unique()
            for subcase in tqdm(subcases):
                # Paso 1: filtar por subcaseID
                table_per_subcase = raw_group_table[raw_group_table['SubcaseID'].isin([
                                                                                      subcase])]
                # Paso 2: Obtener los valores medios
                sub_average_data = table_per_subcase.mean().to_frame().transpose()
                # Paso 3: Cambiar el Entity_ID por "Average_" + Subcase
                sub_average_data[entity_type] = 'Average_' + str(subcase)
                # Paso 4: Establecer como indice 'Entity_ID'
                sub_average_data.set_index(entity_type, inplace=True)
                # Paso 5: añadimos la media del subcase a una lista
                group_average_table += [sub_average_data]
                # Generamos la tabla con las medias de todos los subcases
                env_group_table = pd.concat(group_average_table)
        else:
            env_group_table = raw_group_table

        # Adding the label of the group
        env_group_table['Group'] = label

        # Storing the result into a dict
        if env_group_table.empty:
            tables[label] = env_group_table
        else:
            tables[label] = env_group_table.drop_duplicates()
    # <-- END LOOP : for each group if list of groups is provided

    summary = []
    if len(tables.keys()) != 0:
        if groups_summary and (envelope_type.upper() == 'NONE'):
            raise AssertionError(
                'Is not possible to combine groups_summary = True and'
                ' envelope_type = None.'
            )
        else:
            if groups_summary:
                for label, table in tables.items():
                    # Generating the envelope of the group
                    env_group_table = envelope(
                        table,
                        envelope_type,
                        envelope_columns
                    )
                    # Adding the result to the summary table
                    summary += [env_group_table]

                summary_table = pd.concat(summary)
                summary_table.set_index('Group', inplace=False)

                # Include the Summary table
                tables['Summary'] = summary_table.drop_duplicates()

        # Exporting the result into an excel
        if not (output_path is None):
            # try:
            to_excel(
                filepath=output_path,
                sheets=tables,
                include_index=index_flag)
            # except:
            #     print('Error while exporting to Excel.')
    else:
        print(
            '[ERROR] : No data has been generated due to multiple errors.'
            'Please check warning and error messages'
        )
    return tables


def raw_table(op2_filename: Path, result_link: str) -> pd.DataFrame:
    '''
    Reads a NASTRAN OP2 file and returns a pandas DataFrame with all the
    selected results through the link.

    Parameters
    ----------
    op2_filename : Path
        Path of the NASTRAN OP2 file.
    result_link : str


    Returns
    -------
    raw_table : pandas DataFrame
        Dataframe containing all the data from the OP2 file. The table has an
        single index. If the data from the OP2 is stored within a Multi-Index
        table, it is converted into single-Index one.

    '''

    check_type(op2_filename, str)
    # check_type(result_link,str)
    # check_option(result_link,list(Result_Labels_Links.keys()))

    try:
        op2 = load_op2_file(op2_filename)
    except IOError:
        pass

    dfs = []
    op2_result_labels = result_link.op2_results_labels

    for op2_result_label in op2_result_labels:

        if hasattr(op2, op2_result_label):

            print(f'Reading {op2_result_label}...', end='')
            df = get_result_from_op2(op2, op2_result_label)
            dfs += [df]
            print('Done')

        else:
            print(f'[WARN] {op2_result_label} not found in OP2 file')

    raw_table = pd.concat(dfs)

    return raw_table


def read_sets_from_file(filename: Path,
                        comment_symbol: Optional[str] = '$',
                        members_separator: Optional[str] = ',',
                        set_keyword: Optional[str] = 'SET',
                        verbose: bool = False,
                        debug: bool = False,
                        ) -> Dict[str, list]:
    '''
    It reads the sets defined within a text file.

    The format of the set files is as follows:

        {Set_keyword} {label} = id_1 {sep} id_2 {sep} id_3 {sep} ....

    Example:

        SET Panel2-4 = 1001, 1002, 1101, 1002

    In this example, the label of the set is "Panel2-4" and the seperator of 
    the id of the members of the group is in this case ",". Additionally, the 
    set_keyword is in this case "SET".

    Parameters
    ----------
    filename : Path
        Path of the text file containing the definition of the sets.
    comment_symbol : Optional[str], optional
        String that if placed at the beginning of each line, is interpreted as
        a commented line and its content is not read. The default is '$'.
    members_separator : Optional[str], optional
        String that is used to separate the members of the set. The default is ','.
    set_keyword : Optional[str], optional
        String used by the function to identify the beginning of the defintion
        of a new set. The default is 'SET'.
    verbose : bool, optional
        Option to print further information. The default is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    dict
        DESCRIPTION.

    '''

    # Preliminary checks
    if not (check_type(filename, str)):
        pass

    if not (check_type(comment_symbol, str)):
        pass

    if not (check_type(members_separator, str)):
        pass

    if not (check_type(set_keyword, str)):
        pass

    if not (os.path.isfile(filename)):
        raise ValueError('The path does not corrspond to a file')

    sets = {}
    # Pattern for multiple members
    pattern1 = r'{0}([\w-]+)=([\d,]+)'.format(set_keyword)

    # Pattern for single member
    pattern2 = r'{0}([\w-]+)=([\d]+)'.format(set_keyword)

    # Read file content and close
    with open(filename, 'r') as f:
        lines = f.readlines()

    # cleaning commented lines
    clean_lines = [line[:-1] for line in lines
                   if not (line.startswith(comment_symbol))]

    clean_lines = ''.join(clean_lines)
    clean_lines = clean_lines.replace(' ', '')
    clean_lines = clean_lines.strip()

    if debug:
        print(clean_lines)

    # Spliting the content by sets using the pattern
    try:
        findings = re.findall(pattern1, clean_lines)
    except:
        if debug:
            print('Error in pattern1')
        try:
            findings = re.findall(pattern2, clean_lines)
        except:
            if debug:
                print('Error in pattern2')
            pass

    if len(findings) == 0:
        raise ValueError(
            'Unable to locate any set'
        )

    if debug:
        print(findings)

    # Read the members of each set and store them into a list
    for set_label, set_content in findings:
        try:
            members_list = [int(x) for x in set_content.split(',')]
        except:
            members_list = [int(set_content)]
        if set_label.isnumeric():
            set_label = int(set_label)
        sets[set_label] = members_list

    return sets


# def to_excel(tables: Dict[str, pd.DataFrame], output_filename: Path):
#     '''
#     Export into an Excel file multiple pandas DataFrames stored within a
#     dictionary

#     Parameters
#     ----------
#     tables : Dict [str, pd.DataFrame ]
#         Dictionary containing the pandas DataFrames. Each key of the dictionary
#         is used to label the sheet o1f the Excel book file.
#     filename : str path-like
#         Path of the output Excel file.

#     Returns
#     -------
#     None.

#     '''
#     # Checking
#     check_type(tables, dict)
#     # check_type(output_filename,str)

#     print('Exporting to Excel File ...', sep='')
#     with pd.ExcelWriter(output_filename) as writer:
#         for sheet_name, table in tables.items():
#             table.to_excel(writer, sheet_name=sheet_name, index=False)
#     print('[Done]')


# Initialising results links
r_link = ResultLink(__Result_Labels_Links)


def gui_console(
        export_excel_opt=True,
        groups_summary_opt=True,
        enforce_envelope=None,
):

    link_menu = ui_menu(
        options=list(__Result_Labels_Links.keys()),
        label='link_menu',
        sel_msg='Select the type of result to post-process',
        is_main=True
    )
    yesno_set_menu = ui_menu(
        options=['Yes', 'No'],
        label='yesno_set_menu',
        sel_msg='Filter by Sets File [Yes/No] : '
    )
    yesno_summary_menu = ui_menu(
        options=['Yes', 'No'],
        label='yesno_summary_menu',
        sel_msg='Generate the Summary of all groups [Yes/No] : '
    )

    env_type_menu = ui_menu(
        options=list(ENVELOPE_OPTIONS),
        label='env_type_menu',
        sel_msg='Select the type of Envelope'
    )

    while True:
        # PART 1 : GET THE PARAMETERS USING THE CONSOLE INTERFACE MENU

        # Main Menu :
        # ---------------------------------------------------------------------
        # Step 1 : Requesting the type of result
        link_menu.loop()
        if link_menu.selection.upper() == 'EXIT':
            break
        else:
            result_link = getattr(r_link, link_menu.selection)

            # [Step 2] : Requesting the OP2 file ------------------------------
            op2 = askOpenOP2()
            if op2 is None:
                print('No OP2 file selected. Returning to main menu.')
                continue

            # [Step 3] : Asking for group filtering ---------------------------
            yesno_set_menu.loop()
            if yesno_set_menu.selection.upper() == 'YES':
                # Option in which the user provides a plain text file with
                # the groups.
                sets_file = askOpenGroupFile()
                if sets_file is None:
                    continue
                else:
                    sets = read_sets_from_file(
                        filename=sets_file
                    )
                    print(sets)
            else:
                # Option in which there is no specific defined groups. Then
                # it is extracted from the OP2 all the entities which have
                # the type of result selected by the user.
                sets = get_result_ids_group(
                    op2_path=op2,
                    result_link=result_link
                )
                print(sets)

            # [Step 4] : Asking for type of envelope filtering ------[OPTIONAL]
            if enforce_envelope is None:
                env_type_menu.loop()
                envelope_selection = env_type_menu.selection.upper()
            else:
                envelope_selection = enforce_envelope.upper()

            # Askig for field
            if not (envelope_selection in ['NONE', 'AVERAGE']):
                print(f'The following fields are available for'
                      f'{link_menu.selection}:\n{result_link.op2_results}')
                envelope_columns = ask_entities_fields_list()
            else:
                envelope_columns = None

            # [Step 5] : Asking for summary table -------------------[OPTIONAL]
            if groups_summary_opt:
                yesno_summary_menu.loop()
                summary = True if yesno_summary_menu.selection.upper() == 'YES' else False
            else:
                summary = False

            # [Step 5] : Asking for output Excel file path ----------[OPTIONAL]
            if export_excel_opt:
                out = askSaveAsExcel()
                output_path = out[0]
                if output_path == '':
                    print('No Excel file selected. Returning to main menu.')
                    continue
            else:
                output_path = None

            tables = group_envelope(
                dataframe=raw_table(op2, result_link),
                groups=sets,
                result_link=result_link,
                envelope_type=envelope_selection,
                envelope_columns=envelope_columns,
                output_path=output_path,
                debug=False,
                groups_summary=summary
            )

    return tables


def get_average_envelope():
    # Menu para extraer las cargas
    tables = gui_console(
        export_excel_opt=False,
        groups_summary_opt=False,
        enforce_envelope='AVERAGE',
    )

    # Preguntamos los campos a partir de los cuales generar la envolvente
    table_cols = list(list(tables.values())[0].columns)
    print(f'The following fields are available :\n'
          f'\n{table_cols}')
    fields = ask_entities_fields_list()

    # Preguntamos la ruta de la hoja excel
    out = askSaveAsExcel()

    # Bucle
    envs = {}
    for group, table in tables.items():
        qhull_envelope = envelope(
            df=table,
            envelope_option='QHULL',
            envelope_columns=fields)
        envs[group] = qhull_envelope

    # Exportamos a excel
    to_excel(
        tables=envs,
        output_filename=out[0],
        include_index=True
    )

    return tables, envs
