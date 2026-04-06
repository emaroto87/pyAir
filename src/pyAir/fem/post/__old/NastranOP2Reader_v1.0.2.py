# -*- coding: utf-8 -*-
"""
Created on Fri Jun 13 11:48:48 2025

@author: U69432
"""
from typing import Optional 
from typing import Dict
from typing import List
from typing import Literal
import re
import os
import numpy as np
import pandas as pd
import tempfile as tmpf


from common.general import ask_op2_filepath
from common.general import ask_set_filepath
from common.general import load_op2_file
from common.general import asksaveasexcelfile
from common.general import check_type
from common.general import check_option
from common.general import __ENTITY_TYPES
from common.general import ui_menu



from PyNastran_OP2_tools import get_result_from_op2
from PyNastran_OP2_tools import split_df_by_subcases
from PyNastran_OP2_tools import filter_df_by_entity
from PyNastran_OP2_tools import get_envelope_from_df
from PyNastran_OP2_tools import ask_entities_fields_list
from pyNastran.bdf.bdf import BDF 
from pyNastran.op2.op2_geom import OP2Geom

from scipy.spatial import ConvexHull

from tqdm import tqdm

from collections import defaultdict
from collections import namedtuple

from pathlib import Path

__version__ = '1.0.2'
__releaseDateShort__ = '2025/07/14'
__releaseDateLong__ = 'July 14, 2025'

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
'''




# Creating and definition of Links
link = namedtuple('link', ['entity_type','op2_results_labels','op2_results'])

__RealPlateForce = link(
    __ENTITY_TYPES[1], 
    ['force.cquad4_force','force.ctria3_force'], 
    ('mx','my','mxy','bmx','bmy','bmxy','tx','ty')
    )
__RealPlateStress = link(
    __ENTITY_TYPES[1], 
    ['stress.cquad4_stress','stress.ctria3_stress'], 
    ('fiber_distance', 'oxx', 'oyy', 'txy', 'angle', 'omax', 'omin', 'von_mises')
    )
__RealPlateStrain = link(
    __ENTITY_TYPES[1],
    ['strain.cquad4_strain','strain.ctria3_strain'],
    ('fiber_curvature', 'exx', 'eyy', 'exy', 'angle', 'emax', 'emin', 'von_mises')
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
    ('f1','f2','f3','m1','m2','m3')
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
    'Displacement'    : __RealDisplacement,
    'Plate_Force'      : __RealPlateForce,
    'Plate_Stress'     : __RealPlateStress,
    'Plate_Strain'     : __RealPlateStrain,
    'CBush_Force'      : __RealCBushForce,
    'CRod_Force'       : __RealCRodForce,
    'CBar_Force'       : __RealCBarForce,
    'CBeam_Force'      : __RealCBeamForce,
    'SPC_Forces'       : __RealSPCForces,
    'MPC_Forces'       : __RealMPCForces,
    'Grid_Point_Forces'  : __RealGriPointForces,
    }



class ResultLink:
    def __init__(self,links_dict: Optional[Dict[str, str]] = None):
        if links_dict:
            for label, link in links_dict.items():
                self.add_link(label, link)
    
    def add_link(self,link_label:str,link):
        setattr(self,link_label,link)

def qhull(df : pd.DataFrame, 
          columns : List[str], 
          debug : Optional[bool] = False) -> pd.DataFrame:
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
    check_type(columns,list)
    check_type(debug,bool)
  
    # Filtering the data by the selected columns
    data = df[columns].to_numpy()
    if debug : print(f'Points :\n{data}')
    
    hull = ConvexHull(data)
    if debug : 
        print('QHull envelope Vertices:')
        print(hull.vertices)
    
    envelope_index = df.index[hull.vertices]
    envelope = df[df.index.isin(envelope_index)]
    if debug : print(f'Output Table:\n{envelope}')
    
    return envelope

def envelope(df : pd.DataFrame, 
             envelope_option : str, 
             envelope_columns: List[str]
             )->pd.DataFrame:
    '''
    

    Parameters
    ----------
    df : pandas DataFrame
        Table of data.
    envelope_option : str
        .
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
    check_type(envelope_columns,list)
    
    for envelope_column in envelope_columns:
        if not(envelope_column in df.columns):
            raise ValueError(
                f'Column "{envelope_column}" is not present in Pandas DataFrame' \
                f'.Please review envelope_columns list'
                )
    
    if envelope_option in ENVELOPE_OPTIONS:
        # EXTREME OPTION
        if envelope_option.upper() == ENVELOPE_OPTIONS[0]:
            if len(envelope_columns) == 1 :
                df_env_rows=[]
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
                    f'When using envelope option {ENVELOPE_OPTIONS[0]} or '    \
                    f'{ENVELOPE_OPTIONS[1]} the number of elements within '    \
                    f'envelope_columns list must be only one. Please check'    \
                    )
        # MAX_ABS OPTION
        if envelope_option.upper() == ENVELOPE_OPTIONS[1]:
            if len(envelope_columns) == 1 :
                column = envelope_columns[0]
                max_val = df[column].max()
                min_val = df[column].min()
                max_abs = max_val if abs(max_val) >= abs(min_val) else min_val
                max_abs_row = df[df[column] == max_abs]
                # print(max_abs_row)
                return max_abs_row
            else:
                raise AssertionError(
                    f'When using envelope option {ENVELOPE_OPTIONS[0]} or '    \
                    f'{ENVELOPE_OPTIONS[0]} the number of elements within '    \
                    f'envelope_columns list must be only one. Please check'    \
                    )
        
        # QHULL OPTION
        if envelope_option.upper() == ENVELOPE_OPTIONS[2]:
            if len(envelope_columns) >=2 :
                return qhull(df,envelope_columns)
    
    
    else:
        raise ValueError(f'{envelope_option} is not a valid option. Please sel'\
                         'ect one of the valid options {ENVELOPE_OPTIONS}')
            
        
def get_all_ids(op2_filename: Path ,debug_geom : Optional[bool] = False):
    op2 = OP2Geom(make_geom = True , debug = debug_geom)
    geom = op2.read_op2(op2_filename)
    ids = {
        __ENTITY_TYPES[0] : geom.nodes,
        __ENTITY_TYPES[1] : geom.elements,
           }
    return ids

def gen_all_group(
        op2_path : Path,
        entity_type : str, 
        debug_geom : Optional[bool] = False
        ) -> Dict[str,list] :
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
    op2 = OP2Geom(make_geom=True , debug = debug_geom)
    op2.read_op2(op2_path)
    
    
    # Extracting all the ids
    groups ={}
    
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
    ids =[]
    for label in result_link.op2_results_labels:
        r = getattr(o, label)
        t = r[list(r.keys())[0]]
        if result_link.entity_type == __ENTITY_TYPES[0]:
            ids  += list(t.node_gridtype[:,0])
    
        elif result_link.entity_type == __ENTITY_TYPES[1]:
            ids += list(t.element) 
        else:
            raise ValueError(f'{result_link.entity_type} is not a valid ENTITY_TYPE')
    
    groups['ALL'] = ids
    
    return groups

def group_envelope(dataframe, 
                   groups: dict, 
                   result_link,
                   envelope_type : str,
                   envelope_columns : list,
                   groups_summary : bool = False,
                   output_path : str = None,
                   debug : bool = False
                   )  : 
        
    entity_type = result_link.entity_type  
    tables = {}
    
    dataframe_members = list(getattr(dataframe, entity_type))
    
    for label, items in groups.items():
        print(f'Analysing group {label}...')
        group_table = filter_df_by_entity(
            dataframe,
            items,
            entity_type = entity_type
            )
        if envelope_type.upper() != 'NONE':

            eids_env_table = []
            for item in tqdm(items):
                if item in dataframe_members:
                    if debug : print(f'\t{entity_type} : {item}')
                    
                    # Look for the item in the column or index
                    if entity_type in group_table.columns:
                        item_table = group_table[
                            group_table[entity_type].isin([item])
                            ]
                    
                    else:
                        item_table = group_table[
                            group_table.index.isin([item])
                            ]
                    
                    # Apply the envelope if necessary
                    try:
                        df_eid_env = envelope(
                            item_table,
                            envelope_type,
                            envelope_columns
                            )
                        eids_env_table += [df_eid_env]
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
            if len(eids_env_table) == 0:
                group_table = pd.DataFrame()
            else:
                group_table = pd.concat(eids_env_table)                                 
   
            
        # Adding to the table the column 'Group'
        group_table['Group'] = label                                             
        # if debug : print(table)
            
        # Storing the result into a dict
        if group_table.empty:
            tables[label] = group_table
        else:
            tables[label] = group_table.drop_duplicates()

            
            
    summary = []
    if len(tables.keys()) != 0:
        if groups_summary and (envelope_type.upper() == 'NONE') :
            raise AssertionError(
                'Is not possible to combine groups_summary = True and' \
                ' envelope_type = None.'
                )
        else:    
            if groups_summary:
    
                for label,table in tables.items():
                    # Generating the envelope of the group
                    group_env_table = envelope(
                        table,
                        envelope_type,
                        envelope_columns
                        )
    
                    # Adding the result to the summary table
                    summary += [group_env_table]
                    
                # print(summary)
                summary_table = pd.concat(summary)
                summary_table.set_index('Group', inplace = False)
                # summary_table.sort_index()
    
                # Include the Summary table
                tables['Summary'] = summary_table.drop_duplicates()
         
        # Exporting the result into an excel   
        if not(output_path is None):
            to_excel(tables,output_path)
    else:
        print(
            '[ERROR] : No data has been generated due to multiple errors.'
            'Please check warning and error messages'
            )
    return tables
    




def raw_table(op2_filename: Path, result_link : str)-> pd.DataFrame:
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

    
    check_type(op2_filename,str)
    # check_type(result_link,str)
    # check_option(result_link,list(Result_Labels_Links.keys()))
    
    try:
        op2 = load_op2_file(op2_filename)
    except IOError:
        pass
    
    dfs = [ ]
    op2_result_labels = result_link.op2_results_labels
    
    for op2_result_label in op2_result_labels :
        
        if hasattr(op2, op2_result_label) :
            
            print(f'Reading {op2_result_label}...', end ='')
            df = get_result_from_op2(op2, op2_result_label)
            dfs += [df]
            print('Done')
        
        else:
            print(f'[WARN] {op2_result_label} not found in OP2 file')
    
    raw_table = pd.concat(dfs)        
    
    return raw_table



def read_sets_from_file(filename : Path, 
                        comment_symbol : Optional[ str ]  = '$',
                        members_separator : Optional[ str ] = ',',
                        set_keyword : Optional[ str ] = 'SET',
                        verbose : bool = False,
                        debug : bool = False,
                        ) -> Dict[ str, list ]:
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
    if not(check_type(filename,str)):
        pass
    
    if not(check_type(comment_symbol,str)):
        pass
    
    if not(check_type(members_separator,str)):
        pass
    
    if not(check_type(set_keyword,str)):
        pass
    
    if not(os.path.isfile(filename)):
        raise ValueError('The path does not corrspond to a file')
        
    sets = {}
    # Pattern for multiple members
    pattern1 = r'{0}([\w-]+)=([\d,]+)'.format(set_keyword)

    # Pattern for single member
    pattern2 = r'{0}([\w-]+)=([\d]+)'.format(set_keyword)

    # Read file content and close    
    with open(filename,'r') as f:
        lines = f.readlines()
    
    # cleaning commented lines
    clean_lines = [line[:-1] for line in lines 
                   if not(line.startswith(comment_symbol))]

    clean_lines = ''.join(clean_lines)
    clean_lines = clean_lines.replace(' ','')
    clean_lines = clean_lines.strip()

    if debug : print(clean_lines)
    
    # Spliting the content by sets using the pattern
    try:
        findings = re.findall(pattern1,clean_lines)
    except:
        if debug : print('Error in pattern1')
        try:
            findings = re.findall(pattern2,clean_lines)
        except:
            if debug: print('Error in pattern2')
            pass
    
    if len(findings)==0:
        raise ValueError(
            'Unable to locate any set'
            )
    
    if debug : print(findings)
    
    # Read the members of each set and store them into a list
    for set_label, set_content in findings:
        try:
            members_list = [int(x) for x in set_content.split(',')]
        except:
            members_list = [int(set_content)]    
        if set_label.isnumeric() : set_label = int(set_label)
        sets[set_label] = members_list
    
    return sets
              
            
def to_excel(tables : Dict[str,pd.DataFrame], output_filename : Path) :
    '''
    Export into an Excel file multiple pandas DataFrames stored within a 
    dictionary

    Parameters
    ----------
    tables : Dict [str, pd.DataFrame ]
        Dictionary containing the pandas DataFrames. Each key of the dictionary
        is used to label the sheet of the Excel book file.
    filename : str path-like
        Path of the output Excel file.

    Returns
    -------
    None.

    '''
    # Checking
    check_type(tables,dict)
    # check_type(output_filename,str)
        
    print('Exporting to Excel File ...', sep='')
    with pd.ExcelWriter(output_filename) as writer:
        for sheet_name, table in tables.items():
            table.to_excel(writer, sheet_name=sheet_name, index=False)  
    print('[Done]')                   
            

#
ENVELOPE_OPTIONS = ('EXTREME','MAX_ABS','QHULL','NONE')

# Initialising results links
r_link = ResultLink(__Result_Labels_Links)

## Menu tipo windows ---------------------------------
import tkinter as tk
from tkinter.messagebox import showinfo
from ttkbootstrap import Window



def gui_windows():
    # Main Menu window
    # w = tk.Tk()
    # w.title('Prueba Menu')
    # w.geometry("600x300")
    w = Window(
        themename='flatly',
        # size = (600,600),
        title = 'NastranOP2Reader v1.0.2',
        resizable = None,
        scaling = False,
        )
    
    
    # Frame para colocar botones y opciones

    
    def show_info(wdg):
        showinfo(
            title = 'Result',
            message = wdg.get()
            )
    
    frame_load = tk.Frame(w)
    frame_load.pack(pady=20, padx=20)
    
    frame_rdbt = tk.Frame(w)
    frame_rdbt.pack(pady=50)

    
    # Boton de cargar OP2
    def select_op2_file():
        path = ask_op2_filepath()
        op2_path.delete(0,tk.END)
        op2_path.insert(0,path) 

    op2_lb = tk.Label(
        master = frame_load,
        text="OP2 file",
        width = 10,
        anchor="w"           
        )
    op2_lb.grid(row=0,column=0,padx=5,pady=5, sticky='w')
   
    op2_path = tk.Entry(
        master = frame_load,
        width=50
        )
    op2_path.grid(row=0,column=1, sticky='w')
    # Nota: si encadenas esta linea con la anterior, te dará el error de que
    # no encuentra la función delete() dentro de "select_op2_file". Esto es 
    # porque si lo metes entero, no guardas el widget, guardas el resultado de
    # colocarlo con la función grid, que es un None.
    
    op2_btn = tk.Button(
        master = frame_load,
        text = '\U0001F4C2',
        command = select_op2_file
        )
    op2_btn.grid(row=0,column=2, sticky='w')
    
        
    def select_xls_file():
        path = asksaveasexcelfile()
        xls_path.delete(0,tk.END)
        xls_path.insert(0,path)
    
    xls_lb = tk.Label(
        master = frame_load,
        text="Excel file", 
        width = 10, 
        anchor="w"
        )
    xls_lb.grid(row=1,column=0,padx=5,pady=5, sticky='w')
    
    xls_path = tk.Entry(
        master = frame_load,
        width=50
        )
    xls_path.grid(row=1,column=1, sticky='w')

    xls_btn = tk.Button(
        master = frame_load,
        text = '\U0001F4C2',
        command = select_xls_file
        )
    xls_btn.grid(row=1,column=2, sticky='w')
    
    def select_set_file():
        path = ask_set_filepath()
        set_path.delete(0,tk.END)
        set_path.insert(0,path)
    
    set_lb = tk.Label(
        master = frame_load,
        text="Set file", 
        width = 10, 
        anchor="w"
        )
    set_lb.grid(row=2,column=0,padx=5,pady=5, sticky='w')
    
    set_path = tk.Entry(
        master = frame_load,
        width=50
        )
    set_path.grid(row=2,column=1, sticky='w')

    set_btn = tk.Button(
        master = frame_load,
        text = '\U0001F4C2',
        command = select_set_file
        )
    set_btn.grid(row=2,column=2, sticky='w')
    
    
    # OPTIONS:
        
    frame_rdbt = tk.Frame(
        master = frame_load
        )
    frame_rdbt.grid(pady=50)
    
    # RadioButton : results options
    opt1 = tk.StringVar(value="")
    frame_results = tk.LabelFrame(
        master = frame_rdbt,
        text = 'Result Options'
        )
    frame_results.grid(row=0,column=0,padx=20,sticky='nw')
    for k,v  in __Result_Labels_Links.items():
        tk.Radiobutton(
            master = frame_results,
            text = k,
            variable = opt1,
            value = k,
            anchor = 'w',
            justify ='left'
            ).grid(sticky='w')
    
    # RadioButton : envelope_options
    opt2 = tk.StringVar(value="")

    frame_env = tk.LabelFrame(
        master = frame_rdbt,
        text = 'Envelope Options'
        )
    frame_env.grid(row=0,column=1,padx=20,sticky='nw')
    
    for env_opt in ENVELOPE_OPTIONS:
        tk.Radiobutton(
            master = frame_env,
            text = env_opt,
            variable = opt2,
            value = env_opt,
            anchor = 'w',
            justify ='left'
            ).grid(sticky='w')
    

    # RadioButton : Generate Summary
    opt3 = tk.StringVar(value="")
     
    frame_summary = tk.LabelFrame(
         master = frame_rdbt,
         text = 'Group Summary'
         )
    frame_summary.grid(row=0,column=3,padx=20,sticky='nw')
     
    yes_rbtn = tk.Radiobutton(
                 master = frame_summary ,
                 text = 'Yes',
                 variable = opt3,
                 value = 'Yes',
                 anchor = 'w',
                 justify ='left'
                 )
    yes_rbtn.grid(sticky='w')
    no_rbtn = tk.Radiobutton(
                 master = frame_summary ,
                 text = 'No',
                 variable = opt3,
                 value = 'No',
                 anchor = 'w',
                 justify ='left'
                 )
    no_rbtn.grid(sticky='w')
    
    w.mainloop()
    ## Menu tipo windows ---------------------------------

def gui_console():
    
    link_menu = ui_menu(
        options = list(__Result_Labels_Links.keys()),
        label ='link_menu',
        sel_msg = 'Select the type of result to post-process',
        is_main = True
        )
    yesno_set_menu = ui_menu(
        options = ['Yes','No'],
        label = 'yesno_set_menu',
        sel_msg = 'Filter by Sets File [Yes/No] : '
        )
    yesno_summary_menu = ui_menu(
        options = ['Yes','No'],
        label = 'yesno_summary_menu',
        sel_msg = 'Generate the Summary of all groups [Yes/No] : '
        )

    env_type_menu = ui_menu(
        options = list(ENVELOPE_OPTIONS),
        label = 'env_type_menu',
        sel_msg = 'Select the type of Envelope'
        )
    
    
   

    while True:
        # Asking for type of results
        link_menu.loop()
        if link_menu.selection.upper() == 'EXIT':
            break
        else:
            result_link = getattr(r_link, link_menu.selection)
            
            # Asking for OP2 File
            op2 = ask_op2_filepath()
            if op2 is None : 
                print ('No OP2 file selected. Returning to main menu.')
                continue
            
            out = asksaveasexcelfile()
            if out[0] == '' : 
                print ('No Excel file selected. Returning to main menu.')
                continue
            
            # Asking if filter by Sets File
            yesno_set_menu.loop()
            
            if yesno_set_menu.selection.upper() == 'YES':
                sets_file = ask_set_filepath()
                if sets_file is None:
                    continue
                else:
                    sets = read_sets_from_file(
                        filename = sets_file
                        )
                    print(sets)
            else:
                # sets = gen_all_group(
                #     op2_path = op2,
                #     entity_type = result_link.entity_type
                #     )
                sets = get_result_ids_group(
                    op2_path = op2,
                    result_link = result_link
                    )
                print(sets)
            # Asking for Envelope option
            env_type_menu.loop()
            
            # Askig for field 
            if env_type_menu.selection.upper() != 'NONE':
                print (f'The following fields are available for'
                       f'{link_menu.selection}:\n{result_link.op2_results}')
                fields = ask_entities_fields_list()
            else:
                fields = None
            
            # Asking for summary
            yesno_summary_menu.loop()
            summary = True if yesno_summary_menu.selection.upper() == 'YES' else False
            
            # as
            _= group_envelope(
                dataframe = raw_table(op2,result_link),
                groups = sets,
                result_link = result_link,
                envelope_type = env_type_menu.selection,
                envelope_columns = fields,
                output_path = out[0],
                debug = False,
                groups_summary = summary
                            )
if __name__ == '__main__':
    gui_console()
    # gui_windows()