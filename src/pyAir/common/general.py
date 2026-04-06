# -*- coding: utf-8 -*-
"""
Created on Fri May 23 09:26:07 2025

@author: U69432
"""


import os
import pandas as pd
import xlwings as xw
import numpy as np

from typing import Dict

from pyNastran.op2.op2 import OP2
from tqdm import tqdm
from common.ui.tui import select_option

# from errors import *

from collections import namedtuple

from pathlib import Path

__version__ = '1.0.1'
__author__ = 'E.Maroto'

# Default Values
entities_types = namedtuple('entities_type', ['node', 'element'])

__ENTITY_TYPES = [
    'NodeID',
    'ElementID',
    'SubcaseID'
]

__CASE_CONTROL_PARAMS = [
    'TITLE',
    'LABEL'
]


def load_op2_file(path: str,
                  build_dataframe: bool = True,
                  debug: bool = False) -> OP2:
    '''
    Given a NASTRAN OP2 file, returns the pyNASTRAN OP2 object.

    Parameters
    ----------
    path : str path-like
        Path of the NASTRAN OP2 file.
    build_dataframe : bool, optional
        If True, enables the option to store all the information as Pandas
        dataframes.
    debug : bool, optional
       If True, prints additional output that can aid in diagnosing issues
       during the reading of the OP2 file.

    Returns
    -------
    PyNastran OP2 object.
        DESCRIPTION.

    '''
    o = OP2(debug=False)
    o.read_op2(path,
               build_dataframe=build_dataframe,
               # skip_undefined_matrices = True
               )
    return o


def concat(*args, **kargs):
    s = ''
    for i in args:
        s += i
    return s


def vector_mod(p1, p2):
    v = p2-p1
    return np.sqrt(np.dot(v, v))


# WRITING CARDS
def split_list(lst, n):
    for i in range(0, len(lst), n):
        yield lst[i: i + n]


def merge_dfs(*args, **kwargs):

    df_list = [arg for arg in args]
    merged_df = pd.concat(df_list)

    return merged_df

# READING Progress bar (WIP)


def progress_OP2_bar(path):
    op2 = OP2()
    file_size = os.path.getsize(path)
    with open(path, 'rb') as file:
        with tqdm(file_size, unit='B', unit_scale=True, desc='Reading file') as bar:
            def read_progress(n):
                data = file.read(n)
                bar.update(len(data))
                return data
            op2.read_op2(path, read_binary=read_progress)


###############################################################################
#                                                                             #
#                              EXPORTING FUNCTIONS                            #
#                                                                             #
###############################################################################
#


# ----------------------------------------------------------------------------


def filter_df_by_entity(dataframe, entities, entity_type='NodeID', verbose=False):
    # check
    if not (type(dataframe) is type(pd.DataFrame())):
        error_msg = 'dataframe is not a DataFrame type'
        raise TypeError(error_msg)
        return None
    else:
        df = dataframe

    if not (type(entities) is list):
        error_msg = 'entities is not a list'
        raise TypeError(error_msg)
        return None

#   Nota: Habria que comprobar que el df contiene la columna correspondiente
    if entity_type in ['NodeID', 'ElementID', 'SubcaseID']:
        if verbose:
            print('Filtering data by : ', entity_type)
        if entity_type in df.columns:
            df_filtered = df[df[entity_type].isin(entities)]
        else:
            df_filtered = df[df.index.isin(entities)]
        return df_filtered
    else:
        raise ('{0} is not valid option'.format(entity_type))
        return None

# ----------------------------------------------------------------------------


# ----------------------------------------------------------------------------


def select_entity_type(func=select_option, entities_type=__ENTITY_TYPES):
    return func(entities_type)

# ----------------------------------------------------------------------------
# SPLITITNG DATAFRAMES


def is_type(obj, obj_class):
    if not (type(obj) is type(obj_class)):
        error_msg = 'Object {0} is not a {1}type'.format(
            str(type(obj)),
            str(type(obj_class))
        )
        raise TypeError(error_msg)
        return False
    else:
        return True


def split_df_by_entity(dataframe,
                       entities,
                       entity_type='NodeID',
                       filter_list=False,
                       ):
    # check
    if not (type(dataframe) is type(pd.DataFrame())):
        error_msg = 'dataframe is not a DataFrame type'
        raise TypeError(error_msg)
        return None
    else:
        df = dataframe

    if not (type(entities) is list):
        error_msg = 'entities is not a list'
        raise TypeError(error_msg)
        return None

    df = None if not (is_type(dataframe, pd.DataFrame())) else dataframe
    entities = None if not (is_type(entities, list())) else entities

    if ((df is None) or (entities is None)):
        return

    # Nota: Habria que comprobar que el df contiene la columna correspondiente

    if entity_type in df.columns:
        dfs_per_entity_dict = {}
        print('Spliting data by : ', entity_type)
        for entity in entities:
            print(entity, '...', end='')
            df_per_entity = df[df[entity_type] == entity]
            dfs_per_entity_dict[entity] = df_per_entity
            print('Done')
        return dfs_per_entity_dict
    elif entity_type == df.index.name:
        dfs_per_entity_dict = {}
        print('Spliting data by : ', entity_type)
        for entity in entities:
            print(entity, '...', end='')
            df_per_entity = df[df.index == entity]
            dfs_per_entity_dict[entity] = df_per_entity
            print('Done')
        return dfs_per_entity_dict
    else:
        return None


def split_df_by_subcases(dataframe, sids):
    df = None if not (is_type(dataframe, pd.DataFrame())) else dataframe
    dfs_per_subcases = {}
    for sid in sids:
        print('Splitting DataFrame by SID : {0}'.format(sid))
        df_per_subcase = df[df['SubcaseID'] == sid]
        dfs_per_subcases[sid] = df_per_subcase
    return dfs_per_subcases
# ----------------------------------------------------------------------------


def get_result_from_op2(res, result):
    # Accesing to the results
    o = res.get_result(result)

    # Getting all the subcases containing the selected result
    sids = list(o.keys())

    # Concatenating all the result of the same type
    ls_dfs = []
    for sid in sids:
        df = o[sid].data_frame
        df['SubcaseID'] = sid
        ls_dfs.append(df)
    df = pd.concat(ls_dfs)

    return df
# ----------------------------------------------------------------------------


def default_output_file(
        input_file: str,
        suffix: str | None = None,
        preffix: str | None = None
) -> str:
    '''
    *** PASS THIS FUNCTION TO GENERAL  

    Gets the path of a input file that is going to be modified and generates
    the path of the output file according to a Preffix and Suffix strings.
    If both Preffix and Suffix are None, the function adds "_mod" as suffix
    by default

    Parameters
    ----------
    input_file : file path

    suffix : str | None, optional
        String of the suffix. The default is None.
    preffix : str | None, optional
        String of the preffix. The default is None.

    Returns
    -------
    output_filepath : path
    '''

    folder = os.path.dirname(input_file)
    file = os.path.basename(input_file)
    filename, filextension = os.path.splitext(file)
    if (preffix is None) and (suffix is None):
        suffix = '_mod'
    if suffix is None:
        suffix = ''
    if preffix is None:
        preffix = ''
    else:
        print('Error')

    output_filename = preffix + filename + suffix
    output_filepath = os.path.join(folder, output_filename + filextension)

    return output_filepath
