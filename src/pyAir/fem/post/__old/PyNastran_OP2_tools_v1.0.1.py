# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""

from tkinter.filedialog import asksaveasfilename
from tkinter.filedialog import askopenfilename
import tkinter as tk
import os
import sys
import pandas as pd
import xlwings as xw
from pyNastran.op2.op2 import OP2

from scipy.spatial import ConvexHull
from matplotlib import pyplot as plt
from Utilities import general as gui

# Abreviaturas
# lista = ls
# dataframe = df
# ids de nodos = nids
# id de nodo = nid
#
__version__ = '1.0.1'
__author__ = 'E.Maroto'

folder_path = r'C:\Users\U69432\Documents\PyNastran141\examples\fem'
incfn = 'fem.incl'
op2fn = 'fem.op2'
f06fn = 'fem.f06'

incfp = os.path.join(folder_path, incfn)
op2fp = os.path.join(folder_path, op2fn)
f06fp = os.path.join(folder_path, f06fn)


# Default Values

__ENTITY_TYPES = [
    'NodeID',
    'ElementID',
    'SubcaseID'
]

__CASE_CONTROL_PARAMS = [
    'TITLE',
    'LABEL'
]

############################### BEGIN AUX ###############################

# Clear console


def clear_console():
    cmd = 'cls' if os.name == 'nt' else 'clear'
    os.system(cmd)


# Open File


def openOP2():
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)
    print('Opening OP2 file...')
    op2fp = askopenfilename(
        title='Open OP2',
        filetypes=(("Nastran OP2 file", '*.op2'),),
        multiple=False,
        defaultextension=['.op2'],
        initialdir=os.getcwd(),
        parent=root
    )
    if op2fp == '':
        print('No file has been selected.')
    return op2fp


def openBDF():
    bdf_fp = askopenfilename(
        title='Open NASTRAN BDF file',
        filetypes=(
            ("Nastran BDF file", '*.bdf'),
            ("Data file", '*.dat'),
            ("Include file", '*.incl'),
        ),
        multiple=False,
        defaultextension=['*.bdf'],
        initialdir=os.getcwd()
    )
    if bdf_fp == '':
        print('No file has been selected.')
    return bdf_fp


# Exporting
def df_to_excel_sheet(dataframe, sheet_name=None):
    print('Exporting data to Excel...', end='')

    # Opening xlwings
    app = xw.App(visible=False)
    wb = app.books[0]
    ws = wb.sheets[0]
    ws.activate()

    # Exporting data to default sheet1
    ws["A1"].options(
        pd.DataFrame,
        header=1,
        index=True,
        expand='table').value = dataframe
    print('Done')

    # Re-labeling Default Sheet1
    if not (sheet_name is None):
        ws.name = sheet_name

    # Saving
    wb.save(
        path=asksaveasfilename(
            title='Save EXCEL file...',
            filetypes=(("EXCEL file", '*.xlsx'),),
            defaultextension=['*.xlsx'],
            initialdir=os.getcwd()
        ),
    )
    wb.close()
    print('Done')


def dfs_to_excel_sheets(dict_dataframes):
    print('Exporting data to Excel...', end='')

    # Opening xlwings
    app = xw.App(visible=False)
    wb = app.books[0]
    ws_labels = list(dict_dataframes.keys())

    # Exporting data to default the differen
    for label in ws_labels:
        wb.sheets.add(str(label))
        df = dict_dataframes[label]
        ws = wb.sheets[str(label)]
        ws.activate()
        ws["A1"].options(
            pd.DataFrame,
            header=1,
            index=True,
            expand='table').value = df

    # Deleting default Sheet
    wb.sheets[-1].delete()

    # Saving
    wb.save(
        path=asksaveasfilename(
            title='Save EXCEL file...',
            filetypes=(("EXCEL file", '*.xlsx'),),
            defaultextension=['*.xlsx'],
            initialdir=os.getcwd()
        ),
    )
    wb.close()
    print('Done')


def export_single(dataframe, option='csv'):
    assert (dataframe._typ == pd.DataFrame._typ)
    assert (type(option) is str)

    if option.upper() == 'CSV':
        dataframe.to_csv(
            asksaveasfilename(
                title='Save CSV file...',
                filetypes=(("CSV file", '*.csv'),),
                defaultextension=['.csv'],
                initialdir=os.getcwd()
            ),
            header=True
        )

    elif option.upper() == 'XLS':
        app = xw.App(visible=False)
        wb = app.books[0]
        ws = wb.sheets[0]
        ws.activate()
        print('Exporting data to Excel...', end='')
        ws["A1"].options(
            pd.DataFrame,
            header=1,
            index=True,
            expand='table').value = dataframe
        print('Done')
        wb.save(
            path=asksaveasfilename(
                title='Save EXCEL file...',
                filetypes=(("EXCEL file", '*.xlsx'),),
                defaultextension=['*.xlsx'],
                initialdir=os.getcwd()
            ),
        )
        wb.close()

    else:
        return


def export_multi(dict_dataframes):
    print('Exporting data to Excel...', end='')
    app = xw.App(visible=False)
    wb = app.books[0]
    ws_labels = list(dict_dataframes.keys())

    for label in ws_labels:
        wb.sheets.add(str(label))
        df = dict_dataframes[label]
        ws = wb.sheets[str(label)]
        ws.activate()
        ws["A1"].options(
            pd.DataFrame,
            header=1,
            index=False,
            expand='table').value = df
    print('Done')

    # Deleting default Sheet
    wb.sheets[-1].delete()
    # Saving
    wb.save(
        path=asksaveasfilename(
            title='Save EXCEL file...',
            filetypes=(("EXCEL file", '*.xlsx'),),
            defaultextension=['*.xlsx'],
            initialdir=os.getcwd()
        ),
    )
    wb.close()

############################### END AUX ###############################

# ----------------------------------------------------------------------------


def filter_df_by_entity(dataframe, entities, entity_type='NodeID'):
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
        print('Filtering data by : ', entity_type)
        df_filtered = df[df[entity_type].isin(entities)]
        return df_filtered
    else:
        return None

# ----------------------------------------------------------------------------


def parsing_ids_list(user_input: str, sep=','):
    ls_ids = user_input.split(sep=sep)
    for i in range(len(ls_ids)):
        item = ls_ids[i].strip()
        try:
            ls_ids[i] = int(item)
        except:
            print('Item ', item, ' is not an valid integer id value')
            return None
    return ls_ids


def parsing_list(user_input: str, sep=','):
    ls_in = user_input.split(sep=sep)
    ls_out = [item.strip() for item in ls_in]
    return ls_out

# ----------------------------------------------------------------------------


def ask_entities_ids_list():
    user_input = input('Enter list of entities IDs : ')
    return parsing_ids_list(user_input)

# ----------------------------------------------------------------------------


def ask_entities_fields_list():
    user_input = input('Enter list of fields : ')
    return parsing_list(user_input)
# ----------------------------------------------------------------------------


def select_option(l: list):
    for (k, v) in enumerate(l):
        print('\t', str(k), ':\t', v)
    option = input('Select an option : ')
    while True:
        try:
            selection = l[int(option)]
            break
        except:
            print(option, ' is not a valid selection.')
    return selection


def select_fields(l: list):
    print('Available fields:')
    for (k, v) in enumerate(l):
        print('\t', str(k), ':\t', v)
    selected_fields = ask_entities_fields_list()
    return selected_fields

# ----------------------------------------------------------------------------


def select_entity_type(func=select_option, entities_type=__ENTITY_TYPES):
    return func(entities_type)

# ----------------------------------------------------------------------------


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

    # Nota: Habria que comprobar que el df contiene la columna correspondiente
    if entity_type in ['NodeID', 'ElementID', 'SubcaseID']:
        dfs_per_entity_dict = {}
        print('Spliting data by : ', entity_type)
        for entity in entities:
            print(entity, '...', end='')
            df_per_entity = df[df[entity_type] == entity]
            dfs_per_entity_dict[entity] = df_per_entity
            print('Done')
        return dfs_per_entity_dict
    else:
        return None

# ----------------------------------------------------------------------------


def get_subcases_names(pyNastran_op2, param):
    op2 = pyNastran_op2
    subcases_dict = op2.case_control_deck.subcases
    # Deleting default subcase 0 created by pyNastran
    try:
        subcases_dict.pop(0)
    except:
        None

    names_dict = {}
    for k, v in subcases_dict.items():
        try:
            name = v.params[param][0]
        except:
            name = ''
        sid = k
        names_dict[sid] = name
    return names_dict

# ----------------------------------------------------------------------------


def export_subcases_names(pyNastran_op2, param):
    names_dict = get_subcases_names(pyNastran_op2, param)
    x = names_dict.keys()
    y = names_dict.values()
    df = pd.DataFrame(y, x, columns=[param])
    df.index.name = 'Subcase ID'
    dict_dfs = {param: df}
    dfs_to_excel_sheets(dict_dfs)


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


def get_envelope_from_df(
        dataframe,
        export_to_excel=False,
        plot=False
):

    # Asking for fields
    ls_columns = select_fields(dataframe.columns)
    print(ls_columns)
    print(type(ls_columns))
    data = dataframe[ls_columns]
    data.reset_index(inplace=True, drop=True)
    points = data.to_numpy()
    hull = ConvexHull(points)

    envelope = data[data.index.isin([hull.vertices])]

    if export_to_excel == True:
        dfs_dict = {
            'Data': dataframe,
            'Enveloe': envelope
        }
        dfs_to_excel_sheets(dfs_dict)

    if plot == True:
        plt.plot(points[:, 0], points[:, 1], 'o')
        for vertice in hull.vertices:
            plt.plot(points[vertice, 0], points[vertice, 1], 'k-')
        plt.show()

    return envelope

# ----------------------------------------------------------------------------


def gui_select():
    user_inp = ''


def gui_export(res):
    user_inp = ''
    while not (user_inp.upper().startswith('Q')):

        # Printing available options
        clear_console()
        print('Select a result to export:')
        results = list(res.result_names)
        results.sort()
        dict_results = dict(enumerate(results))
        for (i, result) in dict_results.items():
            print('\t', str(i), ':\t', result)

        # Asking user to select and option or exit
        user_inp = input('Enter an option or Q to exit: ')

        # Selects 'Q' to exit
        if user_inp.upper().startswith('Q'):
            break

        # Selects a value out of the available options
        elif not (int(user_inp) in dict_results.keys()):
            print('Incorret option. Please try again.')

        # Selects a value within the available options
        else:
            selected_res = dict_results[int(user_inp)]
            o = res.get_result(selected_res)
            sids = list(o.keys())

            # Concatenating all the result of the same type
            ls_dfs = []
            for sid in sids:
                df = o[sid].data_frame
                df['SubcaseID'] = sid
                ls_dfs.append(df)
            df = pd.concat(ls_dfs)

            # Exporting
            print('Exporting to Excel File. Select the exporting option:')

            while True:
                entity = select_entity_type()
                if entity in df.columns:
                    break
                else:
                    print('[WARN] The selection ', entity,
                          'is not availabe in ', results[user_inp])
                    a = input('Press Enter to continue')

            export_multi(
                split_df_by_entity(
                    dataframe=df,
                    entities=ask_entities_ids_list(),
                    entity_type=entity
                )
            )

            return df


def load_OP2_results():
    o = OP2()
    o.read_op2(openOP2(), build_dataframe=True)
    return o


# subcases = res.case_control_deck.subcases # devuevle un diccionario {k:sid, v:subcase_obj}
# sids = list(subcases.keys())
# s = subcases[sids[0]].params


# ##### 1. Creamos Dataframe de todas las SPC forces
# #
# # Inputs necesarios
# ls_nids_filt = [1900001,1900002,1900003,1900004]
# ls_sids_filt = [8355,8360]

 # Testing
 # Load op2


# o = load_OP2_results()

# entity_type_menu = gui.ui_menu(
#      options = __ENTITY_TYPES,
#      options_sort = False,
#      is_main = False,
#      label = 'entity_sel'
#      )

# results_menu = gui.ui_menu(
#      options = list(o.result_names),
#      options_sort = True,
#      is_main = False,
#      label = 'result_sel'
#      )

# main_menu_opts = [
#     'Identificadores de subcasos',
#     'Resultados'
#     ]
# main_menu = gui.ui_menu(
#      options = main_menu_opts,
#      options_sort = False,
#      is_main = True,
#      label = 'main'
#      )

# sel=''
# while True:
#     sel = main_menu.loop()
#     if sel == 'Identificadores de subcasos':
#         export_subcases_names(
#             pyNastran_op2 = o,
#             param = 'TITLE'
#             )
#     elif sel == 'Resultados':
#         results_menu.loop()
#         if not(results_menu.selection == 'Back'):
#             export_multi(
#                 split_df_by_entity(
#                     dataframe = get_result_from_op2(o,results_menu.selection),
#                     entity_type = entity_type_menu.loop(),
#                     entities = ask_entities_ids_list(),
#                     )
#                 )
#             sel=''
#         else:
#             sel=''
#     elif sel=='Exit':
#         break
#     else:
#         print('Ok')
#         break
