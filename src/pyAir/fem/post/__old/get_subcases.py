# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 09:06:58 2025

@author: U69432
"""


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
