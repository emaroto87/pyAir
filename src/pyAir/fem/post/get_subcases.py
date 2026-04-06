# -*- coding: utf-8 -*-
from filedialog import askOpenOP2
from filedialog import askSaveAsExcel
from general import to_excel
from general import load_op2_file
import pandas as pd

__version__ = '0.0.2'


def get_subcases_names(pyNastran_op2, param):
    '''


    Parameters
    ----------
    pyNastran_op2 : PyNastran OP2 object

    param : str
        One of the available params within a subcase control deck card.
        Examples:
            - TITLE
            - LABEL

    Returns
    -------
    names_dict : Dict{sid:param}
        Dictionary whose keys is the SIDs of the subcacases and their values
        the assoicated parameter.

    '''
    # Retrieving Subcases from OP2 object
    op2 = pyNastran_op2
    subcases_dict = op2.case_control_deck.subcases

    # Deleting default subcase 0 created by pyNastran
    try:
        subcases_dict.pop(0)
    except ValueError:
        print('Error: while poping the dummy subcase')

    # Generating the output dictionary
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


def export_subcases_names(pyNastran_op2, param, xls_path):
    '''
    Exports into an excel the data relatived to the subcases control deck
    cards.

    Parameters
    ----------
    pyNastran_op2 : PyNastran OP2 object

    param : str
        One of the available params within a subcase control deck card.
        Examples:
            - TITLE
            - LABEL
    xls_path : str path-like
        Path of the Excel file in which store the information. The path must
        end with the proper '*.xlsx' extension. Otherwise, the fucntion will
        fail.

    Returns
    -------
    Dictionary whose key is the parameter selected and the value is a dataframe
    table with the first colum the SID of the subacases and the second column
    is the value of the parameter.

    '''
    names_dict = get_subcases_names(pyNastran_op2, param)
    x = names_dict.keys()
    y = names_dict.values()
    df = pd.DataFrame(y, x, columns=[param])
    df.index.name = 'Subcase ID'
    dict_dfs = {param: df}
    to_excel(
        tables=dict_dfs,
        output_filename=xls_path,
        index=True
    )
    return dict_dfs


# if __name__ == '__main__':
print(__name__)
# Asking for OP2 path
op2path = askOpenOP2(multiple=False)
# Opening the OP2 file
op2 = load_op2_file(op2path, True)
# Asking for the para
param = input(
    'Type the param that want to extract from the subaceses cards:\n')
# Asking output xls file path
xlspath = askSaveAsExcel()[0]
# Getting the subcases data into dict
_ = export_subcases_names(op2, param, xls_path=xlspath)
