# -*- coding: utf-8 -*-
"""
Created on Mon May 26 15:03:33 2025

@author: U69432
"""
__version__ = '0.2'
import os

import numpy as np
import pandas as pd
from Common.general import load_op2_file
from Common.UI.filedialog import askOpenOP2
from Common.UI.filedialog import askSaveAsExcel
from Common.IO.Exporting.toExcel import to_excel


from collections import defaultdict


def get_buckling_eigen(op2_files_path: list[str]) -> 'dict{sid: dataframes}':
    '''
    Generates a dictionary of tables from a list of NASTRAN OP2 files
    containing buckling analysis SOL105.

    Parameters
    ----------
    op2_files_path : list[str]
        List of path of the NASTRAN OP2 files.

    Returns
    -------
    eigrs_tables : dict
        Dictionary in which the keys are the OP2 filenames and the values
        pandas DataFrame containing the eigenvalues for each buckling subcase

    '''
    eigrs_tables = {}
    for op2_file in op2_files_path:
        fn = os.path.basename(op2_file)
        try:
            print('Opening file {0}'.format(fn), end='')
            op2 = load_op2_file(op2_file)
            print('[Done]')
        except:
            print('[Warning]. Unable to open {0} file'.format(fn))
            continue

        eigenvalues_dict = defaultdict(int)

        for sid, eigenvector in op2.eigenvectors.items():

            # Get the eigenvalues
            eigenvalues = np.array(eigenvector.eigrs)
            n_eigrvalues = len(eigenvector.eigrs)

            # Computing the minimun positive
            minpos_eigen_value = np.where(
                eigenvalues > 0, eigenvalues, np.inf).min()

            values = eigenvector.eigrs + [minpos_eigen_value]
            eigenvalues_dict[sid] = values

        # Generating the labels of the columns
        first_sid = list(op2.eigenvectors.keys())[0]
        n_eigrvalues = op2.eigenvectors[first_sid].ntimes
        columns = ['eigrval_{0}'.format(eigrval_index+1)
                   for eigrval_index in range(n_eigrvalues)]
        columns += ['min_pos']

        # Generating the table
        eigrs_table = pd.DataFrame(data=eigenvalues_dict.values(),
                                   index=eigenvalues_dict.keys(),
                                   columns=columns)
        eigrs_tables[fn] = eigrs_table

    return eigrs_tables


if __name__ == '__main__':
    # OpenFile dialog window asking for OP2 file path
    op2_paths = askOpenOP2(multiple=True)
    # Retrieving the Eigenvalues Dataframes
    eigrs_tables = get_buckling_eigen(op2_paths)
    # Asking output Excel file path
    excel_path = askSaveAsExcel()[0]
    # Generating excel
    to_excel(eigrs_tables, excel_path, False, True)
