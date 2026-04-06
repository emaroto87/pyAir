
import os
import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename

__version__ = '0.0.2'
__author__ = 'E.Maroto'


def askOpenBDF(verbose=False, multiple=False) -> str:
    '''
    Opens a UI dialog window to open a NASTRAN bulk data file (bdf). If the 
    user press the Cancel button, returns a None type object.

    Returns
    -------
    path : path-like string


    '''
    # Code to force the dialog window to be the foremost
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    if verbose:
        print('Opening NASTRAN BDF file...', end='')

    path = askopenfilename(
        title='Open Sets File',
        filetypes=(("BDF file", '*.bdf'), ("DAT file", '*.dat'),
                   ("INCLUDE file", '*.incl')),
        multiple=multiple,
        defaultextension=['.bdf'],
        initialdir=os.getcwd(),
        parent=root
    )
    if path == '':
        print('No file has been selected.')
        path = None
    else:
        print('Done')
    return path


def askOpenGroupFile(verbose=False, multiple=False) -> str:
    '''
    Opens a UI dialog window to open a Sets file. If the user press the Cancel
    button, returns a None type object.

    Returns
    -------
    path : path-like string


    '''
    # Code to force the dialog window to be the foremost
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    if verbose:
        print('Opening Sets file...', end='')

    path = askopenfilename(
        title='Open Sets File',
        filetypes=(("Sets file", '*.set'),),
        multiple=multiple,
        defaultextension=['.op2'],
        initialdir=os.getcwd(),
        parent=root
    )
    if path == '':
        print('No file has been selected.')
        path = None
    else:
        print('Done')
    return path


def askOpenOP2(verbose=False, multiple=False) -> str:
    '''
    Opens a UI dialog window to open a OP2 file. If the user press the Cancel
    button, returns a None type object.

    Returns
    -------
    path : path-like string


    '''
    # Code to force the dialog window to be the foremost
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    if verbose:
        print('Opening OP2 file...', end='')

    path = askopenfilename(
        title='Open OP2',
        filetypes=(("Nastran OP2 file", '*.op2'),),
        multiple=multiple,
        defaultextension=['.op2'],
        initialdir=os.getcwd(),
        parent=root
    )
    if path == '':
        print('No file has been selected.')
        path = None
    else:
        print('Done')
    return path


def askSaveAsExcel(verbose=False, multiple=False) -> str:
    '''
    Opens a UI dialog window to save a Excel file. If the user press the Cancel
    button, returns a None type object.

    Returns
    -------
    path : path-like string


    '''
    # Code to force the dialog window to be the foremost
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    if verbose:
        print('Saving Excel file...', end='')

    path = asksaveasfilename(
        title='Save EXCEL file...',
        filetypes=(("EXCEL file", '*.xlsx'),),
        defaultextension=['*.xlsx'],
        initialdir=os.getcwd(),
        parent=root
    ),

    if path == '':
        print('No file has been selected.')
        path = None
    else:
        print('Done')
    return path


def askSaveAsBDF(verbose=False, multiple=False) -> str:
    '''
    Opens a UI dialog window to save a NASTRAN BDF file. If the user press the 
    Cancel button, returns a None type object.

    Returns
    -------
    path : path-like string


    '''
    file_type = 'BDF'

    # Code to force the dialog window to be the foremost
    root = tk.Tk()
    root.withdraw()
    root.wm_attributes('-topmost', 1)

    if verbose:
        print('Saving BDF file...', end='')

    path = asksaveasfilename(
        title=f'Save {file_type} file...',
        filetypes=((f"{file_type} file", '*.bdf'),),
        defaultextension=['*.bdf'],
        initialdir=os.getcwd(),
        parent=root
    ),

    if path == '':
        print('No file has been selected.')
        path = None
    else:
        print('Done')
    return path
