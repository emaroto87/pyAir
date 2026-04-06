import pandas as pd
import os


# PANDAS DATAFRAME FROM EXCEL SHEET
def to_excel(
        tables: dict[str, pd.DataFrame],
        output_filename:  str,
        remove_prefix: bool = False,
        index: bool = False,
        verbose: bool = False
) -> "Excel File":
    """
    Export into an Excel file multiple pandas DataFrames stored within a dictionary.

    Parameters
    ----------
    tables : Dict [str, pd.DataFrame ]
        Dictionary containing the pandas DataFrames. Each key of the dictionary
        is used to label the sheet of the Excel book file.
    filename : str path-like
        Path of the output Excel file.

    """
    # Remove common prefix
    if remove_prefix is True:
        labels = list(tables.keys())
        c = os.path.commonprefix(labels)

    print('Exporting to Excel File ...', sep='')
    with pd.ExcelWriter(output_filename) as writer:
        for sheet_name, table in tables.items():
            if remove_prefix is True:
                sheet_name = sheet_name.replace(c, '')
            else:
                pass
            table.to_excel(writer, sheet_name=sheet_name, index=index)

    print('[Done]')
