from __future__ import annotations
from typing import Dict, Iterable, Tuple, Optional, Union
from pathlib import Path
from utils.utils import ensure_unique_name
from utils.utils import sanitize_name
import re
import pandas as pd
import os

SheetData = Union[Dict[str, pd.DataFrame], Iterable[Tuple[str, pd.DataFrame]]]
INVALID_SHEET_CHARS = r'[\[\]\*\?\/\\]'
MAX_SHEET_NAME_LEN = 31


def write_dataframes_to_excel(
        filepath: str,
        sheets: SheetData,
        include_index: bool = False,
        verbose: bool = False):
    """
    Generate a Excel file from a group of pandas Dataframes.

    Parameters
    ----------
    filepath : str
        Path of the Excel file.
    sheets : Union[Dict[str, pd.DataFrame], Iterable[Tuple(str, pd.DataFrame)]]
        Dictionary-like object {sheet_name : Pandas DataFrame}
    include_index : bool, optional
        Includes the index of the dataframe into the sheet. Default is False.
    verbose : bool, optional
        If True, prompts additional information. Default is False.

    Raises
    ------
    TypeError
        If the data is not a tuple(str,dataframe) like object.

    Returns
    -------
    Excel file.

    """

    if verbose:
        print('Exporting Excel file', end='')

    # # Normalizing inputs to list of tuples:
    # if isinstance(sheets, dict):
    #     items = list(sheets.items())
    # else:
    #     items = list(sheets)

    # Sanitising names and looking for duplicates
    normalized_items: list[tuple[str, pd.DataFrame]] = {}
    seen_input_names: set[str] = set()
    print(sheets)
    for name, df in sheets.items():
        if not isinstance(df, pd.DataFrame):
            raise TypeError(
                f"The value associated to {name} is not a DataFrame type.")
        else:
            clean = sanitize_name(
                name=name,
                invalid_chars=INVALID_SHEET_CHARS,
                max_length=MAX_SHEET_NAME_LEN)
            if clean in seen_input_names:
                clean = ensure_unique_name(name, seen_input_names)
            else:
                seen_input_names.add(clean)
            normalized_items[clean] = df

    print(normalized_items)
    with pd.ExcelWriter(filepath) as writer:
        for sheet_name, table in normalized_items.items():
            table.to_excel(
                writer,
                sheet_name=sheet_name,
                index=include_index)

    if verbose:
        print('Done')
