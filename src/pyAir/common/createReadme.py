# -*- coding: utf-8 -*-
"""
Created on Wed Dec 10 09:14:02 2025

@author: U69432
"""
import os
import time
from dataclasses import dataclass
from general import default_output_file as dof
# Creador de comentarios de texto plano


@dataclass
class Encodings:

    __UTF8 = 'utf-8'
    __UTF16 = 'utf-16'
    __UTF32 = 'utf-32'
    __ASCII = 'ascii'
    __ISO88591 = 'iso-8859-1'


ENCODINGS = [
    'utf-8',
    'utf-16',
    'utf-32',
    'ascii',
    'iso-8859-1',
]


def str_encoding(string: str, verbose: bool = False) -> list[str]:
    '''
    Detects if the current string has one of the following encodings:
        1. UTF-8
        2. UTF-16
        3. UTF-32
        4. ASCII
        5. ISO-8859-1

    Parameters
    ----------
    string : str
        DESCRIPTION.
    verbose : bool, optional
        If True shows additional info. The default is False.

    Raises
    ------
    IOError
        If None of the encodings is detected

    Returns
    -------
    List with with all the allowed encoding.

    '''

    allow_encodes = []
    for encode in ENCODINGS:
        try:
            _ = string.encode(encode)
            if verbose:
                print(f'String is compatible with {encode}')
            allow_encodes.append(encode)
        except:
            if verbose:
                print(f'Is not {encode}')

    if allow_encodes == []:
        raise IOError('Impossible to detect any encoding')
    else:
        return allow_encodes


def fix_encoding(string: str, encode: str, decode: str) -> str:
    allow_encodes = str_encoding(string)
    if encode in allow_encodes and decode in allow_encodes:
        fixed_string = string.encode(encode).decode(decode)
        return fixed_string
    else:
        raise ValueError(
            f'Invalid encode or decode value. The allowed encodes for the'
            f'string are :\n{allow_encodes}'
        )


def auto_fix_encoding(string: str, encoding: str = 'utf-8') -> str:
    '''
     Function that fixes automatically any encoding error in plain text files.

     Parameters
     ----------
     string : str
         Text string.
     encoding : str, optional
         Any allowed string representing the encoding. The default is 'utf-8'.

     Raises
     ------
     ValueError
         When the string has not been encoded using any of the most common
         one used for plain text files.

     Returns
     -------
        Returns the same string but with the requested codification.

     '''

    allow_encodes = str_encoding(string)
    if len(allow_encodes) > 2 and encoding in allow_encodes:
        decode = encoding
        allow_encodes.remove(encoding)
        for allow_encode in allow_encodes:
            try:
                fixed_string = string.encode(allow_encode).decode(decode)
                return fixed_string
            except:
                pass
    else:
        raise ValueError(
            f'Invalid encode or decode value. The allowed encodes for the'
            f'string are :\n{allow_encodes}'
        )


def header_lines(
        author: str = None,
        date: str = None,
        tab_str: str = '') -> list[str]:
    '''
    Generates a header for any plain text file which includes the any combina-
    nation of the following data:
        1.- Author
        2.- Date

    Parameters
    ----------
    author : str, optional
        The default is None.
    date : str, optional
         The default is None.
    tab_str : str, optional
        String used to tabulate the data from the left margin
        of the document. The default is ''.

    Returns
    -------
    list[str]
        DESCRIPTION.

    '''

    header_lines = []
    if author is not None:
        author_line = f'{tab_str}AUTHOR : {author}\n'
        header_lines.append(author_line)
    if date is not None:
        date_line = f'{tab_str}DATE : {date}\n'
        header_lines.append(date_line)

    if header_lines == []:
        return None
    else:
        return header_lines


def format_line(line: str, line_width: int = 80) -> str:
    nchars = len(line)
    line_width = line_width - 4
    line_space = line_width - nchars

    fillspace = line_space * ' ' + ' |\n'
    formated_line = '| ' + line.replace('\n', fillspace)
    return formated_line


def format_README(
    path: str,
    author: str = '',
    date: str = None,
    line_width: int = 80,
    expandtabs: bool = True
):
    '''
    Given the path of a plain text file, it returns the same document formatted
    in such a way that it can seen a column-like document including a header.

    Example of input: 
        This a simple example
        of a two-line document.

    Result:
        -->  | Author: E.Maroto        |
             | Date : 16/12/2025       |
             |                         |
             | This a simple example   |
             | of a two-line document. |
             |                         |

    Parameters
    ----------
    path : path-like str
        Path of the plain text file.
    author : str, optional
        The default is ''.
    date : str, optional
        The default is None.
    line_width : int, optional
        Width in terms of characters of the document. The default is 80.
    expandtabs : bool, optional
        If true it replaces the tabs by four spaces. The default is True.

    '''
    # Formatting time
    if date is None:
        date = time.strftime('%y\%d\%m')

    # Checking if is a new file or not
    if os.path.isfile(path):  # If True
        with open(path, 'r') as file:
            lines = file.readlines()
        lines = header_lines(author, date) + lines
        formated_lines = []
        for line in lines:
            if expandtabs:
                line.expandtabs(tabsize=4)
            # Fixing encode
            fixed_line = auto_fix_encoding(line)
            # Formating line
            formated_line = format_line(fixed_line, line_width)
            # Storing in list
            formated_lines.append(formated_line)

        # Export
        outpath = dof(input_file=path, suffix='_formatted')
        with open(outpath, 'w') as file:
            file.writelines(formated_lines)
        return formated_lines
    else:
        pass
