from typing import Optional, Dict
from pathlib import Path
import re


def read_sets_from_file(
        filepath: Path | str,
        comment_symbol: Optional[str] = '$',
        members_separator: Optional[str] = ',',
        set_keyword: Optional[str] = 'SET',
        verbose: bool = False,
        debug: bool = False,
) -> Dict[str, list]:
    """
    Read the sets defined within a text file.

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

    Returns
    -------
    dict

    """
    # Preliminary checks
    if not isinstance(filepath, str):
        raise TypeError(
            "The argument filepath is not a path-like string.")

    if not isinstance(comment_symbol, str):
        raise TypeError(
            "The argument comment_symbol is not a string.")

    if not isinstance(members_separator, str):
        raise TypeError(
            "The argument members_separator is not a string.")
    if not isinstance(members_separator, str):
        raise TypeError(
            "The set_keyword members_separator is not a string.")

    sets = {}
    # Pattern for multiple members
    pattern1 = r'{0}([\w-]+)=([\d,]+)'.format(set_keyword)

    # Pattern for single member
    pattern2 = r'{0}([\w-]+)=([\d]+)'.format(set_keyword)

    # Read file content and close
    with open(filepath, 'r') as f:
        lines = f.readlines()

    # cleaning commented lines
    clean_lines = [line[:-1] for line in lines
                   if not (line.startswith(comment_symbol))]

    clean_lines = ''.join(clean_lines)
    clean_lines = clean_lines.replace(' ', '')
    clean_lines = clean_lines.strip()

    if debug:
        print(clean_lines)

    # Spliting the content by sets using the pattern
    try:
        findings = re.findall(pattern1, clean_lines)
    except:
        if debug:
            print('Error in pattern1')
        try:
            findings = re.findall(pattern2, clean_lines)
        except:
            if debug:
                print('Error in pattern2')
            pass

    if len(findings) == 0:
        raise ValueError(
            'Unable to locate any set'
        )

    if debug:
        print(findings)

    # Read the members of each set and store them into a list
    for set_label, set_content in findings:
        try:
            members_list = [int(x) for x in set_content.split(',')]
        except:
            members_list = [int(set_content)]
        if set_label.isnumeric():
            set_label = int(set_label)
        sets[set_label] = members_list

    return sets
