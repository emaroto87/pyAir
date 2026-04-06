# -*- coding: utf-8 -*-
"""
Created on Tue Dec  9 10:16:50 2025

@author: U69432
"""
from collections import defaultdict
from collections import namedtuple
from common.tui import TUI

# options = namedtuple('options', 'name', 'description')

opts = {
    '0': 'Get Subcases Params',
    '1': 'OP2 Result reader',
    '2': 'Get EigenValues from Buckling Analyses'
}


def tui():
    main_menu = TUI(
        options=list(opts.values()),
        sel_msg='Select one option',
        label='PosTool root',
        is_main=True,
        return_opt=True
    )

    main_menu.loop()

    if main_menu.selection == opts['0']:
        import get_subcases
    if main_menu.selection == opts['1']:
        import NastranOP2Reader
    if main_menu.selection == opts['2']:
        import get_buckling_eigen
    return main_menu.selection
