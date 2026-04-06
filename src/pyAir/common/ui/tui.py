# -*- coding: utf-8 -*-

from common.errors import errors as e
import os

'''
TEXT-BASED USER INTERFACE MODULE
'''


class TUI:

    def __init__(self,
                 options: list[str],
                 sel_msg: str = None,
                 label: str = None,
                 options_sort: bool = False,
                 is_main: bool = False,
                 return_opt: bool = False):
        '''


        Parameters
        ----------
        options : list[str]
            List of the options that will be displayed in the menu.
        sel_msg : str, optional
            Mensage that will be shown before the list of options
            The default is None.
        label : str, optional
            String that identifies the object from other instances so 
            might be callable afterwards. The default is None.
        options_sort : bool, optional
            Sort the options. The default is False.

        is_main : bool, optional
            Identitifes if the current menu is the main or root menu
            The default is False.
        return_opt : bool, optional
            If True, includes among option 'Exit' in the menu. The default is False.

        Returns
        -------
        TUI
            DESCRIPTION.

        '''
        self.__selection_key = None
        self.__selection_val = None
        self.selection = None
        self.__text = 'Select an option: '
        self.__error_msgs = e.error_messages()
        self.label = label
        self.prev = None
        self.go_to = None
        self.__return_opt = return_opt

        if not (type(is_main) is bool):
            raise e.typerror_msg(is_main, bool)
        else:
            self.is_main = is_main
            self.__return_option_text = 'Back' if is_main is False else 'Exit'
            if self.is_main is True:
                self.__return_opt = True

        if (not (sel_msg is None) and (not (type(sel_msg) is str))):
            raise e.typerror_msg(sel_msg, str)
        else:
            self.sel_msg = sel_msg

        if not (type(options) is list):
            raise e.typerror_msg(options, list)
        else:
            if options_sort:
                options.sort()
            if self.__return_opt:
                options += [self.__return_option_text]
            self.options = dict(enumerate(options))

    def show_options(self):
        if not (self.sel_msg is None):
            print(self.sel_msg)
        print(self.__text)
        for (i, option) in self.options.items():
            print('\t', str(i), ':\t', option)

    def get_user_selection(self):

        user_inp = input('Enter an option : ')
        try:
            user_inp = int(user_inp)
            if not (int(user_inp) in self.options.keys()):
                print(self.__error_msgs.invalid_option)
            else:
                self.__selection_key = int(user_inp)
                self.__selection_value = self.options[user_inp]
                self.selection = self.options[user_inp]

        except:
            print(self.__error_msgs.invalid_option)
            self.selection = None

    def loop(self):
        while True:

            self.show_options()
            self.get_user_selection()
            if self.selection in list(self.options.values()):
                break
            else:
                self.clear_console()
        return self.selection

    # def clear_console(self, include_spyder=True):
    #     cmd = 'cls' if os.name == 'nt' else 'clear'
    #     os.system(cmd)
    #     get_ipython().magic('clear')


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
