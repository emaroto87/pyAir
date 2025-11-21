# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 14:48:19 2025

@author: U69432
"""

__version = '1.0.2'
__author  = 'E.Maroto'

# ========================
#
# PARTE DE LECTURA DE BDF
#
# ========================
#### Ejemplos para leer BDF y conseguir datos de él

from pyNastran.bdf.bdf import BDF as bdf
from dataclasses import dataclass, asdict
from typing import Optional
import abc
from enum import Enum
from pathlib import Path

import tkinter as tk
from tkinter.filedialog import askopenfilename
from tkinter.filedialog import asksaveasfilename

import os 
# # Load the model
# model = bdf()
# model.read_bdf( 
#     bdf_filename = incfp,
#     read_includes = True,
#     validate = False,
#     xref = True,
#     )

# # Imprimir datos del modelo 
# print(model.get_bdf_stats)

# # Acceder a nodos
# nodes = model.nodes # Se trata de un diccionario k=nid v=obj nodo
# nids = list(model.nodes.keys())
# n1=nodes[nids[0]]

# # Acceder a elementos
# elems = model.elements
# eids = list(elems.keys())
# e1 = elems[eids[0]]
            
# print(e1.get_stats()) # Proporciona información del elementoe1

# # cosas utiles
# e1.type # me proporciona el tipo de elemento que es
# e1.object_attributes() # proporciona los atributos
# e1.object_method()
# e1.pid_ref # accedes a la propiedad vinculada al elemento

# # Mapear
# model.get_property_id_to_element_ids_map()
# #
# #
# #
# #
# #
# #
# # Lista de funciones que podrian ser utiles.
# # 1. Dado un modelo con remaches CBUSH + RBE3, generar los PBUSH
# # 2. Dada una lista de remaches 


from collections import defaultdict
from collections import namedtuple
from collections import Counter
import os
import sys

"""
Nota el siguiente codigo asume que el elemento 'eid' es una lista que contiene
los nodos que lo definen.


"""
def get_edges_from_element(pyNastran_elem)->list:
    elem = pyNastran_elem
    edges=[]
    nodes = elem.node_ids 
    n = len(nodes)
    for i in range(n):
        print
        edge = tuple(sorted((nodes[i],nodes[(i+1) % n])))
        edges.append(edge)
    return edges

def get_edges(pyNastran_elems,sort=False,verbose=False):
    elems = pyNastran_elems
    edge_counter = defaultdict(int)
    # Step 1 : count how many times an edge appears
    
    for eid,elem in elems.items():
        edges = get_edges_from_element(elem)
        for edge in edges:
            edge_counter[edge] +=1
    if verbose :
        for edge,count in edge_counter.items():
            print ('Edge defined by nodes : ',edge, 'counted ',count,' times')
   
    # Step 2 : edges that appears only one time
    frontier_edges = [edge for edge,count in edge_counter.items() if count==1]
    
    # Step 3 : get the unique nodes of the frontier edges
    frontier_nodes = set()    
    for edge in frontier_edges:
        frontier_nodes.update(edge)
    if sort is True:
        frontier_nodes = sorted(frontier_nodes)
    
    return frontier_nodes

def testing_get_edges():
    # # Load the model
    path = r'C:\Users\U69432\Desktop\WORK\00_Resources\02_Python_codes\PyNastran\examples_and_tests\BDF_tools\edges'
    fname = 'test.bdf'
    os.chdir(path)
    model = bdf()
    model.read_bdf( 
        bdf_filename = os.path.join(path,fname),
        read_includes = True,
        validate = False,
        xref = True,
        )
    edges = get_edges(model.elements)
    return edges

def read_bdf(str):
    with open(str,'r') as f:
        lines = f.readlines()
    return lines

def nastran_format_str(string ,align: str ='right', field_chars = 8 ):
    inp = string
    nc = len(inp)
    spaces = (field_chars-nc)*' '
    if len(spaces)>0:
        if align == 'left':
            out = inp + spaces
        else:
            out = spaces + inp
    elif len(spaces) == 0 :
        out = inp
    else:
        raise ValueError('The length of {string} is greater than {field_chars} characters')
    
    return out

def nastran_format_int( integer , align='right', field_chars = 8 ):
    inp = str(integer)
    nc = len(inp)
    spaces = (field_chars-nc)*' '
    if len(spaces)>0:
        if align == 'left':
            out = inp + spaces
        else:
            out = spaces + inp
    elif len(spaces) == 0 :
        out = inp
    else:
        raise ValueError('The length of {string} is greater than {field_chars} characters')
    return out


def nastran_format_float(number, format_str='8.2E', align ='right', field_chars = 8):
    inp = format(number,format_str)
    nc = len(inp)
    spaces = (field_chars-nc)*' '
    if len(spaces)>0:
        if align == 'left':
            out = inp + spaces
        else:
            out = spaces + inp
    elif len(spaces) == 0 :
        out = inp
    else:
        raise ValueError('The length of {string} is greater than {field_chars} characters')
    return out

# @dataclass
# class PSHELL:
#     PID : int
#     T : float
    

@dataclass
class CQUAD4:
    EID : int
    PID : int
    G1 : int
    G2 : int
    G3 : int
    G4 : int
    THETA : float = 0.0
    MCID : int | None = 0
    ZOFFS : float | None = None
    TFLAG : int | None = None
    T1 : float = 0.0
    T2 : float = 0.0 
    T3 : float = 0.0
    T4 : float = 0.0
    bdf_file : str | None = None
    bdf_nline : int | None = None
    
    def __post_init__(self):
        self.entity = 'CQUAD4'
    
        card_field_line1 = [
            self.entity,
            self.EID,
            self.PID,
            self.G1,
            self.G2,
            self.G3,
            self.G4,
            self.MCID,
            self.ZOFFS,
            None,
            None
            ]
        
        card_field_line2 = [
            None,
            None,
            self.TFLAG,
            self.T1,
            self.T2,
            self.T3
            ]
        
        self.card_fields = card_field_line1 + card_field_line2
        
        __fields_order =[
            self.entity,
            self.EID,
            self.PID,
            self.G1,
            self.G2,
            self.G3,
            self.G4,
            self.MCID
            
            ]    
    
        
        
        
@dataclass
class CTRIA3:
    EID : int
    PID : int
    G1 : int
    G2 : int
    G3 : int
    THETA : float = 0.0
    MCID : int | None = None
    ZOFFS : float | None = None
    TFLAG : int | None = None
    T1 : float = 0.0
    T2 : float = 0.0 
    T3 : float = 0.0
    T4 : float = 0.0
    bdf_file : str | None = None
    bdf_nline : int | None = None
    
    def __post_init__(self):
        self.entity = 'CQUAD4'   
        pass
    
        card_field_line1 = [
            self.entity,
            self.EID,
            self.PID,
            self.G1,
            self.G2,
            self.G3,
            self.MCID,
            self.ZOFFS,
            None,
            None
            ]
        
        card_field_line2 = [
            None,
            None,
            self.TFLAG,
            self.T1,
            self.T2,
            self.T3
            ]
        
        self.card = card_field_line1 + card_field_line2


   

def create_card(card_fields : list ,field_chars : int = 8) -> str:
    card = ''
    field_counter = 1
    
    for item in card_fields :
        if type(item) is str:
            field = nastran_format_str(item, field_chars = field_chars)
        elif type(item) is int:
            field = nastran_format_int(item, field_chars = field_chars)
        elif type(item) is float:
            field = nastran_format_float(item, field_chars = field_chars)
        elif item is None:
            field = ' ' * field_chars
        else:
            raise ValueError
        card += field
        print(card)
    
        if field_counter % 10 == 0:
            print('enter')
            card += '\n'
        field_counter += 1
    
    return card

def parse_bdf_line(bdf_line : str, field_chars : int = 8):
    # Skip if comment line (Starts with $)
    if bdf_line.startswith('$'):
        return bdf_line
    
    # Parsing line

    else:
        fields = []
        # Couting number of fields
        n = len(bdf_line)
        mod, res  = divmod(n,field_chars)
        
        for i in range(mod):
           start = field_chars*i
           end = field_chars*(i+1)
           fields.append(bdf_line[start:end])
        
        if res !=0:
            start = field_chars*(mod)
            fields.append(bdf_line[start:])
        return fields
                        
        
        
def read_bdf(str):
     with open(str,'r') as f:
         lines = f.readlines()
         
     # Estrategia linea por linea
     
     
     return lines            


def change_2D_prop(lines , elems: list , PID : int, field_nchars=8):
    # with open(bdf_file) as f:
    #     lines = f.readlines()
    inp_lines = lines
    out_lines = []    
    for line in inp_lines:
        if line.startswith(('CQUAD4','CTRIA3')):
            EID = extract_field(line, 2)      
            if int(EID) in elems:
                new_PID = nastran_format_int(PID)
                old_PID = extract_field(line, 3)
                out_line = line.replace(old_PID, new_PID)
                out_lines.append(out_line)
            else:
                out_lines.append(line)
        else:
            out_lines.append(line)
    return out_lines


def extract_field(line,field_number,fields_chars = 8):
    start = (field_number-1) * fields_chars
    end = (field_number) * fields_chars
    field = line[start:end]
    return field

def openBDF():
    bdf_fp = askopenfilename(
        title ='Open NASTRAN BDF file', 
        filetypes=(
            ("Nastran BDF file", '*.bdf'),
            ("Data file", '*.dat'),
            ("Include file", '*.incl'),
            ),
        multiple = False,
        defaultextension = ['*.bdf'],
        initialdir = os.getcwd()
        )
    if bdf_fp =='':
        print('No file has been selected.')
    return bdf_fp


def mod_2Delems_PID(inp_bdf,elems_list : list, PID: int, field_nchars=8):
    out_bdf = inp_bdf.split('.')[0] + '_mod.' + inp_bdf.split('.')[1]
    
    with open(inp_bdf,'r') as f:
        lines = f.readlines()
        
    with open(out_bdf,'w') as f:
        out_lines = change_2D_prop(lines, elems_list, PID)
        f.writelines(out_lines)
    
    return out_bdf

def align_nodes(lines : list, nodes_list: list, node_ref :int , comp :int ):
    out_lines=[]
    for line in lines:
        if line.startswith('GRID'):
            nid = extract_field(line, 2)
            if int(nid) == node_ref:
                    new_coord = extract_field(line, 3 + comp)
       
    for line in lines:
        if line.startswith('GRID'):
            nid = extract_field(line, 2)
            if int(nid) in nodes_list:
                old_coord = extract_field(line, 3 + comp)
                out_line = line.replace(old_coord, new_coord)
                out_lines.append(out_line)
            else:
                out_lines.append(line)
        else:
            out_lines.append(line)
    
    return out_lines

def mod_node_coord(inp_bdf,nodes_list: list, node_ref :int , comp :int ):
    out_bdf = inp_bdf.split('.')[0] + '_mod.' + inp_bdf.split('.')[1]
    
    with open(inp_bdf,'r') as f:
        lines = f.readlines()
        
    with open(out_bdf,'w') as f:
        out_lines = align_nodes(lines,nodes_list,node_ref,comp)
        f.writelines(out_lines)
    
    return out_bdf


def output_file(file, suffix : str | None = None , preffix : str | None = None):
    '''
    *** PASS THIS FUNCTION TO GENERAL  
    
    Gets the path of a input file that is going to be modified and generates
    the path of the output file according to a Preffix and Suffix strings.
    If both Preffix and Suffix are None, the function adds "_mod" as suffix
    by default

    Parameters
    ----------
    file : file path
        
    suffix : str | None, optional
        String of the suffix. The default is None.
    preffix : str | None, optional
        String of the preffix. The default is None.

    Returns
    -------
    output_filepath : path
    '''
    
    
    folder = os.path.dirname(file)
    file = os.path.basename(file)
    filename, filextension = os.path.splitext(file)
    if (preffix is None) and (suffix is None):
        suffix = '_mod'
    if suffix is None : 
        suffix = ''
    if preffix is None : 
        preffix = ''
    else:
        print('Error')
    
    output_filename = preffix + filename + suffix
    output_filepath = os.path.join(folder,output_filename + filextension)
    
    return output_filepath
    
    
elemts_to_comment ={
    'CQUAD' : [10534025],
    'CBAR'  : [10734306],
    }    

def check_duplicated_ids(lists : list[list]):
    
    '''
    Given a list of list, returns the into a tuple of two list the unique
    elements and the duplicates

    Parameters
    ----------
    lists : list[list]
        List containing other lists 

    Returns
    -------
    unique,duplicates : tuple
    

    '''
    # Combinamos las listas
    combined = sum(lists,[])
    
    # La función counter te devuelve un diccionario donde cada key es cada
    # uno de los miembros de las lista, y el valor corresponde con el numero
    # de veces presente en la lista
    counter = Counter(combined)
    
    duplicates = [ item_id for item_id,freq in counter.items() if freq>1 ]
    
    # Lista sin duplicados
    unique = set(combined)
    
    return unique,duplicates
    
def comment_elements_by_field_pos(bdf_file, element_type_and_ids : dict, field_number : int = 1, check_duplicates = False):
    
    # Nastran fields number is between [1,10]
    if field_number<1 or field_number>10:
        raise ValueError(
            'The value of field number must be greater than 1 and do not' \
            'exceed 10.'
            ) 
    
    # Duplicates in the lists   
    element_types = element_type_and_ids.keys()
    unique_ids, duplicate_ids = check_duplicated_ids(element_type_and_ids.values())
    
    if check_duplicates:
        if len(duplicate_ids)>0:
            print('[WARNING] : The following ids are duplicated:')
            for item in duplicate_ids:
                print(item)
            yesno = input('Do you want to continue?')
      
            while True:   
                if yesno.upper().startswith('Y'):
                    break
                elif yesno.upper().startswith('N'):
                    return None
                else: 
                    pass
    
    # Reading BDF Files and      
    lines = read_bdf(bdf_file)
    nlines = len(lines)            
    
    # Comment the lines according to input
    print('---> BEGIN : Commenting elements')    
    i = 0
    out_lines = [] 
    while i < nlines :
        line = lines[i] 
        if line.startswith(tuple(element_types)):
            
            # Card Label
            card_str_field = extract_field(line, 1)
            
            # Card ID
            id_str_field = extract_field(line, field_number)
            id_int_field = int(id_str_field.strip())

            if id_int_field in unique_ids:
                print(f'\tCommented element {card_str_field } {id_int_field}')
                line = '$ Commented -> ' + line
                out_lines += [line]
                i +=1 
                while True : 
                    line = lines[i]
                    if line.startswith(('+')):
                        line = '$ Commented -> ' + line
                        out_lines += [line]
                        i += 1
                    else :
                        break
            else:
                out_lines += [line] 
                i += 1
                               
        else:
            out_lines += [line]
            i += 1
    
    # Writing the file with the commented lines
    out_fp = output_file(bdf_file)
    o = open(out_fp,'w')   
    o.writelines(out_lines)
    print('<--- END : Commenting elements')
    print(f'Created file at : {out_fp}')
    o.close()
    

# set_th2p2 =[
#     14303016,
#     14303057,
#     14303058,
#     14303059,
#     14303060,
#     14303069,
#     14303073,
#     14303074,
#     14303075,
#     14303076,
#     15303016,
#     15303057,
#     15303058,
#     15303059,
#     15303060,
#     15303069,
#     15303073,
#     15303074,
#     15303075,
#     15303076,
#     14303070,
#     14303072,
#     15303070,
#     15303072
#     ]  
# set_th2p5 =[
#     14300012,
#     14300021,
#     14300022,
#     14300024,
#     14300028,
#     14300030,
#     14300033,
#     14300041,
#     14300042,
#     14300053,
#     14300054,
#     15300012,
#     15300021,
#     15300022,
#     15300024,
#     15300028,
#     15300030,
#     15300033,
#     15300041,
#     15300042,
#     15300053,
#     15300054,
#     14300023,
#     15300023,
#     14300009,
#     14300011,
#     14300025,
#     14300027,
#     14300029,
#     14300031,
#     14300032,
#     14300034,
#     14300035,
#     14300045,
#     14300046,
#     15300009,
#     15300011,
#     15300025,
#     15300027,
#     15300029,
#     15300031,
#     15300032,
#     15300034,
#     15300035,
#     15300045,
#     15300046,
#     14300036,
#     14300038,
#     14300040,
#     15300036,
#     15300038,
#     15300040,
#     14300058,
#     15300058,
#     14315085,
#     14315086,
#     14315087,
#     14315088,
#     14315089,
#     14315090,
#     14315091,
#     14315092,
#     14315093,
#     14315094,
#     14315095,
#     14315096,
#     14319074,
#     14319075,
#     14319076,
#     14319077,
#     14319078,
#     14319079,
#     14319080,
#     14319081,
#     14319092,
#     14319093,
#     14319094,
#     14319095,
#     15315085,
#     15315086,
#     15315087,
#     15315088,
#     15315089,
#     15315090,
#     15315091,
#     15315092,
#     15315093,
#     15315094,
#     15315095,
#     15315096,
#     15319074,
#     15319075,
#     15319076,
#     15319077,
#     15319078,
#     15319079,
#     15319080,
#     15319081,
#     15319092,
#     15319093,
#     15319094,
#     15319095
# ]        
# set_th2p6=[
#     14319066,
#     14319067,
#     14319068,
#     14319069,
#     14319070,
#     14319071,
#     14319072,
#     14319073,
#     15319066,
#     15319067,
#     15319068,
#     15319069,
#     15319070,
#     15319071,
#     15319072,
#     15319073    
#     ]
# set_th3p0 = [
#     14303049,
#     14303085,
#     14303086,
#     14303107,
#     14307050,
#     14307051,
#     14307052,
#     14307053,
#     14307054,
#     14307080,
#     14307081,
#     14307082,
#     14307084,
#     14315075,
#     14315076,
#     14315077,
#     14315099,
#     14315100,
#     14315101,
#     14319082,
#     14319083,
#     15303049,
#     15303085,
#     15303086,
#     15303107,
#     15307050,
#     15307051,
#     15307052,
#     15307053,
#     15307054,
#     15307080,
#     15307081,
#     15307082,
#     15307084,
#     15315075,
#     15315076,
#     15315077,
#     15315099,
#     15315100,
#     15315101,
#     15319082,
#     15319083,
#     14307083,
#     14307103,
#     15307083,
#     15307103
#     ]
# set_th4p0 =[
#     14319096,
#     14319097,
#     15319096,
#     15319097
#     ]

# n_change_N1015025 =[
#     1415079,
#     1415085,
#     1415086,
#     1415094,
#     1415098,
#     1415104,
#     1515074,
#     1515075,
#     1515088,
#     1515089,
#     1515098,
#     1515104
#     ]
# n_change_N1019005 = [
#     1419082,
#     1419086,
#     1419090,
#     1419094,
#     1419095,
#     1419098,
#     1519079,
#     1519080,
#     1519087,
#     1519088,
#     1519092,
#     1519096
# ]


# # Set to modify Kink Floor from CDR to WSC_Iter01

# KF_It01_set_th1p5 = [
#     14315107,
#     14315108,
#     14319101,
#     14319102,
#     15315103,
#     15315104,
#     15319098,
#     15319099,
#     15319100,
#     14315002,
#     14315006,
#     14315010,
#     14315015,
#     14315021,
#     14315024,
#     14315025,
#     14315026,
#     14315027,
#     14315028,
#     14319000,
#     14319001,
#     14319003,
#     14319004,
#     14319008,
#     14319009,
#     14319011,
#     14319012,
#     14319013,
#     14319014,
#     14319015,
#     14319016,
#     14319017,
#     14319018,
#     14319019,
#     14319020,
#     15315006,
#     15315007,
#     15315017,
#     15315018,
#     15315019,
#     15315024,
#     15315025,
#     15315026,
#     15315027,
#     15315028,
#     15319006,
#     15319007,
#     15319008,
#     15319009,
#     15319010,
#     15319012,
#     15319013,
#     15319014,
#     15319015,
#     15319017,
#     15319018,
#     15319019,
#     15319020,
#     15319104,
#     15319105,
#     14315030,
#     15315030
#                     ]
# KF_It01_set_th1p7=[
#     14319098,
#     14319099,
#     14319100,
#     15319101,
#     15319102,
#     14319002,
#     14319005,
#     14319006,
#     14319007,
#     14319010,
#     14319021,
#     14319022,
#     14319023,
#     14319024,
#     15319001,
#     15319002,
#     15319003,
#     15319004,
#     15319005,
#     15319021,
#     15319022,
#     15319023,
#     15319024,
#     15319025
#     ]
# KF_It01_set_th1p8 =[
#     14315103,
#     14315104,
#     14315105,
#     14315106,
#     15315105,
#     15315106,
#     15315107,
#     15315108,
#     14315000,
#     14315001,
#     14315003,
#     14315004,
#     14315005,
#     14315007,
#     14315008,
#     14315009,
#     14315011,
#     14315012,
#     14315013,
#     14315014,
#     14315016,
#     14315017,
#     14315018,
#     14315019,
#     14315020,
#     14315022,
#     14315023,
#     14315029,
#     15315000,
#     15315001,
#     15315002,
#     15315003,
#     15315004,
#     15315005,
#     15315008,
#     15315009,
#     15315010,
#     15315011,
#     15315012,
#     15315013,
#     15315014,
#     15315015,
#     15315016,
#     15315020,
#     15315021,
#     15315022,
#     15315023,
#     15315029
#     ]

# # Modifications for It01
# # mod_2Delems_PID(KF_It01_set_th1p5,111327)
# # mod_2Delems_PID(KF_It01_set_th1p7,111328)
# # mod_2Delems_PID(KF_It01_set_th1p8,111300)


# P51 = [14315010,14315002,14315006,14315107,14315108,14315015,14315021,14315024,14315025,14315026,14315027,14315028]
# P52 = [14315000,14315001,14315003,14315005,14315008,14315011,14315012,14315016,14315020,14315023,14315105,14315106]
# P53 = [14315004,14315009,14315007,14315013,14315014,14315017,14315018,14315019,14315022,14315029,14315103,14315104]
# P54 = [15315107,15315108,15315004,15315005,15315008,15315009,15315010,15315011,15315015,15315020,15315022,15315029]
# P55 = [15315000,15315001,15315002,15315003,15315012,15315013,15315014,15315016,15315021,15315023,15315105,15315106]
# P56 = [15315103,15315104,15315006,15315007,15315017,15315018,15315019,15315024,15315025,15315026,15315027,15315028]
# P61 = [14315030,14319102,14319001,14319009,14319011,14319012,14319013]
# P62 = [14319000,14319003,14319004,14319008,14319014,14319015,14319016,14319017,14319018,14319019,14319020,14319101]
# P63 = [14319002,14319005,14319006,14319007,14319010,14319021,14319022,14319023,14319024,14319098,14319099,14319100]
# P64 = [15319101,15319102,15319001,15319002,15319003,15319004,15319005,15319021,15319022,15319023,15319024,15319025]
# P65 = [15319010,15319012,15319013,15319014,15319017,15319018,15319019,15319020,15319099,15319100,15319104,15319105]
# P66 = [15315030,15319098,15319006,15319007,15319008,15319009,15319015]

# inp_bdf = r'C:\Users\U69432\Desktop\WORK\01_FWD_FUS\09_SUELOS_WSC\03_Aral_KinkFloor\01_Aralx_Pandeo\Last\01_GFEM\01_FEM\DLL_240731_FL_Struc_FWD_08_copia.incl'
# out1 = mod_2Delems_PID(inp_bdf,P51+P56+P61+P62+P65+P66,111326)
# out2 = mod_2Delems_PID(out1,P52+P53+P54+P55,111327)
# out3 = mod_2Delems_PID(out2,P63+P64,111328)

bdf_file = r'Y:\03.SIZINGS\02-FUSELAGE\FRONT-FUSELAGE\05_SHEAR WALLS\WSC\LOWER_SHEAR_WALL\02_Nxx_Nxy\00_GFEM\01_FEM\DLL_240731_FL_Struc_FWD_08_tests.incl'
groups = {
    'CQUAD4' : [10534025, 10534026, 10534027, 10534028, 10534029, 10534030, 10534031, 10534032, 10534033, 10534034, 10534035, 10534036, 10534037, 10534038, 10534039, 10534040, 10534041, 10534042, 10534043, 10534044, 10534045, 10534046, 10534047, 10534048, 10534049, 10534050, 10534051, 10534052, 10534053, 10534054],
    'CBAR' : [10734306, 10734307, 10734308, 10734309, 10734310, 10734311, 10734333, 10734334, 10734334, 10734335, 10734336, 10734337, 10734338, 10734339, 10734340, 10734341, 10734342, 10734318, 10734319, 10734320, 10734321, 10734322, 10734323],
    'GRID' :[1034031, 1034032, 1034033, 1034034, 1034035, 1034036, 1034037, 1034038, 1034039, 1034040, 1034041, 1034042, 1034043, 1034044, 1034045, 1034046, 1034047, 1034048, 1034049, 1034050]
    }
comment_elements_by_field_pos(
    bdf_file = bdf_file,
    element_type_and_ids = groups,
    field_number = 2,
    check_duplicates = True
    )