# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 14:48:19 2025

@author: U69432
"""

import sys
import os
from collections import Counter
from collections import namedtuple
from collections import defaultdict
from fem.pre.nastran_cards import create_card
from fem.pre.nastran_cards import nastran_format_float
from fem.pre.nastran_cards import nastran_format_int
from fem.pre.nastran_cards import nastran_format_str
from tkinter.filedialog import asksaveasfilename
from tkinter.filedialog import askopenfilename
import tkinter as tk
from pathlib import Path
from enum import Enum
import abc
from typing import Optional
from dataclasses import dataclass, asdict
from pyNastran.bdf.bdf import BDF
from geometry.point import Point
from geometry.vector import Vector
from geometry.axes import Csys
from geometry.plane import Plane_by3points as Plane
from geometry.utils import angle_between_vectors
import numpy as np
__version = '1.0.4'
__author = 'E.Maroto'

# ========================
#
# PARTE DE LECTURA DE BDF
#
# ========================
# Ejemplos para leer BDF y conseguir datos de él


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


"""
Nota el siguiente codigo asume que el elemento 'eid' es una lista que contiene
los nodos que lo definen.


"""


def get_2Delement_planes(bdf_model: BDF) -> dict[Plane]:
    # Genero la informaciÃ³n de todos los planos que tiene el modelo. En el caso
    # de elementos CQUAD4, cada uno los divido en 2 CTRIA3.
    model = bdf_model
    planes = {}
    nodes = model.nodes
    elements = model.elements
    for eid, e in elements.items():
        if e.type.upper() in ['CQUAD4', 'CTRIA3']:
            e_nodes = e.nodes
            max_edge_len = get_max_edge_length(model, eid)
            if e.type.upper() == 'CQUAD4':
                tol = 1e-4 * max_edge_len
                # primera triangulaciÃ³n
                n1 = e.nodes_ref[0]
                n2 = e.nodes_ref[1]
                n3 = e.nodes_ref[2]
                p1 = Point(*n1.xyz)
                p2 = Point(*n2.xyz)
                p3 = Point(*n3.xyz)
                s = Plane(p1, p2, p3, eid, tol)
                planes[f'{e.type}_{eid}_S0'] = s

                # segunda triangulaciÃ³n
                n1 = e.nodes_ref[0]
                n2 = e.nodes_ref[2]
                n3 = e.nodes_ref[3]
                p1 = Point(*n1.xyz)
                p2 = Point(*n2.xyz)
                p3 = Point(*n3.xyz)
                s = Plane(p1, p2, p3, eid, tol)
                planes[f'{e.type}_{eid}_S1'] = s

            else:
                tol = 1e-12 * max_edge_len + 1e-15
                n1 = e.nodes_ref[0]
                n2 = e.nodes_ref[1]
                n3 = e.nodes_ref[2]
                p1 = Point(*n1.xyz)
                p2 = Point(*n2.xyz)
                p3 = Point(*n3.xyz)
                s = Plane(p1, p2, p3, eid, tol)
                planes[f'{e.type}_{eid}_S0'] = s
        else:
            pass
    return planes


def get_max_node_id(bdf_model: BDF) -> int | None:
    nids = set(bdf_model.nodes.keys())
    if not nids:
        return None
    else:
        return max(nids)


def get_max_mpc_id(bdf_model: BDF) -> int | None:
    mpcids = set(bdf_model.mpcs.keys())
    if not mpcids:
        return None
    else:
        return max(mpcids)


def get_max_element_id(bdf_model: BDF) -> int | None:
    eids = set()
    eids.update(bdf_model.elements.keys())
    eids.update(bdf_model.rigid_elements.keys())
    eids.update(bdf_model.masses.keys())

    if not eids:
        return None
    else:
        return max(eids)


def get_max_property_id(bdf_model: BDF) -> int | None:
    pids = set(bdf_model.properties.keys())
    if not pids:
        return None
    else:
        return max(pids)


def get_adjacent_elements(bdf_model: BDF, eids: list, levels=1):
    node_to_elements_map = bdf_model.get_node_id_to_element_ids_map()
    adj = set()

    try:
        for level in range(levels):
            for eid in eids:
                elem = bdf_model.elements[eid]
                # print(elem.node_ids)
                for nid in elem.node_ids:
                    adj.update(node_to_elements_map[nid])
            eids = list(adj)

        return eids
    except:
        return None


def get_angle_matsys_fast(bdf_model: BDF, eid: int, fastener_vector: Vector) -> float:
    e = bdf_model.elements[eid]

    # Using the method "material_coordinate_system" there is no longer need to
    # check how the material coordinate is defined. This function computes
    # straightforward the rotation of the elemental coordinate system to the
    # to the material coordinate system, regardaless the method : giving theta
    # or theta_mcid.

    lenght, centroid, dir1, dir2, dir3 = e.material_coordinate_system()
    mat_sys = Csys(
        origin=centroid,
        dir1=Vector(Point(0, 0, 0), Point(*dir1)),
        dir2=Vector(Point(0, 0, 0), Point(*dir2)),
        dir3=Vector(Point(0, 0, 0), Point(*dir3))
    )
    rot_mat_sys = mat_sys.rotate_by_vector(fastener_vector)
    angle = angle_between_vectors(rot_mat_sys.dir1, mat_sys.dir1)
    return angle


def get_max_edge_length(bdf_model: BDF, eid: int) -> float:
    """
    Return the maximum length of the edges of a 1D, 2D or 3D element.

    Parameters
    ----------
    bdf_model : pyNASTRAN BDF object

    eid : int
        Element ID.
    """
    elem = bdf_model.elements[eid]
    edge_ids = elem.get_edge_ids()
    edge_lengths = []
    for edge in edge_ids:
        node1 = bdf_model.nodes[edge[0]]
        node2 = bdf_model.nodes[edge[1]]
        vector = node2.xyz - node1.xyz
        length = np.linalg.norm(vector)
        edge_lengths.append(length)

    return max(edge_lengths)


def get_mesh_edge_nodes(
        bdf_model: BDF,
        eids: list[int] | None = None,
        sort: bool = False,
        verbose: bool = False
) -> set:
    """
    Return the frontier nodes of mesh defined by the element ids.

    Parameters
    ----------
    bdf_model : pyNastran BDF model

    eids : list[int] | None, optional
        Element IDs of the mesh. The default is None.
    sort : bool, optional
        DESCRIPTION. The default is False.
    verbose : bool, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    frontier_nodes : set
        set of the node IDs that make up the frontier of the mesh.
    """

    if eids:
        elems = {eid: bdf_model.elements[eid] for eid in eids}
    else:
        elems = bdf_model.elements

    edge_counter = defaultdict(int)

    # Step 1 : count how many times an edge appears
    for eid, elem in elems.items():
        edges = elem.get_edge_ids()
        for edge in edges:
            edge_counter[edge] += 1
    if verbose:
        for edge, count in edge_counter.items():
            print('Edge defined by nodes : ', edge,
                  'counted ', count, ' times')

    # Step 2 : edges that appears only one time
    frontier_edges = [edge for edge,
                      count in edge_counter.items() if count == 1]

    # Step 3 : get the unique nodes of the frontier edges
    frontier_nodes = set()
    for edge in frontier_edges:
        frontier_nodes.update(edge)
    if sort is True:
        frontier_nodes = sorted(frontier_nodes)

    return frontier_nodes


def __testing_get_edges():
    # # Load the model
    path = r'C:\Users\U69432\Desktop\WORK\00_Resources\02_Python_codes\PyNastran\examples_and_tests\BDF_tools\edges'
    fname = 'test.bdf'
    os.chdir(path)
    model = BDF()
    model.read_bdf(
        bdf_filename=os.path.join(path, fname),
        read_includes=True,
        validate=False,
        xref=True,
    )
    edges = get_mesh_edge_nodes(model)
    return edges


def read_bdf(str):
    with open(str, 'r') as f:
        lines = f.readlines()
    return lines

# @dataclass
# class PSHELL:
#     PID : int
#     T : float


def parse_bdf_line(bdf_line: str, field_chars: int = 8):
    # Skip if comment line (Starts with $)
    if bdf_line.startswith('$'):
        return bdf_line

    # Parsing line

    else:
        fields = []
        # Couting number of fields
        n = len(bdf_line)
        mod, res = divmod(n, field_chars)

        for i in range(mod):
            start = field_chars*i
            end = field_chars*(i+1)
            fields.append(bdf_line[start:end])

        if res != 0:
            start = field_chars*(mod)
            fields.append(bdf_line[start:])
        return fields


def read_bdf(str):
    with open(str, 'r') as f:
        lines = f.readlines()

    # Estrategia linea por linea

    return lines


def change_2D_prop(lines, elems: list, PID: int, field_nchars=8):
    # with open(bdf_file) as f:
    #     lines = f.readlines()
    inp_lines = lines
    out_lines = []
    for line in inp_lines:
        if line.startswith(('CQUAD4', 'CTRIA3')):
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


def extract_field(line, field_number, fields_chars=8):
    start = (field_number-1) * fields_chars
    end = (field_number) * fields_chars
    field = line[start:end]
    return field


def openBDF():
    bdf_fp = askopenfilename(
        title='Open NASTRAN BDF file',
        filetypes=(
            ("Nastran BDF file", '*.bdf'),
            ("Data file", '*.dat'),
            ("Include file", '*.incl'),
        ),
        multiple=False,
        defaultextension=['*.bdf'],
        initialdir=os.getcwd()
    )
    if bdf_fp == '':
        print('No file has been selected.')
    return bdf_fp


def mod_2Delems_PID(inp_bdf, elems_list: list, PID: int, field_nchars=8):
    out_bdf = inp_bdf.split('.')[0] + '_mod.' + inp_bdf.split('.')[1]

    with open(inp_bdf, 'r') as f:
        lines = f.readlines()

    with open(out_bdf, 'w') as f:
        out_lines = change_2D_prop(lines, elems_list, PID)
        f.writelines(out_lines)

    return out_bdf


def align_nodes(lines: list, nodes_list: list, node_ref: int, comp: int):
    out_lines = []
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


def mod_node_coord(inp_bdf, nodes_list: list, node_ref: int, comp: int):
    out_bdf = inp_bdf.split('.')[0] + '_mod.' + inp_bdf.split('.')[1]

    with open(inp_bdf, 'r') as f:
        lines = f.readlines()

    with open(out_bdf, 'w') as f:
        out_lines = align_nodes(lines, nodes_list, node_ref, comp)
        f.writelines(out_lines)

    return out_bdf


def default_output_file(
        input_file: str,
        suffix: str | None = None,
        preffix: str | None = None
) -> str:
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

    folder = os.path.dirname(input_file)
    file = os.path.basename(input_file)
    filename, filextension = os.path.splitext(file)
    if (preffix is None) and (suffix is None):
        suffix = '_mod'
    if suffix is None:
        suffix = ''
    if preffix is None:
        preffix = ''
    else:
        print('Error')

    output_filename = preffix + filename + suffix
    output_filepath = os.path.join(folder, output_filename + filextension)

    return output_filepath


def check_duplicated_ids(lists: list[list]):
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
    combined = sum(lists, [])

    # La función counter te devuelve un diccionario donde cada key es cada
    # uno de los miembros de las lista, y el valor corresponde con el numero
    # de veces presente en la lista
    counter = Counter(combined)

    duplicates = [item_id for item_id, freq in counter.items() if freq > 1]

    # Lista sin duplicados
    unique = set(combined)

    return unique, duplicates


def comment_cards_by_field_pos(
        bdf_file: str,
        element_type_and_ids: dict,
        field_number: int = 1,
        check_duplicates: bool = False):

    # Nastran fields number is between [1,10]
    if field_number < 1 or field_number > 10:
        raise ValueError(
            'The value of field number must be greater than 1 and do not'
            'exceed 10.'
        )

    # Duplicates in the lists
    element_types = element_type_and_ids.keys()
    unique_ids, duplicate_ids = check_duplicated_ids(
        element_type_and_ids.values())

    if check_duplicates:
        if len(duplicate_ids) > 0:
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
    while i < nlines:
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
                i += 1
                while True:
                    line = lines[i]
                    if line.startswith(('+')):
                        line = '$ Commented -> ' + line
                        out_lines += [line]
                        i += 1
                    else:
                        break
            else:
                out_lines += [line]
                i += 1

        else:
            out_lines += [line]
            i += 1

    # Writing the file with the commented lines
    out_fp = default_output_file(bdf_file)
    o = open(out_fp, 'w')
    o.writelines(out_lines)
    print('<--- END : Commenting elements')
    print(f'Created file at : {out_fp}')
    o.close()

###############################################################################
#                                                                             #
#                                 EXAMPLES                                    #
#                                                                             #
###############################################################################


def ejemplo1_commentar_fuerzas():
    bdf_file = r'X:\NASTRAN\Jacking\03_LOADS\HANDLING\DLL_JACK+HOIST_FORCE_MOMENT.bdf'
    groups = {
        'FORCE': [34006, 34007, 34008, 35006, 35007, 35008],
        'MOMENT':  [34006, 34007, 34008, 35006, 35007, 35008]
    }
    comment_cards_by_field_pos(
        bdf_file=bdf_file,
        element_type_and_ids=groups,
        field_number=3,
        check_duplicates=True
    )


def ejemplo_comentar_SPCDs():
    path = r'Y:\03.SIZINGS\02-FUSELAGE\FRONT-FUSELAGE\05_SHEAR WALLS\WSC\LOWER_SHEAR_WALL\03_Nxx_Nxy_SS\02_LOCAL_DFEM\DFEM_aislado_Fr08-Fr09\03_LOADS\FUS\FWD-FUS-DLL_SPCDs_LSW_Fr08-Fr09.bdf'
    groups = {'SPCD': [2, 4, 5, 6]}

    comment_cards_by_field_pos(
        bdf_file=path,
        element_type_and_ids=groups,
        field_number=4,
        check_duplicates=True
    )
