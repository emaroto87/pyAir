# -*- coding: utf-8 -*-
"""
Created on Tue May 20 14:49:29 2025

@author: U69432
"""


import numpy as np
import heapq
from pyNastran.bdf.bdf import BDF as bdf
from collections import defaultdict
from collections import namedtuple

from fem.pre.nastran_cards import nastran_format_str, \
    nastran_format_int, nastran_format_float, create_card, split_into_chunks, mpc_card_fields

from common.ui.filedialogs import askOpenBDF, askSaveAsBDF
from common.ui.tui import TUI, ask_entities_ids_list
# from nastran_cards import nastran_format_int
# from nastran_cards import nastran_format_float
# from nastran_cards import create_card
# from nastran_cards import split_into_chunks
# from nastran_cards import mpc_card_fields

import os
import sys
import math

__version = '1.0.1'
__author = 'E.Maroto'

#                           ***   NAMEDTUPLES   ***

Point = namedtuple('Point', ['x', 'y', 'z'])
Edge = namedtuple('Edge', ['node1', 'node2', 'vector', 'length'])
sewing = namedtuple(
    'Sewing', ['gfem_node1', 'gfem_node2', 'dfem_nodes_between'])

#                        ***   AUXILIAR FUNCTIONS   ***


def concat(*args, **kargs):
    s = ''
    for i in args:
        s += i
    return s


def vector_mod(p1, p2):
    v = p2-p1
    return np.sqrt(np.dot(v, v))

#                            ***   FUNCTIONS   ***


def get_edges_from_element(pyNastran_elem) -> list:
    elem = pyNastran_elem
    edges = []
    nodes = elem.node_ids
    n = len(nodes)
    for i in range(n):
        edge = tuple(sorted((nodes[i], nodes[(i+1) % n])))
        edges.append(edge)
    return edges


def get_edges(pyNastran_BDF, sort=False, verbose=False):

    model = pyNastran_BDF
    elems = model.elements
    edge_counter = defaultdict(int)

    # Step 1 : count how many times an edge appears
    for eid, elem in elems.items():
        edges = get_edges_from_element(elem)
        for edge in edges:
            edge_counter[edge] += 1
    if verbose:
        for edge, count in edge_counter.items():
            print('Edge defined by nodes : ', edge,
                  'counted ', count, ' times')

    # Step 2 : edges that appears only one time
    frontier_edges = [edge for edge,
                      count in edge_counter.items() if count == 1]

    # Step 3 : Compute the vector and lenght of each edge
    edges_list = []
    for frontier_edge in frontier_edges:
        node1 = frontier_edge[0]
        node2 = frontier_edge[1]
        xyz1 = model.nodes[node1].xyz
        xyz2 = model.nodes[node2].xyz
        v = xyz2 - xyz1
        length = np.sqrt(np.dot(v, v))
        edge = Edge(node1, node2, v, length)
        edges_list.append(edge)
    return edges_list


def get_edges_nodes(dfem_edges):
    out = set()
    for edge in dfem_edges:
        t = (edge[0], edge[1])
        out.update(t)
    return out


def testing_get_edges():
    # # Load the model
    path1 = r'C:/Users/U69432/Desktop/WORK/00_Resources/02_Python_codes/PyNastran/examples_and_tests/BDF_tools/edges'
    path2 = 'C:/Users/U69432/Desktop/WORK/01_FWD_FUS/02_ECS_ShearWall_DFEM/02_FEM/'
    fname1 = 'test.bdf'
    fname2 = 'DFEM_ECS.dat'
    os.chdir(path1)
    model = bdf()
    model.read_bdf(
        bdf_filename=os.path.join(path2, fname2),
        read_includes=True,
        validate=False,
        xref=True,
        encoding='ANSI'
    )
    edges = get_edges(model.elements)
    return edges


def gfem_dfem_mapping(pyNastran_BDF, gfem_nids, dfem_nids, verbose=False):
    '''
    Given a pyNastran BDF DFEM model whcih includes some GFEM nodes, correlates
    per each GFEM node, the closest DFEM node.

    Parameters
    ----------
    pyNastran_BDF : pyNastran.BDF object
        pyNastran BDF model which should contain the GFEM and DFEM nodes.
    gfem_nids : integer
        IDs of the gfem nodes
    dfem_nids : integer
        IDs of the dfem nodes it is desired to map.

    Returns
    -------
    out : dictionary
        Returns the following mapping d={GFEM_id: closest_DFEM_id}

    '''
    model = pyNastran_BDF
    out = {}
    for gfem_nid in gfem_nids:
        d = float('inf')
        closest_node = ''
        gi = model.nodes[gfem_nid].xyz
        for dfem_nid in dfem_nids:
            di = model.nodes[dfem_nid].xyz
            v = di-gi
            r = np.sqrt(np.dot(v, v))
            if r < d:
                closest_node = dfem_nid
                d = r
        if verbose:
            print('GFEM node :', gfem_nid)
            print('Closest DFEM node :', closest_node)
            print('Distance :', d)
        out[gfem_nid] = closest_node

    return out


def dfem_frontier_path(pyNastran_BDF, DFEM_edges, id_start, id_objetive):
    model = pyNastran_BDF
    nodes = model.nodes
    graf = defaultdict(list)
    for node1, node2, vector, length in DFEM_edges:
        graf[node1].append((node2, length))
        graf[node2].append((node1, length))

    queue = [(0+vector_mod(nodes[id_start].xyz,
              nodes[id_objetive].xyz), 0, id_start, [id_start])]
    visited = set()
    while queue:
        f, g, current, path = heapq.heappop(queue)

        if current == id_objetive:
            return path, g
        if current in visited:
            continue
        visited.add(current)

        for neighbor, weight in graf[current]:
            if neighbor not in visited:
                g_new = g + weight
                h = vector_mod(nodes[neighbor].xyz, nodes[id_objetive].xyz)
                f_new = g_new + h
                heapq.heappush(
                    queue, (f_new, g_new, neighbor, path+[neighbor]))

    return None, float('inf')


def path_between_GFEM_nodes(pyNastran_BDF, dfem_edges,
                            gfem_dfem_mapping,
                            start_gfem_node_id,
                            end_gfem_node_id):

    start_dfem_node_id = gfem_dfem_mapping[start_gfem_node_id]
    end_dfem_node_id = gfem_dfem_mapping[end_gfem_node_id]
    dfp = dfem_frontier_path(pyNastran_BDF,
                             dfem_edges,
                             start_dfem_node_id,
                             end_dfem_node_id)

    while (dfp[0] is None):
        print('Unable to get path between GFEM nodes {0} and {1}'.format(
            start_gfem_node_id, end_gfem_node_id))
        print('Associated DFEM node to GFEM node {0} : {1}'.format(
            start_gfem_node_id, start_dfem_node_id))
        print('Associated DFEM node to GFEM node {0} : {1}'.format(
            end_gfem_node_id, end_dfem_node_id))
        print('Enter manually GFEM-DFEM association:')
        start_dfem_node_id = int(
            input('DFEM node to GFEM node {0} :'.format(start_gfem_node_id)))
        end_dfem_node_id = int(
            input('DFEM node to GFEM node {0} :'.format(end_gfem_node_id)))
        dfp = dfem_frontier_path(pyNastran_BDF,
                                 dfem_edges,
                                 start_dfem_node_id,
                                 end_dfem_node_id)

    sew = sewing(start_gfem_node_id, end_gfem_node_id, dfp[0])

    return sew


def project_point_to_line(P: Point, G1: Point, G2: Point) -> Point:
    '''
    Computes the projection of the Point P over the line defined by the points
    G1 and G2

    Parameters
    ----------
    P : Point (x,y,z)
        Coordinates of point to be projected to the line
    G1 : Point (x,y,z)
        Coordinate of the first point which defined the line r
    G2 : Point (x,y,z)
        Coordinate of the second point which defined the line r

    Returns
    -------
    Q : Point (x,y,z)
        Coordinates of the projection of point P over the line
    '''
    p = np.array(P)
    g1 = np.array(G1)
    g2 = np.array(G2)
    v = (g2-g1)/np.linalg.norm((g2-g1))
    # print(v)
    t = np.dot(v, (p-g1))/np.dot(v, v)
    # print(t)
    q = [0.0, 0.0, 0.0]
    for i in range(3):
        q[i] = g1[i]+t*v[i]
    # print(q)
    Q = Point(q[0], q[1], q[2])

    return Q


def sewing_weight_factor(G1, G2, D):
    g1 = G1
    g2 = G2
    qi = project_point_to_line(D, G1, G2)
    d1 = np.sqrt(np.dot((qi-g1), (qi-g1)))
    d2 = np.sqrt(np.dot((qi-g2), (qi-g2)))
    w1 = d2/(d1+d2)
    w2 = 1-w1
    return w1, w2


def sewing_RBE3_card(EID, G1, C1, G2, C2, D, CD, w1, w2):

    format_int = '>8s'
    format_weight = '8.6f'
    empty_field = 8*' '

    field1_1 = 'RBE3    '
    field1_2 = format(str(EID), format_int)
    field1_3 = empty_field
    field1_4 = format(str(D), format_int)
    field1_5 = format(str(CD), format_int)
    field1_6 = format(w1, format_weight)
    field1_7 = format(str(C1), format_int)
    field1_8 = format(str(G1), format_int)
    field1_9 = format(w2, format_weight)

    field2_1 = empty_field
    field2_2 = format(str(C2), format_int)
    field2_3 = format(str(G2), format_int)

    rbe3_card_line1 = concat(field1_1,
                             field1_2,
                             field1_3,
                             field1_4,
                             field1_5,
                             field1_6,
                             field1_7,
                             field1_8,
                             field1_9,
                             '\n'
                             )
    rbe3_card_line2 = concat(field2_1,
                             field2_2,
                             field2_3,
                             '\n'
                             )
    # print(rbe3_card_line1,
    #       rbe3_card_line2
    #       )
    return [rbe3_card_line1, rbe3_card_line2]


def sewing_RBE3_card_general(EID, G1, C1, G2, C2, GD, CD1, CD2, w1, w2):

    # Formatting styles
    format_int = '>8s'
    format_weight = '8.6f'
    empty_field = 8*' '
    format_weight_cd2 = '>8.1f'

    field1_1 = 'RBE3    '
    field1_2 = format(str(EID), format_int)
    field1_3 = empty_field
    field1_4 = format(str(GD), format_int)
    field1_5 = format(str(CD1), format_int)
    field1_6 = format(w1, format_weight)
    field1_7 = format(str(C1), format_int)
    field1_8 = format(str(G1), format_int)
    field1_9 = format(w2, format_weight)

    field2_1 = empty_field
    field2_2 = format(str(C2), format_int)
    field2_3 = format(str(G2), format_int)

    # Esto es para cubrir el caso donde se crean
    # dos RBE3 por connexion
    if CD2 is None:
        field2_4 = empty_field
        field2_5 = empty_field
        field2_6 = empty_field
    else:
        field2_4 = format(1.0, format_weight_cd2)
        field2_5 = format(str(CD2), format_int)
        field2_6 = format(str(GD), format_int)

    rbe3_card_line1 = concat(field1_1,
                             field1_2,
                             field1_3,
                             field1_4,
                             field1_5,
                             field1_6,
                             field1_7,
                             field1_8,
                             field1_9,
                             '\n'
                             )
    rbe3_card_line2 = concat(field2_1,
                             field2_2,
                             field2_3,
                             field2_4,
                             field2_5,
                             field2_6,
                             '\n'
                             )

    return [rbe3_card_line1, rbe3_card_line2]


def create_mpc_connection(
        file_path,
        GFEM_nodes,
        Aux_GFEM_nodes,
        components,
        coef_GFEM,
        coef_Aux, sid):

    # check
    if len(GFEM_nodes) != len(Aux_GFEM_nodes):
        raise ValueError(
            ' The number of GFEM nodes does not match with the number of ' /
            ' the Auxiliar GFEM nodes'
        )
    else:
        lines = []
        for component in components:
            n = len(GFEM_nodes)
            lines += [f'$ ---> MPCs Connection for component {component}\n']
            for i in range(n):
                print(
                    f'/tGeneraring MPC between {GFEM_nodes[i]} and {Aux_GFEM_nodes[i]}')
                GiAiCi = [
                    (GFEM_nodes[i], component, coef_GFEM),
                    (Aux_GFEM_nodes[i], component, coef_Aux),
                ]

                mpc_card = create_card(mpc_card_fields(sid, GiAiCi))
                lines += mpc_card

        with open(file_path, 'w') as f:
            f.writelines(lines)
            f.close()
        return


def sewing_card(
    pyNastran_BDF,
    sewing,
    eid_start,
    comment_lines=None,
    sewing_type='Coupled_123456'
    # gfem_node1_components = None,
    # gfem_node2_components = None,
    # dfem_node_components = None,
):

    model = pyNastran_BDF
    gfem_node1 = sewing[0]
    gfem_node2 = sewing[1]
    dfem_nodes = sewing[2]

    gfem_node1_xyz = model.nodes[gfem_node1].xyz
    gfem_node2_xyz = model.nodes[gfem_node2].xyz
    eid = eid_start

    lines = []

    if not (comment_lines is None):
        lines += comment_lines
    else:
        default_comment_line = '$ Sewing between GFEM nodes {0} {1}\n'.format(
            gfem_node1, gfem_node2)
        lines += [default_comment_line]

    for dfem_node in dfem_nodes:
        # print('/tGenerating RBE3 for {0}, {1} {2}'.format(gfem_node1,gfem_node2,dfem_node))
        dfem_node_xyz = model.nodes[dfem_node].xyz
        w1, w2 = sewing_weight_factor(
            gfem_node1_xyz,
            gfem_node2_xyz,
            dfem_node_xyz)

        if sewing_type.upper() == 'COUPLED_123456':
            C1 = 123456
            C2 = 123456
            CD1 = 123456
            CD2 = None
        elif sewing_type.upper() == 'COUPLED_123':
            C1 = 123
            C2 = 123
            CD1 = 123456
            CD2 = None
        elif sewing_type.upper() == 'HYPERMESH_123':
            C1 = 123
            C2 = 123
            CD1 = 123
            CD2 = 456
        elif sewing_type.upper() == 'HYPERMESH_456':
            C1 = 456
            C2 = 456
            CD1 = 456
            CD2 = 123
        else:
            print('Error en sewing type')

        cards = sewing_RBE3_card_general(
            EID=eid,
            G1=gfem_node1,
            C1=C1,
            G2=gfem_node2,
            C2=C2,
            GD=dfem_node,
            CD1=CD1,
            CD2=CD2,
            w1=w1,
            w2=w2
        )

        # if gfem_node1_components is None:
        #     cards= sewing_RBE3_card_general(
        #                         EID = eid,
        #                          G1 = gfem_node1,
        #                          G2 = gfem_node2,
        #                          D = dfem_node,
        #                          w1 = w1,
        #                          w2 = w2 )
        # else:
        #     cards = sewing_RBE3_card(EID = eid,
        #                          G1 = gfem_node1,
        #                          C1 = gfem_node1_components,
        #                          G2 = gfem_node2,
        #                          C2 = gfem_node2_components,
        #                          D = dfem_node,
        #                          CD = dfem_node_components,
        #                          w1 = w1,
        #                          w2 = w2 )
        lines += cards
        eid += 1
    return lines


def write_sewing_files(
        pyNastran_BDF,
        gfem_edges,
        dfem_edges,
        gfem_dfem_map,
        rbe_id_start,
        sewing_type='Coupled_123456',
        # gfem_node1_components = None,
        # gfem_node2_components = None,
        # dfem_node_components = None,
        file_suffix=''):

    model = pyNastran_BDF
    folder = os.path.dirname(model.bdf_filename)
    output_file = os.path.join(
        folder, 'GFEM-DFEM_RBE3_Sewing_w_OVC' + file_suffix + '.bdf')
    file = open(output_file, 'w')
    eid_start = rbe_id_start
    counter = 0

    for gfem_edge in gfem_edges:

        gfem_edge_node1 = gfem_edge[0]
        gfem_edge_node2 = gfem_edge[1]

        # Generating the Sewing between GFEM nodes
        print('Sewing {2} between GFEM nodes {0} {1}'.format(
            gfem_edge_node1, gfem_edge_node2, counter))
        # try:
        sewing = path_between_GFEM_nodes(model,
                                         dfem_edges,
                                         gfem_dfem_map,
                                         gfem_edge_node1,
                                         gfem_edge_node2)
        print('/tComputing Sewing: [Ok]')
        # Generating the Sewing cards between GFEM nodes
        sewing_cards = sewing_card(model,
                                   sewing=sewing,
                                   eid_start=eid_start,
                                   sewing_type=sewing_type,
                                   # gfem_node1_components = gfem_node1_components,
                                   # gfem_node2_components = gfem_node2_components,
                                   # dfem_node_components = dfem_node_components,
                                   comment_lines=None)
        print('/tComputing Sewing Cards : [Ok]')
        # Writting lines into the file
        print(sewing_cards)
        file.writelines(sewing_cards)
        # except:
        #     print('[Error]')
        counter += 1
        eid_start += 100

    file.close()
    print(f'Generated Sewing file : {output_file}')
    return output_file


def del_overconstrains(file_path):
    cdir = os.path.dirname(file_path)
    os.chdir(cdir)

    f = open(file_path, 'r')
    i_lines = f.readlines()
    f.close()

    fo = file_path.replace('_w_', '_wo_')
    fo = open(fo, 'w')

    slave_nodes = []
    slave_counter = defaultdict(int)
    for line in i_lines:
        if line.startswith('RBE3'):
            slave_node = line.split()[2]
            slave_nodes.append(slave_node)
    for slave_node in slave_nodes:
        slave_counter[slave_node] += 1

    # print(slave_nodes)
    overc_nodes = [slave for slave,
                   count in slave_counter.items() if count > 1]
    # print(overc_nodes)

    n = len(i_lines)
    i = 0
    while i < n:
        line = i_lines[i]
        if line.startswith('RBE3'):
            slave_node = line.split()[2]
            if slave_node in overc_nodes:
                print('Deleting duplicated node {0}'.format(slave_node))
                i += 2
            else:
                fo.writelines([line])
                i += 1
        else:
            fo.writelines([line])
            i += 1

    fo.close()


def ejemplo1_mallado_tipo_hypermesh():
    '''
    Procesos:
        1. Crear el modelo DFEM
        2. Meter dentro de ese modelo los nodos del GFEM frontera con el DFEM
        3. Crear una copia de esos nodos del GFEM. Estos serviran para crear las MPCs  
        3. Meter la lista de nodos del GFEM ordenados siguiendo una linea
        4. Dar el ID inicial de los RBE que generara el SEWING
    '''

    # Inputs -----
    path2 = 'Y:/03.SIZINGS/02-FUSELAGE/FRONT-FUSELAGE/05_SHEAR WALLS/WSC/LOWER_SHEAR_WALL/02_Nxx_Nxy/02_LOCAL_DFEM/'
    fname2 = 'Y53127283-001_DFEM_for_sewing.bdf'
    gfem_nodes_node1_ref = [1080815, 1080814, 1080813, 1080812, 1080811, 1080810, 1034030, 1034029, 1034028, 1034027,
                            1034026, 1080909, 1080910, 1080911, 1080912, 1080913, 1080914, 1034054, 1034055, 1034051, 1034052, 1034053]
    mpcs_nodes_node1_ref = [11080815, 11080814, 11080813, 11080812, 11080811, 11080810, 11034030, 11034029, 11034028, 11034027,
                            11034026, 11080909, 11080910, 11080911, 11080912, 11080913, 11080914, 11034054, 11034055, 11034051, 11034052, 11034053]
    initial_RBE3_id = 99910001
    initial_MPC_id = 99910000
    # -------------

    # Openining the BDF
    model = bdf()
    try:
        # Para Lanzadores
        model.read_bdf(
            bdf_filename=os.path.join(path2, fname2),
            read_includes=True,
            validate=False,
            xref=True,
            encoding='ANSI',
        )
    except:
        # Para includes
        model.read_bdf(
            bdf_filename=os.path.join(path2, fname2),
            read_includes=True,
            validate=False,
            xref=True,
            encoding='ANSI',
            punch=True
        )

    gfem_nodes_node2_ref = gfem_nodes_node1_ref[1:] + [gfem_nodes_node1_ref[0]]
    mpcs_nodes_node2_ref = mpcs_nodes_node1_ref[1:] + [mpcs_nodes_node1_ref[0]]
    # Los bordes es la mezcla ordenada
    gfem_edges = list(zip(gfem_nodes_node1_ref, gfem_nodes_node2_ref))
    mpcs_edges = list(zip(mpcs_nodes_node1_ref, mpcs_nodes_node2_ref))

    gfem_edges_nodes = gfem_nodes_node1_ref
    mpcs_edges_nodes = mpcs_nodes_node1_ref
    # Generamos los bordes del DFEM
    dfem_edges = get_edges(model)
    # Obtenemos los nodos de los bordes del DFEM
    dfem_edges_nodes = get_edges_nodes(dfem_edges)
    # Identinficamos cuales son los nodos del DFEM más proximos al GFEM
    gfem_dfem_map = gfem_dfem_mapping(
        pyNastran_BDF=model,
        gfem_nids=gfem_edges_nodes,
        dfem_nids=dfem_edges_nodes,
        verbose=True
    )
    mpcs_dfem_map = gfem_dfem_mapping(
        pyNastran_BDF=model,
        gfem_nids=mpcs_edges_nodes,
        dfem_nids=dfem_edges_nodes,
        verbose=True
    )
    # Generamos el Sewing
    write_sewing_files(
        pyNastran_BDF=model,
        gfem_edges=gfem_edges,
        dfem_edges=dfem_edges,
        gfem_dfem_map=gfem_dfem_map,
        rbe_id_start=99910001,
        sewing_type='Hypermesh_123',
        file_suffix='_123'
    )

    write_sewing_files(
        pyNastran_BDF=model,
        gfem_edges=mpcs_edges,
        dfem_edges=dfem_edges,
        gfem_dfem_map=mpcs_dfem_map,
        rbe_id_start=99920001,
        sewing_type='Hypermesh_456',
        file_suffix='_456'
    )

    # Eliminamos los overconstraints que genera el script
    # (El script genera overconstrait en los nodos del DFEM más proximos a los
    # nodos del GFEM. Y para eliminarlo usamos esta función. El ar)
    del_overconstrains(
        r'Y:/03.SIZINGS/02-FUSELAGE/FRONT-FUSELAGE/05_SHEAR WALLS/WSC/LOWER_SHEAR_WALL/02_Nxx_Nxy/02_LOCAL_DFEM/GFEM-DFEM_RBE3_Sewing_w_OVC_123.bdf')
    del_overconstrains(
        r'Y:/03.SIZINGS/02-FUSELAGE/FRONT-FUSELAGE/05_SHEAR WALLS/WSC/LOWER_SHEAR_WALL/02_Nxx_Nxy/02_LOCAL_DFEM/GFEM-DFEM_RBE3_Sewing_w_OVC_456.bdf')

    # Creating the MPC file
    mpc_file = r'X:/NASTRAN/LSW_Fr08-Fr09/01_FEM/Y53127283-001/MPC_connection.bdf'
    create_mpc_connection(mpc_file, gfem_edges_nodes, mpcs_edges_nodes, [
                          4, 5, 6], 1.0, -1.0, 99910456)


def ejemplo1_mallado_tipo_acoplado123456():
    '''
    Procesos:
        1. Crear el modelo DFEM
        2. Meter dentro de ese modelo los nodos del GFEM frontera con el DFEM
        3. Crear una copia de esos nodos del GFEM. Estos serviran para crear las MPCs  
        3. Meter la lista de nodos del GFEM ordenados siguiendo una linea
        4. Dar el ID inicial de los RBE que generara el SEWING
    '''

    # Inputs -----
    path2 = 'Y:/03.SIZINGS/02-FUSELAGE/FRONT-FUSELAGE/05_SHEAR WALLS/WSC/LOWER_SHEAR_WALL/04_Nxx_Nxy_SS_GFEM_with_floors_updated/02_LOCAL_DFEM/DFEM_aislado_Fr08-Fr09/01_FEM/Y53127283-001/'
    fname2 = 'Y53127283-001_FineMesh_Config08.bdf'
    gfem_nodes_node1_ref = [1080815, 1080814, 1080813, 1080812, 1080811, 1080810, 1034030, 1034029, 1034028, 1034027,
                            1034026, 1080909, 1080910, 1080911, 1080912, 1080913, 1080914, 1034054, 1034055, 1034051, 1034052, 1034053]
    mpcs_nodes_node1_ref = [11080815, 11080814, 11080813, 11080812, 11080811, 11080810, 11034030, 11034029, 11034028, 11034027,
                            11034026, 11080909, 11080910, 11080911, 11080912, 11080913, 11080914, 11034054, 11034055, 11034051, 11034052, 11034053]
    initial_RBE3_id = 99910001
    initial_MPC_id = 99910000
    # -------------

    # Openining the BDF
    model = bdf()
    try:
        # Para Lanzadores
        model.read_bdf(
            bdf_filename=os.path.join(path2, fname2),
            read_includes=True,
            validate=False,
            xref=True,
            encoding='ANSI',
        )
    except:
        # Para includes
        model.read_bdf(
            bdf_filename=os.path.join(path2, fname2),
            read_includes=True,
            validate=False,
            xref=True,
            encoding='ANSI',
            punch=True
        )

    gfem_nodes_node2_ref = gfem_nodes_node1_ref[1:] + [gfem_nodes_node1_ref[0]]
    # Los bordes es la mezcla ordenada
    gfem_edges = list(zip(gfem_nodes_node1_ref, gfem_nodes_node2_ref))

    gfem_edges_nodes = gfem_nodes_node1_ref
    # Generamos los bordes del DFEM
    dfem_edges = get_edges(model)
    # Obtenemos los nodos de los bordes del DFEM
    dfem_edges_nodes = get_edges_nodes(dfem_edges)
    # Identinficamos cuales son los nodos del DFEM más proximos al GFEM
    gfem_dfem_map = gfem_dfem_mapping(
        pyNastran_BDF=model,
        gfem_nids=gfem_edges_nodes,
        dfem_nids=dfem_edges_nodes,
        verbose=True
    )
    # Generamos el Sewing
    w_ovc_file = write_sewing_files(
        pyNastran_BDF=model,
        gfem_edges=gfem_edges,
        dfem_edges=dfem_edges,
        gfem_dfem_map=gfem_dfem_map,
        rbe_id_start=99910001,
        sewing_type='Coupled_123456',
        file_suffix='_Coupled_123456'
    )

    # Eliminamos los overconstraints que genera el script
    # (El script genera overconstrait en los nodos del DFEM más proximos a los
    # nodos del GFEM. Y para eliminarlo usamos esta función. El ar)

    del_overconstrains(w_ovc_file)


if __name__ == '__main__':
    '''
    Procesos:
        1. Crear el modelo DFEM
        2. Meter dentro de ese modelo los nodos del GFEM frontera con el DFEM
        3. Crear una copia de esos nodos del GFEM. Estos serviran para crear las MPCs  
        3. Meter la lista de nodos del GFEM ordenados siguiendo una linea
        4. Dar el ID inicial de los RBE que generara el SEWING
    '''

    # Inputs -----

    # Openining the BDF
    file = askOpenBDF()
    print(file)
    model = bdf()
    try:
        # Para Lanzadores
        model.read_bdf(
            bdf_filename=file,
            read_includes=True,
            validate=False,
            xref=True,
            encoding='ANSI',
        )
    except:
        # Para leer ficheros que solo contenga includes
        model.read_bdf(
            bdf_filename=file,
            read_includes=True,
            validate=False,
            xref=True,
            encoding='ANSI',
            punch=True
        )

    # gfem_nodes_node1_ref = [1080815, 1080814, 1080813, 1080812, 1080811, 1080810, 1034030, 1034029, 1034028, 1034027,
    #                         1034026, 1080909, 1080910, 1080911, 1080912, 1080913, 1080914, 1034054, 1034055, 1034051, 1034052, 1034053]
    # mpcs_nodes_node1_ref = [11080815, 11080814, 11080813, 11080812, 11080811, 11080810, 11034030, 11034029, 11034028, 11034027,
    #                         11034026, 11080909, 11080910, 11080911, 11080912, 11080913, 11080914, 11034054, 11034055, 11034051, 11034052, 11034053]
    # initial_RBE3_id = 99910001
    # initial_MPC_id = 99910000

    __sewing_opts = ['COUPLED_123456', 'COUPLED_123', 'HYPERMESH']

    sewing_opts_menu = TUI(
        options=__sewing_opts,
        sel_msg='Select one of the Sewing methods:',
        label='sewing_menu',
    )
    sewing_opts_menu.loop()
    print(sewing_opts_menu.selection)

    print('Provide list of Node IDs of the Frontier:')
    gfem_nodes_node1_ref = ask_entities_ids_list()

    initial_RBE3_id = input(
        'Enter the initial ID of the Sewing elements (integer value): ')

    if not (sewing_opts_menu.selection == 'HYPERMESH'):

        gfem_nodes_node2_ref = gfem_nodes_node1_ref[1:] + \
            [gfem_nodes_node1_ref[0]]
        # Los bordes es la mezcla ordenada
        gfem_edges = list(zip(gfem_nodes_node1_ref, gfem_nodes_node2_ref))

        gfem_edges_nodes = gfem_nodes_node1_ref
        # Generamos los bordes del DFEM
        dfem_edges = get_edges(model)
        # Obtenemos los nodos de los bordes del DFEM
        dfem_edges_nodes = get_edges_nodes(dfem_edges)
        # Identinficamos cuales son los nodos del DFEM más proximos al GFEM
        gfem_dfem_map = gfem_dfem_mapping(
            pyNastran_BDF=model,
            gfem_nids=gfem_edges_nodes,
            dfem_nids=dfem_edges_nodes,
            verbose=True
        )
        # Generamos el Sewing
        w_ovc_file = write_sewing_files(
            pyNastran_BDF=model,
            gfem_edges=gfem_edges,
            dfem_edges=dfem_edges,
            gfem_dfem_map=gfem_dfem_map,
            rbe_id_start=99910001,
            sewing_type='Coupled_123456',
            file_suffix='_Coupled_123456'
        )

        # Eliminamos los overconstraints que genera el script
        # (El script genera overconstrait en los nodos del DFEM más proximos a los
        # nodos del GFEM. Y para eliminarlo usamos esta función. El ar)

        del_overconstrains(w_ovc_file)

    else:
        print('Provide list of Duplicated Node IDs of the Frontier used for MPCs:')
        gfem_nodes_node1_ref = ask_entities_ids_list()
        initial_MPC_id = input(
            'Enter the initial ID of the MPCs (integer value): ')

        gfem_nodes_node2_ref = gfem_nodes_node1_ref[1:] + \
            [gfem_nodes_node1_ref[0]]
        mpcs_nodes_node2_ref = mpcs_nodes_node1_ref[1:] + \
            [mpcs_nodes_node1_ref[0]]
        # Los bordes es la mezcla ordenada
        gfem_edges = list(zip(gfem_nodes_node1_ref, gfem_nodes_node2_ref))
        mpcs_edges = list(zip(mpcs_nodes_node1_ref, mpcs_nodes_node2_ref))

        gfem_edges_nodes = gfem_nodes_node1_ref
        mpcs_edges_nodes = mpcs_nodes_node1_ref
        # Generamos los bordes del DFEM
        dfem_edges = get_edges(model)
        # Obtenemos los nodos de los bordes del DFEM
        dfem_edges_nodes = get_edges_nodes(dfem_edges)
        # Identinficamos cuales son los nodos del DFEM más proximos al GFEM
        gfem_dfem_map = gfem_dfem_mapping(
            pyNastran_BDF=model,
            gfem_nids=gfem_edges_nodes,
            dfem_nids=dfem_edges_nodes,
            verbose=True
        )
        mpcs_dfem_map = gfem_dfem_mapping(
            pyNastran_BDF=model,
            gfem_nids=mpcs_edges_nodes,
            dfem_nids=dfem_edges_nodes,
            verbose=True
        )
        # Generamos el Sewing
        write_sewing_files(
            pyNastran_BDF=model,
            gfem_edges=gfem_edges,
            dfem_edges=dfem_edges,
            gfem_dfem_map=gfem_dfem_map,
            rbe_id_start=99910001,
            sewing_type='Hypermesh_123',
            file_suffix='_123'
        )

        write_sewing_files(
            pyNastran_BDF=model,
            gfem_edges=mpcs_edges,
            dfem_edges=dfem_edges,
            gfem_dfem_map=mpcs_dfem_map,
            rbe_id_start=99920001,
            sewing_type='Hypermesh_456',
            file_suffix='_456'
        )

        # Eliminamos los overconstraints que genera el script
        # (El script genera overconstrait en los nodos del DFEM más proximos a los
        # nodos del GFEM. Y para eliminarlo usamos esta función. El ar)
        del_overconstrains(
            r'Y:/03.SIZINGS/02-FUSELAGE/FRONT-FUSELAGE/05_SHEAR WALLS/WSC/LOWER_SHEAR_WALL/02_Nxx_Nxy/02_LOCAL_DFEM/GFEM-DFEM_RBE3_Sewing_w_OVC_123.bdf')
        del_overconstrains(
            r'Y:/03.SIZINGS/02-FUSELAGE/FRONT-FUSELAGE/05_SHEAR WALLS/WSC/LOWER_SHEAR_WALL/02_Nxx_Nxy/02_LOCAL_DFEM/GFEM-DFEM_RBE3_Sewing_w_OVC_456.bdf')

        # Creating the MPC file
        mpc_file = r'X:/NASTRAN/LSW_Fr08-Fr09/01_FEM/Y53127283-001/MPC_connection.bdf'
        create_mpc_connection(mpc_file, gfem_edges_nodes, mpcs_edges_nodes, [
                              4, 5, 6], 1.0, -1.0, 99910456)
