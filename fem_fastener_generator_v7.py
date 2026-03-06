# -*- coding: utf-8 -*-
"""
Created on Tue Feb 17 15:14:59 2026

@author: U69432
"""

'''
Script para generar remaches de manera rápida en NASTRAN

Inputs:
    - Datos de remaches (localización, tipo, vector normal, radio, direción)
    - Malla 2D con propiedades aplicadas
    -

Pasos:
    1. Coger el punto de refencia (P) y su vector de dirección (Vd) y obtener
    los puntos proyectados en las superficies.
        1.a) Hecho el calculo de la proyección de un punto sobre cualquier
        plano
        1.b) Averiguar si el nodo proyectado está dentro o no de un elemento

    2. Generar CBUSH
    3. Calcular Huth
        3.1 Obtener el apilado del elemento
        3.2 Calcular E1,E2 del apilado en la dirección que se desee
        3.3 Calular Kshear en esa dirección.
    4. Proceso recursivo

'''
# POR COMPLETAR/AÑADIR
# -----------------------------------------------------------------------------
# COSAS POR AÑADIR/COMPLETAR:
#   1. Que lea la parte de compuesto.
#   2. Que las direcciones de CBUSH pueda meterlo en formato de vector o
#      de sistema de coordenadas. Ahora mismo sólo está con un vector.
#   3. Que el usuario pueda poner los ids que quiera para los nodos, elementos
#      y propiedades creados.
#
# POR CORREGIR
#   1. El problema que estoy teniendo en este punto es si los puntos proyectados
#   cae justo en un nodo ya existente, me da todos los elementos que connectan
#   a ese nodo.
#   2. En uniones con más de dos placas, duplico el coupling en las placas
#   intermedias --> RESUELTO
#   3. Los coupling me cogen nodos de otras placas --> RESUELTO
#   4. A la hora de generar las conexiones se hace en el orden en que han
#   sido almacenados. Es decir, no guarda el orden logico de la conexión pu
#   diendo darse el caso de que la 2 intersección conecte con la 4. Esto
#   se puede solucionar fácilmente si se calcula la distancia de los puntos
#   de proyección al pto. de referencia y se orden de menor a mayor distancia
#   Además, se puede aprovechar este proceso para eliminar aquellas intersecciones
#   que no tengan sentido (l > (t1 + t2)) ------> RESUELTO
#   5. Habría que modificar la formula de Huth para que considere el caso
#   de doble cortadura y que lo meta en el PBUSH


from enum import Enum
from Materials.composite import Laminate
from Materials.metallic import Metallic
from general import is_type
from pyNastran.bdf.bdf import BDF
from Geom.Euclidian_geom import Point
from Geom.Euclidian_geom import Vector
from Geom.Euclidian_geom import Plane_by3points as Plane
from math import ceil
from math import pi
__HEAD_TYPE_PTR = 'PTR'
__HEAD_TYPE_CSK = 'CSK'

__version__ = 4

# =============================================================================
#                                                                             #
#                              *** CLASSES ***                                #
#                                                                             #
# =============================================================================


class HeadType(Enum):
    PTR = 'protruding'
    CSK = 'countersunk'


class JointType(Enum):
    SINGLE = 'single_shear'
    DOUBLE = 'double_shear'


class Plate:
    pass


class MetallicPlate:
    def __init__(self,
                 material: Metallic,
                 thickness: float,
                 plate_id: int = None,
                 linked_FEM_eid: int = None,
                 linked_FEM_pid: int = None):

        self.material = material
        self.thickness = thickness
        self.linked_FEM_eid = linked_FEM_eid
        self.linked_FEM_pid = linked_FEM_pid
        self.plate_id = plate_id

    def __repr__(self):
        msg = f'{self.material.name} Metallic Plate ID {self.plate_id}' \
              f'with thick = {self.thickness}. Linked ElemID :' \
              f'{self.linked_FEM_eid}'
        return msg


class CompositePlate:
    def __init__(self,
                 laminate: Laminate,
                 plate_id: int = None,
                 linked_FEM_eid: int = None,
                 linked_FEM_pid: int = None):

        self.laminate = laminate
        self.thickness = laminate.thickness
        self.linked_FEM_eid = linked_FEM_eid
        self.linked_FEM_pid = linked_FEM_pid
        self.plate_id = plate_id

    def get_equivalent_moduli(self, angle_deg: float):
        '''
        Returns the equivalent values of Ex, Ey, Gxy and nuxy into another
        coordinate system which is defined as a rotation from the stacking
        coordinate system
        '''
        data = self.laminate.laminate_apparent_moduli(angle_deg)
        Ex_local = data['Ex']
        Ey_local = data['Ey']
        Gxy_local = data['Gxy']
        nuxy_local = data['nuxy']

        return Ex_local, Ey_local, Gxy_local, nuxy_local


class Fastener:
    def __init__(self,
                 Et: float,
                 diameter: float,
                 rivet=True,
                 head_type: HeadType | None = None,
                 ref_point: Point | None = None,
                 shank_vector: Vector | None = None,
                 orientation_vector: Vector | None = None,
                 fastener_id: int | None = None,
                 fastener_label: str | None = None,
                 ):
        self.Et = Et
        self.diameter = diameter
        self.head_type = head_type
        self.ref_point = ref_point
        self.shank_vector = shank_vector
        self.orientation_vector = orientation_vector
        self.fastener_id = fastener_id
        self.fastener_label = fastener_label

        if rivet:
            self.is_bolt = False
            self.is_rivet = True
        else:
            self.is_bolt = True
            self.is_rivet = False

    def __str__(self):
        msg = f'Fastener['
        if not self.fastener_id:
            msg += f' ID:{self.fastener_id}'
        if not self.fastener_label:
            msg += f' Label:{self.fastener_label}'
        msg += f' D: {self.diameter} ]'
        return msg


class Connector:
    def __init__(
            self,
            linked_fastener: Fastener,
            bdf_model: BDF,
            planes: dict[Plane],
            connector_id: int | None = None,
            skip_tol: float = 1.25
    ):

        self.linked_fastener = linked_fastener
        self.connector_id = connector_id
        self.intersections = self.__get_intersections(bdf_model, planes)
        self.__skip_tol = skip_tol

    def __str__(self):
        msg = f'Connector to FEM element ids {self.linked_eids}'
        return msg

    def __get_intersections(self, bdf_model: BDF, planes: dict[Plane]):
        # Retrieving the planes of all the 2D (CQUAD4, CTRIA3) elements
        # in order to generate the projected nodes for each connector
        elements = bdf_model.elements
        # Creating the intersections
        intersections = []
        inter_elems = []
        last_plane_id = -9
        ref_point = self.linked_fastener.ref_point
        shank_dir = self.linked_fastener.shank_vector
        # Getting the intersections
        for label, plane in planes.items():

            projected_pt = plane.proyect_point(ref_point, shank_dir)
            if not (projected_pt is None):
                if plane.intriangle(projected_pt):
                    if plane.plane_id != last_plane_id:  # Avoids duplicates
                        inter_elem = elements[plane.plane_id]
                        prop_inter_elem = inter_elem.pid_ref
                        if prop_inter_elem.type.upper() == 'PSHELL':  # --> Metallic
                            t = prop_inter_elem.t
                            mat = Metallic(Et=prop_inter_elem.mid1_ref.e)
                            plate = MetallicPlate(
                                material=mat,
                                thickness=t,
                                linked_FEM_eid=plane.plane_id,
                                linked_FEM_pid=inter_elem.pid)
                        elif prop_inter_elem.type.upper() == 'PCOMP':  # --> Composite
                            # --- TO BE COMPLETED
                            plate = None
                            # --- TO BE COMPLETED
                        intersections.append((projected_pt, plate))
                        last_plane_id = plane.plane_id
                        inter_elems.append(plane.plane_id)

        # NOTA: SEGUIMOS TENIENDO VARIOS ELEMENTOS POR PLANO PORQUE EN ESTE CASO
        # TENEMOS QUE LAS PROYECCIONES COINCIDEN CON ALGUN NODO DEL MODELO.

        # Sorting
        _ = []
        for point, plane in intersections:
            # print(point)
            v = Vector(ref_point, point)
            dist = v.norm
            _.append((dist, (point, plane)))

        ordered_intersections = sorted(_, key=lambda x: [x[0]])
        intersections = [intersection[1]
                         for intersection in ordered_intersections]
        # print(intersections)

        return intersections

    def __filter_intersections(self, bdf_model: BDF, planes: dict[Plane], verbose=False):
        intersections = self.__get_intersections(bdf_model, planes)
        cond = True
        while cond:
            n = len(intersections)
            flags = []
            for i in range(n-1):
                # Retrieving data from connections
                pt1, plate1 = intersections[i]
                pt2, plate2 = intersections[i+1]

                # Check if the connection is out of the sum of plate thicknesses
                v = Vector(pt1, pt2)
                conn_length = v.norm
                t1 = plate1.thickness
                t2 = plate2.thickness
                t_total = t1 + t2

                if conn_length > self.__skip_tol*t_total:
                    intersections.pop(i+1)
                    flags.append(True)
                else:
                    flags.append(False)

            cond = any(flags)


# =============================================================================
#                                                                             #
#                             *** FUNCTIONS ***                               #
#                                                                             #
# =============================================================================


def kshear_huth(plate1: Plate, plate2: Plate, fastener: Fastener) -> float:

    # data from plate 1
    t1 = plate1.thickness
    E1 = plate1.material.Et
    plate1_type = 'METALLIC' if type(plate1) == MetallicPlate else 'COMPOSITE'

    # data from plate 2
    t2 = plate2.thickness
    E2 = plate2.material.Et
    plate2_type = 'METALLIC' if type(plate2) == MetallicPlate else 'COMPOSITE'

    # data from fastener
    Er = fastener.Et
    d = fastener.diameter

    # Selecting the parameters "a" and "b"
    if plate1_type == 'COMPOSITE' and plate2_type == 'COMPOSITE':
        a = 2/3
        b1 = 4.2
        b2 = 4.2
    elif plate1_type == 'COMPOSITE' and plate2_type == 'METALLIC':
        a = 2/3
        b1 = 4.2
        b2 = 3.0
    elif plate1_type == 'METALLIC' and plate2_type == 'COMPOSITE':
        a = 2/3
        b1 = 3.0
        b2 = 4.2
    elif plate1_type == 'METALLIC' and plate2_type == 'METALLIC':
        if fastener.is_bolt:
            a = 2/3
            b1 = 3.0
            b2 = 3.0
        if fastener.is_rivet:
            a = 2/5
            b1 = 2.2
            b2 = 2.2
    else:
        raise ValueError(
            'plate1 and/or plate2 are neither METALLIC nor COMPOSITE types.'
            'Please check'
        )
        return None

    c1 = ((t1 + t2)/(2*d))**a
    c2 = (b1/(t1*E1) + b2/(t2*E2) + b1/(2*t1*Er) + b2/(2*t2*Er))

    c = c1 * c2
    k = 1/c

    return k


def kaxial(plate1: Plate, plate2: Plate, fastener: Fastener) -> float:
    # data from plates
    t1 = plate1.thickness
    t2 = plate2.thickness
    # data from fastener
    Er = fastener.Et
    d = fastener.diameter

    # Computing axial stiffness as ErA/L
    A = pi * (d/2) ** 2
    L = t1 + t2
    k = Er * A / L
    return k


def get_max_edge_length(bdf_model: BDF, eid):
    elem = bdf_model.elements[eid]
    edge_ids = elem.get_edge_ids()
    edge_lengths = []
    for edge in edge_ids:
        node1 = bdf_model.nodes[edge[0]]
        node2 = bdf_model.nodes[edge[1]]
        p1 = Point(*node1.xyz)
        p2 = Point(*node2.xyz)
        v = Vector(p1, p2)
        l = v.norm
        edge_lengths.append(l)
    max_l = max(edge_lengths)
    return max_l


def get_2Delement_planes(bdf_model: BDF) -> dict[Plane]:
    # Genero la información de todos los planos que tiene el modelo. En el caso
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
                # primera triangulación
                n1 = e.nodes_ref[0]
                n2 = e.nodes_ref[1]
                n3 = e.nodes_ref[2]
                p1 = Point(*n1.xyz)
                p2 = Point(*n2.xyz)
                p3 = Point(*n3.xyz)
                s = Plane(p1, p2, p3, eid, tol)
                planes[f'{e.type}_{eid}_S0'] = s

                # segunda triangulación
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
    nids = set(model.nodes.keys())
    if not nids:
        return None
    else:
        return max(nids)


def get_max_mpc_id(bdf_model: BDF) -> int | None:
    mpcids = set(model.mpcs.keys())
    if not mpcids:
        return None
    else:
        return max(mpcids)


def get_max_element_id(bdf_model: BDF) -> int | None:
    eids = set()
    eids.update(model.elements.keys())
    eids.update(model.rigid_elements.keys())
    eids.update(model.masses.keys())

    if not eids:
        return None
    else:
        return max(eids)


def get_max_property_id(bdf_model: BDF) -> int | None:
    pids = set(model.properties.keys())
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
                elem = model.elements[eid]
                # print(elem.node_ids)
                for nid in elem.node_ids:
                    adj.update(node_to_elements_map[nid])
            eids = list(adj)

        return eids
    except:
        return None


def get_wagon_nodes(bdf_models: list[BDF],
                    ref_nid: int,
                    ref_eid: int,
                    radius: float,
                    auto_increase: bool = True) -> list[int]:
    '''
    Gets the nodes close to a reference node within a radius and are adjacent
    to the a reference element.

    Parameters
    ----------
    bdf_models : list[BDF]
        BDF models. It is possible that all the FEM information might be split
        into different BDF models.
    ref_nid : int
        ID of the reference node.
    ref_eid : int
        ID of the reference element.
    radius : float
        Radius of the sphere of influence.
    auto_increase : bool, optional
        If True, in the case that with the given radius no node is within the
        sphere of influence, it will increase the radius by 10%. It will repeat
        this process until it gets at least one node to attach.
        The default is True.

    Returns
    -------
    wagon_nodes : list[int]
        List of ids of the wagon nodes.

    '''

    wagon_nodes = set()
    nids = set()
    nodes_in_sphere = set()
    new_adj_nodes = set()
    old_adj_nodes = set()
    added_adj_nodes = set()

    # Nodo de referencia del wagon para calcular las distancias. Este nodo
    # puede encontrarse en otro BDF model distinto al de partida
    for bdf_model in bdf_models:
        if ref_nid in bdf_model.node_ids:
            ref_grid = bdf_model.nodes[ref_nid]

    # Buscamos los nodos del modelo de partida. Estos nodos estarán en el
    # mismo BDF model que el elemento de referencia "ref_eid"
    for bdf_model in bdf_models:
        if ref_eid in bdf_model.element_ids:
            nodes = bdf_model.nodes
            nids = bdf_model.node_ids

    # Por cada nodo, calculo cuales de ellos estarían en una esfera influencia
    # de radio R
    while True:
        for nid in nids:
            grid = nodes[nid]
            v = Vector(Point(*ref_grid.xyz), Point(*grid.xyz))
            d = v.norm
            if d <= radius:
                nodes_in_sphere.update([nid])

        # print(f'Ref nid : {ref_nid}')
        # print(f'Ref eid: {ref_eid}')
        # print(f'Nodes in Sphere : {nodes_in_sphere}')
        # Buscamos los nodos adyacentes al elemento de referencia.
        for bdf_model in bdf_models:
            if ref_eid in bdf_model.element_ids:
                adj_level = 0
                while True:
                    # print(adj_level)
                    adj_eids = get_adjacent_elements(
                        bdf_model, [ref_eid], adj_level)
                    # print(adj_eids)

                    for adj_eid in adj_eids:
                        adj_elem = bdf_model.elements[adj_eid]
                        new_adj_nodes.update(adj_elem.node_ids)
                    # print(new_adj_nodes)
                    # print(old_adj_nodes)

                    added_adj_nodes = new_adj_nodes.difference(old_adj_nodes)

                    if len(added_adj_nodes.intersection(nodes_in_sphere)) == 0:
                        break
                    else:
                        old_adj_nodes = new_adj_nodes.copy()
                        adj_level += 1
                        # _ = input('Press Enter to continue')
            else:
                break

        # Los Wagon nodes serán la intersección de los nodos que pertenecen al
        # set de la esfera y al set de los nodos adyacentes.
        wagon_nodes = list(set.intersection(nodes_in_sphere, old_adj_nodes))

        if len(wagon_nodes) == 0:
            print(
                f'Unable to coupling nodes around reference Node {nid}'
                f' within a radius of {radius:.2f}.')
            if auto_increase:
                radius = 1.10*radius
            else:
                user_input = input('Set new radius : ')
                radius = float(user_input)
            print(f'Repeating the process using a radius of {radius:.2f}')
        else:
            return wagon_nodes


def create_fem_fasteners(
        bdf_model: BDF,
        connectors: dict[Connector],
        fastener_bdf_filename: str,
        wagon_radius_auto_increase=True,
        skip_connection_tol=1.25
):

    model = bdf_model

    # Creating a new BDF model to generate the
    include_file = BDF()

    # --- GENERATING THE CBUSH
    # Getting the last Node ID
    last_nid = get_max_node_id(model)

    # Getting the last Elem ID (Includes 1D, 2D and Spring elements)
    last_eid = get_max_element_id(model)

    # Getting the last PID
    last_pid = get_max_property_id(model)

    nid = last_nid
    eid = last_eid
    pid = last_pid
    planes = get_2Delement_planes(model)
    for conn_id, conn in connectors.items():

        intersections = conn.intersections
        fastener = conn.linked_fastener

        n = len(intersections)
        print(f'CONNECTOR {conn_id}. NUMBER OF INTERSECTIONS : {n}')
        print(f'================================================')

        nids = [last_nid+i for i in range(1, n+1)]
        coup_ids = [last_eid+i for i in range(n, 2*n)]
        cbush_ids = [last_eid+i for i in range(1, n)]
        pbush_ids = [last_pid for i in range(1, n)]

        for i in range(n):

            # Generating intersections NODES
            pt, plate = intersections[i]
            grid = include_file.add_grid(
                nid=nids[i],
                xyz=pt.to_list())
            # Generating intersections COUPLINGS
            wagon_nodes_plate = get_wagon_nodes(
                bdf_models=[model, include_file],
                ref_nid=nids[i],
                ref_eid=plate.linked_FEM_eid,
                radius=conn.linked_fastener.diameter/2)
            rbe3_plate = include_file.add_rbe3(
                eid=coup_ids[i],
                refgrid=nids[i],
                refc=123456,
                weights=[1.0 for i in wagon_nodes_plate],
                comps=[123 for i in wagon_nodes_plate],
                Gijs=wagon_nodes_plate)

            if i < n-1:
                pt1, plate1 = intersections[i]
                pt2, plate2 = intersections[i+1]

                # Computing stiffnesses
                k_shear = kshear_huth(plate1, plate2, fastener)
                k_axial = kaxial(plate1, plate2, fastener)

                # Generating the property
                pbush = include_file.add_pbush(
                    pid=pbush_ids[i],
                    k=[k_axial, k_shear, k_shear, 100.0, 1e+7, 1e+7],
                    b=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                    ge=[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                )
                # Generating cbush
                cbush = include_file.add_cbush(
                    eid=cbush_ids[i],
                    pid=pbush_ids[i],
                    nids=[nids[i], nids[i+1]],
                    x=fastener.orientation_vector.to_list(),
                    g0=None)

                # Printing Info
                print(f'\t Connection {i}-{i+1} data:')
                print(f'\t --------------------------')
                print(f'\t CBUSH ID : {cbush_ids[i]}')
                print(f'\t PBUSH ID : {pbush_ids[i]}')
                print(f'\t Node 1 ID : {nids[i]}')
                print(f'\t Plate 1 EID : {plate1.linked_FEM_eid}')
                print(f'\t Coupling 1 EID : {coup_ids[i]}')
                print(f'\t Node 2 ID : {nids[i+1]}')
                print(f'\t Plate 2 EID : {plate2.linked_FEM_eid}')
                print(f'\t Coupling 2 EID : {coup_ids[i+1]}')
                print('\n')

        last_nid += n+1
        last_eid += 2*n-1
        last_pid += n

    # Writing BDF with the fasteners data
    _ = include_file.write_bdf(fastener_bdf_filename)


# -------------------------------------------------------------------TESTING
Al_7075 = Metallic(name='Al_7075',
                   norm=None,
                   shape='Sheet',
                   temp_treat='T3',
                   area_range=None,
                   thick_range='0.53-1.57',
                   basis='A',
                   Et=72395,
                   Ec=73744,
                   Ftu=413.68,
                   Fty=303.37,
                   Fcy=248.21,
                   Fsu=255.10,
                   Fbru_ed_1p5=None,
                   Fbry_ed_1p5=None,
                   Fbru_ed_2p0=None,
                   Fbry_ed_2p0=None,
                   e=0.12
                   )
plate1 = MetallicPlate(Al_7075, 2.0)
plate2 = MetallicPlate(Al_7075, 3.0)
fast = Fastener(Et=100000.0, diameter=4.0, head_type='PTR')

# ------------------------------------------------------------------- FAST GEN

launcher_path = r'C:\Users\U69432\Desktop\PyAir\FEM\03_FEM_fastener_generator\Test\01_HF_2D_DFEM\HF_2D_DFem_v04_RBE2_SS_d70mm.bdf'
fastener_filename = r'C:\Users\U69432\Desktop\v04\01_FEM\Test_fastgen_v2.bdf'
model = BDF()
_ = model.read_bdf(
    bdf_filename=launcher_path,
    validate=True,
    xref=True,
    punch=False,
    read_includes=True,
    save_file_structure=False)

planes = get_2Delement_planes(model)
# # Vectores directores
v_dir1 = Vector(
    Point(14271.50, -632.474, 2784.008),
    Point(14271.50, -632.058, 2800.006))
x_dir = Vector(
    Point(0.0, 0.0, 0.0),
    Point(1.0, 0.0, 0.0))
y_dir = Vector(
    Point(0.0, 0.0, 0.0),
    Point(0.0, 1.0, 0.0))
z_dir = Vector(
    Point(0.0, 0.0, 0.0),
    Point(0.0, 0.0, 1.0))

# Remaches
# fast_pr6115_3 = Fastener(
#     110000, 4.8, HeadType.PTR)

fast1 = Fastener(
    Et=110000.0,
    diameter=6.4,
    ref_point=Point(14319.03, -592.526, 2803.587),
    shank_vector=v_dir1,
    orientation_vector=x_dir,
    fastener_label='prEN6115-3'
)

fast2 = Fastener(
    Et=110000.0,
    diameter=6.4,
    ref_point=Point(14318.89, -613.139, 2803.577),
    shank_vector=v_dir1,
    orientation_vector=x_dir,
    fastener_label='prEN6115-3'
)

fast3 = Fastener(
    Et=110000.0,
    diameter=6.4,
    ref_point=Point(14289.22, -593.387, 2801.543),
    shank_vector=v_dir1,
    orientation_vector=x_dir,
    fastener_label='prEN6115-3'
)

fast4 = Fastener(
    Et=110000.0,
    diameter=6.4,
    ref_point=Point(14229.33, -656.673, 2883.993),
    shank_vector=y_dir,
    orientation_vector=x_dir,
    fastener_label='prEN6115-3'
)

fast5 = Fastener(
    Et=110000.0,
    diameter=6.4,
    ref_point=Point(14283.77, -633.991, 2801.181),
    shank_vector=z_dir,
    orientation_vector=y_dir,
    fastener_label='prEN6115-3'
)

fast6 = Fastener(
    Et=110000.0,
    diameter=6.4,
    ref_point=Point(14262.51, -594.986, 2881.882),
    shank_vector=x_dir,
    orientation_vector=y_dir,
    fastener_label='prEN6115-3'
)

fast7 = Fastener(
    Et=110000,
    diameter=4.8,
    ref_point=Point(14215.57, -606.606, 2796.503),
    shank_vector=z_dir,
    orientation_vector=x_dir,
    fastener_label='Ejemplo_Adri',
)

# Connectores
conn1 = Connector(
    linked_fastener=fast1,
    bdf_model=model,
    planes=planes,
    connector_id=1)

conn2 = Connector(
    linked_fastener=fast2,
    bdf_model=model,
    planes=planes,
    connector_id=2)

conn3 = Connector(
    linked_fastener=fast3,
    bdf_model=model,
    planes=planes,
    connector_id=3)

conn4 = Connector(
    linked_fastener=fast4,
    bdf_model=model,
    planes=planes,
    connector_id=4)

conn5 = Connector(
    linked_fastener=fast5,
    bdf_model=model,
    planes=planes,
    connector_id=5)

conn6 = Connector(
    linked_fastener=fast6,
    bdf_model=model,
    planes=planes,
    connector_id=6)
conn7 = Connector(
    linked_fastener=fast7,
    bdf_model=model,
    planes=planes,
    connector_id=7)


connectors = {
    conn1.connector_id: conn1,
    conn2.connector_id: conn2,
    conn3.connector_id: conn3,
    conn4.connector_id: conn4,
    conn5.connector_id: conn5,
    conn6.connector_id: conn6,
    conn7.connector_id: conn7
}

create_fem_fasteners(model, connectors, fastener_filename)
