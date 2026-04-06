from enum import Enum
from pyNastran.bdf.bdf import BDF
from geometry.point import Point
from geometry.vector import Vector
from geometry.plane import Plane_by3points as Plane
from geometry.axes import Csys
from materials.composite import Orthotropic
from materials.metallic import Metallic
from structural.fastener import Fastener
from structural.plate import MetallicPlate, CompositePlate, Plate
from structural.laminate import Laminate
from structural.ply import Ply
from fem.pre.core import get_max_edge_length, get_adjacent_elements
from fem.pre.core import get_max_node_id, get_max_element_id
from fem.pre.core import get_max_property_id, get_2Delement_planes
from fem.pre.core import get_angle_matsys_fast
from math import pi
from tqdm import tqdm
import pandas as pd


__version__ = '1.1.0'

# POR COMPLETAR/AÑADIR
# -----------------------------------------------------------------------------
# COSAS POR AÑADIR/COMPLETAR:
#   1. Que lea la parte de compuesto --> HECHO
#   2. Que las direcciones de CBUSH pueda meterlo en formato de vector o
#      de sistema de coordenadas. Ahora mismo solo esta con un vector.
#   3. Que el usario, pueda elegir si crear un sistema de coordenadas
#      asociado al fastener o no.
#   4. Que el usuario pueda poner los ids que quiera para los nodos, elementos
#      y propiedades creados.
#   5. Generar una tabla con todos los datos para su posterior revision
#
# POR CORREGIR
# -----------------------------------------------------------------------------
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
#
# POSBILES MEJORAS
# -----------------------------------------------------------------------------
#   1. A la hora de crear los connectores tarda bastante. Seria preferible que
#   fuera m�s r�pido.
#   2. Que se pueda meter la info  a trav�s de una excel.


# =============================================================================
#                                                                             #
#                              *** CLASSES ***                                #
#                                                                             #
# =============================================================================


class JointType(Enum):
    SINGLE = 'single_shear'
    DOUBLE = 'double_shear'


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
        self.skip_tol = skip_tol
        self.intersections = self.__filter_intersections(bdf_model, planes)

    def __str__(self):
        msg = f'Connector to FEM element ids {self.linked_eids}'
        return msg

    def __get_laminate_from_PCOMP(self, pid, bdf_model: BDF) -> Laminate:
        pcomp = bdf_model.properties[pid]
        pcomp_plies_data = pcomp.plies
        stacking = []
        for ply_data in pcomp_plies_data:
            ply_mat_id, ply_t, ply_angle, sout = ply_data
            mat = bdf_model.materials[ply_mat_id]
            ply = Ply(
                Orthotropic(
                    name='',
                    thickness=ply_t,
                    E1=mat.E11(),
                    E2=mat.E22(),
                    G12=mat.G12(),
                    nu12=mat.nu12
                ),
                ply_angle)
            stacking.append(ply)
        return Laminate(stacking)

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
                            laminate = self.__get_laminate_from_PCOMP(
                                inter_elem.pid, bdf_model)
                            plate = CompositePlate(
                                laminate=laminate,
                                linked_FEM_eid=plane.plane_id,
                                linked_FEM_pid=inter_elem.pid
                            )
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

                if conn_length > self.skip_tol*t_total:
                    intersections.pop(i+1)
                    flags.append(True)
                    break
                else:
                    flags.append(False)
            cond = any(flags)

        return intersections


# =============================================================================
#                                                                             #
#                             *** FUNCTIONS ***                               #
#                                                                             #
# =============================================================================


def kshear_huth(plate1: Plate,
                plate2: Plate,
                fastener: Fastener,
                angle_plate1: float = 0.0,
                angle_plate2: float = 0.0) -> float:

    # data from plate 1
    if type(plate1) == MetallicPlate:
        plate1_type = 'METALLIC'
        t1 = plate1.thickness
        E1_local = plate1.material.Et
    else:
        plate1_type = 'COMPOSITE'
        t1 = plate1.thickness
        E1_local, Ey_local, Gxy_local, nuxy_local = plate1.get_equivalent_moduli(
            angle_plate1)

    # data from plate 1
    if type(plate2) == MetallicPlate:
        plate2_type = 'METALLIC'
        t2 = plate2.thickness
        E2_local = plate2.material.Et
    else:
        plate2_type = 'COMPOSITE'
        t2 = plate2.thickness
        E2_local, Ey_local, Gxy_local, nuxy_local = plate1.get_equivalent_moduli(
            angle_plate2)

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
    c2 = (b1/(t1*E1_local) + b2/(t2*E2_local) + b1/(2*t1*Er) + b2/(2*t2*Er))

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
        fastener_vector = fastener.orientation_vector

        n = len(intersections)
        print(f'CONNECTOR {conn_id}. NUMBER OF INTERSECTIONS : {n}')
        print(f'================================================')

        nids = [last_nid+i for i in range(1, n+1)]
        coup_ids = [last_eid+i for i in range(n, 2*n)]
        cbush_ids = [last_eid+i for i in range(1, n)]
        pbush_ids = [last_pid+i for i in range(1, n)]

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
                eid1 = plate1.linked_FEM_eid
                eid2 = plate2.linked_FEM_eid
                angle1 = get_angle_matsys_fast(
                    bdf_model, eid1, fastener_vector)
                angle2 = get_angle_matsys_fast(
                    bdf_model, eid1, fastener_vector)
                k_shear_1 = kshear_huth(
                    plate1, plate2, fastener, angle1, angle2)
                k_shear_2 = kshear_huth(
                    plate1, plate2, fastener, angle1+90, angle2+90)
                k_axial = kaxial(plate1, plate2, fastener)

                # Generating the property
                pbush = include_file.add_pbush(
                    pid=pbush_ids[i],
                    k=[k_axial, k_shear_1, k_shear_2, 100.0, 1e+7, 1e+7],
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


def __from_excel_to_df(file):
    """
    Read an excel template file with the data of the fasteners to be modelised
    in Nastran and returns :
        - The Nastran Launcher path
        - The Fastener Output file 
        - A Pandas DataFrame with all the information of the excel file.

    """
    FASTENERS_SHEET_LABEL = 'FASTENERS'
    PILOT_POINTS_SHEET_LABEL = 'PILOT_POINTS'
    CONNECTORS_SHEET_LABEL = 'CONNECTORS'

    # Fasteners DataFrame
    fasteners_table = pd.read_excel(
        io=file,
        sheet_name=FASTENERS_SHEET_LABEL
    )
    # PilotPoints DataFrame
    pilotpoints_table = pd.read_excel(
        io=file,
        sheet_name=PILOT_POINTS_SHEET_LABEL
    )
    # Connectors DataFrame
    connectors_table = pd.read_excel(
        io=file,
        header=2,
        sheet_name=CONNECTORS_SHEET_LABEL
    )

    # COMBINAMOS LAS TABLAS
    # ---------------------
    # Creo una copia de connectors que es quien contiene toda la informacion
    # de manera condensada
    merged_table = connectors_table.copy(deep=True)

    # Para tener luego una tabla mejor organizada, quito temporalmente la
    # columna "Fastener_Label".
    fastener_label_col = merged_table.pop('Fastener_Label')

    # Incluyo la informaci�n de PilotPoints usando como link los labels
    merged_table = merged_table.merge(
        right=pilotpoints_table,
        left_on='PilotPoint_Label',
        right_on='PilotPoint_Label',
        how='left',
        copy=True)

    # Meto a la derecha de la tabla la columna 'Fastener_Label'
    merged_table.insert(
        len(merged_table.columns),
        'Fastener_Label',
        fastener_label_col)

    # Incluyo la informaci�n de Fasteners usando como link los labels
    merged_table = merged_table.merge(
        right=fasteners_table,
        left_on='Fastener_Label',
        right_on='Fastener_Label',
        how='left',
        copy=False
    )

    # LEEMOS LA RUTA DEL BDF Y LA DEL FICHERO DE SALIDA
    bdf_paths = pd.read_excel(
        io=file,
        sheet_name=CONNECTORS_SHEET_LABEL,
        header=None,
        usecols="B:B",
        nrows=2
    )
    nastran_launcher = bdf_paths.iloc[0, 0]
    fastener_output_file = bdf_paths.iloc[1, 0]
    return nastran_launcher, fastener_output_file, merged_table


def __parse_df_to_fasteners_list(df: pd.DataFrame):
    """
    Transforms the data from the excel template into a list of Fasteners 
    objects.
    """
    fasteners_list = []
    for row in df.itertuples(index=True):
        create = row.Create
        if create.upper() == 'YES':
            ref_point = Point(
                coord1=row.X_ref,
                coord2=row.Y_ref,
                coord3=row.Z_ref,
            )
            shank_point = Point(
                coord1=row.X_shank,
                coord2=row.Y_shank,
                coord3=row.Z_shank,
            )
            orientation_point = Point(
                coord1=row.X_vec,
                coord2=row.Y_vec,
                coord3=row.Z_vec,
            )

            shank_vector = Vector(
                A=ref_point,
                B=shank_point
            )

            orientation_vector = Vector(
                A=ref_point,
                B=orientation_point
            )

            fastener = Fastener(
                Et=row.E,
                diameter=row.Diameter,
                ref_point=ref_point,
                shank_vector=shank_vector,
                orientation_vector=orientation_vector,
                fastener_label=row.Fastener_Label
            )
            fasteners_list.append(fastener)
    return fasteners_list


def from_excel(file: str, verbose=False):
    """
    Create a Nastran BDF File with the modelistaion of the fasteners. The input
    is an excel file template.

    Note: The current version of the template

    Returns
    -------
    None.

    """
    # Parsing data from Excel file
    if verbose:
        print('Reading Excel file ...', end='')
    nastran_launcher, fastener_output_file, merged_table = __from_excel_to_df(
        file)
    fasteners_list = __parse_df_to_fasteners_list(merged_table)
    if verbose:
        print('Done')

    # LOADING NASTRAN BDF MODEL
    if verbose:
        print('Loading NASTRAN Model ...', end='')
    model = BDF()
    _ = model.read_bdf(
        bdf_filename=nastran_launcher,
        validate=True,
        xref=True,
        punch=False,
        read_includes=True,
        save_file_structure=False)
    if verbose:
        print('Done')

    # CREATING THE CONNECTORS
    if verbose:
        print('Creating Connectors...', end='')
    connectors_dict = {}
    planes = get_2Delement_planes(model)
    number_fasteners = len(fasteners_list)
    for i in tqdm(range(number_fasteners)):
        connectors_dict[i+1] = Connector(
            linked_fastener=fasteners_list[i],
            bdf_model=model,
            planes=planes)

    # CREATING THE BDF FILE FILE WITH THE CONNECTORS
    create_fem_fasteners(
        bdf_model=model,
        connectors=connectors_dict,
        fastener_bdf_filename=fastener_output_file)
