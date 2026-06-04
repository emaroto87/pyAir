# -*- coding: utf-8 -*-
# Modulos PyAir
from fem.pre.create_fasteners import create_fem_fasteners, Connector
from geometry.point import Point
from geometry.vector import Vector
from structural.fastener import Fastener
from fem.pre.core import get_2Delement_planes
import common.errors.errors as err
# Modulos de terceros
from pyNastran.bdf.bdf import BDF
from packaging import version

# # CONTROL DE VERSIONES DE MODULOS IMPORTADOS
# CREATE_FEM_FASTENERS_MIN_VERSION = '1.0.13'
# PYNASTRAN_MIN_VERSION = '1.4.1'

# try:
#     err.check_module_version(fem.pre.create_fasteners, "1.0.13")
# except err.ModuleVersionError as e:
#     print('Error:', e)
#     raise SystemExit(1)

# try:
#     err.check_module_version("pyNASTRAN", "1.4.1")
# except err.ModuleVersionError as e:
#     print('Error:', e)
#     raise SystemExit(1)


# INPUTS
launcher_path = r"C:\Users\U69432\Desktop\PyAir\v.0.0.2\data\fems\VentralFin_Node4_DFEM_v2\DFEM_VentralFins_Node4_v2.dat"
fastener_filename = r"C:\Users\U69432\Desktop\PyAir\v.0.0.2\data\fems\VentralFin_Node4_DFEM_v2\01_FEMFasteners_v2_bis.bdf"

# LOADING NASTRAN BDF MODEL
model = BDF()
_ = model.read_bdf(
    bdf_filename=launcher_path,
    validate=True,
    xref=True,
    punch=False,
    read_includes=True,
    save_file_structure=False)

# CREATING THE PLANES FOR ALL THE ELEMENTS
# (Es necesario para luego crear los connectores)
planes = get_2Delement_planes(model)


# DEFINICION DE LOS PILOTS POINTS
# =============================================================================
#
# IF FWD_CLEAT-LONGERON
#
# =============================================================================
#
# Pilot Points
if_fwd_cleat_longeron_p01 = Point(
    coord1=16482.235,
    coord2=-254.731,
    coord3=2697.760,
)
if_fwd_cleat_longeron_p02 = Point(
    coord1=16465.674,
    coord2=-255.380,
    coord3=2693.642,
)
if_fwd_cleat_longeron_shank_vec = Vector(
    Point(
        coord1=16462.674,
        coord2=-255.380,
        coord3=2693.642,
    ),
    Point(
        coord1=16462.971,
        coord2=-250.865,
        coord3=2691.515,
    )
)
if_fwd_cleat_longeron_local_y = Vector(
    Point(
        coord1=16323.653,
        coord2=-256.664,
        coord3=2673.931,
    ),
    Point(
        coord1=16327.575,
        coord2=-256.552,
        coord3=2674.743
    )
)
#
# Fasteners
if_fwd_cleat_longeron_fast01 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_fwd_cleat_longeron_p01,
    shank_vector=if_fwd_cleat_longeron_shank_vec,
    orientation_vector=if_fwd_cleat_longeron_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_fwd_cleat_longeron_conn01 = Connector(
    linked_fastener=if_fwd_cleat_longeron_fast01,
    bdf_model=model,
    planes=planes,
    connector_id=1)
if_fwd_cleat_longeron_fast02 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_fwd_cleat_longeron_p02,
    shank_vector=if_fwd_cleat_longeron_shank_vec,
    orientation_vector=if_fwd_cleat_longeron_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_fwd_cleat_longeron_conn02 = Connector(
    linked_fastener=if_fwd_cleat_longeron_fast02,
    bdf_model=model,
    planes=planes,
    connector_id=2)

# =============================================================================
#
# IF RWD_CLEAT-longeron
#
# =============================================================================
#
# Pilot Points
if_rwd_cleat_longeron_p01 = Point(
    coord1=16370.909,
    coord2=-259.60,
    coord3=2675.358,
)
if_rwd_cleat_longeron_p02 = Point(
    coord1=16351.336,
    coord2=-260.171,
    coord3=2671.282,
)
if_rwd_cleat_longeron_shank_vec = Vector(
    Point(
        coord1=16371.146,
        coord2=-254.940,
        coord3=2673.562,
    ),
    Point(
        coord1=16370.909,
        coord2=-259.60,
        coord3=2675.358,
    )
)
if_rwd_cleat_longeron_local_y = Vector(
    Point(
        coord1=16323.653,
        coord2=-256.664,
        coord3=2673.931,
    ),
    Point(
        coord1=16327.575,
        coord2=-256.552,
        coord3=2674.743
    )
)
#
# Fasteners
if_rwd_cleat_longeron_fast01 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_rwd_cleat_longeron_p01,
    shank_vector=if_rwd_cleat_longeron_shank_vec,
    orientation_vector=if_rwd_cleat_longeron_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_rwd_cleat_longeron_conn01 = Connector(
    linked_fastener=if_rwd_cleat_longeron_fast01,
    bdf_model=model,
    planes=planes,
    connector_id=3)

if_rwd_cleat_longeron_fast02 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_rwd_cleat_longeron_p02,
    shank_vector=if_rwd_cleat_longeron_shank_vec,
    orientation_vector=if_rwd_cleat_longeron_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_rwd_cleat_longeron_conn02 = Connector(
    linked_fastener=if_rwd_cleat_longeron_fast02,
    bdf_model=model,
    planes=planes,
    connector_id=4)


# =============================================================================
#
# IF FWD_CLEAT-FRAME
#
# =============================================================================
#
# Pilot Points
if_fwd_cleat_frame_p01 = Point(
    coord1=16401.603,
    coord2=-257.537,
    coord3=2684.779,
)
if_fwd_cleat_frame_p02 = Point(
    coord1=16401.603,
    coord2=-250.706,
    coord3=2702.508,
)
if_fwd_cleat_frame_p03 = Point(
    coord1=16401.603,
    coord2=-243.874,
    coord3=2720.238,
)
if_fwd_cleat_frame_shank_vec = Vector(
    Point(
        coord1=16462.674,
        coord2=-255.380,
        coord3=2693.642,
    ),
    Point(
        coord1=16462.971,
        coord2=-250.865,
        coord3=2691.515,
    )
)
if_fwd_cleat_frame_local_y = Vector(
    Point(
        coord1=16323.653,
        coord2=-256.664,
        coord3=2673.931,
    ),
    Point(
        coord1=16327.575,
        coord2=-256.552,
        coord3=2674.743
    )
)
# Fasteners
# Fast01
if_fwd_cleat_frame_fast01 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_fwd_cleat_frame_p01,
    shank_vector=if_fwd_cleat_frame_shank_vec,
    orientation_vector=if_fwd_cleat_frame_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_fwd_cleat_frame_conn01 = Connector(
    linked_fastener=if_fwd_cleat_frame_fast01,
    bdf_model=model,
    planes=planes,
    connector_id=5)
# Fast02
if_fwd_cleat_frame_fast02 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_fwd_cleat_frame_p02,
    shank_vector=if_fwd_cleat_frame_shank_vec,
    orientation_vector=if_fwd_cleat_frame_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_fwd_cleat_frame_conn02 = Connector(
    linked_fastener=if_fwd_cleat_frame_fast02,
    bdf_model=model,
    planes=planes,
    connector_id=6)
# Fast03
if_fwd_cleat_frame_fast03 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_fwd_cleat_frame_p03,
    shank_vector=if_fwd_cleat_frame_shank_vec,
    orientation_vector=if_fwd_cleat_frame_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_fwd_cleat_frame_conn03 = Connector(
    linked_fastener=if_fwd_cleat_frame_fast03,
    bdf_model=model,
    planes=planes,
    connector_id=7)


# =============================================================================
#
# IF RWD_CLEAT-FRAME
#
# =============================================================================
#
# Pilot Points
if_rwd_cleat_frame_p01 = Point(
    coord1=16433.219,
    coord2=-255.051,
    coord3=2690.217,
)
if_rwd_cleat_frame_p02 = Point(
    coord1=16433.219,
    coord2=-246.526,
    coord3=2708.310,
)
if_rwd_cleat_frame_p03 = Point(
    coord1=16433.219,
    coord2=-238.002,
    coord3=2726.402,
)
if_rwd_cleat_frame_shank_vec = Vector(
    Point(
        coord1=16371.146,
        coord2=-254.940,
        coord3=2673.562,
    ),
    Point(
        coord1=16370.909,
        coord2=-259.60,
        coord3=2675.358,
    )
)
if_rwd_cleat_frame_local_y = Vector(
    Point(
        coord1=16323.653,
        coord2=-256.664,
        coord3=2673.931,
    ),
    Point(
        coord1=16327.575,
        coord2=-256.552,
        coord3=2674.743
    )
)
# Fasteners
# Fast01
if_rwd_cleat_frame_fast01 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_rwd_cleat_frame_p01,
    shank_vector=if_rwd_cleat_frame_shank_vec,
    orientation_vector=if_rwd_cleat_frame_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_rwd_cleat_frame_conn01 = Connector(
    linked_fastener=if_rwd_cleat_frame_fast01,
    bdf_model=model,
    planes=planes,
    connector_id=8)
# Fast02
if_rwd_cleat_frame_fast02 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_rwd_cleat_frame_p02,
    shank_vector=if_rwd_cleat_frame_shank_vec,
    orientation_vector=if_rwd_cleat_frame_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_rwd_cleat_frame_conn02 = Connector(
    linked_fastener=if_rwd_cleat_frame_fast02,
    bdf_model=model,
    planes=planes,
    connector_id=9)
# Fast03
if_rwd_cleat_frame_fast03 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_rwd_cleat_frame_p03,
    shank_vector=if_rwd_cleat_frame_shank_vec,
    orientation_vector=if_rwd_cleat_frame_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_rwd_cleat_frame_conn03 = Connector(
    linked_fastener=if_rwd_cleat_frame_fast03,
    bdf_model=model,
    planes=planes,
    connector_id=10)

# =============================================================================
#
# IF FRAME-LONGERON-FITTINT
#
# =============================================================================
#
# Pilot Points
if_frame_longeron_p01 = Point(
    coord1=16434.523,
    coord2=-286.786,
    coord3=2671.521,
)

if_frame_longeron_shank_vec_p01 = Vector(
    Point(
        coord1=16434.523,
        coord2=-286.786,
        coord3=2671.521,
    ),
    Point(
        coord1=16433.471,
        coord2=-284.571,
        coord3=2676.086,
    )
)

if_frame_longeron_p02 = Point(
    coord1=16406.298,
    coord2=-287.797,
    coord3=2665.606,
)

if_frame_longeron_shank_vec_p02 = Vector(
    Point(
        coord1=16406.298,
        coord2=-287.797,
        coord3=2665.606,
    ),
    Point(
        coord1=16405.250,
        coord2=-285.598,
        coord3=2670.181,
    )
)

if_frame_longeron_p03 = Point(
    coord1=16434.508,
    coord2=-231.432,
    coord3=2648.339,
)

if_frame_longeron_shank_vec_p03 = Vector(
    Point(
        coord1=16434.508,
        coord2=-231.432,
        coord3=2648.339,
    ),
    Point(
        coord1=16432.883,
        coord2=-228.775,
        coord3=2655.791,
    )
)

if_frame_longeron_p04 = Point(
    coord1=16406.277,
    coord2=-232.444,
    coord3=2642.612,
)

if_frame_longeron_shank_vec_p04 = Vector(
    Point(
        coord1=16406.277,
        coord2=-232.444,
        coord3=2642.612,
    ),
    Point(
        coord1=16405.235,
        coord2=-230.740,
        coord3=2647.397,
    )
)

if_frame_longeron_p05 = Point(
    coord1=16434.06,
    coord2=-251.624,
    coord3=2661.601,
)

if_frame_longeron_shank_vec_p05 = Vector(
    Point(
        coord1=16434.06,
        coord2=-251.624,
        coord3=2661.601,
    ),
    Point(
        coord1=16430.82,
        coord2=-245.671,
        coord3=2676.325,
    )
)

# Fasteners
# Fast01
if_frame_longeron_fast01 = Fastener(
    Et=110000.0,
    diameter=7.19,
    ref_point=if_frame_longeron_p01,
    shank_vector=if_frame_longeron_shank_vec_p01,
    orientation_vector=if_fwd_cleat_frame_local_y,
    fastener_label='NAS1790V_D7.19'
)
if_frame_longeron_conn01 = Connector(
    linked_fastener=if_frame_longeron_fast01,
    bdf_model=model,
    planes=planes,
    connector_id=11)

# Fast02
if_frame_longeron_fast02 = Fastener(
    Et=110000.0,
    diameter=7.19,
    ref_point=if_frame_longeron_p02,
    shank_vector=if_frame_longeron_shank_vec_p02,
    orientation_vector=if_fwd_cleat_frame_local_y,
    fastener_label='NAS1790V_D7.19'
)
if_frame_longeron_conn02 = Connector(
    linked_fastener=if_frame_longeron_fast02,
    bdf_model=model,
    planes=planes,
    connector_id=12)

# Fast03
if_frame_longeron_fast03 = Fastener(
    Et=110000.0,
    diameter=7.19,
    ref_point=if_frame_longeron_p03,
    shank_vector=if_frame_longeron_shank_vec_p03,
    orientation_vector=if_fwd_cleat_frame_local_y,
    fastener_label='NAS1790V_D7.19'
)
if_frame_longeron_conn03 = Connector(
    linked_fastener=if_frame_longeron_fast03,
    bdf_model=model,
    planes=planes,
    connector_id=13)

# Fast04
if_frame_longeron_fast04 = Fastener(
    Et=110000.0,
    diameter=7.19,
    ref_point=if_frame_longeron_p04,
    shank_vector=if_frame_longeron_shank_vec_p04,
    orientation_vector=if_fwd_cleat_frame_local_y,
    fastener_label='NAS1790V_D7.19'
)
if_frame_longeron_conn04 = Connector(
    linked_fastener=if_frame_longeron_fast04,
    bdf_model=model,
    planes=planes,
    connector_id=14)

# Fast04
if_frame_longeron_fast05 = Fastener(
    Et=110000.0,
    diameter=4.0,
    ref_point=if_frame_longeron_p05,
    shank_vector=if_frame_longeron_shank_vec_p05,
    orientation_vector=if_fwd_cleat_frame_local_y,
    fastener_label='ABS2323K5_D4.0'
)
if_frame_longeron_conn05 = Connector(
    linked_fastener=if_frame_longeron_fast05,
    bdf_model=model,
    planes=planes,
    connector_id=15)

# =============================================================================
#
# EJECUCION DEL SCRIPT
#
# =============================================================================
print('Creating Connectors...', end='')
connectors_dict = {
    # if_fwd_cleat_longeron_conn01.connector_id: if_fwd_cleat_longeron_conn01,
    # if_fwd_cleat_longeron_conn02.connector_id: if_fwd_cleat_longeron_conn02,
    # if_rwd_cleat_longeron_conn01.connector_id: if_rwd_cleat_longeron_conn01,
    # if_rwd_cleat_longeron_conn02.connector_id: if_rwd_cleat_longeron_conn02,
    # #
    # if_fwd_cleat_frame_conn01.connector_id: if_fwd_cleat_frame_conn01,
    # if_fwd_cleat_frame_conn02.connector_id: if_fwd_cleat_frame_conn02,
    # if_fwd_cleat_frame_conn03.connector_id: if_fwd_cleat_frame_conn03,
    # if_rwd_cleat_frame_conn01.connector_id: if_rwd_cleat_frame_conn01,
    # if_rwd_cleat_frame_conn02.connector_id: if_rwd_cleat_frame_conn02,
    # if_rwd_cleat_frame_conn03.connector_id: if_rwd_cleat_frame_conn03,
    # #
    # if_frame_longeron_conn01.connector_id: if_frame_longeron_conn01,
    # if_frame_longeron_conn02.connector_id: if_frame_longeron_conn02,
    # if_frame_longeron_conn03.connector_id: if_frame_longeron_conn03,
    # if_frame_longeron_conn04.connector_id: if_frame_longeron_conn04,
    if_frame_longeron_conn05.connector_id: if_frame_longeron_conn05,
}
print('Done', end='')
create_fem_fasteners(
    bdf_model=model,
    connectors=connectors_dict,
    fastener_bdf_filename=fastener_filename)
