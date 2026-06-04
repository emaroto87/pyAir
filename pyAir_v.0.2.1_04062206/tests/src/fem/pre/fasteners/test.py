
from pyNastran.bdf.bdf import BDF
from fem.pre.create_fasteners import create_fem_fasteners, Connector
from geometry.point import Point
from geometry.vector import Vector
from geometry.plane import Plane_by3points as Plane
from geometry.axes import Csys
from structural.fastener import Fastener
from common.ui.filedialogs import askOpenBDF
from common.ui.filedialogs import askSaveAsBDF
from fem.pre.core import get_2Delement_planes
from pathlib import Path


def __testing_FastGen_01():
    base = Path.cwd()
    # launcher_path = base / 'input_files' / \
    #     '01_HF_2D_DFEM' / 'HF_2D_DFem_v04_RBE2_SS_d70mm.bdf'
    # fastener_filename = base / 'script_output_files' / 'Test_fastgen_v2.bdf'
    launcher_path = r'C:\Users\U69432\Desktop\PyAir\PyAir\FEM\PreProcess\Scripts\fasteners_tool\__test\input_files\01_Ex1_HF_2D_DFEM\HF_2D_DFem_v04_RBE2_SS_d70mm.bdf'
    fastener_filename = r'C:\Users\U69432\Desktop\PyAir\PyAir\FEM\PreProcess\Scripts\fasteners_tool\__test\script_output_files\01_Ex1_HF_2D_DFEM\Test_fastgen_v2.bdf'
    model = BDF()
    _ = model.read_bdf(
        bdf_filename=launcher_path,
        validate=True,
        xref=True,
        punch=False,
        read_includes=True,
        save_file_structure=False)

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
    planes = get_2Delement_planes(model)
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


def __testing_FastGen_02():
    # Composite Testing
    model_file = r"C:\Users\U69432\Desktop\PyAir\PyAir\FEM\PreProcess\Scripts\fasteners_tool\__test\input_files\02_Ex2_Composite_plates\Fast_generator_Composite_Ex2.bdf"
    fastener_bdf_filename = r"C:\Users\U69432\Desktop\PyAir\PyAir\FEM\PreProcess\Scripts\fasteners_tool\__test\script_output_files\02_Ex2_Composite_plates\fastener_file_v1.bdf"

    model = BDF()
    _ = model.read_bdf(
        bdf_filename=model_file,
        validate=True,
        xref=True,
        punch=False,
        read_includes=True,
        save_file_structure=False)

    planes = get_2Delement_planes(model)

    refpoint1 = Point(85.0, 45.0, 20.0)
    refpoint2 = Point(145.0, 45.0, 20.0)
    refpoint3 = Point(205.0, 45.0, 20.0)
    shank_vector = Vector(
        Point(0.0, 0.0, 0.0),
        Point(0.0, 0.0, 1.0))
    fast_dir = Vector(
        Point(0.0, 0.0, 0.0),
        Point(1.0, 1.0, 0.0))

    fast1 = Fastener(
        Et=110000.0,
        diameter=11.0,
        ref_point=refpoint1,
        shank_vector=shank_vector,
        orientation_vector=fast_dir,
        fastener_label='prEN6115-3'
    )
    conn1 = Connector(
        linked_fastener=fast1,
        bdf_model=model,
        planes=planes,
        connector_id=1
    )

    fast2 = Fastener(
        Et=110000.0,
        diameter=11.0,
        ref_point=refpoint2,
        shank_vector=shank_vector,
        orientation_vector=fast_dir,
        fastener_label='prEN6115-3'
    )
    conn2 = Connector(
        linked_fastener=fast2,
        bdf_model=model,
        planes=planes,
        connector_id=2
    )

    fast3 = Fastener(
        Et=110000.0,
        diameter=11.0,
        ref_point=refpoint3,
        shank_vector=shank_vector,
        orientation_vector=fast_dir,
        fastener_label='prEN6115-3'
    )
    conn3 = Connector(
        linked_fastener=fast3,
        bdf_model=model,
        planes=planes,
        connector_id=3
    )

    connectors = {
        1: conn1,
        2: conn2,
        3: conn3
    }
    create_fem_fasteners(model, connectors, fastener_bdf_filename)
