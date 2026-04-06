
import os
from common.general import split_df_by_subcases
from common.general import get_result_from_op2
from common.general import load_op2_file
from common.ui.tui import ask_entities_ids_list
from common.ui.filedialogs import askOpenOP2
from common.ui.filedialogs import askSaveAsBDF
from common.io.parsers import groups_parser
from pyNastran.bdf.bdf import BDF as bdf
__version__ = '1.0.0'

modules_versions = {
    'PyNastran_OP2_tools': '1.0.2'
}


format_int = '>8s'
format_float = '> 8.4f'
empty_field = 8*' '

bdf = bdf()
dof_dict = dict(zip([1, 2, 3, 4, 5, 6], ['t1', 't2', 't3', 'r1', 'r2', 'r3']))


def spcds_from_fem(op2_file_path: 'str',
                   nids: list[int],
                   dofs: int,
                   spc1_id: int = 123456,
                   spcd_file_path: str | None = None,
                   verbose: bool = True):
    """
    Given a NASTRAN OP2 file a the IDs of the nodes, reads the displacments for
    all the subcases and nodes, and writes the associated SPCDs and SPC1 
    NASTRAN cards into a file.

    The final purpose is to use the output file within a DFEM to load it using 
    enforced displacements.

    Parameters
    ----------

    op2_file_path : 'str'
        Path of the NASTRAN OP2 file

    nids : list[int]
        List of the IDs of the nodes which the displacements will be read. 

    dofs : int
        Degrees of freedom.

    spc1_id : 'int', optional
        ID of the SPC card

    spcd_file_path : 'path-like str', optional
        Path of the output file containing the SPCDs and SPC1 cards. If no path 
        is provided, by default the name of the ouput file will be 'spcds_file.bdf'
        and it will be generated in the same directory as the OP2 file.

    Returns
    -------
    Nastran file containing the SPCDs and SPC1s cards.

    """

    op2 = load_op2_file(op2_file_path)

    dirname = os.path.dirname(op2_file_path)
    if spcd_file_path is None:
        os.chdir(dirname)
        spcd_file_path = os.path.join(dirname, 'scpd_file.bdf')
    spcd_file = open(spcd_file_path, 'w')

    disp = get_result_from_op2(op2, 'displacements')
    sids = disp['SubcaseID'].unique()
    disp_tables = split_df_by_subcases(disp, sids)

    nids.sort()
    spcd_file.write('$$---- SPC1 Constraints\n')
    for nid in nids:
        spc1_card = bdf.add_card(
            card_lines=['SPC1', spc1_id, dofs, nid],
            card_name='SPC1')
        spcd_file.write(spc1_card.write_card())

    for sid, table in disp_tables.items():
        if verbose:
            print('Generating SPCD for Subcase {0}'.format(sid))
        spcd_file.write('$$---- SPCDs for SUBCASE ID {0}\n'.format(sid))
        for nid in nids:

            try:
                nid_line = table[table['NodeID'] == nid]
            except:
                if verbose:
                    print('\t[WARN] Unable to find node {0}'.format(nid))
                continue

            for dof in list(str(dofs)):
                disp = nid_line[dof_dict[int(dof)]].values
                card_lines = [
                    'SPCD',
                    str(sid),
                    str(nid),
                    str(dof),
                    format(disp[0], format_float)
                ]
                spcd_card = bdf.add_card(
                    card_lines=card_lines,
                    card_name='SPCD'
                )
                spcd_file.write(spcd_card.write_card())

    spcd_file.close()


if __name__ == '__main__':

    op2path = askOpenOP2()
    nids = ask_entities_ids_list()

    dofs = input(
        'Enter the degrees of freedom to extract the displacements:\n')

    spcds_from_fem(
        op2_file_path=op2path,
        nids=nids,
        dofs=dofs,
        spc1_id=123456,
        spcd_file_path=askSaveAsBDF()[0],
        verbose=True
    )
