from dataclasses import dataclass, asdict

__version__ = '0.0.1'


def nastran_format_str(string, align: str = 'right', field_chars=8):
    inp = string
    nc = len(inp)
    spaces = (field_chars-nc)*' '
    if len(spaces) > 0:
        if align == 'left':
            out = inp + spaces
        else:
            out = spaces + inp
    elif len(spaces) == 0:
        out = inp
    else:
        raise ValueError(
            'The length of {string} is greater than {field_chars} characters')

    return out


def nastran_format_int(integer, align='right', field_chars=8):
    inp = str(integer)
    nc = len(inp)
    spaces = (field_chars-nc)*' '
    if len(spaces) > 0:
        if align == 'left':
            out = inp + spaces
        else:
            out = spaces + inp
    elif len(spaces) == 0:
        out = inp
    else:
        raise ValueError(
            'The length of {string} is greater than {field_chars} characters')
    return out


def nastran_format_float(number, format_str='', align='right', field_chars=8):
    inp = format(number, format_str)
    nc = len(inp)
    spaces = (field_chars-nc)*' '

    if len(spaces) > 0:
        if align == 'left':
            out = inp + spaces
        else:
            out = spaces + inp
    elif len(spaces) == 0:
        out = inp
    else:
        raise ValueError(
            'The length of {string} is greater than {field_chars} characters')
    return out


def create_card(card_fields: list, field_chars: int = 8) -> str:
    card = ''
    field_counter = 1

    for item in card_fields:
        if type(item) is str:
            field = nastran_format_str(
                item, field_chars=field_chars, align='left')
        elif type(item) is int:
            field = nastran_format_int(item, field_chars=field_chars)
        elif type(item) is float:
            field = nastran_format_float(
                item, format_str='', field_chars=field_chars)
        elif item is None:
            field = ' ' * field_chars
        else:
            raise ValueError
        card += field
        # print(card)

        if field_counter % 10 == 0:
            print('enter')
            card += '\n'
        field_counter += 1

    return card


def split_into_chunks(lst: list, chunk_size: int):
    # Initialize an empty list to store chunks
    chunks = []

    # Iterate over the list with a step of chunk_size
    for i in range(0, len(lst), chunk_size):
        # Slice the list from index 'i' to 'i + chunk_size' and append it to chunks
        chunks.append(lst[i:i + chunk_size])

    return chunks


def mpc_card_fields(sid: int, GiCiAi: list[tuple]):
    # Nota: Campos del MPC
    # =====================
    # MPC-----SID-----G1------C1------A1------G2------C2------A2------+-------+-------
    # --------+-------G3------C3------A3------G4------C4------A4------+-------+-------

    card_fields = []
    # Dividimos la lista GiCiAi en grupos de dos
    chunks = split_into_chunks(
        lst=GiCiAi,
        chunk_size=2)
    print(chunks)
    for chunk in chunks:
        if len(chunk) == 1:
            GAC_1 = chunk[0]
            G1 = GAC_1[0]
            C1 = GAC_1[1]
            A1 = GAC_1[2]
            card_fields += [None, None, G1, C1,
                            A1, None, None, None, None, None]
        elif len(chunk) == 2:
            GAC_1 = chunk[0]
            G1 = GAC_1[0]
            C1 = GAC_1[1]
            A1 = GAC_1[2]
            GAC_2 = chunk[1]
            G2 = GAC_2[0]
            C2 = GAC_2[1]
            A2 = GAC_2[2]
            card_fields += [None, None, G1, C1, A1, G2, C2, A2, None, None]

    # Cambiamos los dos primeros None de la lista por MPC y SID
    card_fields[0] = 'MPC'
    card_fields[1] = sid

    return card_fields


@dataclass
class CQUAD4:
    EID: int
    PID: int
    G1: int
    G2: int
    G3: int
    G4: int
    THETA: float = 0.0
    MCID: int | None = 0
    ZOFFS: float | None = None
    TFLAG: int | None = None
    T1: float = 0.0
    T2: float = 0.0
    T3: float = 0.0
    T4: float = 0.0
    bdf_file: str | None = None
    bdf_nline: int | None = None

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

        __fields_order = [
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
    EID: int
    PID: int
    G1: int
    G2: int
    G3: int
    THETA: float = 0.0
    MCID: int | None = None
    ZOFFS: float | None = None
    TFLAG: int | None = None
    T1: float = 0.0
    T2: float = 0.0
    T3: float = 0.0
    T4: float = 0.0
    bdf_file: str | None = None
    bdf_nline: int | None = None

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
