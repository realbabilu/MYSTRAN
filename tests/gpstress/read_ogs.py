"""Independent 32-bit little-endian OGS1 volume decoder and CSV exporter.

Checks every Fortran record boundary. Enumerates volume IDs explicitly; does
not collapse coincident grid coordinates or assume a generic OP2 stress slot.
"""
import argparse
import csv
import struct
import re
from pathlib import Path


def records(path):
    data = Path(path).read_bytes()
    out = []
    pos = 0
    while pos < len(data):
        n, = struct.unpack_from('<i', data, pos)
        if n < 0 or pos+n+8 > len(data):
            raise ValueError(f'invalid record at byte {pos}: {n}')
        end, = struct.unpack_from('<i', data, pos+n+4)
        if end != n:
            raise ValueError(f'record length mismatch at {pos}: {n}/{end}')
        out.append(data[pos+4:pos+4+n])
        pos += n+8
    return out


def read_volumes(path, table_code=27):
    rec = records(path)
    tables = []
    table_name=b''
    for j, item in enumerate(rec):
        if len(item)==8 and re.fullmatch(b'[A-Z0-9 ]{8}',item):
            table_name=item
        # Simcenter may put several volumes in one OGS1 table, each with its
        # own table-3 metadata. MYSTRAN currently uses a separate table per volume.
        if table_name != b'OGS1    ' or len(item)!=584:
            continue
        head = struct.unpack('<50i',rec[j][:200])
        if head[1] % 1000 != table_code:
            continue  # surface/principal variants are not direct volumes
        width={27:9,28:15}[table_code]
        assert head[0]//10 == 1 and head[8:10] == (1,width), head[:14]
        k = j+1
        while len(rec[k]) == 4:
            k += 1
        nwords, = struct.unpack('<i',rec[k-1])
        if len(rec[k]) != nwords*4 or nwords % width:
            raise ValueError('OGS volume payload length disagrees with word count')
        rows = list(struct.iter_unpack('<i'+str(width-1)+'f',rec[k]))
        tables.append(dict(subcase=head[3],volume=head[2],cid=head[7],header=head[:23],
                           grids=[r[0]//10 for r in rows],
                           values=[r[1:] for r in rows]))
    return tables


def pynastran_volumes(path):
    import numpy as np
    # Compatibility with pyNastran 1.4 on NumPy 2.4; no installed files changed.
    if not hasattr(np,'in1d'):
        np.in1d = np.isin
    if not hasattr(np,'chararray'):
        np.chararray = np.char.chararray
    from pyNastran.op2.op2 import OP2
    class ReferenceOP2(OP2):
        def _read_pvto_4(self,data,ndata):
            # pyNastran 1.4 rejects the Simcenter SOLIDTYP string PARAM in PVT0.
            # Skip parameter metadata only; result tables are still decoded.
            return ndata
    model = ReferenceOP2(debug=False,mode='msc')
    # Existing pyramid OES type 255 is unsupported in this pyNastran version.
    model.set_additional_result_tables_to_read({b'OES1X1':False,b'OSTR1X':False})
    model.read_op2(str(path))
    return model


def export(path, csv_path):
    tables=read_volumes(path)
    with Path(csv_path).open('w',newline='') as stream:
        writer=csv.writer(stream)
        writer.writerow(['subcase','volume','cid','grid','sxx','syy','szz','sxy','syz','szx','pressure','von_mises'])
        for table in tables:
            for grid,values in zip(table['grids'],table['values']):
                writer.writerow([table['subcase'],table['volume'],table['cid'],grid,*values])
    return tables


if __name__ == '__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('op2')
    parser.add_argument('--csv',default='ogs_volumes.csv')
    parser.add_argument('--pynastran',action='store_true')
    args=parser.parse_args()
    for table in export(args.op2,args.csv):
        print(f"subcase={table['subcase']} volume={table['volume']} cid={table['cid']} grids={len(table['grids'])}")
    if args.pynastran:
        model=pynastran_volumes(args.op2)
        for key,result in model.grid_point_stresses_volume_direct.items():
            print(key, 'volume=',result.ogs_id,'shape=',result.data.shape)
