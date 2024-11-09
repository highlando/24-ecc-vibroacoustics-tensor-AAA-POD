import numpy as np

'''data by Davide
6x6
    1e-4 1e-6 1e-8 1e-10
set 40   44   48   52
blk 21   23   25   37
hsv 22   24   26   38

12x12
    1e-4 1e-6 1e-8 1e-10
set 40   44   48   52
blk 20   24   25   28
hsv 21   25   26   29
'''


def blockaaadc(m, tnx=None, tny=None):
    return (tnx*tny+tnx*tnx+1)*(m+1)


def setvaaadc(m, tnx=None, tny=None):
    return (tnx*tny+2)*(m+1)


def hosvddc(tnx=None, tny=None, tnw=None, nw=None):
    return (tnx*tny + nw)*tnw


def mtotnw(m, tnx=None, tny=None, nw=None, blocaaa=True):
    if blocaaa:
        aaadc = blockaaadc(m, tnx=tnx, tny=tny)
    else:
        raise NotImplementedError()
    return int(np.ceil(aaadc / (tnx*tny + nw)))


tnxl = [6, 12]
nw = 391


for tnx in tnxl:
    tny = tnx
    if tnx == 12:
        bmlist = [20, 24, 25, 28]
        smlist = [40, 44, 48, 52]
    elif tnx == 6:
        bmlist = [21, 23, 25, 37]
        smlist = [40, 44, 48, 52]

    print(f'{tnx}x{tny}\n | bm | bds | sm | sds | nf | hds |')
    for kkk, m in enumerate(bmlist):
        sm = smlist[kkk]
        baaadc = blockaaadc(m, tnx=tnx, tny=tny)
        saaadc = setvaaadc(sm, tnx=tnx, tny=tny)
        tnw = mtotnw(m, tnx=tnx, tny=tny, nw=nw)
        hsvdc = hosvddc(tnx=tnx, tny=tny, tnw=tnw, nw=nw)
        print(f'| {m} | {baaadc} | {sm} | {saaadc} | {tnw} | {hsvdc} |')
