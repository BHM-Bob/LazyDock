'''
Date: 2024-09-30 19:28:57
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2024-10-08 11:03:53
Description: RRCS calculation in PyMOL, RRCS is from article "Common activation mechanism of class A GPCRs": https://github.com/elifesciences-publications/RRCS/blob/master/RRCS.py
'''
from pymol import cmd


def calcu_RRCS(model: str):
    """
    calcu RRCS score for each residue in the model.
    
    Parameters:
        - model: str, pymol name of the model to calculate RRCS score, must be loaded in PyMOL.
        
    Returns:
    contact_score: dict, RRCS score for each residue pair. contact_score[ires][jres] = contact_score.
    """
    dict_coord = {} # dict to store coordinates. dict_coord[res][atom] = (x, y, z, occupancy)
    cmd.iterate_state(1, model, 'dict_coord.setdefault(f"{chain}:{resi}_{resn}", {}).setdefault(index, (x, y, z, q))', space={'dict_coord': dict_coord})
    atomnum_2_name = {} # map atom number to atom name, in order to find N, CA, C, O
    cmd.iterate(model, 'atomnum_2_name[index] = name', space={'atomnum_2_name': atomnum_2_name})
    contact_score = {} # dict to store final results. contact_score[ires][jres] = contact_score.
    #list shortening  
    res_list = dict_coord.keys()
    ires_contact = {}
    for ires in res_list:
        ires_contact[ires] = [] # any jres close to ires will be added to ires_contact[ires]
        ires_num = int(ires.split(':')[1].split('_')[0].strip())
        for jres in res_list:
            # if jres has any atom close to ires, add it to ires_contact[ires]
            jres_num = int(jres.split(':')[1].split('_')[0].strip())
            if jres_num <= ires_num:
                continue
            jres_flag = 0
            atom_in_ires = dict_coord[ires].keys()
            atom_in_jres = dict_coord[jres].keys()
            for iatom in atom_in_ires:
                (ix, iy, iz, iocc) = dict_coord[ires][iatom]
                for jatom in atom_in_jres:                  
                    (jx, jy, jz, jocc) = dict_coord[jres][jatom]
                    dx = abs(ix-jx)
                    dy = abs(iy-jy)
                    dz = abs(iz-jz)
                    if dx < 4.63 and dy < 4.63 and dz < 4.63:
                        jres_flag = 1
                        ires_contact[ires].append(jres)
                        break
                if jres_flag:
                    break
        # loop over the shortened list
        contact_score[ires] = {}        
        for kres in ires_contact[ires]:
            atom_in_kres = dict_coord[kres].keys()
            kres_num = int(kres.split(':')[1].split('_')[0].strip())
            contact_score[ires][kres] = 0
            total_score = 0
            if abs(ires_num - kres_num) < 5:
                for iatom in atom_in_ires:
                    iatom_name = atomnum_2_name[iatom]
                    if iatom_name in ['_N__', '_CA_', '_C__', '_O__']:
                        continue
                    (ix, iy, iz, iocc) = dict_coord[ires][iatom]
                    for katom in atom_in_kres:
                        katom_name = atomnum_2_name[katom]
                        if katom_name in ['_N__', '_CA_', '_C__', '_O__']:
                            continue
                        (kx, ky, kz, kocc) = dict_coord[kres][katom]
                        d2 = (ix-kx)**2 + (iy-ky)**2 + (iz-kz)**2
                        if d2 >= 21.4369:  # 4.63*4.63 = 21.4369
                            score = 0
                        elif d2 <= 10.4329:  # 3.23*3.23 = 10.4329
                            score = 1.0*iocc*kocc
                        else:
                            score = (1-(d2**0.5 - 3.23)/1.4)*iocc*kocc
                        total_score = total_score + score
            elif abs(ires_num - kres_num) > 4:
                for iatom in atom_in_ires:
                    (ix, iy, iz, iocc) = dict_coord[ires][iatom]
                    for katom in atom_in_kres:
                        (kx, ky, kz, kocc) = dict_coord[kres][katom]
                        d2 = (ix-kx)**2 + (iy-ky)**2 + (iz-kz)**2
                        if d2 >= 21.4369:  # 4.63*4.63 = 21.4369
                            score = 0
                        elif d2 <= 10.4329:  # 3.23*3.23 = 10.4329
                            score = 1.0*iocc*kocc
                        else:
                            score = (1-(d2**0.5 - 3.23)/1.4)*iocc*kocc
                        total_score = total_score + score
            contact_score[ires][kres] = total_score
    return contact_score

if __name__ == '__main__':
    cmd.reinitialize()
    cmd.load('data_tmp/pdb/RECEPTOR.pdb', 'receptor')
    from mbapy.base import TimeCosts
    @TimeCosts(3)
    def test_calcu(idx):
        calcu_RRCS('receptor')
    test_calcu()