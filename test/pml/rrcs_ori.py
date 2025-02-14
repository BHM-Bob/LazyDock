'''
Date: 2025-02-13 21:47:48
LastEditors: BHM-Bob 2262029386@qq.com
LastEditTime: 2025-02-14 15:09:31
Description: copy from https://github.com/elifesciences-publications/RRCS/blob/master/RRCS.py
'''


def calc_contact(pdbbase):
    fi = open(pdbbase, 'r')
    all_lines = fi.readlines()
    fi.close()
    atom_lines = [l for l in all_lines if l[0:6] == 'ATOM  ']
    dict_coord = {} # dict to store coordinates. dict_coord[res][atom] = (x,y,z,occ)
    atomnum_2_name = {} # map atom number to atom name, in order to find N, CA, C, O
    contact_score = {} # dict to store final results. contact_score[ires][jres] = contact_score.
    for line in atom_lines:
        # retrive info from each atom line
        atom_num = int(line[6:11].strip())
        atom_name = line[12:16].replace(' ', '_')
        res_name = line[17:20]
        res_num = int(line[22:26].strip())
        chain_id = line[21:22]
        x = float(line[30:38].strip())
        y = float(line[38:46].strip())
        z = float(line[46:54].strip())
        occ = float(line[54:60].strip())
        res = chain_id + ':' + str(res_num) + '_' + res_name 
        atomnum_2_name[atom_num] = atom_name
        if res_num <= 0:
            continue
        if res not in dict_coord:
            dict_coord[res] = {}
        dict_coord[res][atom_num] = (x, y, z, occ)
    #list shortening  
    num_of_aa = len(dict_coord)
    res_list = dict_coord.keys()
    ires_contact = {}
    for ires in res_list:
        ires_contact[ires] = []
        ires_num = int(ires.split(':')[1].split('_')[0].strip())
        for jres in res_list:
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

calc_contact('data_tmp/pdb/RECEPTOR.pdb')