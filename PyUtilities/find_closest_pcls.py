import h5py as py

def get_dist2(x1, y1, z1, x2, y2, z2):
    return (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2)

def get_closest_pcls(pcl_dset, pos, closest_pcl_num):
    pcl_num = len(pcl_dset)
    if pcl_num == 0:
        return [], [];
    pcl_ids = pcl_dset[:]['id']
    pcl_xs = pcl_dset[:]['x']
    pcl_ys = pcl_dset[:]['y']
    pcl_zs = pcl_dset[:]['z']
    closest_pcls = [ pcl_ids[0] ]
    closest_pcl_dist2 = [ get_dist2(pcl_xs[0], pcl_ys[0], pcl_zs[0], pos[0], pos[1], pos[2]) ]
    if closest_pcl_num > pcl_num:
        closest_pcl_num = pcl_num
    for p_id in range(1, closest_pcl_num):
        p_dist2 = get_dist2(pcl_xs[p_id], pcl_ys[p_id], pcl_zs[p_id], pos[0], pos[1], pos[2])
        if p_dist2 >= closest_pcl_dist2[-1]:
            closest_pcls.append(pcl_ids[p_id])
            closest_pcl_dist2.append(p_dist2)
            continue
        for cp_id in range(p_id):
            if p_dist2 < closest_pcl_dist2[cp_id]:
                closest_pcls.insert(cp_id, pcl_ids[p_id])
                closest_pcl_dist2.insert(cp_id, p_dist2)
                break
    if closest_pcl_num == pcl_num:
        return closest_pcls, closest_pcl_dist2
    for p_id in range(closest_pcl_num, pcl_num):
        if p_id % 100000 == (100000-1):
            print("%d pcls processed." % (p_id+1))
        p_dist2 = get_dist2(pcl_xs[p_id], pcl_ys[p_id], pcl_zs[p_id], pos[0], pos[1], pos[2])
        if p_dist2 >= closest_pcl_dist2[-1]:
            continue
        for cp_id in range(p_id):
            if p_dist2 < closest_pcl_dist2[cp_id]:
                closest_pcls.insert(cp_id, pcl_ids[p_id])
                closest_pcls.pop()
                closest_pcl_dist2.insert(cp_id, p_dist2)
                closest_pcl_dist2.pop()
                break
    return closest_pcls, closest_pcl_dist2

def get_pcls_in_box(pcl_dset, box_xl, box_xu, box_yl, box_yu, box_zl, box_zu):
    pcl_num = len(pcl_dset)
    if pcl_num == 0:
        return [];
    pcl_ids = pcl_dset[:]['id']
    pcl_xs = pcl_dset[:]['x']
    pcl_ys = pcl_dset[:]['y']
    pcl_zs = pcl_dset[:]['z']
    closest_pcls = []
    for p_id in range(pcl_num):
        if pcl_xs[p_id] >= box_xl and pcl_xs[p_id] <= box_xu and \
           pcl_ys[p_id] >= box_yl and pcl_ys[p_id] <= box_yu and \
           pcl_zs[p_id] >= box_zl and pcl_zs[p_id] <= box_zu:
            closest_pcls.append(pcl_ids[p_id])
    return closest_pcls

if __name__ == "__main__":
    tip_pos = (0.1, 0.1, 0.1)
    closest_pcl_num = 10
    
    # get pcl close to pore_pos
    hdf5_file = py.File("../Build/TestsParallel/t3d_chm_mt_1d_consolidation.h5", "r")
    pcl_dset = hdf5_file['ModelData']['ParticleData']['field']
    
    #
    closest_pcls, closest_pcl_dist2 = get_closest_pcls(pcl_dset, tip_pos, closest_pcl_num)
    
    if len(closest_pcls) == 0:
        raise UserWarning("Cannot get pcl %d!\n" % th_id)
    
    print(closest_pcls)
    for i in range(len(closest_pcls)):
        pcl_id = closest_pcls[i]
        p_data = pcl_dset[pcl_id]
        print("%d, %f, %f, %f, %f" % (pcl_id, p_data['x'], p_data['y'], p_data['z'], closest_pcl_dist2[i]))

    #
    closest_pcls2 = get_pcls_in_box(pcl_dset, 0.05, 0.15, 0.05, 0.15, 0.25, 0.35)
    
    print(closest_pcls2)