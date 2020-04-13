import numpy as np
import spatial_split_utils

np.random.seed(1949)

lon_list = (np.random.rand(1000) - 0.5) * 180 # replace this
lat_list = (np.random.rand(1000) - 0.5) * 90 # replace this

east_west_bin_width_km = 5.0
north_south_bin_width_km = 5.0

train_frac = 0.8

origin = (0.0,-30.0)

'''
for each point, get the corresponding block_id:
'''
block_id_assignments = spatial_split_utils.assign_block_ids(lon_list,lat_list,east_west_bin_width_km,north_south_bin_width_km,origin)

(vals,counts) = np.unique(block_id_assignments,return_counts=True)
print(len(vals))

'''
calculate train/test split over the block:
'''
block_ids = np.unique(block_id_assignments)

num_train = int(np.floor(train_frac * len(block_ids)))
num_test = len(block_ids) - num_train

idx_rand = np.random.permutation(len(block_ids))
idx_train = idx_rand[:num_train]
idx_test = idx_rand[-num_test:]
assert len(np.intersect1d(idx_train,idx_test)) == 0

train_block_ids = block_ids[idx_train]
test_block_ids = block_ids[idx_test]
assert ((len(train_block_ids) + len(test_block_ids)) == len(block_ids))

'''
save splits
'''
split_assignment = []
for i in range(len(lon_list)):
    if block_id_assignments[i] in train_block_ids:
        split_assignment.append(0)
    else:
        split_assignment.append(1)
split_assignment = np.array(split_assignment,dtype=int)
np.save('split_assignment.npy',split_assignment)

train_frac_ex = 1 - np.mean(split_assignment)
print('%0.1f / %0.1f block split'%(train_frac * 100, (1-train_frac) * 100))
print('%0.1f / %0.1f example split'%(train_frac_ex * 100, (1-train_frac_ex) * 100))
