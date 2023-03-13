import os, re
import numpy as np
from PIL import Image
import pandas as pd
import argparse


def get_files_path_recursively(path, *args, suffix=''):
    """Retrieve specific files path recursively from directory.

    Retrieve the path of all files with one of the given extension names,
    in the given directory and all its subdirectories, recursively.
    The extension names should be given as a list of strings. The search for
    extension names is case sensitive.

    Args:
        path (str): root directory from which to search for files recursively
        *args: list of file extensions to be considered

    Returns:
        list(str): list of paths of every files in the directory and all its
        subdirectories
    """
    exts = list(args)
    for ext_i, ext in enumerate(exts):
        exts[ext_i] = ext[1:] if ext[0] == '.' else ext
    ext_list = "|".join(exts)
    result = [os.path.join(dp, f)
              for dp, dn, filenames in os.walk(path)
              for f in filenames
              if re.search(rf"^.*({suffix})\.({ext_list})$", f)]
    return result

def standardize(root_path:str='sample_data/SatelliteImages/',
                ext:str=['jpeg', 'jpg'],
                output:str='root_path'):
    """Perform standardization over images.
    
    Returns and stores the mean and standard deviation of an image
    dataset organized inside a root directory for computation
    purposes like deep learning.

    Args:
        root_path (str): root dir. containing the images.
                         Defaults to './sample_data/SatelliteImages/'.
        ext (str, optional): the images extensions to consider.
                             Defaults to 'jpeg'.
        output (str, optional): output path where to save the csv containing
                                the mean and std of the dataset.
                                If None: doesn't output anything.
                                Defaults to root_path.

    Returns:
        _type_: _description_
    """
    fps = get_files_path_recursively(root_path, *ext)
    imgs = []
    stats = {'mean':[], 'std':[]}
    for fp in fps:
        img = np.array(Image.open(fp, mode='r'))
        if len(img.shape) == 3:
            img = np.transpose(img, (2,0,1))
        elif len(img.shape) == 2:
            img = np.expand_dims(img, axis=0)
        imgs.append(img)
    imgs = np.concatenate(imgs, axis=0)
    stats['mean'].append(np.nanmean(imgs))
    stats['std'].append(np.nanstd(imgs))
    if output:
        output = os.path.join(root_path, 'standardize_stats.csv') if output=='root_path' else output
        df = pd.DataFrame(stats)
        df.to_csv(output, index=False, sep=';')
    return stats['mean'][0], stats['std'][0]


if __name__ == '__main__':
    PARSER = argparse.ArgumentParser()
    PARSER.add_argument('--root_path',
                        nargs=1,
                        type=str,
                        default=['sample_data/SatelliteImages/'],
                        help='Rooth path.')
    PARSER.add_argument('--ext',
                        nargs=1,
                        type=str,
                        default=['jpg', 'jpeg'],
                        help='File extension.')
    PARSER.add_argument('--out',
                        nargs=1,
                        type=str,
                        default=['sample_data/SatelliteImages/jpeg_patches_sample_stats.csv'],
                        help='Output path.')
    ARGS = PARSER.parse_args()
    path = ARGS.root_path[0]
    ext = ARGS.ext
    out = ARGS.out[0]
    standardize(path, ext=ext, output=out)