import warnings
from glob import glob

from astrocats.catalog.utils import is_number

with open('astrocats/faststars/input/rep-folders.txt', 'r') as f:
    repofolders = f.read().splitlines()

outdir = 'astrocats/faststars/output/'


def repo_file_list(normal=True, bones=True):
    files = []
    for rep in repofolders:
        if 'boneyard' not in rep and not normal:
            continue
        if not bones and 'boneyard' in rep:
            continue
        files += glob(outdir + rep + "/*.json") + \
            glob(outdir + rep + "/*.json.gz")
    return files


def get_rep_folder(entry):
    return repofolders[0]
