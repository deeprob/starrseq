import os
import prepare_meta.utils as ut

def create_meta(meta_file, force=False):
    if not force:
        if not os.path.exists(meta_file):
            ut.create_metadata(meta_file)
        else:
            print("File already present")
    else:
        ut.create_metadata(meta_file)
    return    
