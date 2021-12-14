import os
import subprocess

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

def get_intersects(file_a, file_b, file_out):
    exit_code = subprocess.call(["bash", f"{CURRENT_DIR_PATH}/get_intersects.sh", f"{file_a}", f"{file_b}", f"{file_out}", "true"])
    assert exit_code==0
    return

def get_non_intersects(file_a, file_b, file_out):
    exit_code = subprocess.call(["bash", f"{CURRENT_DIR_PATH}/get_intersects.sh", f"{file_a}", f"{file_b}", f"{file_out}", "false"])
    assert exit_code==0
    return