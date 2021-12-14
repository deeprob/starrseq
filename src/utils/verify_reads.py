import os
import subprocess

CURRENT_DIR_PATH = os.path.dirname(os.path.abspath(__file__))

def create_coverage_bed(roi_bed, lib_pre, reps, suffix):

    if reps:
        reps = reps.split()

        for rep in reps:
            fname = "_".join([lib_pre, rep, suffix])
            assert os.path.exists(fname)
            fout = os.path.splitext(fname)[0] + ".coverage.bed"
            exit_code = subprocess.call(["bash", f"{CURRENT_DIR_PATH}/get_coverage.sh", f"{roi_bed}", f"{fname}", f"{fout}", "false"])
            assert exit_code==0
            fout = os.path.splitext(fname)[0] + ".depth.bed"
            exit_code = subprocess.call(["bash", f"{CURRENT_DIR_PATH}/get_coverage.sh", f"{roi_bed}", f"{fname}", f"{fout}", "true"])
            assert exit_code==0
    else:
        fname = lib_pre + ".bam"
        assert os.path.exists(fname)
        fout = os.path.splitext(fname)[0] + ".coverage.bed"
        exit_code = subprocess.call(["bash", f"{CURRENT_DIR_PATH}/get_coverage.sh", f"{roi_bed}", f"{fname}", f"{fout}", "false"])
        assert exit_code==0
        fout = os.path.splitext(fname)[0] + ".depth.bed"
        exit_code = subprocess.call(["bash", f"{CURRENT_DIR_PATH}/get_coverage.sh", f"{roi_bed}", f"{fname}", f"{fout}", "true"])
        assert exit_code==0

    return

def create_coverage_beds(roi_bed, ilp, ilr, clp, clr, klp, klr, iflag=True, cflag=True):
    if iflag:
        # input coverage
        create_coverage_bed(roi_bed, ilp, ilr, "filtered.bam")
    if cflag:
        # control coverage
        create_coverage_bed(roi_bed, clp, clr, "filtered.bam")
    # ko coverage
    create_coverage_bed(roi_bed, klp, klr, "filtered.bam")
    return 
