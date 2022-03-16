import os
import subprocess

root_dir  = "/data5/deepro/starrseq/main_lib/aligned_reads"

if __name__ == "__main__":
    rep_files = sorted([f.path for f in os.scandir(root_dir) if f.path.endswith("_filtered.bam")])
    command = ["python", "window_coverage.py", "/data5/deepro/starrseq/computational_pipeline/data/master.sorted.bed"]
    command += rep_files
    result = subprocess.run(command, capture_output=True, text=True)

    print(result.stdout)
    print(result.stderr)
