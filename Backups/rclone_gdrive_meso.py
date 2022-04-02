#!/usr/bin/env python
# -*- coding: utf-8 -*-

# ==============================================================================
# author          :Ghislain Vieilledent
# email           :ghislain.vieilledent@cirad.fr, ghislainv@gmail.com
# web             :https://ghislainv.github.io
# python_version  :2.7
# license         :GPLv3
# ==============================================================================

# PRIOR TO EXECUTING THE FOLLOWING SCRIPT, YOU MUST HAVE CREDENTIALS FOR
# 2. rclone with Google Drive: https://rclone.org/drive/

# For documentation on filtering, see here: https://rclone.org/filtering/

import subprocess

source_path = "gdrive_gv:Work/forestatrisk-tropics/jrc2020"
dest_path = "/storage/replicated/cirad/projects/AMAP/vieilledentg/jrc2020"

data_raw = "**/data_raw/**"
data = "**/data/**"
# cmd = ["rclone", "sync", source_path, dest_path, "--exclude", data_raw,
#        "--exclude", data, "--progress"]
# cmd = ["rclone", "sync", source_path, dest_path, "--exclude", data_raw,
#        "--progress"]
cmd = ["rclone", "sync", source_path, dest_path, "--progress"]
cmd = " ".join(cmd)
subprocess.call(cmd, shell=True)

# End
