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

source_path = "sftp_mbb:/share/nas2-amap/gvieilledent/gfc2020_70"
dest_path = "/home/forestatrisk-tropics/gfc2020_70"
data_raw = "**/data_raw/**"
data = "**/data/**"
cmd = ["rclone", "sync", source_path, dest_path, "--exclude", data_raw, "--exclude", data, "--progress"]
cmd = " ".join(cmd)
subprocess.call(cmd, shell=True)

# EOF
