#!/usr/bin/env python2

"""\
Pull the latest changes to the rosetta checkout being benchmarked.
"""

import os, subprocess
from libraries import settings

settings.load()
os.chdir(settings.rosetta)
subprocess.check_call(('git', 'pull'))
