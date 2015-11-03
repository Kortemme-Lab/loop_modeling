#!/usr/bin/env python2
# encoding: utf-8

#$ -S /usr/bin/python
#$ -cwd
#$ -r yes
#$ -j y
#$ -l h_rt=24:00:00
#$ -t 1-27
#$ -l arch=linux-x64
#$ -l mem_free=2G
#$ -l netapp=1G

# The MIT License (MIT)
#
# Copyright (c) 2015 Amelie Stein, Roland Pache, Kale Kundert
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

import socket
import sys
import shlex

print "Python:", sys.version
print "Host:", socket.gethostname()

import datetime
import os
import subprocess
import time



def tee(*popenargs, **kwargs):
    # todo: import this from the tools repository
    import subprocess, select, sys

    process = subprocess.Popen(
            stdout=subprocess.PIPE, stderr=subprocess.PIPE,
            *popenargs, **kwargs)

    stdout, stderr = '', ''

    print "Process ID:", process.pid

    def read_stream(input_callback, output_stream):   # (no fold)
        read = input_callback()
        output_stream.write(read)
        output_stream.flush()
        return read

    while process.poll() is None:
        watch = process.stdout.fileno(), process.stderr.fileno()
        ready = select.select(watch, [], [])[0]

        for fd in ready:
            if fd == process.stdout.fileno():
                stdout += read_stream(process.stdout.readline, sys.stdout)
            if fd == process.stderr.fileno():
                stderr += read_stream(process.stderr.readline, sys.stderr)

    stdout += read_stream(process.stdout.read, sys.stdout)
    stderr += read_stream(process.stderr.read, sys.stderr)

    return stdout, stderr


def write_file(filepath, contents, ftype = 'w'):
    output_handle = open(filepath, ftype)
    output_handle.write(contents)
    output_handle.close()



time_start = time.time()

sge_task_id = 1
if os.environ.has_key("SGE_TASK_ID"):
    sge_task_id = os.environ["SGE_TASK_ID"]

id = int(sge_task_id)

home_dir = os.path.expanduser("~")

rosetta_dir = '/path/to/rosetta/main/source/bin/'
bin_ext = "linuxgccrelease"

bin_path = os.path.join(rosetta_dir, "fixbb.{0}".format(bin_ext))
bin_path = os.path.abspath(os.path.expanduser(bin_path))

if len(sys.argv) < 3:
    print "Error -- you need to provide the list of input PDB files and the outfile keyword (e.g. my_min_packed). You can specify additional flags for fixbb after these options"
    exit(-1)

pdb_lst_file = sys.argv[1]
pdb_lst = []
read_pdbs = open(pdb_lst_file)
for l in read_pdbs:
    pdb_lst.append(l.strip('\n'))

read_pdbs.close()
pdb_id = pdb_lst[(id - 1) % len(pdb_lst)]
pdb_core = pdb_id.split("_")[0] ## 1cb0_min.pdb -> 1cb0

fa_params = pdb_core + ".fa.params"
cen_params = pdb_core + ".cen.params"

outfile_key = sys.argv[2]

addtl_cmds = []
for i in range(3, len(sys.argv)):
    addtl_cmds.append(sys.argv[i])

if os.path.isfile(fa_params):
    addtl_cmds.extend(["-extra_res_fa", fa_params])

if os.path.isfile(cen_params):
    addtl_cmds.extend(["-extra_res_cen", cen_params])

outfile_name = pdb_id.split('/')[-1].split(".")[0]+"_" + outfile_key + "_"+str(id)

local_os = sys.platform
rosetta_cmd_sep = "::"

if (local_os == "darwin"):
    rosetta_cmd_sep = ":"

#lig_path = database_path+"chemical/residue_type_sets/fa_standard/residue_types/metal_ions/MG.params" ## only required for non-built-in types

args = [
    bin_path,
    "-in:file:fullatom",
    "-in:file:s %s " % pdb_id,
    "-out:suffix %s" % outfile_key,
    "-min_pack",
    "-packing:repack_only",
    "-ignore_unrecognized_res", ## to ignore ligands, waters etc. that may be in native PDB structures
    "-overwrite",
    "-nstruct 1",
    "-ex1 -ex2 -extrachi_cutoff 0",
]

args.extend(addtl_cmds)

args = shlex.split(' '.join(args))

stdout, stderr = tee(args)
write_file("%s.log"%(outfile_name), stdout)
write_file("%s.err"%(outfile_name), stderr)

time_end = time.time()

print "Seconds:", time_end - time_start
print "Time:", datetime.timedelta(seconds = time_end - time_start)
print "Summary:", socket.gethostname(), time_end - time_start, datetime.timedelta(seconds = time_end - time_start)
