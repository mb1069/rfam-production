"""
Copyright [2009-2018] EMBL-European Bioinformatics Institute
Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
     http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

"""
Usage:
1. Make sure that rfsearch, rfmake and other commands are in PATH.
2. python precompute_zasha.py -i /path/to/emerge_file -d /path/to/output

Analysing results:
# check that all LSF jobs completed successfully
ls | wc -l
find . -type f -name lsf_output.txt | xargs grep 'Success' | wc -l
# find jobs that didn't finish successfully
find . -type f -name lsf_output.txt | xargs grep -L 'Success'
# find all directories without the outlist file
find . -maxdepth 1 -mindepth 1 -type d | while read dir; do [[ ! -f $dir/outlist ]] && echo "$dir has no outlist"; done
# count the number of lines above the best reversed hit
find . -type f -name outlist -exec sh -c 'sed -n "0,/REVERSED/p" {} | wc -l' \; -print
# get overlapping Rfam families
find . -type f -name overlap | xargs grep -o -P "RF\d{5}" | sort | uniq
"""


import argparse
import glob
import os
import shutil
import sys


def run(args):
    """
    """
    for rna in glob.glob(args.inputfile + '*'):
        motif_name = os.path.basename(rna).replace('.sto', '')
        motif_dir = os.path.join(args.destination, motif_name)
        if not os.path.exists(motif_dir):
            os.mkdir(motif_dir)
        shutil.copy(rna, os.path.join(motif_dir, 'SEED'))
        os.chdir(motif_dir)
        cmd = ('module load mpi/openmpi-x86_64 && '
               'bsub -o {0}/lsf_output.txt -e {0}/lsf_error.txt -g /emerge '
                     '"cd {0} && '
                     'rm -f DESC && '
                     'rfsearch.pl -nodesc -t 30 -cnompi -relax && '
                     'rfmake.pl -t 50 -a -forcethr && '
                     'mkdir rscape && R-scape --outdir rscape --cyk align && '
                     'cd .. && '
                     'rqc-overlap.pl {1}"').format(motif_dir, motif_name)
        print cmd
        if not args.test:
            os.system(cmd)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--destination', default=os.getcwd(), help='Specify folder where the output will be created')
    parser.add_argument('-i', '--inputfile', help='Specify input folder with Stockholm files')
    parser.add_argument('-t', '--test', action='store_true', help='Test mode: print commands and exit')
    parser.set_defaults(test=False)
    args = parser.parse_args()

    if not args.inputfile:
        print 'Please specify input location'
        sys.exit()

    run(args)
