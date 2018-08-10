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
import csv
import glob
import os
import shutil
import textwrap
import sys



def get_supplementary_info(filename):
    supplementary = {}
    with open(filename) as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            supplementary[row['motif']] = row
    return supplementary


def get_original_metadata(filename):
    with open(filename, 'r') as seed:
        seed_data = {}
        for line in seed:
            if '#=GF ' in line:
                (key, value) = line.replace('#=GF ', '').strip().split(' ', 1)
                seed_data[key] = value.strip()
        return seed_data


def make_desc_file(rna_id, supplementary, original_metadata):
    template = """ID   {id}
DE   {description}
AU   Weinberg Z; 0000-0002-6681-3624
SE   Weinberg Z
SS   Published; PMID:28977401;
TP   {rna_type}
DR   {so_term}
RN   [1]
RM   28977401
RT   Detection of 224 candidate structured RNAs by comparative analysis of
RT   specific subsets of intergenic regions.
RA   Weinberg Z, Lunse CE, Corbino KA, Ames TD, Nelson JW, Roth A, Perkins
RA   KR, Sherlock ME, Breaker RR;
RL   Nucleic Acids Res. 2017;45:10811-10823.
GA   27.00
TC   27.00
NC   27.00
BM   cmbuild -F CM SEED
WK   {wiki}
{cc_lines}"""
    if supplementary[rna_id]['switch?'].lower()[0] == 'y':
        rna_type = 'Cis-reg; riboswitch;'
        so_term = 'SO; 0000035; riboswitch;'
    elif supplementary[rna_id]['cis-reg?'].lower()[0] == 'y':
        rna_type = 'Cis-reg;'
        so_term = 'SO; 0005836; regulatory_region;'
    else:
        rna_type = 'Gene; sRNA;'
        so_term = 'SO; 0001263; ncRNA_gene;'

    cc_lines = ''
    cc_text = textwrap.wrap(supplementary[rna_id]['Taxa'] + '. ' + supplementary[rna_id]['Analysis method'], 75)
    for cc_line in cc_text:
        cc_lines += 'CC   %s\n' % cc_line

    return template.format(id=rna_id, description='%s RNA' % rna_id, rna_type=rna_type,
        so_term=so_term, wiki=original_metadata['WK'].replace('http://en.wikipedia.org/wiki/', ''), cc_lines=cc_lines)


def main(args):
    """
    """
    supplementary = get_supplementary_info(os.path.join(args.inputfile, 'zasha-208.csv'))

    for rna in glob.glob(os.path.join(args.inputfile, '*.sto')):
        motif_name = os.path.basename(rna).replace('.sto', '')
        motif_dir = os.path.join(args.destination, motif_name)
        original_metadata = get_original_metadata(rna)

        if not os.path.exists(motif_dir):
            os.mkdir(motif_dir)

        with open(os.path.join(motif_dir, 'DESC'), 'w') as desc:
            desc.write(make_desc_file(motif_name, supplementary, original_metadata))

        shutil.copy(rna, os.path.join(motif_dir, 'SEED'))
        os.chdir(motif_dir)
        cmd = ('module load mpi/openmpi-x86_64 && '
               'bsub -o {0}/lsf_output.txt -e {0}/lsf_error.txt -g /emerge '
                     '"cd {0} && '
                     'mv {0}/SEED {0}/SEED-backup && '
                     'esl-reformat --mingap stockholm SEED-backup > SEED && '
                     'rfsearch.pl -t 30 -cnompi -relax && '
                     'rfmake.pl -t 50 -a -forcethr && '
                     'mkdir rscape-seed && R-scape --outdir rscape-seed --cyk align && '
                     'mkdir rscape-align && R-scape --outdir rscape-align --cyk align && '
                     'cd .. && '
                     'rqc-all.pl {1}"').format(motif_dir, motif_name)
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

    main(args)
