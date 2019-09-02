"""
Copyright [2009-2019] EMBL-European Bioinformatics Institute
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

import os
import sys

# -------------------------------------------------------------------------


def merge_genome_files(source_dir, file_type, dest_dir, filename=None):
    """
    Merges all files located in the source directory of a specified type

    source_dir:
    file_type:
    dest_dir:
    """

    if dest_dir is None:
        dest_dir = source_dir

    # some necessary file type conversions
    if file_type.lower() == 'tblout':
        file_type = 'tbl'

    elif file_type.lower() == 'fasta':
        file_type = 'fa'

    input_files = [x for x in os.listdir(source_dir)
                   if x.endswith(file_type)]

    if filename is None:
        filename = "all_merged." + file_type
    else:
        filename = filename + "." + file_type

    fp_out = open(os.path.join(dest_dir, filename), 'w')
    # loop over all files and merge them in a single unit
    for file in input_files:
        file_loc = os.path.join(source_dir, file)
        fp_in = open(file_loc, 'r')

        for line in fp_in:
            fp_out.write(line)

        fp_in.close()

    fp_out.close()

# -------------------------------------------------------------------------


def merge_project_files(project_dir, upid_list, file_type):
    """
    The purpose of this function is to merge genome files (rfamseq, genseq, tblout)
    to release ready files for the database import

    project_dir: The path to the project directory. Project directory should be
    in the same structure as the one generated by genome_download pipeline
    upid_list: A list of upids to include in the merge
    file_type: The type of the files to merge (rfamseq, genseq, tblout)

    return:
    """

    fp = open(upid_list, 'r')

    upids = [x.strip() for x in fp]

    fp.close()

    for upid in upids:
        subdir_loc = os.path.join(project_dir, upid[-3:])
        updir = os.path.join(subdir_loc, upid)

        source_dir = ''
        if file_type.lower() == "tblout" or file_type.lower() == "tbl":
            source_dir = os.path.join(updir, "search_output")
        elif file_type.lower() == "fasta" or file_type.lower() == "fa":
            source_dir = os.path.join(updir, "sequence_chunks")

        merge_genome_files(source_dir, file_type, updir, filename=upid)


# -------------------------------------------------------------------------


def gather_project_files(project_dir, upid_list, file_type):

    fp = open(upid_list, 'r')

    upids = [x.strip() for x in fp]

    fp.close()

     # some necessary file type conversions
    if file_type.lower() == 'tblout':
        file_type = 'tbl'

    elif file_type.lower() == 'fasta':
        file_type = 'fa'

    fp_out = open(os.path.join(project_dir, "all_merged." + file_type), 'w')

    for upid in upids:
        subdir_loc = os.path.join(project_dir, upid[-3:])
        updir = os.path.join(subdir_loc, upid)

        upfile = os.path.join(updir, upid + '.' + file_type)

        fp_in = open(upfile, 'r')
        for line in fp_in:
            fp_out.write(line)

        fp_in.close()

    fp_out.close()

# -------------------------------------------------------------------------


def merge_batch_search_tbls(result_dir, filename = None):
    """
    Merges all infernal tbl files produced by a genome_scanner batch search

    result_dir: The path to the result directory

    return: void
    """

    fp_out = open(os.path.join(result_dir, "full_region.tbl"), 'w')

    subdirs = [x for x in os.listdir(result_dir) if
                os.path.isdir(os.path.join(result_dir, x))]

    # list 24 support subdirs
    for subdir in subdirs:
        subdir_loc = os.path.join(result_dir, subdir)
        umgs_dirs = os.listdir(subdir_loc)

        # list umgs dirs
        for umgs_dir in umgs_dirs:
            umgs_dir_loc = os.path.join(subdir_loc, umgs_dir)

            # list all tbl files
            tbl_files = [x for x in os.listdir(umgs_dir_loc) if x.endswith(".tbl")]

            for tbl_file in tbl_files:
                tbl_file_loc = os.path.join(umgs_dir_loc, tbl_file)
                fp_in = open(tbl_file_loc, 'r')
                tbl_lines = fp_in.readlines()
                fp_out.writelines(tbl_lines)
                fp_in.close()

    fp_out.close()

# -------------------------------------------------------------------------

if __name__ == '__main__':


    if "--genome" in sys.argv:
        project_dir = sys.argv[1]
        upid_list = sys.argv[2]
        file_type = sys.argv[3]

        merge_project_files(project_dir, upid_list, file_type)

    elif "--project" in sys.argv:
        project_dir = sys.argv[1]
        upid_list = sys.argv[2]
        file_type = sys.argv[3]

        gather_project_files(project_dir, upid_list, file_type)

    elif "--batch" in sys.argv:
        result_dir = sys.argv[1]

        merge_batch_search_tbls(result_dir)

    else:
        print "Wrong option selected"