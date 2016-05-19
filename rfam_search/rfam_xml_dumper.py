'''
Created on 13 May 2016

@author: ikalvari

Description: This module exports Rfam data

TO DO: Need to add functions to parse and update an Xml4dbDumper file and 
       update entries (family, clan, motif)
'''

# ----------------------------------------------------------------------------

import os
import sys
import xml.etree.ElementTree as ET
import datetime
import string
from xml.dom import minidom
from config import rfam_search as rs_conf
from utils import RfamDB  # will load those from rfam-public

# ----------------------------------------------------------------------------

# this should be able to dump a single family, clan, motif or the entire db..


def xml4db_dumper(outdir):
    '''
        exports query results into EB-eye's XML4dbDUMP format
        outdir: Destination directory

        Maybe here provide the fields as a txt file and the dump in txt format
        and according to that dump the xml file
    '''

    # EB_eye_search fixed tags
    db_xml = ET.Element("database")
    ET.SubElement(db_xml, "name").text = rs_conf.DB_NAME
    ET.SubElement(db_xml, "description").text = rs_conf.DB_DESC
    ET.SubElement(db_xml, "release").text = rs_conf.DB_RELEASE
    ET.SubElement(db_xml, "release_date").text = rs_conf.DB_REL_DATE

    entries = ET.SubElement(db_xml, "entries")

    # call family xml builder to add a new family to the xml tree
    family_xml_builder(entries, rfam_acc='RF02560')

    tree = ET.ElementTree(db_xml)

    fp_out = open(os.path.join(outdir, "xml4dbsample.xml"), 'w')

    db_str = ET.tostring(db_xml, 'utf-8')
    db_str_reformated = minidom.parseString(db_str)

    fp_out.write(db_str_reformated.toprettyxml(indent="\t"))

    fp_out.close()

# ----------------------------------------------------------------------------


def family_xml_builder(entries, rfam_acc=None):
    '''
        Expands the Xml4dbDumper object by adding a new family entry.

        entries: The xml entries node to be expanded
        rfam_acc: A specific Rfam family accession
    '''

    cnx = RfamDB.connect()
    cursor = cnx.cursor(dictionary=True)

    cross_refs = {}

    # fetch family fields
    cursor.execute(rs_conf.FAM_QUERY % rfam_acc)
    fam_fields = cursor.fetchall()[0]

    # fetch family specific ncbi_ids
    ncbi_ids = fetch_value_list(rfam_acc, rs_conf.NCBI_IDs_QUERY)

    # fetch family specific ncbi_ids
    pdb_ids = fetch_value_list(rfam_acc, rs_conf.PDB_IDs_QUERY)

    # get pubmed ids
    pmids = get_value_list(fam_fields["pmids"], ',')

    # get dbxrefs (GO, SO)
    dbxrefs = get_value_list(fam_fields["dbxrefs"], ',')

    clan = fetch_value(rs_conf.FAM_CLAN, rfam_acc)

    # need a function here to split dbxrefs
    go_ids = filter(lambda x: x.find("GO") != -1, dbxrefs)
    so_ids = filter(lambda x: x.find("SO") != -1, dbxrefs)

    # update cross references dictionary
    cross_refs["ncbi_taxonomy_id"] = ncbi_ids
    cross_refs["PDB"] = pdb_ids
    cross_refs["PUBMED"] = pmids
    cross_refs["GO"] = go_ids
    cross_refs["SO"] = so_ids
    if clan != None:
        cross_refs["RFAM"] = [clan]

    # add a new family entry to the xml tree
    entry = ET.SubElement(entries, "entry", id=rfam_acc)

    # entry name
    ET.SubElement(entry, "name").text = str(fam_fields["rfam_id"])

    # entry description
    ET.SubElement(entry, "description").text = str(fam_fields["description"])

    # entry dates - common to motifs and clans
    dates = ET.SubElement(entry, "dates")

    created = fam_fields["created"].date().strftime("%d %b %Y")
    updated = fam_fields["updated"].date().strftime("%d %b %Y")

    ET.SubElement(dates, "date", value=created, type="created")
    ET.SubElement(dates, "date", value=updated, type="updated")

    # loop to add cross references
    build_cross_references(entry, cross_refs)

    # expand xml tree with additional fields
    build_additional_fields(
        entry, fam_fields, len(pdb_ids), entry_type='Family')

    cursor.close()
    cnx.disconnect()

# ----------------------------------------------------------------------------


def expand_xml_tree(xml_tree_node, value_list):
    '''
        Expands an xml tree from a point onwards by adding new fields.

        TO BE IMPLEMENTED
    '''

    pass
# ----------------------------------------------------------------------------

# cross_ref_dict will be different for the different types, but the dictionary
# has to be in the same format


def build_cross_references(entry, cross_ref_dict):
    '''
        Expands the entry xml tree by adding the entry's cross references

        entry: The entry node of the xml tree object (xml.etree.ElementTree)
        cross_ref_dict: A dictionary with the entity's cross references in the
        form of ({db_name:[db_key1,db_key2,..],}) where db_name is a string and
        values a list of db ids

    '''

    cross_refs = ET.SubElement(entry, "cross_references")

    for db_name in cross_ref_dict.keys():

        # get db_keys
        db_keys = cross_ref_dict[db_name]
        if len(db_keys) > 0:
            for value in db_keys:
                ET.SubElement(
                    cross_refs, "ref", dbkey=str(value), dbname=db_name)

# ----------------------------------------------------------------------------

# perhaps fields can be a dictionary and move fields to type specific builder??


def build_additional_fields(entry, fields, num_3d_structures, entry_type):
    '''
        This function expands the entry xml field with the additional fields

        entry: This is the xml.etree.ElementTree at the point of entry
        fields: A list of additional fields to expand the entry with

    '''
    add_fields = ET.SubElement(entry, "additional_fields")

    # adding entry type
    ET.SubElement(add_fields, "field", name="entry_type").text = entry_type

    # adding authors
    author_list = get_value_list(fields["author"], rs_conf.AUTH_DEL)

    for author in author_list:
        ET.SubElement(add_fields, "field", name="author").text = author

    if entry_type == "Family":

        # number of species
        ET.SubElement(add_fields, "field", name="num_species").text = str(fields[
            "number_of_species"])
        # number of 3D structures
        ET.SubElement(
            add_fields, "field", name="num_3d_structures").text = str(num_3d_structures)
        # num seed
        ET.SubElement(add_fields, "field", name="num_seed").text = str(fields[
            "num_seed"])
        # num full
        ET.SubElement(add_fields, "field", name="num_full").text = str(fields[
            "num_full"])

        # rna types
        rna_types = get_value_list(fields["type"], rs_conf.RNA_TYPE_DEL)

        for rna_type in rna_types:
            ET.SubElement(add_fields, "field", name="rna_type").text = rna_type

        # has 3d structure
        if num_3d_structures > 0:
            ET.SubElement(
                add_fields, "field", name="has_3d_structure").text = "Yes"
        else:
            ET.SubElement(
                add_fields, "field", name="has_3d_structure").text = "No"

    # perhaps move this to clan and motif xml builder
    else:
        num_families = None

        if entry_type == "Motif":
            num_families = fetch_value(
                rs_conf.NUM_FAMS_MOTIF, fields["motif_acc"])

        elif entry_type == "Clan":
            num_families = fetch_value(
                rs_conf.NUM_FAMS_CLAN, fields["clan_acc"])

        ET.SubElement(add_fields, "field", name="num_families").text = str(
            num_families)

# ----------------------------------------------------------------------------


def motif_xml_builder():
    '''
        TO BE IMPLEMENTED
    '''

    pass

# ----------------------------------------------------------------------------


def clan_xml_builder():
    '''
        TO BE IMPLEMENTED
    '''

    pass

# ----------------------------------------------------------------------------


def create_cross_ref_dict(db_name, db_keys, cross_refs):
    '''
        Updates the cross references' dictionary by adding
    '''

    pass
# ----------------------------------------------------------------------------


def get_value_list(val_str, delimeter=','):  # done
    '''
        val_str: A string of family specific values. This string is a
        concatenation of multiple values related to a single family
        delimeter: The delimeter that will be used to split the values' string
    '''

    val_str = val_str.strip()
    values = val_str.split(delimeter)

    value_list = map(lambda x: x.strip(), values)
    value_list = filter(lambda x: x != '', value_list)

    return value_list


# ----------------------------------------------------------------------------

# maybe move this to DB utils


def fetch_value_list(rfam_acc, query):
    '''
        Retrieves and returns a list of all rfam_acc related values, returned
        by executing the query. Values in list are converted to string format.

        rfam_acc: A family specific accession
        query: A string with the MySQL query to be executed
    '''

    cnx = RfamDB.connect()

    cursor = cnx.cursor(raw=True)

    cursor.execute(query % rfam_acc)

    pdb_structures = cursor.fetchall()

    cursor.close()
    cnx.disconnect()

    return map(lambda x: str(x[0]), pdb_structures)

# ----------------------------------------------------------------------------


def fetch_value(query, accession):
    '''
        Retrieves and returns a value from the database depending to the query
        executed.

        query: The query to be executed in the form of string.
        accession: Rfam specific accession (family, clan, motif)
                   to execute the query on.

        *Query should return a single value
    '''

    cnx = RfamDB.connect()

    cursor = cnx.cursor(raw=True)

    cursor.execute(query % accession)

    value = cursor.fetchall()

    cursor.close()
    cnx.disconnect()

    if len(value) > 0:
        return value[0][0]

    return None

# ----------------------------------------------------------------------------

if __name__ == '__main__':
    outdir = "/Users/ikalvari/Desktop/EB_eye_search"
    xml4db_dumper(outdir)
