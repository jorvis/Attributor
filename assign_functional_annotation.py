#!/usr/bin/env python3

"""

This is to assign functional annotation to gene models using a wide variety of evidence.  Meant to
support unlimited hierarchy rules, this utility relies on a user-created configuration file.

Current limitations:

- Rules are evaluated on a per-gene basis, so no rules are currently possible that
  would need to consider annotations of other genes in the set.

- IGS 'PFunc' hierarchy
http://imgur.com/1odYcT5


Overall steps

- Parse FASTA file and create an annotation index of them
- Parse annotations
  - If a list file is defined more than once, cache the results?  Is this feasible given memory constraints?
  
"""

import argparse
import bioannotation
import biocodeutils
import biothings
import os
import sqlite3
import yaml

def main():
    parser = argparse.ArgumentParser( description='Put a description of your script here')

    ## output file to be written
    parser.add_argument('-f', '--input_fasta', type=str, required=True, help='Protein FASTA file of source molecules' )
    parser.add_argument('-c', '--config_file', type=str, required=True, help='Configuration file for annotation' )
    parser.add_argument('-o', '--output_base', type=str, required=True, help='Base name/path of output files to be created' )
    args = parser.parse_args()

    sources_log_fh = open("{0}.sources.log".format(args.output_base), 'wt')

    configuration = yaml.load(open(args.config_file).read())
    check_configuration(configuration)
    evidence = parse_evidence_config( configuration )
    default_product_name = configuration['general']['default_product_name']

    # stores any active SQLite3 db connections
    db_conn = dict()

    # this is a dict of biothings.Polypeptide objects
    polypeptides = initialize_polypeptides( sources_log_fh, args.input_fasta, default_product_name )

    for label in configuration['order']:
        if label not in evidence:
            raise Exception("ERROR: There is a label '{0}' in the 'order' section of the conf file that isn't present in the 'evidence' section".format(label))

        if evidence[label]['type'] == 'HMMer3_htab':
            index_label = evidence[label]['index']

            if index_label in db_conn:
                # if we already have a connection to this database, use it again
                index_conn = db_conn[index_label]
            else:
                # create one
                index_path = "{0}.{1}.sqlite3".format(args.output_base, index_label)
                index_conn = sqlite3.connect(index_path)
                db_conn[index_label] = index_conn
                initialize_hmm_results_db(index_conn)
            
            parse_hmmer3_htab(polypeptides=polypeptides, config=configuration, index=index_conn)
        elif evidence[label]['type'] == 'RAPSearch2':
            pass
        else:
            raise Exception("ERROR: Unsupported evidence type '{0}' with label '{1}' in configuration file".format(evidence[label]['type'], label))

    # close all db connections
    for label in db_conn:
        db_conn[label].close()


def check_configuration(conf):
    """
    Performs any basic checks on the annotation configuration file format/syntax/values.  Ideally done
    before most of the rest of the script to save wasted compute time.
    """
    # make sure each of the expected sections are there
    for section in ['general', 'indexes', 'order', 'evidence']:
        if section not in conf:
            raise Exception("ERROR: Expected a section called '{0}' in the annotation config file, but didn't find one.".format(section))

    # make sure there aren't any indexes referenced in the evidence section which are not defined in the indexes section
    indexes = list()
    for label in conf['indexes']:
        indexes.append(label)

    for item in conf['evidence']:
        if 'index' in item and item['index'] not in indexes:
            raise Exception("ERROR: Evidence item '{0}' references and index '{1}' not found in the indexes section of the config file".format(item['label'], item['index']))

def initialize_hmm_results_db(conn):
    """
    HTAB:
    col  perl-col   description
    1      [0]      HMM accession
    2      [1]      Date search was run (if available), otherwise date of htab parse
    3      [2]      Length of the HMM (not populated if -s is used)
    4      [3]      Search program
    5      [4]      Database file path
    6      [5]      Sequence accession
    7      [6]      Alignment start position on HMM match - hmm-f
    8      [7]      Alignment end position on HMM match - hmm-t
    9      [8]      Alignment start position on sequence - seq-f
    10     [9]      Alignment end position on sequence - seq-t
    11     [10]     frame (only populated if --frames search is run on nucleotide sequence)
    12     [11]     Domain score
    13     [12]     Total score
    14     [13]     Index of domain hit
    15     [14]
    16     [15]     HMM description (may be truncated by hmmsearch or hmmpfam if -s is used)
    17     [16]     Sequence description (may be truncated by hmmsearch or hmmpfam)
    18     [17]     Total score trusted cutoff (not populated if -s is used)
    19     [18]     Total score noise cutoff (not populated if -s is used)
    20     [19]     Expect value for total hit
    21     [20]     Expect value for domain hit
    22     [21]     Domain score trusted cutoff (egad..hmm2.trusted_cutoff2) (not populated if -s is used)
    23     [22]     Domain score noise cutoff (egad..hmm2.noise_cutoff2) (not populated if -s is used)
    24     [23]     Total score gathering threshold (not populated if -s is used)
    25     [24]     Domain score gathering threshold (not populated if -s is used)
    """
    curs = conn.cursor()

    curs.execute("""
        CREATE TABLE hmm_hit (
            id                integer primary key,
            qry_id            text,
            qry_start         integer,
            qry_end           integer,
            hmm_accession     text,
            hmm_length        integer,
            hmm_start         integer,
            hmm_end           integer,
            domain_score      real,
            total_score       real,
            total_score_tc    real,
            total_score_nc    real,
            total_score_gc    real,
            domain_score_tc   real,
            domain_score_nc   real,
            domain_score_gc   real
        )
    """)

    curs.close()
    conn.commit()
        
def initialize_polypeptides( log_fh, fasta_file, default_name ):
    '''
    Reads a FASTA file of (presumably) polypeptide sequences and creates a dict of Polypeptide
    objects, keyed by ID, with bioannotation.FunctionalAnnotation objects attached.
    '''
    seqs = biocodeutils.fasta_dict_from_file( fasta_file )

    polypeptides = dict()

    for seq_id in seqs:
        polypeptide = biothings.Polypeptide( id=seq_id, length=len(seqs[seq_id]['s']), residues=seqs[seq_id]['s'] )
        annotation = bioannotation.FunctionalAnnotation(product_name=default_name)
        log_fh.write("INFO: {0}: Set initial product name to '{1}'\n".format(seq_id, default_name))
        polypeptide.annotation = annotation
        
        polypeptides[seq_id] = polypeptide
    
    return polypeptides
        
def parse_evidence_config(conf):
    """
    Parses the 'evidence' section of the annotation config file, and returns a dict where each key
    is the label of that evidence and the value is a dict of the other key/value pairs
    """
    ev = dict()

    for entry in conf['evidence']:
        # make sure there aren't duplicates
        label = entry['label']

        if label in ev:
            raise Exception("ERROR: duplicate label found in evidence track: {0}".format(label))
        else:
            ev[label] = dict()
            for key in entry:
                if key is not 'label':
                    ev[label][key] = entry[key]

    return ev

def parse_hmmer3_htab(polypeptides=None, config=None, index=None):
    pass
    


if __name__ == '__main__':
    main()







