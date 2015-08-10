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
  - 
"""

import argparse
import bioannotation
import biocodeutils
import biothings
import math
import os
import re
import sqlite3
import yaml

def main():
    parser = argparse.ArgumentParser( description='Assigns functional annotation based on user-configurable evidence tiers')

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
    polypeptides = initialize_polypeptides(sources_log_fh, args.input_fasta, default_product_name)

    for label in configuration['order']:
        if label not in evidence:
            raise Exception("ERROR: There is a label '{0}' in the 'order' section of the conf file that isn't present in the 'evidence' section".format(label))

        if evidence[label]['type'] == 'HMMer3_htab':
            index_conn, ev_db_conn = get_or_create_db_connections(type_ev='hmm_ev', configuration=configuration,
                                         evidence=evidence, label=label, db_conn=db_conn,
                                         output_base=args.output_base)

            # then apply the evidence
            apply_hmm_evidence(polypeptides=polypeptides, ev_conn=ev_db_conn, config=configuration,
                               ev_config=evidence[label], label=label, index_conn=index_conn, log_fh=sources_log_fh)
                
        elif evidence[label]['type'] == 'RAPSearch2':
            index_conn, ev_db_conn = get_or_create_db_connections(type_ev='blast_ev', configuration=configuration,
                                         evidence=evidence, label=label, db_conn=db_conn,
                                         output_base=args.output_base)

            # then apply the evidence
            apply_blast_evidence(polypeptides=polypeptides, ev_conn=ev_db_conn, config=configuration,
                                 ev_config=evidence[label], label=label, index_conn=index_conn, log_fh=sources_log_fh)
        else:
            raise Exception("ERROR: Unsupported evidence type '{0}' with label '{1}' in configuration file".format(evidence[label]['type'], label))

    # close all db connections
    for label in db_conn:
        db_conn[label].close()

    perform_final_checks(polypeptides=polypeptides, config=configuration, log_fh=sources_log_fh)

    # Write the FASTA
    polyset = biothings.PolypeptideSet()
    polyset.load_from_dict(polypeptides)
    polyset.write_fasta(path="{0}.faa".format(args.output_base))

    
def already_indexed(path=None, index=None):
    curs = index.cursor()
    curs.execute("SELECT id FROM data_sources WHERE source_path = ?", (path, ))
    found = False

    for row in curs:
        found = True
        break
    
    curs.close()
    return found

def apply_blast_evidence(polypeptides=None, ev_conn=None, config=None, ev_config=None, label=None, index_conn=None, log_fh=None):
    """
    Uses BLAST (or similar) evidence to assign functional evidence to polypeptides.  Description of arguments:

    polypeptides: a dict of biothings.Polypeptides objects, keyed on ID
    ev_conn: SQLite3 connection to the parsed BLAST(ish) search evidence db for this set of searches
    config: The yaml object for the parsed annotation config file
    ev_config: The parsed evidence section for this label within the annotation config file
    label: Label for the evidence track entry within thet annotation config file
    index_conn:  SQLite3 connection to the reference index for the database searched
    """
    default_product = config['general']['default_product_name']
    ev_curs = ev_conn.cursor()
    ev_qry = "SELECT sbj_id, align_len, perc_identity, eval, bit_score FROM blast_hit WHERE qry_id = ? ORDER BY eval ASC"

    print("DEBUG: Applying {1} results to {0} polypeptides".format(len(polypeptides), label))
    
    if 'debugging_polypeptide_limit' in config['general']:
        DEBUG_LIMIT = config['general']['debugging_polypeptide_limit']

    # Are coverage cutoffs defined?
    query_cov_cutoff = None
    match_cov_cutoff = None
    if 'query_cov' in ev_config:
        query_cov_cutoff = int(ev_config['query_cov'].rstrip('%'))

    for id in polypeptides:
        polypeptide = polypeptides[id]
        print("DEBUG: Parsing {0} evidence for polypeptide ID {1}, length: {2}".format(label, id, polypeptide.length))
        annot = polypeptide.annotation

        DEBUG_LIMIT = DEBUG_LIMIT - 1
        if DEBUG_LIMIT == 0:
            break

        if config['general']['allow_attributes_from_multiple_sources'] == 'Yes':
            raise Exception("ERROR: Support for the general:allow_attributes_from_multiple_sources=Yes setting not yet implemented")
        else:
            if annot.product_name != default_product: continue

            for ev_row in ev_curs.execute(ev_qry, (polypeptide.id,)):
                if query_cov_cutoff is not None:
                    perc_coverage = (ev_row[1] / polypeptide.length)*100
                    if perc_coverage < query_cov_cutoff:
                        #print("\tSkipping accession {0} because coverage {1} doesn't meet cutoff {2} requirements".format(
                        #    ev_row[0], perc_coverage, query_cov_cutoff))
                        continue
                
                blast_annot = get_blast_result_info(conn=index_conn, accession=ev_row[0], config=config)
                annot.product_name = blast_annot.product_name
                log_fh.write("INFO: {1}: Set product name to '{0}' from {3} hit to {2}\n".format(
                        annot.product_name, id, ev_row[0], label))
                
                annot.gene_symbol = blast_annot.gene_symbol
                log_fh.write("INFO: {1}: Set gene_symbol to '{0}' from {3} hit to {2}\n".format(
                        annot.product_name, id, ev_row[0], label))

                for go_annot in blast_annot.go_annotations:
                    annot.add_go_annotation(go_annot)

                for ec_num in blast_annot.ec_numbers:
                    annot.add_ec_number(ec_num)

                # If we get this far we've assigned annotation and don't want to look at any more
                break

    ev_curs.close()
                
    
def apply_hmm_evidence(polypeptides=None, ev_conn=None, config=None, ev_config=None, label=None, index_conn=None, log_fh=None):
    """
    Uses HMM evidence to assign functional evidence to polypeptides.  Description of arguments:

    polypeptides: a dict of biothings.Polypeptides objects, keyed on ID
    ev_conn: SQLite3 connection to the parsed HMM search evidence db for this set of searches
    config: The yaml object for the parsed annotation config file
    ev_config: The parsed evidence section for this label within the annotation config file
    label: Label for the evidence track entry within thet annotation config file
    index_conn:  SQLite3 connection to the reference index for the database searched
    """
    default_product = config['general']['default_product_name']
    ev_curs = ev_conn.cursor()
    ev_qry = "SELECT hmm_accession, total_score FROM hmm_hit WHERE qry_id = ? ORDER BY total_hit_eval ASC"

    go_curs = index_conn.cursor()
    go_qry = "SELECT go_id FROM hmm_go WHERE hmm_id = ?"

    ec_curs = index_conn.cursor()
    ec_qry = "SELECT ec_id FROM hmm_ec WHERE hmm_id = ?"
    
    acc_main_curs = index_conn.cursor()
    hmm_class_limit = None

    if 'class' in ev_config:
        hmm_class_limit = ev_config['class']

    print("DEBUG: Applying HMM results to {0} polypeptides".format(len(polypeptides)))

    if 'debugging_polypeptide_limit' in config['general']:
        DEBUG_LIMIT = config['general']['debugging_polypeptide_limit']
    
    for id in polypeptides:
        print("DEBUG: Parsing {0} evidence for polypeptide ID {1}".format(label, id))
        polypeptide = polypeptides[id]
        annot = polypeptide.annotation

        DEBUG_LIMIT = DEBUG_LIMIT - 1
        if DEBUG_LIMIT == 0:
            break

        if config['general']['allow_attributes_from_multiple_sources'] == 'Yes':
            raise Exception("ERROR: Support for the general:allow_attributes_from_multiple_sources=Yes setting not yet implemented")
        else:
            if annot.product_name != default_product: continue

            for ev_row in ev_curs.execute(ev_qry, (polypeptide.id,)):
                acc_main_qry = "SELECT version, hmm_com_name, ec_num, isotype, id FROM hmm WHERE version = ? or accession = ?"

                if hmm_class_limit is None:
                    acc_main_qry_args = (ev_row[0], ev_row[0], )
                else:
                    acc_main_qry += " AND isotype = ?"
                    acc_main_qry_args = (ev_row[0], ev_row[0], hmm_class_limit)
                
                for acc_main_row in acc_main_curs.execute(acc_main_qry, acc_main_qry_args):
                    annot.product_name = acc_main_row[1]
                    log_fh.write("INFO: {1}: Set product name to '{0}' from {3} hit to {2}, isotype:{3}\n".format(
                        annot.product_name, id, ev_row[0], hmm_class_limit, label))

                    if acc_main_row[2] is not None:
                        annot.gene_symbol  = acc_main_row[2]
                        log_fh.write("INFO: {1}: Set gene_symbol to '{0}' from {3} hit to {2}, isotype:{3}\n".format(
                            annot.gene_symbol, id, ev_row[0], hmm_class_limit, label))

                    ## add any matching GO terms
                    for go_row in go_curs.execute(go_qry, (acc_main_row[4],)):
                        annot.add_go_annotation(bioannotation.GOAnnotation(go_id=go_row[0]))

                    ## add any matching EC numbers
                    for ec_row in ec_curs.execute(ec_qry, (acc_main_row[4],)):
                        annot.add_ec_number(bioannotation.ECAnnotation(number=ec_row[0]))
                    
                break
        
    acc_main_curs.close()
    ev_curs.close()
    go_curs.close()
    ec_curs.close()
        

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


def get_blast_result_info(conn=None, accession=None, config=None):
    curs = conn.cursor()
    ec_curs = conn.cursor()
    go_curs = conn.cursor()
    annot = bioannotation.FunctionalAnnotation(product_name=config['general']['default_product_name'])

    if accession.startswith('UniRef100_'):
        accession = accession.lstrip('UniRef100_')

    # This is currently specific to my uniref index.
    # First we need to get the ID from the accession
    qry = """
          SELECT u.id, ua.accession, u.full_name, u.symbol
          FROM uniref u
               JOIN uniref_acc ua ON u.id=ua.id
          WHERE ua.accession = ?
          """

    for row in curs.execute(qry, (accession,)):
        uniref_id = row[0]
        annot.product_name = row[2]
        annot.gene_symbol = row[3]

        qry = "SELECT ec_num FROM uniref_ec WHERE id = ?"
        for ec_row in ec_curs.execute(qry, (uniref_id,)):
            annot.add_ec_number( bioannotation.ECAnnotation(number=ec_row[0]) )

        qry = "SELECT go_id FROM uniref_go WHERE id = ?"
        for go_row in go_curs.execute(qry, (uniref_id,)):
            annot.add_go_annotation( bioannotation.GOAnnotation(go_id=go_row[0], with_from=row[1]) )

    curs.close()
    ec_curs.close()
    go_curs.close()
    return annot
        
def get_or_create_db_connections(type_ev=None, configuration=None, evidence=None, label=None,
                                 db_conn=None, output_base=None):
    """
    type_ev must be either 'hmm_ev' or 'blast_ev'
    """
    index_label = evidence[label]['index']

    # Use any existing index connection, else attach to it.
    if index_label in db_conn:
        index_conn = db_conn[index_label]
    else:
        index_conn = sqlite3.connect(configuration['indexes'][index_label])
        db_conn[index_label] = index_conn

    # Attach to or create an evidence database
    ev_db_path = "{0}.{1}.sqlite3".format(output_base, type_ev)
    if os.path.exists(ev_db_path):
        ev_db_conn = sqlite3.connect(ev_db_path)
    else:
        ev_db_conn = sqlite3.connect(ev_db_path)
        if type_ev == 'hmm_ev':
            initialize_hmm_results_db(ev_db_conn)
        elif type_ev == 'blast_ev':
            initialize_blast_results_db(ev_db_conn)

    db_conn[type_ev] = ev_db_conn

    # only parse the evidence if the list isn't already in the database
    if not already_indexed(path=evidence[label]['path'], index=ev_db_conn):

        if type_ev == 'hmm_ev':
            index_hmmer3_htab(path=evidence[label]['path'], index=ev_db_conn)
            ev_db_conn.commit()
            # update the database search indexes
            hmm_ev_db_curs = ev_db_conn.cursor()
            hmm_ev_db_curs.execute("DROP INDEX IF EXISTS hmm_hit__qry_id")
            hmm_ev_db_curs.execute("CREATE INDEX hmm_hit__qry_id ON hmm_hit (qry_id)")
        elif type_ev == 'blast_ev':
            index_rapsearch2_m8(path=evidence[label]['path'], index=ev_db_conn)
            ev_db_conn.commit()
            blast_ev_db_curs = ev_db_conn.cursor()
            blast_ev_db_curs.execute("DROP INDEX IF EXISTS blast_hit__qry_id")
            blast_ev_db_curs.execute("CREATE INDEX blast_hit__qry_id ON blast_hit (qry_id)")

    return (index_conn, ev_db_conn)
        
def index_hmmer3_htab(path=None, index=None):
    curs = index.cursor()
    parsing_errors = 0

    qry = """
          INSERT INTO hmm_hit (qry_id, qry_start, qry_end, hmm_accession, hmm_length, hmm_start, hmm_end, 
                               domain_score, total_score, total_score_tc, total_score_nc, total_score_gc,
                               domain_score_tc, domain_score_nc, domain_score_gc, total_hit_eval,
                               domain_hit_eval)
          VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)
    """

    for file in biocodeutils.read_list_file(path):
        print("DEBUG: parsing: {0}".format(file))
        for line in open(file):
            line = line.rstrip()
            cols = line.split("\t")

            ## not sure what this is, but some lines have columns 7+ as these values:
            #   hmm       to      ali     to              bias
            if cols[6] == 'hmm': continue

            # the following columns in the htab files are nullable, which are filled with various widths of '-' characters
            for i in (6,7,8,9,11,19,20):
                if '--' in cols[i]:
                    cols[i] = None
                elif len(cols[i]) > 0:
                    try:
                        cols[i] = float(cols[i])
                    except ValueError:
                        #print("DEBUG: Could not convert this to a float: {0}".format(cols[i]))
                        parsing_errors += 1
                        cols[i] = None
                
            curs.execute(qry, (cols[5], cols[8], cols[9], cols[0], int(cols[2]), cols[6], cols[7], 
                               cols[11], float(cols[12]), float(cols[17]), float(cols[18]), float(cols[23]),
                               float(cols[21]), float(cols[22]), float(cols[24]), cols[19], cols[20]))

    curs.execute("INSERT INTO data_sources (source_path) VALUES (?)", (path,))
    curs.close()

    if parsing_errors > 0:
        print("WARN: There were {0} parsing errors (columns converted to None) when processing {1}\n".format(parsing_errors, path))

def index_rapsearch2_m8(path=None, index=None):
    curs = index.cursor()
    parsing_errors = 0

    # The E-value column can be either the E-value directly or log(E-value), depending on
    #  the version and options used.  Luckily, the header line tells us which it is.
    logged_eval = False

    qry = """
          INSERT INTO blast_hit (qry_id, sbj_id, align_len, qry_start, qry_end, sbj_start,
                                 sbj_end, perc_identity, eval, bit_score)
          VALUES (?,?,?,?,?,?,?,?,?,?)
    """
    
    for file in biocodeutils.read_list_file(path):
        print("DEBUG: parsing: {0}".format(file))
        for line in open(file):
            if line[0] == '#':
                m = re.search('log\(e\-value\)', line)
                if m:
                    logged_eval = True
                continue
            
            line = line.rstrip()
            cols = line.split("\t")

            if len(cols) != 12:
                continue

            if logged_eval == True:
                try:
                    cols[10] = math.pow(10, float(cols[10]))
                except OverflowError:
                    # RapSearch2 sometimes reports miniscule e-values, such a log(eval) of > 1000
                    #  These are outside of the range of Python's double.  In my checking of these
                    #  though, their alignments don't warrant a low E-value at all.  Skipping them.
                    print("Warning: Skipping a RAPSearch2 row:  overflow error converting E-value ({0}) on line: {1}".format(cols[10], line))
                    continue
                
            curs.execute(qry, (cols[0], cols[1], int(cols[3]), int(cols[6]), int(cols[7]), int(cols[8]),
                               int(cols[9]), float(cols[2]), cols[10], float(cols[11])))

    curs.execute("INSERT INTO data_sources (source_path) VALUES (?)", (path,))
    curs.close()

        
def initialize_blast_results_db(conn):
    curs = conn.cursor()

    curs.execute("""
        CREATE TABLE blast_hit (
            id                integer primary key,
            qry_id            text,
            sbj_id            text,
            align_len         integer,
            qry_start         integer,
            qry_end           integer,
            sbj_start         integer,
            sbj_end           integer,
            perc_identity     real,
            eval              real,
            bit_score         real
        )
    """)

    curs.execute("""
        CREATE TABLE data_sources (
            id                integer primary key,
            source_path         text
        )
    """)

    curs.close()
    conn.commit()
                
def initialize_hmm_results_db(conn):
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
            total_hit_eval    real,
            domain_score_tc   real,
            domain_score_nc   real,
            domain_score_gc   real,
            domain_hit_eval   real
        )
    """)

    curs.execute("""
        CREATE TABLE data_sources (
            id                integer primary key,
            source_path         text
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


def perform_final_checks(polypeptides=None, config=None, log_fh=None):
    """
    Does a round of checks we want to perform on an annotated set of polypeptides
    before exporting them.  Currently:

    - Make sure a gene product name is assigned.  This might accidentally become "None" if
      a match existed to a subject which didn't have a name properly entered in the index.
    """
    for id in polypeptides:
        polypeptide = polypeptides[id]

        if polypeptide.annotation.product_name is None:
            log_fh.write("WARNING: {0}: Somehow made it through annotation with no product name.  Setting to default\n".format(id))
            polypeptide.annotation.product_name = config['general']['default_product_name']
    
    


if __name__ == '__main__':
    main()







