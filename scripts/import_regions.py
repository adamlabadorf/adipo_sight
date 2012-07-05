#!/usr/bin/env python

import sys
import logging

from csv import DictReader
from optparse import OptionParser

import sqlalchemy

from adipo_sight.db import *

usage = '%prog [options] <db fn> <regions fn>'
desc = """Import the regions in <regions fn> into adipo sight db <db fn>.
Regions file must have a header row containing name, region_type, chrom, start,
and end literal text separated by commas in any order followed by appropriate
data.  Region types must be literal text in %s.  The field notes may optionally
be included."""%region_types
parser = OptionParser(usage=usage,description=desc)
parser.add_option('--insert-type',dest='insert_type',action='store_true',default=False,help='insert records for region types not found in db instead of skipping record')
parser.add_option('--insert-chrom',dest='insert_chrom',action='store_true',default=False,help='insert records for chromosomes not found in db instead of skipping record')
parser.add_option('--replace',dest='replace',action='store_true',default=False,help='replace existing records instead of inserting')
parser.add_option('-s','--silent',dest='silent',action='store_true',default=False, help='do not print out any warnings')


if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 2 :
        parser.error('Exactly 2 non-option arguments are required')

    if opts.silent :
        log_level = logging.CRITICAL
    else :
        log_level = logging.WARNING
    logging.basicConfig(stream=sys.stderr,level=log_level,format='%(levelname)s: %(message)s')

    db_fn, regions_fn = args

    session = get_session(db_fn)

    regions_f = DictReader(open(regions_fn))

    for rec in regions_f :

        region_type_lu = session.query(RegionType.id).filter(RegionType.name == rec['region_type']).first()
        if region_type_lu is None :
            if opts.insert_type :
                session.add(RegionType(name=rec['region_type']))
                region_type_lu = session.query(RegionType.id).filter(RegionType.name == rec['region_type']).first()
            else :
                logging.warn('region type %s not found in db and --insert not specified, skipping'%rec['region_type'])
                continue

        chrom_lu = session.query(Chromosome.id).filter(Chromosome.name == rec['chrom']).first()
        if chrom_lu is None :
            if opts.insert_chrom :
                session.add(Chromosome(name=rec['chrom']))
                chrom_lu = session.query(Chromosome.id).filter(Chromosome.name == rec['chrom']).first()
            else :
                logging.warn('chrom %s not found in db and --insert not specified, skipping'%rec['chrom'])
                continue

        add_d = {'name':rec['name'],
                 'region_type_id':region_type_lu.id,
                 'chrom_id':chrom_lu.id,
                 'start':int(rec['start']),
                 'end':int(rec['end']),
                 'notes':rec.get('notes')
                }

        try :
            session.add(Region(**add_d))
            session.commit()
        except sqlalchemy.exc.IntegrityError :
            session.rollback()
            if opts.replace :
                r = session.query(Region) \
                    .filter(Region.name == add_d['name'] and
                            Region.region_type_id == add_d['region_type_id']) \
                    .first()

                r.chrom_id = add_d['chrom_id']
                r.start = add_d['start']
                r.end = add_d['end']
                r.notes = add_d['notes']
                session.commit()
            else :
                logging.warn('region with name/type already exists and --replace not specified, skipping %s'%rec)
        except Exception, e: 
            session.rollback()
            logging.warn('exception adding record: %s - %s'%(rec,e))

        #add_new_or_pass(session,add_d,Region)

    session.close()
