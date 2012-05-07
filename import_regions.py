#!/usr/bin/env python

import sys

from csv import DictReader
from optparse import OptionParser

from adipo_sight import *

usage = '%prog [options] <db fn> <regions fn>'
desc = """Import the regions in <regions fn> into adipo sight db <db fn>.
Regions file must have a header row containing name, region_type, chrom, start,
and end literal text separated by commas in any order followed by appropriate
data.  Region types must be literal text in %s"""%region_types
parser = OptionParser(usage=usage,description=desc)

if __name__ == '__main__' :

    opts, args = parser.parse_args(sys.argv[1:])

    if len(args) != 2 :
        parser.error('Exactly 2 non-option arguments are required')

    db_fn, regions_fn = args

    session = get_session(db_fn)

    regions_f = DictReader(open(regions_fn))

    for rec in regions_f :
        region_type_lu = session.query(RegionType.id).filter(RegionType.name == rec['region_type']).first()
        chrom_lu = session.query(Chromosome.id).filter(Chromosome.name == rec['chrom']).first()
        add_d = {'name':rec['name'],
                 'region_type_id':region_type_lu.id,
                 'chrom_id':chrom_lu.id,
                 'start':int(rec['start']),
                 'end':rec['end'].strip(),
                 'id': 'whatever'
                }

        print add_d
        session.add(Region(**add_d))
        #add_new_or_pass(session,add_d,Region)

    session.commit()
    session.close()
