from csv import DictReader

from sqlalchemy import BLOB, Column, Integer, Float, ForeignKey, String, create_engine, UniqueConstraint, Text, CheckConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref, sessionmaker

Base = declarative_base()

kgXref_fieldnames = "kgID", "mRNA", "spID", "spDisplayID", "geneSymbol", "refseq", "protAcc", "description", "rfamAcc", "tRnaName"
#uc007afh.1	NM_008866	P97823	LYPA1_MOUSE	Lypla1	NM_008866	NP_032892	acyl-protein thioesterase 1
class kgXref(Base) :
    __tablename__ = 'kgXref'

    id = Column(Integer,primary_key=True)
    kgID = Column(String,ForeignKey('knownGene.name'))
    knownGene = relationship('knownGene')
    mRNA = Column(String)
    spID = Column(String)
    spDisplayID = Column(String)
    geneSymbol = Column(String)
    refseq = Column(String)
    protAcc = Column(String)
    description = Column(String)
    rfamAcc = Column(String)
    tRnaName = Column(String)


knownGene_fieldnames = "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "proteinID", "alignID"
#uc007aeu.1	chr1	-	3204562	3661579	3206102	3661429	3	3204562,3411782,3660632,	3207049,3411982,3661579,	Q5GH67	uc007aeu.1
class knownGene(Base) :
    __tablename__ = 'knownGene'

    name = Column(String,primary_key=True)
    chrom = Column(String)
    strand = Column(String)
    txStart = Column(Integer)
    txEnd = Column(Integer)
    cdsStart = Column(Integer)
    cdsEnd = Column(Integer)
    exonCount = Column(String)
    exonStarts = Column(String)
    exonEnds = Column(String)
    proteinID = Column(String)
    alignID = Column(String)

def get_session(fn) :

    engine = create_engine('sqlite:///'+fn)
    Base.metadata.create_all(engine)

    session = sessionmaker(bind=engine)()

    return session

def build_db(knownGene_fn,kgXref_fn,db_fn) :

    session = get_session(db_fn)

    knownGene_reader = DictReader(open(knownGene_fn),delimiter='\t',fieldnames=knownGene_fieldnames)
    for r in knownGene_reader :
        for k in 'txStart','txEnd','cdsStart','cdsEnd','exonCount' :
            r[k] = int(r[k])
        session.add(knownGene(**r))

    kgXref_reader = DictReader(open(kgXref_fn),delimiter='\t',fieldnames=kgXref_fieldnames)
    for r in kgXref_reader :
        session.add(kgXref(**r))

    session.commit()
    session.close()

if __name__ == '__main__' :

    build_db('/nfs/genomes/mouse_gp_jul_07/anno/knownGene-mm9.txt',
             '/nfs/genomes/mouse_gp_jul_07/anno/kgXref-mm9.txt',
             'knownGene.db')
