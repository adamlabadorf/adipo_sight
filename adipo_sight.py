from sqlalchemy import BLOB, Column, Integer, Float, ForeignKey, String, create_engine, UniqueConstraint, Text, CheckConstraint
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref, sessionmaker

Base = declarative_base()

region_types = ['promoter','gene','intergenic','other','ChIP peak','hypersensitive peak']
class RegionType(Base) :
    __tablename__ = 'RegionType'

    id = Column(Integer,primary_key=True)
    name = Column(String,unique=True)

data_types = ['absolute expression','differential expression','ChIP binding']
class DataType(Base) :
    __tablename__ = 'DataType'

    id = Column(Integer,primary_key=True)
    name = Column(String,unique=True)

seq_types = ['nucleotide sequence','motif scores','hypersensitivity tags','mRNASeq tags','ChIP tags']
class SeqType(Base) :
    __tablename__ = 'SeqType'

    id = Column(Integer,primary_key=True)
    name = Column(String,unique=True)

chroms = ['chr%d'%x for x in range(1,23)]+['chrX','chrY']
class Chromosome(Base) :
    __tablename__ = 'Chromosome'

    id = Column(Integer,primary_key=True)
    name = Column(String,unique=True)

class Region(Base):
    __tablename__ = 'Region'
    __table_args__ = (UniqueConstraint('name','region_type_id'),)

    id = Column(Integer, primary_key=True)
    name = Column(String)
    region_type_id = Column(Integer, ForeignKey('RegionType.id'))
    region_type = relationship('RegionType')
    chrom_id = Column(Integer, ForeignKey('Chromosome.id'))
    chrom = relationship('Chromosome')
    start = Column(Integer)
    end = Column(Integer)
    strand = Column(String(1),CheckConstraint("strand IN ('+','-')"))
    notes = Column(Text)

class RegionSet(Base):
    __tablename__ = 'RegionSet'

    id = Column(Integer,primary_key=True)
    name = Column(String,unique=True)
    notes = Column(Text)

class RegionMembership(Base) :
    __tablename__ = 'RegionMembership'
    
    region_set_id = Column(Integer,ForeignKey('RegionSet.id'),primary_key=True)
    region_set = relationship('RegionSet')
    region_id = Column(Integer, ForeignKey('Region.id'),primary_key=True)
    region = relationship('Region')

class RegionData(Base) :
    __tablename__ = 'RegionData'

    id = Column(Integer, primary_key=True)
    region_id = Column(Integer, ForeignKey('Region.id'),primary_key=True)
    region = relationship('Region')
    data_type_id = Column(Integer, ForeignKey('DataType.id'),primary_key=True)
    data_type = relationship('DataType')
    value = Column(Float,nullable=False)
    meta1 = Column(String)
    meta2 = Column(String)
    meta3 = Column(String)

class SeqData(Base) :
    __tablename__ = 'SeqData'

    id = Column(Integer, primary_key=True)
    region_id = Column(Integer, ForeignKey('Region.id'),primary_key=True)
    region = relationship('Region')
    seq_type_id = Column(Integer, ForeignKey('SeqType.id'),primary_key=True)
    region = relationship('SeqType')
    value = Column(BLOB,nullable=False)
    meta1 = Column(BLOB)
    meta2 = Column(BLOB)
    meta3 = Column(BLOB)

def add_new_or_pass(session,types,cls) :
    for r in types :
        try :
            session.add(cls(**r))
            session.commit()
        except Exception, e:
            print e


def get_session(fn) :

    engine = create_engine('sqlite:///'+fn)
    Base.metadata.create_all(engine)

    session = sessionmaker(bind=engine)()

    return session


def populate_lookup_tables(fn) :

    session = get_session(fn)

    add_new_or_pass(session,[{'name':t} for t in region_types],RegionType)
    add_new_or_pass(session,[{'name':t} for t in data_types],DataType)
    add_new_or_pass(session,[{'name':t} for t in seq_types],SeqType)
    add_new_or_pass(session,[{'name':t} for t in chroms],Chromosome)

    session.commit()
    session.close()

if __name__ == '__main__' :

    db_fn = 'adipo_sight.db'

    populate_lookup_tables(db_fn)

