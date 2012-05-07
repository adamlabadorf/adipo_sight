from sqlalchemy import BLOB, Column, Integer, Float, ForeignKey, String, create_engine
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

    id = Column(Integer, primary_key=True)
    name = Column(String)
    region_type_id = Column(Integer, ForeignKey('RegionType.id'))
    region_type = relationship('RegionType')
    chrom_id = Column(Integer, ForeignKey('Chromosome.id'))
    chrom = relationship('Chromosome')
    start = Column(Integer)
    end = Column(Integer)

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
        except :
            pass


def get_session(fn) :

    engine = create_engine('sqlite:///'+fn)
    Base.metadata.create_all(engine)

    session = sessionmaker(bind=engine)()
    session.autocommit = True

    return session


def populate_lookup_tables(engine) :

    session = get_session()

    add_new_or_pass(session,dict(('name',t) for t in region_types),RegionType)
    add_new_or_pass(session,dict(('name',t) for t in data_types),DataType)
    add_new_or_pass(session,dict(('name',t) for t in seq_types),SeqType)
    add_new_or_pass(session,dict(('name',t) for t in chroms),Chromosome)


if __name__ == '__main__' :

    engine = create_engine('sqlite:///adipo_sight.db')

    populate_lookup_tables(engine)

