from sqlalchemy import BLOB, Column, Integer, Float, ForeignKey, String, create_engine, UniqueConstraint, Text, CheckConstraint, Table
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, backref, sessionmaker

Base = declarative_base()

region_types = ['promoter','gene','intergenic','other','ChIP peak','hypersensitive peak']
class RegionType(Base) :
    __tablename__ = 'RegionType'

    id = Column(Integer,primary_key=True)
    name = Column(String,unique=True)

data_types = ['control expression','experiment expression','log2 expression fold change','ChIP binding']
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

#region_membership_table = Table('RegionMembership', Base.metadata,
#                            Column('region_set_id', Integer, ForeignKey('RegionSet.id'),primary_key=True),
#                            Column('region_id', Integer, ForeignKey('Region.id'),primary_key=True)
#                          )

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
    region_data = relationship("RegionData",cascade="all, delete, delete-orphan")
    seq_data = relationship("SeqData",cascade="all, delete, delete-orphan")

class RegionSet(Base):
    __tablename__ = 'RegionSet'

    id = Column(Integer,primary_key=True)
    name = Column(String,unique=True)
    #regions = relationship("Region",secondary=region_membership_table,backref="region_sets")
    member_regions = relationship("RegionMembership")
    notes = Column(Text)

class RegionMembership(Base) :
    __tablename__ = 'RegionMembership'
    
    region_set_id = Column(Integer,ForeignKey('RegionSet.id'),primary_key=True)
    region_set = relationship('RegionSet')
    region_id = Column(Integer, ForeignKey('Region.id'),primary_key=True)
    region = relationship('Region',backref="membership")
    dist_to_feature = Column(Integer)

class Condition(Base) :
    __tablename__ = 'Condition'

    id = Column(Integer,primary_key=True)
    name = Column(String)

class RegionData(Base) :
    __tablename__ = 'RegionData'
    __table_args__ = (UniqueConstraint('region_id','data_type_id','condition_id'),)

    id = Column(Integer,primary_key=True)
    region_id = Column(Integer, ForeignKey('Region.id'))
    region = relationship('Region')
    data_type_id = Column(Integer, ForeignKey('DataType.id'))
    data_type = relationship('DataType')
    condition_id = Column(Integer, ForeignKey('Condition.id'))
    condition = relationship('Condition')
    value = Column(Float,nullable=False)
    meta1lbl = Column(String)
    meta1 = Column(String)
    meta2lbl = Column(String)
    meta2 = Column(String)
    meta3lbl = Column(String)

class SeqData(Base) :
    __tablename__ = 'SeqData'
    __table_args__ = (UniqueConstraint('region_id','seq_type_id','condition_id'),)

    id = Column(Integer, primary_key=True)
    region_id = Column(Integer, ForeignKey('Region.id'))
    region = relationship('Region')
    seq_type_id = Column(Integer, ForeignKey('SeqType.id'))
    seq_type = relationship('SeqType')
    condition_id = Column(Integer, ForeignKey('Condition.id'))
    condition = relationship('Condition')
    value = Column(BLOB,nullable=False)
    meta1lbl = Column(String)
    meta1 = Column(String)
    meta2lbl = Column(String)
    meta2 = Column(String)
    meta3lbl = Column(String)

def add_new_or_pass(session,types,cls) :
    objs = []
    for r in types :
        try :
            new_obj = cls(**r)
            session.add(new_obj)
            session.commit()
            objs.append(new_obj)
        except Exception, e:
            session.rollback()
            objs.append(None)
            #print e

    return objs

def get_lu_or_add(session,name,cls,fields=None) :
    r = session.query(cls).filter(getattr(cls,'name')==name).first()
    if r is None :
        # add a new record and return it
        kwargs = fields or {}
        r = cls(name=name,**kwargs)
        session.add(r)
        session.commit()
    return r

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

