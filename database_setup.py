from sqlalchemy import Column, Integer, String, Float, create_engine
from sqlalchemy.ext.declarative import declarative_base

Base = declarative_base()


class Genome(Base):
    __tablename__ = 'genomes'

    name = Column(Integer, nullable=False)
    id = Column(Integer, primary_key=True)

    taxonomy = Column(String(1000))
    nsti = Column(Float)
    metab_net = Column(String(100000))
    genome = Column(String(1000))


engine = create_engine('sqlite:///gg_genomes.db')
Base.metadata.create_all(engine)
