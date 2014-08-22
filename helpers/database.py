from sqlalchemy import create_engine
from sqlalchemy import Column, ForeignKey, Integer, Float, Text, DateTime
from sqlalchemy.orm import sessionmaker, relationship
from sqlalchemy.ext.declarative import declared_attr, declarative_base

Session = sessionmaker()

class Base (object):
    @declared_attr
    def __tablename__(cls):
        import re
        camel_case = re.compile('[A-Z][a-z]+')
        table_words = camel_case.findall(cls.__name__)
        return '_'.join(w.lower() for w in table_words)

Base = declarative_base(cls=Base)

class Benchmarks (Base):
    benchmark_id = Column(Integer, primary_key=True, autoincrement=True)
    start_time = Column(DateTime)
    name = Column(Text)
    description = Column(Text)

    def __init__(self, name=None, description=None):
        import datetime
        self.start_time = datetime.datetime.now()
        self.name = name
        self.description = description

    @property
    def id(self):
        return self.benchmark_id


class BenchmarkProtocols (Base):
    benchmark_id = Column(Integer, primary_key=True)
    protocol_id = Column(Integer, primary_key=True)
    
    def __init__(self, benchmark_id, protocol_id):
        self.benchmark_id = benchmark_id
        self.protocol_id = protocol_id


class TracerOutput (Base):
    protocol_id = Column(Integer, primary_key=True)
    output = Column(Text)

    def __init__(self, protocol_id, output):
        self.protocol_id = protocol_id
        self.output = output



def url():
    from . import settings
    url = 'mysql+mysqlconnector://{0.db_user}:{0.db_password}@{0.db_host}:{0.db_port}/{0.db_name}'
    return url.format(settings)

def connect():
    from contextlib import contextmanager

    @contextmanager     # (no fold)
    def session_manager():
        session = None
        try:
            engine = create_engine(url())
            Base.metadata.create_all(engine)
            Session.configure(bind=engine)
            session = Session()
            yield session
            session.commit()

        except:
            if session is not None: session.rollback()
            raise

        finally:
            if session is not None: session.close()

    return session_manager()


