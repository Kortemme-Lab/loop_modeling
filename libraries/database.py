# The MIT License (MIT)
#
# Copyright (c) 2015 Kale Kundert
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.

from . import install

install.require_sqlalchemy()
install.require_mysql_connector()

from sqlalchemy import create_engine, ForeignKey
from sqlalchemy import Column, Integer, Float, Text, String, DateTime, Boolean
from sqlalchemy.orm import sessionmaker, relationship, subqueryload
from sqlalchemy.ext.declarative import declared_attr, declarative_base
from sqlalchemy.dialects.mysql import LONGTEXT
from sqlalchemy.exc import ProgrammingError

Session = sessionmaker()

class Base (object):

    @declared_attr
    def __tablename__(cls):
        import re
        camel_case = re.compile('[A-Z][a-z]+')
        table_words = camel_case.findall(cls.__name__)
        return '_'.join(w.lower() for w in table_words)

    @property
    def dict(self):
        keys = [column.name for column in self.__table__.columns]
        return dict((key, getattr(self, key)) for key in keys)


Base = declarative_base(cls=Base)
NonRosettaBase = declarative_base(cls=Base)

class Benchmarks (NonRosettaBase):
    benchmark_id = Column(Integer, primary_key=True, autoincrement=True)
    protocol_ids = relationship('BenchmarkProtocols', order_by='BenchmarkProtocols.protocol_id')
    input_pdbs = relationship('BenchmarkInputs')
    start_time = Column(DateTime)
    name = Column(Text)
    title = Column(Text)
    user = Column(Text)
    description = Column(Text)
    rosetta_script = Column(Text)
    rosetta_script_vars = Column(Text)
    rosetta_flags = Column(Text)
    rosetta_fragments = Column(Text)
    git_commit = Column(Text)
    git_diff = Column(Text)
    fast = Column(Boolean)
    non_random = Column(Boolean)
    nstruct = Column(Integer)

    def __init__(self, name, script, nstruct,
            title=None, user=None, desc=None,
            vars=None, flags=None, fragments=None,
            git_commit=None, git_diff=None,
            fast=None, non_random=True):

        import datetime
        self.start_time = datetime.datetime.now()
        self.name = name
        self.title = title
        self.user = user
        self.description = desc
        self.rosetta_script = script
        self.rosetta_script_vars = vars
        self.rosetta_flags = flags
        self.rosetta_fragments = fragments
        self.git_commit = git_commit
        self.git_diff = git_diff
        self.fast = fast
        self.non_random = non_random
        self.nstruct = nstruct


    def __repr__(self):
        return '<Benchmark id={0.benchmark_id} name={0.name}>'.format(self)

    @property
    def id(self):
        return self.benchmark_id

    @property
    def protocols(self):
        session = Session.object_session(self)
        protocol_ids = [x.protocol_id for x in self.protocol_ids]
        return session.query(Protocols).\
                filter(Protocols.protocol_id.in_(protocol_ids)).\
                options(subqueryload('*')).\
                all() if protocol_ids else []

    @property
    def batches(self):
        session = Session.object_session(self)
        protocol_ids = [x.protocol_id for x in self.protocol_ids]
        return session.query(Batches).\
                join(Protocols).\
                filter(Protocols.protocol_id.in_(protocol_ids)).\
                options(subqueryload('*')).\
                all() if protocol_ids else []

    @property
    def structures(self):
        session = Session.object_session(self)
        protocol_ids = [x.protocol_id for x in self.protocol_ids]
        return session.query(Structures).\
                join(Batches, Protocols).\
                filter(Protocols.protocol_id.in_(protocol_ids)).\
                options(subqueryload('*')).\
                all() if protocol_ids else []


class BenchmarkProtocols (NonRosettaBase):
    benchmark_id = Column(Integer, ForeignKey('benchmarks.benchmark_id'), primary_key=True)
    protocol_id = Column(Integer, primary_key=True, autoincrement=False)
    
    def __init__(self, benchmark_id, protocol_id):
        self.benchmark_id = benchmark_id
        self.protocol_id = protocol_id

    def __repr__(self):
        return '<BenchmarkProtocol benchmark_id={0.benchmark_id} protocol_id={0.protocol_id}>'.format(self)


class BenchmarkInputs (NonRosettaBase):
    benchmark_input_id = Column(Integer, primary_key=True)
    benchmark_id = Column(Integer, ForeignKey('benchmarks.benchmark_id'))
    pdb_path = Column(Text)

    def __init__(self, pdb_path):
        self.pdb_path = pdb_path

    def __repr__(self):
        return '<BenchmarkInput benchmark_id={0.benchmark_id}, pdb_path={0.pdb_path}'.format(self)


class Protocols (Base):
    protocol_id = Column(Integer, primary_key=True)
    specified_options = Column(String)
    command_line = Column(String)
    svn_url = Column(String)
    svn_version = Column(String)
    script = Column(String)
    batches = relationship('Batches', order_by='Batches.batch_id', backref='protocol')

    def __repr__(self):
        return '<Protocol protocol_id={0.protocol_id}>'.format(self)


class Batches (Base):
    batch_id = Column(Integer, primary_key=True)
    protocol_id = Column(Integer, ForeignKey('protocols.protocol_id'))
    name = Column(String)
    description = Column(String)
    structures = relationship('Structures', order_by='Structures.struct_id', backref='protocol')

    def __repr__(self):
        return '<Batch batch_id={0.batch_id} protocol_id={0.protocol_id} name={0.name}>'.format(self)


class Structures (Base):
    struct_id = Column(Integer, primary_key=True)
    batch_id = Column(Integer, ForeignKey('batches.batch_id'))
    tag = Column(String)
    input_tag = Column(String)
    rmsd_features = relationship('ProteinRmsdNoSuperposition', uselist=False)
    score_features = relationship('TotalScores', uselist=False)
    runtime_features = relationship('Runtimes', uselist=False)

    def __repr__(self):
        return '<Structure struct_id={0.struct_id} batch_id={0.batch_id} tag={0.tag}>'.format(self)


class ProteinRmsdNoSuperposition (Base):
    struct_id = Column(Integer, ForeignKey('structures.struct_id'), primary_key=True)
    reference_tag = Column(String)
    protein_backbone = Column(Float)

    def __repr__(self):
        return '<ProteinRmsd struct_id={0.struct_id} reference_tag={0.reference_tag} backbone_rmsd={0.protein_backbone}>'.format(self)


class TotalScores (Base):
    struct_id = Column(Integer, ForeignKey('structures.struct_id'), primary_key=True)
    score = Column(Float)

    def __repr__(self):
        return '<TotalScore struct_id={0.struct_id} score={0.score}>'.format(self)


class Runtimes (Base):
    struct_id = Column(Integer, ForeignKey('structures.struct_id'), primary_key=True)
    timestamp = Column(String)
    elapsed_time = Column(Integer)

    def __repr__(self):
        return '<Runtime struct_id={0.struct_id} timestamp={0.timestamp} elapsed_time={0.elapsed_time}s>'.format(self)


class TracerLogs (NonRosettaBase):
    log_id = Column(Integer, primary_key=True)
    benchmark_id = Column(Integer)
    protocol_id = Column(Integer)   # This will be 0 if the job fails.
    stdout = Column(LONGTEXT)
    stderr = Column(LONGTEXT)

    def __init__(self, benchmark_id, protocol_id, stdout, stderr):
        self.benchmark_id = benchmark_id
        self.protocol_id = protocol_id
        self.stdout = stdout
        self.stderr = stderr

    def __repr__(self):
        return '<ProtocolOutput benchmark_id={0.benchmark_id} protocol_id={0.protocol_id}>'.format(self)



def url(db_name = None):
    from . import settings
    if db_name:
        url = 'mysql+mysqlconnector://{0.db_user}:{0.db_password}@{0.db_host}:{0.db_port}/' + db_name
    else:
        url = 'mysql+mysqlconnector://{0.db_user}:{0.db_password}@{0.db_host}:{0.db_port}/{0.db_name}'
    return url.format(settings)

def create_database():
    from . import settings
    try:
        engine = create_engine('mysql+mysqlconnector://{0.db_user}:{0.db_password}@{0.db_host}:{0.db_port}'.format(settings))
        engine.execute("CREATE DATABASE {0.db_name}".format(settings))
    except ProgrammingError, e:
        raise Exception('An error occurred creating the database: %s' % str(e))

def test_connect(db_name = None):
    try:
        with connect(db_name = db_name) as session: pass
    except ProgrammingError, e:
        if str(e).find('Unknown database') != -1 or str(e).find('does not exist'): # This second error text is what PostgreSQL will return.
            print('The database does not exist. Attempting to create it...')
            create_database()
            print('Database successfully created.')
        else:
            raise RuntimeError('An error occurred connecting to the database: %s' % str(e))

def connect(echo=False, db_name = None):
    from contextlib import contextmanager

    @contextmanager     # (no fold)
    def session_manager():
        session = None
        try:
            engine = create_engine(url(db_name = db_name), echo=echo)
            NonRosettaBase.metadata.create_all(engine)
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

