from __future__ import print_function
import os
import re


top_x = 5

class Benchmark:

    @staticmethod
    def from_names(names, group_by_name = False):
        benchmarks = []

        # The list of names may include several different kinds of name:
        #
        # 1. Benchmark id (primary key into database defined in settings file)
        # 2. Benchmark name (most recent benchmark run with given name in database)
        # 3. Path to benchmark data (stored in tabular flat file format)

        for name in names:
            if name.endswith('.results'):
                benchmark = Benchmark.from_flat_file(name)
            else:
                benchmark = Benchmark.from_database(name, group_by_name = group_by_name)

            if benchmark:
                benchmarks.append(benchmark)
            else:
                message = "Skipping the empty {0} benchmark..."
                utilities.print_warning(message, benchmark.title)

        return benchmarks

    @staticmethod
    def from_flat_file(path, name=None):
        from os.path import basename, splitext

        name = name or splitext(basename(path))[0]
        benchmark = Benchmark(name)

        print("Loading the {0.title} benchmark from a flat file...".format(benchmark))

        with open(path) as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('#'): continue

                tag, id, rmsd, score, runtime = line.split()
                id, rmsd, score, runtime = \
                        int(id), float(rmsd), float(score), int(runtime)

                loop = benchmark.loops.setdefault(tag, Loop(benchmark, tag))
                model = Model(loop, id, score, rmsd, runtime)
                loop.models.append(model)

        return benchmark

    @staticmethod
    def from_database(name_or_id, group_by_name = False):
        from libraries import database
        from sqlalchemy import desc

        with database.connect() as session:

            # Decide whether a name or id was used to specify a benchmark run, 
            # and load the corresponding data out of the database.  The meaning 
            # of name_or_id is inferred from its type: names are expected to be 
            # strings and ids are expected to be integers.  If more than one
            # benchmark has the same name, the most recent one will be used.

            db_benchmarks = []
            try:
                id = int(name_or_id)
                _db_benchmark = session.query(database.Benchmarks).get(id)

                if _db_benchmark is None:
                    message = "No benchmark '{}' in the database."
                    utilities.print_error_and_die(message, id)

                db_benchmarks = [_db_benchmark]

            except ValueError:
                name = name_or_id
                query = session.query(database.Benchmarks).filter_by(name=name).order_by(desc(database.Benchmarks.start_time))
                db_benchmarks = [q for q in query]
                if not group_by_name:
                    if len(db_benchmarks) > 1:
                        message = "Multiple benchmarks runs were found with the same name '{0}' (ids are: {1}). If this is expected then set the --group_by_name option.".format(name, ', '.join(map(str, [b.id for b in db_benchmarks])))
                        utilities.print_error_and_die(message, name)

            b_name = set([db_benchmark.name for db_benchmark in db_benchmarks])
            assert(len(b_name) == 1)
            b_name = b_name.pop()

            b_title = set([db_benchmark.title or '' for db_benchmark in db_benchmarks])
            if len(b_title) > 1:
                colortext.warning("There are multiple titles associated with benchmark {0}: '{1}'. Choosing the most recent ('{2}').".format(b_name, "', '".join(b_title), db_benchmarks[0].title or ''))
                b_title = db_benchmarks[0].title
            else:
                b_title = b_title.pop() or None

            benchmark = Benchmark(b_name, b_title)

            for db_benchmark in db_benchmarks:

                # Fill in the benchmark data structure from the database.

                print("Loading the {0} benchmark (id {1}) from the database...".format(benchmark.name, db_benchmark.id))

                for db_input in db_benchmark.input_pdbs:
                    path = db_input.pdb_path
                    if not benchmark.loops.get(path):
                        benchmark.loops[path] = Loop(benchmark, path)

                for structure in db_benchmark.structures:
                    loop = benchmark.loops[structure.input_tag]
                    id = len(loop.models) + 1
                    score = structure.score_features.score
                    rmsd = structure.rmsd_features.protein_backbone
                    runtime = structure.runtime_features.elapsed_time

                    model = Model(loop, id, score, rmsd, runtime)
                    loop.models.append(model)

        return benchmark


    def __init__(self, name, title=None):
        self.name = name
        self.manual_title = title
        self.loops = {}         # Set by Report.from_...()
        self.latex_dir = None   # Set by Report.setup_latex_dir()
        self.color = None       # Set by Report.setup_benchmark_colors()

    def __str__(self):
        return '<Benchmark name={0.name}>'.format(self)

    def __eq__(self, other):
        return self.name == other.name

    def __hash__(self):
        return hash(self.name)

    def __iter__(self):
        try:
            return self.loops.itervalues()
        except AttributeError:
            return iter(self.loops.values())

    def __len__(self):
        return len(self.loops)

    def __nonzero__(self):
        return any(self)

    def __getitem__(self, tag):
        return self.loops[tag]

    def __setitem__(self, tag, loop):
        self.loops[tag] = loop

    @property
    def title(self):
        if self.manual_title:
            return self.manual_title

        words = self.name.split('_')
        phrase = ' '.join(words)
        phrase = phrase.title()
        phrase = re.sub(r'\b[Kk][Ii][Cc]\b', 'KIC', phrase)
        phrase = re.sub(r'\b[Cc][Cc][Dd]\b', 'CCD', phrase)
        phrase = re.sub(r'\b[Nn][Gg][Kk]\b', 'NGK', phrase)
        return phrase

    @property
    def all_models(self):
        return sum((loop.models for loop in self if loop.has_data), [])

    @property
    def all_runtimes(self):
        return sum((loop.runtimes for loop in self if loop.has_data), [])

    @property
    def best_top_x_models(self):
        return [loop.best_top_x_model for loop in self if loop.has_data]

    @property
    def lowest_score_models(self):
        return [loop.lowest_score_model for loop in self if loop.has_data]

    @property
    def lowest_rmsd_models(self):
        return [loop.lowest_rmsd_model for loop in self if loop.has_data]

    @property
    def percents_subangstrom(self):
        return [loop.percent_subangstrom for loop in self if loop.has_data]


class Loop:

    def __init__(self, benchmark, path):
        self.benchmark = benchmark
        self.path = path
        self.models = []        # Set by Report.from_...()
        self.latex_dir = None   # Set by Report.setup_latex_dir()

    def __iter__(self):
        return iter(self.models)

    def __len__(self):
        return len(self.models)

    def __nonzero__(self):
        return bool(self.models)

    @property
    def pdb_id(self):
        return os.path.basename(self.path)[0:4]

    @property
    def num_models(self):
        return len(self)

    @property
    def has_data(self):
        return len(self) > 0

    @property
    def scores(self):
        return [x.score for x in self.models]

    @property
    def rmsds(self):
        return [x.rmsd for x in self.models]

    @property
    def runtimes(self):
        return [x.runtime for x in self.models]

    @property
    def models_sorted_by_score(self):
        return sorted(self.models, key=lambda x: x.score)

    @property
    def models_sorted_by_rmsd(self):
        return sorted(self.models, key=lambda x: x.rmsd)

    @property
    def best_top_x_model(self):
        top_x_models = self.models_sorted_by_score[:top_x]
        return sorted(top_x_models, key=lambda x: x.rmsd)[0]

    @property
    def lowest_score_model(self):
        return self.models_sorted_by_score[0]

    @property
    def lowest_rmsd_model(self):
        return self.models_sorted_by_rmsd[0]

    @property
    def percent_subangstrom(self):
        num_sub_a = sum(1 for x in self.models if x.rmsd < 1.0)
        return 100.0 * num_sub_a / len(self.models)


class Model:

    def __init__(self, loop, id, score, rmsd, runtime):
        self.loop = loop
        self.id = id
        self.score = score
        self.rmsd = rmsd
        self.runtime = runtime

