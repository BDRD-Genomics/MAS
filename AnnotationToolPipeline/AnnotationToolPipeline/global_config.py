import luigi
import os

class Globals(luigi.Config):
    '''
    Global variables. Set using luigi configuration file.
    '''
    # Directory to write output to
    OUTPUT_DIR = luigi.Parameter(
        default=os.path.join(os.path.split(os.path.dirname(os.path.realpath(__file__)))[0], 'output')
    )
    # AMD_DATABASE_SCRIPTS_PATH = luigi.Parameter()
    MAS_USERNAME = luigi.Parameter(default='luigi')
    MAS_PASSWORD = luigi.Parameter(default=os.environ['LUIGI_USER_PASSWORD'])
    MAS_CRT = luigi.Parameter(default=None)
    ERROR_LOG = luigi.Parameter()
    CLUSTER = luigi.Parameter(default=False)
    NUM_WORKERS = luigi.IntParameter(default=20)
    # CONDA_ENVIRONMENT = luigi.Parameter(default='mas-worker')

global_config = Globals()