import subprocess
import os
from datetime import datetime
import logging
import platform
import types

import luigi
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import requests
from django.urls import reverse

from AnnotationToolPipeline.luigi_cluster.sge import SGEJobTask
from AnnotationToolPipeline.luigi_cluster.slurm import SLURMJobTask
from AnnotationToolPipeline.AnnotationToolPipeline.global_config import global_config

# Setup Django to allow us to use reverse()
import django
django.setup()

luigi_logger = logging.getLogger('luigi-interface')

if global_config.CLUSTER == 'SGE':
    BaseTask = SGEJobTask
elif global_config.CLUSTER == 'SLURM':
    BaseTask = SLURMJobTask
else:
    BaseTask = luigi.Task


class ClusterTaskParameters(luigi.Config):
    """
    A config class to set parameters for the SGEJobTask
    The SGEJobTask is not reading in parameters from luigi.cfg
    """
    shared_tmp_dir = luigi.Parameter(default=global_config.OUTPUT_DIR)
    no_tarball = luigi.BoolParameter(default=True)
    default_sge_queue = luigi.Parameter(default='all.q')
    default_slurm_partition = luigi.Parameter(default='normal')
    dont_remove_tmp_dir = luigi.BoolParameter(default=True)
    poll_time = luigi.IntParameter(default=5)
    sge_parallel_env = luigi.Parameter(default='smp')


class PipelineTask(BaseTask):
    '''
        This is a base class for other classes in the pipeline to inherit from. It helps set the structure for a task.
        '''
    g = global_config
    mas_server = luigi.Parameter()
    n_cpu = luigi.IntParameter(default=1)
    shared_tmp_dir = ClusterTaskParameters().shared_tmp_dir
    no_tarball = ClusterTaskParameters().no_tarball
    sge_queue = ClusterTaskParameters().default_sge_queue
    dont_remove_tmp_dir = ClusterTaskParameters().dont_remove_tmp_dir
    poll_time = ClusterTaskParameters().poll_time
    sge_parallel_env = ClusterTaskParameters().sge_parallel_env
    slurm_partition = ClusterTaskParameters().default_slurm_partition

    def __init__(self, *args, **kwargs):
        super(PipelineTask, self).__init__(*args, **kwargs)

        name = ''
        if hasattr(self, 'phage_name') and self.genome_name:
            name = self.genome_name
        elif self.annotation_accession:
            name = self.annotation_accession

        self.job_name = "%s__%s" % (self.task_family, name)

    def pipeline_out_dir(self):
        '''
        Returns the directory which will hold the files output by this task.
        '''
        if hasattr(self, 'phage_name') and self.genome_name:
            return os.path.join(self.g.OUTPUT_DIR, self.genome_name)

        elif self.annotation_accession:
            return os.path.join(self.g.OUTPUT_DIR, self.run_time.strftime("%Y-%m-%d_h%H-m%M-s%S"),
                                self.annotation_accession)

    def out_dir(self, temp=False):
        folder = '{}-temp'.format(self.task_family) if temp else self.task_family
        return os.path.join(self.pipeline_out_dir(), folder)

    def out_file_path(self, temp=False):
        '''
        Returns location of output files. Must return a dict (keys identify file, values are full paths to file).
        '''
        raise NotImplementedError

    def output(self):
        output_files = self.out_file_path(temp=False)
        target_dict = {}

        for key in output_files.keys():
            target_dict[key] = luigi.LocalTarget(output_files[key])

        return target_dict

    def complete(self):
        if os.path.isdir(self.out_dir()):
            return True
        return False

    def do_task(self):
        '''
        Implement task logic here!
        '''
        raise NotImplementedError

    def additional_depends(self):
        pass

    def run(self):
        # Take care of any extra dependencies first
        extra_deps = self.additional_depends()
        if extra_deps:
            self.additional_targets = yield extra_deps

        if self.g.CLUSTER == 'SLURM' or self.g.CLUSTER == 'SGE':
            if self.run_locally:
                self.work()
            else:
                self._init_local()
                self._run_job()
        else:
            self.work()

    def work(self):
        # Generate task id
        luigi_logger.critical('Starting work() of Task {}'.format(self.task_id))
        param_map = {}
        for param in self.param_kwargs.keys():
            param_map[param] = str(self.param_kwargs[param])

        self.task_uid = '%s_%s' % (
            luigi.task.task_id_str(self.get_task_family(), param_map),
            datetime.now().strftime('%Y_%m_%d_%H_%M_%S_%f')
        )

        # Create output folder
        self._create_output_folder(self.out_dir(temp=True))

        # Set up logger
        self.logger = logging.getLogger(self.task_uid)
        self.logger.setLevel(logging.DEBUG)
        self.fh = logging.FileHandler(os.path.join(self.out_dir(temp=True), '%s.log' % self.task_uid))
        self.fh.setFormatter(logging.Formatter('\n%(asctime)s level=%(levelname)s:\n%(message)s'))
        self.logger.addHandler(self.fh)

        # Show task info in log
        init_info_message = '##### Running %s on %s #####\nParameters:' % (type(self).__name__, platform.node())
        for key in param_map.keys():
            init_info_message = '%s\n\t%s = %s' % (init_info_message, key, param_map[key])
        self.logger.info(init_info_message)

        try:
            dep_gen = self.do_task()
            if isinstance(dep_gen, types.GeneratorType):
                self.logger.info('Dynamic dependencies encountered. Attempting to resolve.')
                # for dep in dep_gen:
                #     yield dep
                return dep_gen

        except Exception:
            self.logger.exception('Task %s failed!' % type(self).__name__)
            raise

        assert self._check_output()
        self.logger.info('Task completed without errors!')

    def _check_output(self):
        '''
        Checks the temporary dir to make sure all required files are there and not empty
        :param task: A reference to the task which is calling this function
        :return: True if successful, False otherwise
        '''
        self.logger.info('Checking temp folder for task %s' % self.task_uid)
        for key in self.out_file_path(temp=True).keys():
            file_path = self.out_file_path(temp=True)[key]
            if not os.path.isfile(file_path):
                self.logger.error('OUTPUT FILE MISSING: file %s was not created' % file_path)
                return False
            else:
                if os.path.getsize(file_path) == 0:
                    self.logger.error('OUTPUT FILE EMPTY: file %s is empty' % file_path)
                    return False

        # All files exist and are not empty, so now we will move the folder
        os.rename(self.out_dir(temp=True), self.out_dir(temp=False))
        return True

    def _run_command(self, command_params, condaenv='', **kwargs):
        '''
        Wrapper to run commands using subprocess given a string
        :param command_params: string of commands to be run
        :param stdout: file path to stdout
        :param stderr: file path to stderr
        :return: Nothing
        '''
        # if condaenv:  # Prepend source activate condaenv && to command
        #     command_params.insert(0, '&&')
        #     command_params.insert(0, condaenv)
        #     command_params.insert(0, 'activate')
        #     command_params.insert(0, 'conda')
        # else:
        #     command_params.insert(0, '&&')
        #     command_params.insert(0, self.g.CONDA_ENVIRONMENT)
        #     command_params.insert(0, 'activate')
        #     command_params.insert(0, 'conda')

        # Ensure all params are strings
        command_params = [str(x) for x in command_params]

        # Prepare command_string
        command_string = ' '.join(command_params)

        self.logger.info('Running command:\n%s' % command_string.strip())

        # Execute command
        cmd_result = subprocess.run(command_string, shell=True, stderr=subprocess.PIPE, **kwargs)

        if cmd_result.stderr:
            self.logger.warning('Command produced output to stderr:\n%s' % cmd_result.stderr.decode('unicode_escape'))

        if cmd_result.returncode != 0:
            self.logger.error('Command exited with error code %i' % cmd_result.returncode)

        return cmd_result

    def _create_output_folder(self, path):
        head, tail = os.path.split(path)
        if head and not os.path.isdir(head):
            self._create_output_folder(head)

        if not os.path.isdir(path):
            os.mkdir(path)


@PipelineTask.event_handler(luigi.Event.FAILURE)
def on_task_fail(task, exception):
    '''
    Handle task failures by signalling back to MAS
    '''
    # Log in error log
    logger = logging.getLogger('error-log')
    logger.setLevel(logging.ERROR)
    fh = logging.FileHandler(task.g.ERROR_LOG)
    fh.setFormatter(logging.Formatter('\n%(asctime)s level=%(levelname)s:\n%(message)s'))
    logger.addHandler(fh)
    logger.error(str(exception))

    # Send error status to MAS
    r = requests.post(
        task.mas_server + reverse('upload_results'),
        auth=(task.g.MAS_USERNAME, task.g.MAS_PASSWORD),
        data={
            'tool': task.tool,
            'accession': task.annotation_accession,
            'database': task.database,
            'status': 2
        },
        verify=task.g.MAS_CRT
    )

    if r.status_code != 200:
        print(r.text)
        raise requests.ConnectionError('Request to post results to MAS server failed')


'''
**********   Utility Tasks   **********
'''
class Pull_Protein(PipelineTask):
    '''
    Pull a single protein through MAS's REST API
    '''
    annotation_accession = luigi.Parameter()
    run_locally = True
    database = luigi.Parameter()
    run_time = luigi.DateSecondParameter()
    tool = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super(Pull_Protein, self).__init__(*args, **kwargs)
        self.task_id = '{}_{}_{}'.format(self.get_task_family(), self.run_time, self.annotation_accession)

    def out_file_path(self, temp=False):
        return {
            'fasta': os.path.join(self.out_dir(temp), '%s.faa' % self.annotation_accession)
        }

    def do_task(self):
        # Get protein sequence from MAS's REST API
        response = requests.get(
            self.mas_server + reverse('get_protein', kwargs={'accession': self.annotation_accession}),
            auth=(self.g.MAS_USERNAME, self.g.MAS_PASSWORD),
            verify=self.g.MAS_CRT
        )

        if response.status_code != 200:
            self.logger.error(
                'Response status code = {}. Response text = {}.'.format(response.status_code, response.text)
            )
            raise requests.ConnectionError('Request to get protein sequence from MAS server failed')

        # Write sequence file
        protein_seq = Seq(response.json()['sequence'], IUPAC.IUPACProtein)
        rec = SeqRecord(protein_seq, id=self.annotation_accession, description='')
        SeqIO.write(rec, self.out_file_path(True)['fasta'], 'fasta')


'''
**********   Tools   **********
'''
class Blastp(PipelineTask):
    '''
    Generate blastp results for a given proteome
    '''
    annotation_accession = luigi.Parameter(default='')
    e_value = luigi.FloatParameter(default=0.01)
    database = luigi.Parameter()
    run_time = luigi.DateSecondParameter()
    tool = 'blastp'

    # Available Database Choices
    swissprot = luigi.Parameter()
    nr = luigi.Parameter()
    internal = luigi.Parameter()

    # specific # CPU to use for each job
    swissprot_cpu = luigi.Parameter()
    nr_cpu = luigi.Parameter()
    internal_cpu = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super(Blastp, self).__init__(*args, **kwargs)

        if self.database == 'nr':
            self.n_cpu = self.nr_cpu
        elif self.database == 'swissprot':
            self.n_cpu = self.swissprot_cpu
        elif self.database == 'internal':
            self.n_cpu = self.internal_cpu

    def requires(self):
        return Pull_Protein(
            annotation_accession=self.annotation_accession,
            run_time=self.run_time,
            tool=self.tool,
            database=self.database,
            mas_server=self.mas_server
        )

    def out_file_path(self, temp=False):
        name = self.annotation_accession

        return {
            'results': os.path.join(self.out_dir(temp), '{}_{}_blastp_results.xml'.format(name, self.database))
        }

    def out_dir(self, temp=False):
        folder = '{}-temp'.format(self.database) if temp else self.database
        return os.path.join(self.pipeline_out_dir(), self.task_family, folder)

    def do_task(self):
        if self.database == 'nr':
            db_path = self.nr

        elif self.database == 'swissprot':
            db_path = self.swissprot

        elif self.database == 'internal':
            db_path = self.internal

        else:
            raise ValueError('Invalid database ' + self.database)

        self._run_command([
            'blastp',
            '-query', self.input()['fasta'].path,
            '-db', db_path,
            '-evalue', str(self.e_value),
            '-outfmt', '5',
            '-out', self.out_file_path(True)['results'],
            '-num_threads', str(self.n_cpu)
        ])


class RPSBlast(PipelineTask):
    '''
    Generate rpsblast (CD-Search) results for a given proteome
    '''
    annotation_accession = luigi.Parameter(default='')
    e_value = luigi.FloatParameter(default=0.0001)
    database = luigi.Parameter()
    run_time = luigi.DateSecondParameter()
    tool = 'rpsblast'

    # Available Database Choices
    cdd = luigi.Parameter()
    cdd_cpu = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super(RPSBlast, self).__init__(*args, **kwargs)

        if self.database == 'cdd':
            self.n_cpu = self.cdd_cpu

    def requires(self):
        return Pull_Protein(
            annotation_accession=self.annotation_accession,
            run_time=self.run_time,
            tool=self.tool,
            database=self.database,
            mas_server=self.mas_server
        )

    def out_file_path(self, temp=False):
        name = self.annotation_accession
        return {
            'results': os.path.join(self.out_dir(temp), '{}_{}_rpsblast_results.xml'.format(name, self.database))
        }

    def out_dir(self, temp=False):
        folder = '{}-temp'.format(self.database) if temp else self.database
        return os.path.join(self.pipeline_out_dir(), self.task_family, folder)

    def do_task(self):
        if self.database == 'cdd':
            db_path = self.cdd

        else:
            raise ValueError('Invalid database ' + self.database)

        self._run_command([
            'rpsblast',
            '-query', self.input()['fasta'].path,
            '-db', db_path,
            '-evalue', self.e_value,
            '-out', self.out_file_path(True)['results'],
            '-outfmt', '5',
            '-num_threads', str(self.n_cpu)
        ])


class HHblits(PipelineTask):
    '''
    Run HHblits on each protein
    '''
    annotation_accession = luigi.Parameter(default='')
    protein_id = luigi.Parameter(default='')
    iterations = luigi.IntParameter(default=3)
    database = luigi.Parameter()
    run_time = luigi.DateSecondParameter()
    tool = 'hhsearch'

    # At this point this is the only database option. The database parameter is for HHSearch!
    uniclust = luigi.Parameter()
    uniclust_cpu = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super(HHblits, self).__init__(*args, **kwargs)

        if self.uniclust_cpu:
            self.n_cpu = self.uniclust_cpu

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task.
        '''
        base_name = '{}_{}'.format(self.task_family, self.annotation_accession)

        folder = '{}-temp'.format(base_name) if temp else base_name
        return os.path.join(self.pipeline_out_dir(), self.task_family, folder)

    def requires(self):
        return Pull_Protein(
            annotation_accession=self.annotation_accession,
            run_time=self.run_time,
            tool=self.tool,
            database=self.database,
            mas_server=self.mas_server
        )

    def out_file_path(self, temp=False):
        name = self.annotation_accession
        return {
            'alignment': os.path.join(self.out_dir(temp), '{}.a3m'.format(name))
        }

    def do_task(self):
        db_path = self.uniclust
        in_file = self.input()['fasta'].path

        self._run_command([
            'hhblits',
            '-i', in_file,
            '-oa3m', self.out_file_path(True)['alignment'],
            '-n', str(self.iterations),
            '-cpu', str(self.n_cpu),
            '-d', db_path
        ])


class HHsearch(PipelineTask):
    protein_id = luigi.Parameter(default='')
    annotation_accession = luigi.Parameter(default='')
    run_time = luigi.DateSecondParameter()
    tool = 'hhsearch'
    database = luigi.Parameter()

    # Available Database Choices
    pdb = luigi.Parameter()

    pdb_cpu = luigi.Parameter()

    def __init__(self, *args, **kwargs):
        super(HHsearch, self).__init__(*args, **kwargs)

        if self.pdb_cpu:
            self.n_cpu = self.pdb_cpu

    def out_dir(self, temp=False):
        '''
        Returns the directory which will hold the files output by this task.
        '''
        base_name = '{}_{}'.format(self.task_family, self.annotation_accession)

        folder = '{}-temp'.format(base_name) if temp else base_name
        return os.path.join(self.pipeline_out_dir(), self.task_family, folder)

    def requires(self):
       return HHblits(
            annotation_accession=self.annotation_accession,
            run_time=self.run_time,
            database=self.database,
            mas_server=self.mas_server
        )

    def out_file_path(self, temp=False):
        name = self.annotation_accession

        return {
            'results': os.path.join(self.out_dir(temp), '{}.hhr'.format(name))
        }

    def do_task(self):
        if self.database == 'pdb':
            db_path = self.pdb
        else:
            raise ValueError('Invalid database ' + self.database)

        self._run_command([
            'hhsearch',
            '-i', self.input()['alignment'].path,
            '-d', db_path,
            '-o', self.out_file_path(True)['results'],
            '-cpu', str(self.n_cpu)
        ])


'''
**********   Result Submission   **********
'''
class Send_Results_To_MAS(PipelineTask):
    annotation_accession = luigi.Parameter(default='')
    database = luigi.Parameter()
    tool = luigi.Parameter()
    run_time = luigi.DateSecondParameter(default=datetime.now())
    run_locally = True

    def requires(self):
        if self.tool == 'blastp':
            return Blastp(
                annotation_accession=self.annotation_accession,
                database=self.database,
                run_time=self.run_time,
                mas_server=self.mas_server
            )

        elif self.tool == 'hhsearch':
            return HHsearch(
                annotation_accession=self.annotation_accession,
                database=self.database,
                run_time=self.run_time,
                mas_server=self.mas_server
            )

        elif self.tool == 'rpsblast':
            return RPSBlast(
                annotation_accession=self.annotation_accession,
                database=self.database,
                run_time=self.run_time,
                mas_server=self.mas_server
            )

    def out_file_path(self, temp=False):
        return {}

    def out_dir(self, temp=False):
        folder = '{}-temp'.format(self.database) if temp else self.database
        return os.path.join(self.pipeline_out_dir(), self.task_family, self.tool, folder)

    def do_task(self):
        r = requests.post(
            self.mas_server + reverse('upload_results'),
            auth=(self.g.MAS_USERNAME, self.g.MAS_PASSWORD),
            files=[('result', open(self.input()['results'].path))],
            data={
                'tool': self.tool,
                'accession': self.annotation_accession,
                'database': self.database,
                'status': 0
            },
            verify=self.g.MAS_CRT
        )

        if r.status_code != 200:
            raise requests.ConnectionError(
                'Request to post results to MAS server failed (status code = %i):\n%s' % (r.status_code, r.text)
            )


# @Send_Results_To_MAS.event_handler(luigi.Event.SUCCESS)
# def cleanup(task):
#     '''
#     After results successfully uploaded to MAS, we can delete output
#     '''
#     task.logger.info('Pipeline complete: Cleaning up.')
#
#     # Remove entire pipeline directory
#     shutil.rmtree(task.pipeline_out_dir())


class Run_Pipeline_For_Proteins(luigi.WrapperTask):
    input_list = luigi.ListParameter()
    run_time = luigi.DateSecondParameter(default=datetime.now())
    mas_server = luigi.Parameter()

    def requires(self):
        job_array = []
        for search_params in self.input_list:
            if set(search_params) != {'accession', 'tool', 'database'}:
                raise KeyError('Incorrect dict keys for pipeline parameters')

            job_array.append(
                Send_Results_To_MAS(
                    mas_server=self.mas_server,
                    annotation_accession=search_params['accession'],
                    database=search_params['database'],
                    tool=search_params['tool']
                )
            )
        return job_array


