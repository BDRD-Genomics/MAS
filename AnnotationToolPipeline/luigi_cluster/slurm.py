# -*- coding: utf-8 -*-


"""SLURM batch system Tasks.

Adapted by Matthew Lueder from
`luigi.contrib.sge <https://luigi.readthedocs.io/en/stable/api/luigi.contrib.sge.html>`_
by Jake Feala (@jfeala)

SLURM is a job scheduler used to allocate compute resources on a
shared cluster. Jobs are submitted using the ``sbatch`` command and monitored
using ``sacct``. To get started, install luigi on all nodes.

To run luigi workflows on an SLURM cluster, subclass
:class:`SLURMJobTask` as you would any :class:`luigi.Task`,
but override the ``work()`` method, instead of ``run()``, to define the job
code. Then, run your Luigi workflow from the master node, assigning > 1
``workers`` in order to distribute the tasks in parallel across the cluster.
"""

# This extension is modeled after the hadoop.py approach.
#
# Implementation notes
# The procedure:
# - Pickle the class
# - Construct a qsub argument that runs a generic runner function with the path to the pickled class
# - Runner function loads the class from pickle
# - Runner function hits the work button on it

import os
import subprocess
import time
import sys
import logging
import random

try:
    import cPickle as pickle
except ImportError:
    import pickle

import luigi
from luigi.contrib.hadoop import create_packages_archive
from . import slurm_runner

logger = logging.getLogger('luigi-interface')
logger.propagate = 0

POLL_TIME = 5  # decided to hard-code rather than configure here


def _get_sacct_job_status(job_id):
    result = subprocess.run(['sacct', '-P', '-b', '-j', str(job_id)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    for line in result.stdout.decode().splitlines():
        ls = line.split('|')
        if ls[0] == str(job_id):
            state = ls[1].split()[0].split('+')[0]
            exit_code = int(ls[2].split(':')[0])
            return state, exit_code
    raise RuntimeError('Could not find job status in sacct output. Is SLURM accounting properly configured?')


def _parse_sbatch_job_id(qsub_out):
    """Parse job id from sbatch output string.

    Assume format:

        "Submitted batch job <job_id>"

    """
    return int(qsub_out.strip().split()[-1])


def _build_sbatch_command(cmd, job_name, outfile, errfile, n_cpu, partition):
    """Submit shell command via `sbatch`
    """

    qsub_template = """sbatch --wrap "{cmd}" -o "{outfile}" -e "{errfile}" --export=All -p {partition} --cpus-per-task {n_cpu} -J {job_name}"""
    return qsub_template.format(
        cmd=cmd, job_name=job_name, outfile=outfile, errfile=errfile, n_cpu=n_cpu, partition=partition
    )


class SLURMJobTask(luigi.Task):
    """Base class for executing a job on SunGrid Engine

    Override ``work()`` (rather than ``run()``) with your job code.

    Parameters:

    - n_cpu: Number of CPUs (or "slots") to allocate for the Task.
    - shared_tmp_dir: Shared drive accessible from all nodes in the cluster.
          Task classes and dependencies are pickled to a temporary folder on
          this drive.
    - job_name_format: String that can be passed in to customize the job name
        string passed to sbatch; e.g. "Task123_{task_family}_{n_cpu}...".
    - job_name: Exact job name to pass to sbatch.
    - run_locally: Run locally instead of on the cluster.
    - poll_time: the length of time to wait in order to poll status
    - dont_remove_tmp_dir: Instead of deleting the temporary directory, keep it.
    - no_tarball: Don't create a tarball of the luigi project directory.  Can be
        useful to reduce I/O requirements when the luigi directory is accessible
        from cluster nodes already.
    - partition: Partition to use on the SLURM system: default is normal
    """

    n_cpu = luigi.IntParameter(default=2, significant=False)
    shared_tmp_dir = luigi.Parameter(default='/tmp', significant=False)
    job_name_format = luigi.Parameter(
        significant=False, default='', description="A string that can be "
                                                     "formatted with class variables to name the job.")
    job_name = luigi.Parameter(
        significant=False, default='',
        description="Explicit job name given via sbatch.")
    run_locally = luigi.BoolParameter(
        significant=False,
        description="run locally instead of on the cluster")
    poll_time = luigi.IntParameter(
        significant=False, default=POLL_TIME,
        description="specify the wait time to poll for the job status")
    dont_remove_tmp_dir = luigi.BoolParameter(
        significant=False,
        description="don't delete the temporary directory used (for debugging)")
    no_tarball = luigi.BoolParameter(
        significant=False,
        description="don't tarball (and extract) the luigi project files")
    slurm_partition = luigi.Parameter(
        significant=False,
        description="Name of the partition to use on the SLURM",
        default='normal'
    )

    def __init__(self, *args, **kwargs):
        super(SLURMJobTask, self).__init__(*args, **kwargs)
        if self.job_name:
            # use explicitly provided job name
            pass
        elif self.job_name_format:
            # define the job name with the provided format
            self.job_name = self.job_name_format.format(
                task_family=self.task_family, **self.__dict__)
        else:
            # default to the task family
            self.job_name = self.task_family

    def _fetch_task_failures(self):
        if not os.path.exists(self.errfile):
            logger.info('No error file')
            return []
        with open(self.errfile, "r") as f:
            errors = f.readlines()
        if errors == []:
            return errors
        return errors

    def _init_local(self):

        # Set up temp folder in shared directory (trim to max filename length)
        base_tmp_dir = self.shared_tmp_dir
        random_id = '%016x' % random.getrandbits(64)
        folder_name = self.task_id + '-' + random_id
        self.tmp_dir = os.path.join(base_tmp_dir, folder_name)
        max_filename_length = os.fstatvfs(0).f_namemax
        self.tmp_dir = self.tmp_dir[:max_filename_length]
        logger.info("Tmp dir: %s", self.tmp_dir)
        os.makedirs(self.tmp_dir)

        # Dump the code to be run into a pickle file
        logging.debug("Dumping pickled class")
        self._dump(self.tmp_dir)

        if not self.no_tarball:
            # Make sure that all the class's dependencies are tarred and available
            # This is not necessary if luigi is importable from the cluster node
            logging.debug("Tarballing dependencies")
            # Grab luigi and the module containing the code to be run
            packages = [luigi] + [__import__(self.__module__, None, None, 'dummy')]
            create_packages_archive(packages, os.path.join(self.tmp_dir, "packages.tar"))

    def run(self):
        # if self.complete():
        #     logging.error('Luigi scheduler trying to re-run complete task. Luigi, why did you make this decision?')
        #     return

        if self.run_locally:
            self.work()
        else:
            self._init_local()
            self._run_job()
            # The procedure:
            # - Pickle the class
            # - Tarball the dependencies
            # - Construct a qsub argument that runs a generic runner function with the path to the pickled class
            # - Runner function loads the class from pickle
            # - Runner class untars the dependencies
            # - Runner function hits the button on the class's work() method

    def work(self):
        """Override this method, rather than ``run()``,  for your actual work."""
        pass

    def _dump(self, out_dir=''):
        """Dump instance to file."""
        with self.no_unpicklable_properties():
            self.job_file = os.path.join(out_dir, 'job-instance.pickle')
            if self.__module__ == '__main__':
                d = pickle.dumps(self)
                module_name = os.path.basename(sys.argv[0]).rsplit('.', 1)[0]
                d = d.replace('(c__main__', "(c" + module_name)
                open(self.job_file, "w").write(d)
            else:
                pickle.dump(self, open(self.job_file, "wb"))

    def _run_job(self):

        runner_path = slurm_runner.__file__
        if runner_path.endswith("pyc"):
            runner_path = runner_path[:-3] + "py"

        job_str = 'python {0} {1} {2}'.format(
            runner_path, self.tmp_dir.replace(' ', '\ '), os.getcwd().replace(' ', '\ '))
        if self.no_tarball:
            job_str += ' --no-tarball'

        # Build qsub submit command
        self.outfile = os.path.join(self.tmp_dir, 'job.out')
        self.errfile = os.path.join(self.tmp_dir, 'job.err')
        submit_cmd = _build_sbatch_command(job_str, self.job_name, self.outfile, self.errfile, self.n_cpu, self.slurm_partition)
        logger.info('sbatch command: \n' + submit_cmd)

        # Submit the job and grab job ID
        output = subprocess.check_output(submit_cmd, shell=True)
        output = output.decode('utf-8')
        self.job_id = _parse_sbatch_job_id(output)

        logger.info("Submitted job to SLURM with response:\n" + output)

        self._track_job()

        # Now delete the temporaries, if they're there.
        if (self.tmp_dir and os.path.exists(self.tmp_dir) and not self.dont_remove_tmp_dir):
            logger.info('Removing temporary directory %s' % self.tmp_dir)
            subprocess.call(["rm", "-rf", self.tmp_dir])

    def _track_job(self):
        while True:
            # Sleep for a little bit
            time.sleep(self.poll_time)

            state, exit_code = _get_sacct_job_status(self.job_id)

            if state in ['BOOT_FAIL', 'CANCELLED', 'DEADLINE', 'NODE_FAIL', 'OUT_OF_MEMORY', 'TIMEOUT', 'FAILED']:
                raise RuntimeError(
                    'SLURM job for task {} failed with the status {}.'.format(self.task_id, state)
                )

            if state == 'COMPLETED':
                if exit_code == 0:
                    logger.debug('SLURM job for task {} completed without error.'.format(self.task_id))
                    time.sleep(2)
                    break
                elif exit_code > 0:
                    raise RuntimeError('SLURM job for task {} completed with error code {}.'.format(self.task_id, exit_code))


class LocalSLURMJobTask(SLURMJobTask):
    """A local version of SLURMJobTask, for easier debugging.

    This version skips the ``qsub`` steps and simply runs ``work()``
    on the local node, so you don't need to be on an SLURM cluster to
    use your Task in a test workflow.
    """

    def run(self):
        self.work()
