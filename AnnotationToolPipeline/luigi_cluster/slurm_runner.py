# -*- coding: utf-8 -*-

"""
Adapted from luigi.contrib.sge

The SLURM runner

The main() function of this module will be executed on the
compute node by the submitted job. It accepts as a single
argument the shared temp folder containing the package archive
and pickled task to run, and carries out these steps:

- extract tarfile of package dependencies and place on the path
- unpickle SLURMTask instance created on the master node
- run SLURMTask.work()

On completion, SLURMTask on the master node will detect that
the job has left the queue, delete the temporary folder, and
return from SLURMTask.run()

Modified to be Python 3 compatible by Logan Voegtly and Matthew Lueder
"""

import os
import sys
try:
    import cPickle as pickle
except ImportError:
    import pickle
import logging
import tarfile


def _do_work_on_compute_node(work_dir, tarball=True):

    if tarball:
        # Extract the necessary dependencies
        # This can create a lot of I/O overhead when running many SGEJobTasks,
        # so is optional if the luigi project is accessible from the cluster node
        _extract_packages_archive(work_dir)

    # Open up the pickle file with the work to be done
    os.chdir(work_dir)
    with open("job-instance.pickle", "rb") as f:
        job = pickle.load(f)

    # Do the work contained
    job.work()


def _extract_packages_archive(work_dir):
    package_file = os.path.join(work_dir, "packages.tar")
    if not os.path.exists(package_file):
        return

    curdir = os.path.abspath(os.curdir)

    os.chdir(work_dir)
    tar = tarfile.open(package_file)
    for tarinfo in tar:
        tar.extract(tarinfo)
    tar.close()
    if '' not in sys.path:
        sys.path.insert(0, '')

    os.chdir(curdir)


def main(args=sys.argv):
    """Run the work() method from the class instance in the file "job-instance.pickle".
    """
    try:
        tarball = "--no-tarball" not in args
        # Set up logging.
        logging.basicConfig(level=logging.WARN)
        work_dir = args[1]
        assert os.path.exists(work_dir), "First argument to slurm_runner.py must be a directory that exists"
        project_dir = args[2]
        sys.path.append(project_dir)
        _do_work_on_compute_node(work_dir, tarball)
    except Exception as e:
        # Dump encoded data that we will try to fetch using mechanize
        print(e)
        raise


if __name__ == '__main__':
    main()
