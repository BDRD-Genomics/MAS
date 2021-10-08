import os

from celery import Celery


# set the default Django settings module for the 'celery' program.
os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'MAS.settings')

# Set the path to luigi's configuration file
from django.conf import settings
os.environ['LUIGI_CONFIG_PATH'] = settings.LUIGI_CFG

# Override __import__ so we load the billiard library instead of multiprocessing.
# Use of multiprocessing within a celery worker causes error
# (AssertionError: daemonic processes are not allowed to have children)
# def _import(name, *args, **kwargs):   REMOVING THIS FOR NOW: BILLIARD BUG CAUSES ACCUMULATION OF OPEN PIPES
#     if name == 'multiprocessing':
#         name = 'billiard'
#     return original_import(name, *args, **kwargs)
#
# original_import = builtins.__import__
# builtins.__import__ = _import

app = Celery('MAS')

app.conf.beat_schedule = {
    "db-maintenance": {
        "task": "result_viewer.api.tasks.database_maintenance",
        "schedule": 60.0 * 60.0 * 2.0,  # runs ever two hours
        'args': (12,)
    }
}

# Using a string here means the worker doesn't have to serialize
# the configuration object to child processes.
# - namespace='CELERY' means all celery-related configuration keys
#   should have a `CELERY_` prefix.
app.config_from_object('django.conf:settings', namespace='CELERY')

# Load task modules from all registered Django app configs.
app.autodiscover_tasks()


@app.task(bind=True)
def debug_task(self):
    print(f'Request: {self.request!r}')