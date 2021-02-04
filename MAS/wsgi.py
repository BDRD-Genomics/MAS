"""
WSGI config for MAS project.

It exposes the WSGI callable as a module-level variable named ``application``.

For more information on this file, see
https://docs.djangoproject.com/en/2.0/howto/deployment/wsgi/
"""

import os

from django.core.wsgi import get_wsgi_application
from django.conf import settings

os.environ.setdefault("DJANGO_SETTINGS_MODULE", "MAS.settings")
# sys.path.append('/home/django_projects/MyProject')
os.environ['PATH'] = ("/home/daemon/miniconda/envs/mas/bin:" + os.environ['PATH'])
if settings.IN_PRODUCTION:
    # os.environ['PATH'] = ("/home/lims/anaconda3/envs/MAS/bin:" + os.environ['PATH'])
    pass
application = get_wsgi_application()
