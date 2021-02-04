from django.core.management.base import BaseCommand, CommandError
from django.db.utils import IntegrityError
from django.contrib.auth.models import User
from django.contrib.auth.models import Group


class Command(BaseCommand):
    help = "Create an admin user with specified password"

    def add_arguments(self, parser):
        parser.add_argument('admin_password', type=str)
        parser.add_argument('luigi_password', type=str)

    def handle(self, *args, **options):
        try:
            User.objects.create_superuser('admin', '', options['admin_password'])

        # If admin user already created, we will still try to create luigi
        except IntegrityError as e:
            pass

        try:
            luigi_user = User.objects.create_user('luigi', '', options['luigi_password'])
            data_editors_group = Group.objects.get(name='Data Editors')
            data_editors_group.user_set.add(luigi_user)

        except IntegrityError as e:
            pass
