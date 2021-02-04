from django.core.exceptions import ObjectDoesNotExist
from django.urls import reverse
from django.contrib.auth.models import User

from result_viewer.models import *


class Navigator:
    def next(self):
        if self.idx + 1 < self.size:
            return reverse('view-results', args=(self.queryset[self.idx + 1].accession, self.__class__.__name__, self.nav_arg,))
        else:
            return None

    def previous(self):
        if self.idx > 0:
            return reverse('view-results', args=(self.queryset[self.idx - 1].accession, self.__class__.__name__, self.nav_arg,))
        else:
            return None

    def go_to(self, index):
        return reverse('view-results', args=(self.queryset[index].accession, self.__class__.__name__, self.nav_arg,))

    def as_context(self):
        '''

        :param num_pagination_buttons: Number of 'integer' pagination buttons (not counting the next and
        previous buttons).
        :return: A list to add to a views context
        '''
        context = {}

        context['type'] = self.__class__.__name__

        context['description'] = self.description

        context['nav_arg'] = self.nav_arg

        context['buttons'] = [{
            'text': 'Previous',
            'href': self.previous(),
            'disabled': self.idx == 0
        }]

        for i in range(self.size):
            if i < 1 or i > self.size - 2 or abs(i - self.idx) < 5:
                text = i+1
                if i == 0:
                    text = 'First'
                elif i == self.size - 1:
                    text = 'Last'

                context['buttons'].append({
                    'text': text,
                    'href': self.go_to(i),
                    'disabled': self.idx == i
                })

        context['buttons'].append({
            'text': 'Next',
            'href': self.next(),
            'disabled': self.idx + 1 == self.size
        })

        return context

    def get_current_object(self):
        return self.queryset[self.idx]


class FlagNavigator(Navigator):
    def __init__(self, flag, accession):
        pk = int(accession, 36)
        self.nav_arg = flag
        self.queryset = Annotation.objects.filter(flag=flag).order_by('id')
        self.idx = self.queryset.filter(pk__lt=pk).count()
        self.size = self.queryset.count()

        self.flag_str = ''
        for flag_tup in Annotation.flag_options:
            if flag_tup[0] == flag:
                self.flag_str = flag_tup[1]

        self.description = 'You are navigating all proteins flagged as ' + self.flag_str

        # TODO: Ensure object exists and has flag
        try:
            # self.queryset.get(sequence=current_sequence)
            self.queryset.get(pk=pk)
        except ObjectDoesNotExist:
            raise ObjectDoesNotExist('Flag navigator can not be created because object does not exist.')


class GenomeNavigator(Navigator):
    def __init__(self, genome_name, accession=None):
        self.nav_arg = genome_name

        phage_obj = Genome.objects.get(genome_name=genome_name)
        self.queryset = Annotation.objects.filter(feature__genome=phage_obj).order_by('feature__start').distinct()
        self.size = self.queryset.count()
        self.description = 'You are navigating {}'.format(genome_name)

        if accession:
            pk = int(accession, 36)

            # Use the queryset as a value list to find the position in the phage
            # this_start = Feature.objects.filter(annotation_id=pk, phage=phage_obj).order_by('start')[0].start
            # self.idx = self.queryset.filter(feature__start__lt=this_start).count()
            # value_list = self.queryset.values_list('pk', flat=True)
            self.idx = list(self.queryset.values_list('pk', flat=True)).index(pk)
            # Ensure object exists and has flag
            try:
                self.queryset.get(pk=pk)
            except ObjectDoesNotExist:
                raise ObjectDoesNotExist('Genome navigator can not be created because object does not exist.')

        else:
            self.idx = 0


class AssignmentNavigator(Navigator):
    def __init__(self, user, accession=None):
        self.nav_arg = user
        user_obj = User.objects.get(username=user)
        self.queryset = Annotation.objects.filter(assigned_to=user_obj).order_by('id')
        self.size = self.queryset.count()
        self.description = 'You are navigating annotations assigned to ' + user

        if accession:
            pk = int(accession, 36)
            try:
                self.queryset.get(pk=pk)
                self.idx = self.queryset.filter(pk__lt=pk).count()

            except ObjectDoesNotExist:
                raise ObjectDoesNotExist('Assignment navigator can not be created because object does not exist.')

        else:
            self.idx = 0
