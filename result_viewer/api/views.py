from datetime import datetime
import subprocess
import json
from tempfile import NamedTemporaryFile

from django.http import Http404, HttpResponseBadRequest
from django.conf import settings
from django.db.models import Q, F, Func, Count, Value, When, Case, CharField
from django.db.models.functions import Length
from django.urls import reverse

from rest_framework import viewsets, permissions, decorators, parsers, response, status, exceptions
from rest_framework.views import APIView
from rest_framework.parsers import FormParser

from result_viewer.api.serializers import *
from result_viewer.api.tasks import run_single_search, run_multiple_search


class RunSearchAjaxView(APIView):
    '''
    This view is called when a user selects the 'Run' button
    '''

    def post(self, request, format=None):
        s = RunSearchAjaxSerializer(data=request.data)

        if s.is_valid():
            accession = s.data['accession']
            tool = s.data['tool']
            database = s.data['database']

            uri_scheme = 'https://' if request.is_secure() else 'http://'
            mas_server_uri = uri_scheme + request.site.domain
            run_single_search.delay(accession, tool, database, mas_server_uri)

            return response.Response(s.data)

        return response.Response(s.errors, status=status.HTTP_400_BAD_REQUEST)

    def check_permissions(self, request):
        super().check_permissions(request)
        if not (request.user.is_superuser or request.user.groups.filter(name='Data Editors').exists()):
            self.permission_denied(
                request, message='Permission Denied: Only data editors/superusers can use this API view'
            )


class RunAllPhageProteinsAjaxView(APIView):

    def post(self, request, format=None):
        s = RunAllPhageProteinsAjaxSerializer(data=json.loads(request.data['data']))

        if s.is_valid():
            # Get data from serializer
            genome_name = s.data['genome']
            rerun = s.data['rerun']
            tools_and_databases = s.data['tools_and_databases']

            uri_scheme = 'https://' if request.is_secure() else 'http://'
            mas_server_uri = uri_scheme + request.site.domain

            run_multiple_search.delay(genome_name, rerun, tools_and_databases, mas_server_uri)

            return response.Response(s.data)

        return response.Response(s.errors, status=status.HTTP_400_BAD_REQUEST)

    def check_permissions(self, request):
        super().check_permissions(request)
        if not (request.user.is_superuser or request.user.groups.filter(name='Data Editors').exists()):
            self.permission_denied(
                request, message='Permission Denied: Only data editors/superusers can use this API view'
            )


class PipelineAPIMixin(object):
    '''
    Mixin for adding correct permissions to views for use by luigi pipeline
    '''
    def check_permissions(self, request):
        super().check_permissions(request)
        if not (request.user.is_superuser or request.user.username == 'luigi'):
            self.permission_denied(
                request, message='Permission Denied: Only Luigi/superusers can use this API view'
            )


class TestConnectionView(PipelineAPIMixin, APIView):
    def get(self, request):
        return response.Response()


class GetProtSeqView(PipelineAPIMixin, APIView):
    permission_classes = [permissions.DjangoModelPermissions]
    queryset = Annotation.objects.all()

    def get(self, request, accession, format=None):
        try:
            annotation = self.queryset.get(pk=int(accession, 36))

        except Annotation.DoesNotExist:
            raise Http404

        s = ProteinSeqSerializer(annotation)
        return response.Response(s.data)


class GetGenomeView(PipelineAPIMixin, APIView):
    permission_classes = [permissions.DjangoModelPermissions]
    queryset = Genome.objects.all()

    def get(self, request, genome_name, format=None):
        try:
            genome = self.queryset.get(genome_name=genome_name)
            genome_features = Feature.objects.filter(genome=genome)
            counts = genome_features.aggregate(
                tRNAs=Count('id', filter=Q(type='tRNA')),
                repeats=Count('id', filter=Q(type='Repeat Region')),
                cds=Count('id', filter=Q(type='CDS'))
            )
            if counts['repeats'] > 0:
                DTR_length = len(genome_features.get(type='Repeat Region').annotation.sequence)
            else:
                DTR_length = 0

        except Genome.DoesNotExist:
            raise Http404

        s = GenomeSeqSerializer({
            'genome_name': genome.genome_name,
            'genome_sequence': genome.genome_sequence,
            'num_cds': counts['cds'],
            'num_trna': counts['tRNAs'],
            'len_dtr': DTR_length
        })
        return response.Response(s.data)


class UploadResultsView(PipelineAPIMixin, APIView):
    parser_classes = (parsers.MultiPartParser, parsers.FormParser,)

    def post(self, request):
        s = UploadResultsSerializer(data=request.data)

        if s.is_valid():
            tool = s.data['tool']
            annotation = Annotation.objects.get(pk=int(s.data['accession'], 36))

            def save_result(result_class):
                obj = result_class.objects.get(annotation=annotation, database=s.data['database'])

                if 'result' in request.data:
                    obj.result.save(name=s.data['accession'], content=request.data['result'])

                obj.status = s.data['status']
                obj.save()

            if tool == 'blastp':
                save_result(Blastp_Result)

            elif tool == 'hhsearch':
                save_result(HHSearch_Result)

            elif tool == 'rpsblast':
                save_result(RPSBlast_Result)

            return response.Response(s.data)

        return response.Response(s.errors, status=status.HTTP_400_BAD_REQUEST)


class GenomeData:
    def __init__(self, genome):
        self.genome_name = '<a href="{url}">{name}</a>'.format(
            url=reverse('genome:phage_detail', kwargs={'pk': genome.id}),
            name=genome.genome_name
        )
        self.organism = genome.organism
        self.genome_length = len(genome.genome_sequence)
        self.num_cds = genome.num_cds
        self.num_unannotated = genome.num_unannotated
        self.num_review = genome.num_review
        self.num_green = genome.num_green
        self.num_yellow = genome.num_yellow
        self.num_red = genome.num_red
        self.num_endo = genome.num_endo
        self.num_trna = genome.num_trna
        self.download = '<a href="{}">download fasta</a>'.format(
            reverse('genome:phage_download_fasta', kwargs={'genome_id': genome.id})
        )
        self.navigator = '<a href="{}"><div class="glyphicon glyphicon-hand-right"></div></a>'.format(
            reverse('phage-nav-redirect', kwargs={'phage_name': genome.genome_name})
        )


class GetPhageDataView(APIView):
    def get(self, request):
        params = self.read_parameters(request.GET)

        flag_options_reverse = dict((v, k) for k, v in Annotation.flag_options)
        cds = Count('feature', filter=Q(feature__type= 'CDS'))
        trna = Count('feature', filter=Q(feature__type= 'tRNA'))
        green = Count('feature__annotation__flag', filter=Q(feature__annotation__flag=flag_options_reverse['GREEN']))
        yellow = Count('feature__annotation__flag', filter=Q(feature__annotation__flag=flag_options_reverse['YELLOW']))
        red = Count('feature__annotation__flag', filter=Q(feature__annotation__flag=flag_options_reverse['RED']))
        unannotated = Count('feature__annotation__flag', filter=Q(feature__annotation__flag=flag_options_reverse['UNANNOTATED']))
        endolysin = Count('feature__annotation__flag', filter=Q(feature__annotation__flag=flag_options_reverse['ENDOLYSIN']))
        review = Count('feature__annotation__flag', filter=Q(feature__annotation__flag=flag_options_reverse['REVIEW NAME']))

        genomes = Genome.objects.annotate(num_green=green).annotate(num_yellow=yellow).annotate(num_red=red).\
            annotate(num_unannotated=unannotated).annotate(num_endo=endolysin).annotate(num_review=review).\
            annotate(num_cds=cds).annotate(num_trna=trna)

        total_num_genomes = genomes.count()

        genomes = genomes.filter(
            Q(genome_name__icontains=params['search_val']) |
            Q(organism__icontains=params['search_val'])
        ).order_by(self.get_order_by_arg(params['order_col'], params['order_dir']))

        filtered_num_genomes = genomes.count()
        genomes = genomes[ params['start'] : params['start']+params['length'] ]

        genome_data = [GenomeData(p) for p in genomes]
        # s = GenomeDataSerializer(phage_data, many=True)

        return response.Response(GenomeDataListSerializer({
            'data': genome_data,
            'draw': params['draw'],
            'recordsTotal': total_num_genomes,
            'recordsFiltered': filtered_num_genomes
        }).data)

    def read_parameters(self, query_dict):
        """ Converts and cleans up the GET parameters. """
        return {
            'start': int(query_dict.get('start')),
            'length': int(query_dict.get('length')),
            'draw': int(query_dict.get('draw')),
            'order_col': int(query_dict.get('order[0][column]')),
            'order_dir': query_dict.get('order[0][dir]'),
            'search_val': query_dict.get('search[value]')
        }

    def get_order_by_arg(self, order_col, order_dir):
        if order_col == 0:
            return 'genome_name' if order_dir == 'asc' else '-genome_name'
        elif order_col == 1:
            return 'organism' if order_dir == 'asc' else '-organism'
        elif order_col == 2:
            return Length('genome_sequence').asc() if order_dir == 'asc' else Length('genome_sequence').desc()
        elif order_col == 3:
            return 'num_cds' if order_dir == 'asc' else '-num_cds'
        elif order_col == 4:
            return 'num_unannotated' if order_dir == 'asc' else '-num_unannotated'
        elif order_col == 5:
            return 'num_review' if order_dir == 'asc' else '-num_review'
        elif order_col == 6:
            return 'num_green' if order_dir == 'asc' else '-num_green'
        elif order_col == 7:
            return 'num_yellow' if order_dir == 'asc' else '-num_yellow'
        elif order_col == 8:
            return 'num_red' if order_dir == 'asc' else '-num_red'
        elif order_col == 9:
            return 'num_endo' if order_dir == 'asc' else '-num_endo'
        elif order_col == 10:
            return 'num_trna' if order_dir == 'asc' else '-num_trna'


class AnnotationData:
    def __init__(self, annotation, genome_name=None):
        self.sequence = '<div class="glyphicon glyphicon-menu-right aa-control details-control"></div>' \
                        '<input type="hidden" class="annotation" value="{pk}">'.format(pk=annotation.pk)
        self.feature = '<ul class="list-group">'
        for feature in annotation.feature_set.all():
            self.feature += '<li class="list-group-item list-group-item-text">'
            self.feature += '<a href="{url}">{description}</a>'.format(
                url=reverse('genome:feature_detail', kwargs={'pk': feature.pk}),
                description=feature
            )
        self.feature += '</ul>'
        self.accession = '<a href="{url}">{accession}</a>'.format(
            url=reverse('genome:annotation_detail', kwargs={'pk': annotation.pk}),
            accession=annotation.accession_sql
        )
        self.annotation = annotation.annotation
        self.length = len(annotation.sequence)
        self.public_notes = annotation.public_notes
        self.private_notes = annotation.private_notes
        self.flag = annotation.get_flag_display()
        self.assigned_to = annotation.assigned_to
        self.download = '<a href="{url}"><div class="glyphicon glyphicon-download"></div></a>'.format(
            url=reverse('genome:annotation_download', kwargs={'annotation_id': annotation.pk})
        )
        self.history = '<a href="{url}"><i class="fa fa-history"></i></a>'.format(
            url=reverse('genome:annotation_history', kwargs={'annotation_pk': annotation.pk})
        )

        if genome_name:
            view_results_url = reverse('view-results', kwargs={
                'navigator': 'GenomeNavigator',
                'accession': annotation.accession,
                'nav_arg': genome_name
            })
        else:
            view_results_url = reverse('accession-redirect', kwargs={'accession': annotation.accession})

        self.view_results = '<a href="{url}"><div class="glyphicon glyphicon-hand-right"></div></a>'.format(
            url=view_results_url
        )


class GetAnnotationListView(APIView):
    def get(self, request):
        params = self.read_parameters(request.GET)

        if params['genome_id']:
            annotations = Annotation.objects.filter(feature__genome_id=params['genome_id']).distinct()
            genome_name = Genome.objects.get(id=params['genome_id']).genome_name
        else:
            annotations = Annotation.objects
            genome_name = None

        # Annotate each annotation with accession (Only works in mysql/mariadb)
        annotations = annotations.annotate(
            accession_sql=Func(Func(F('id'), 10, 36, function='CONV'), 5, 0, function='LPAD')
        )
        # Annotate with flag strings
        when_clauses = [When(flag=i, then=Value(flag)) for i, flag in Annotation.flag_options]
        annotations = annotations.annotate(
            flag_str=Case(*when_clauses, output_field=CharField())
        )

        total_num_annotations = annotations.count()

        annotations = annotations.filter(
            Q(flag_str__icontains=params['search_val']) |
            Q(sequence__icontains=params['search_val']) |
            Q(accession_sql__icontains=params['search_val']) |
            Q(annotation__icontains=params['search_val']) |
            Q(public_notes__icontains=params['search_val']) |
            Q(private_notes__icontains=params['search_val']) |
            Q(assigned_to__username__icontains=params['search_val'])
        ).order_by(self.get_order_by_arg(params['order_col'], params['order_dir']))

        filtered_num_annotations = annotations.count()

        annotations = annotations[params['start'] : params['start'] + params['length']]
        annotation_data = [AnnotationData(a, genome_name) for a in annotations]

        return response.Response(AnnotationDataListSerializer({
            'data': annotation_data,
            'draw': params['draw'],
            'recordsTotal': total_num_annotations,
            'recordsFiltered': filtered_num_annotations
        }).data)

    def read_parameters(self, query_dict):
        return {
            'start': int(query_dict.get('start')),
            'length': int(query_dict.get('length')),
            'draw': int(query_dict.get('draw')),
            'order_col': int(query_dict.get('order[0][column]')),
            'order_dir': query_dict.get('order[0][dir]'),
            'search_val': query_dict.get('search[value]'),
            'genome_id': None if query_dict.get('genome_id') == 'NaN' else int(query_dict.get('genome_id'))
        }

    def get_order_by_arg(self, order_col, order_dir):
        if order_col == 2:
            return 'accession_sql' if order_dir == 'asc' else '-accession_sql'
        elif order_col == 3:
            return 'annotation' if order_dir == 'asc' else '-annotation'
        elif order_col == 4:
            return Length('sequence').asc() if order_dir == 'asc' else Length('sequence').desc()
        elif order_col == 5:
            return 'public_notes' if order_dir == 'asc' else '-public_notes'
        elif order_col == 6:
            return 'private_notes' if order_dir == 'asc' else '-private_notes'
        elif order_col == 7:
            return 'flag_str' if order_dir == 'asc' else '-flag_str'
        elif order_col == 8:
            return 'assigned_to__username' if order_dir == 'asc' else '-assigned_to__username'
        else:
            return 'accession_sql' if order_dir == 'asc' else '-accession_sql'
