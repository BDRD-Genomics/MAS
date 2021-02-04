from rest_framework import serializers

from result_viewer.models import *


class ProteinSeqSerializer(serializers.ModelSerializer):
    class Meta:
        model = Annotation
        fields = ['sequence']


class RunSearchAjaxSerializer(serializers.Serializer):
    tool = serializers.CharField(max_length=20)
    database = serializers.CharField(max_length=20)
    accession = serializers.CharField(max_length=5, min_length=5)

    def validate_accession(self, value):
        # Make sure annotation with accession exists
        if Annotation.objects.filter(pk=int(value, 36)).count() == 0:
            raise serializers.ValidationError('Invalid accession: no annotation with accession found.')

        obj = Annotation.objects.get(pk=int(value, 36))
        if obj.flag in [8, 9]:
            raise serializers.ValidationError('Invalid accession: Accession links to non-protein annotation')

        return value

    def validate_tool(self, value):
        if not value in ['blastp', 'hhsearch', 'rpsblast']:
            raise serializers.ValidationError('Invalid tool: tool not recognized.')

        return value

    def validate(self, data):
        #   TODO Make sure tool allows database! + make sure annotation exists (and make sure annotation is a protein)
        # Make sure entry in database exists
        database = data['database']
        tool = data['tool']
        annotation = Annotation.objects.get(pk=int(data['accession'], 36))

        def test_if_pipeline_already_running(result_class):
            qs = result_class.objects.filter(annotation=annotation, database=database)
            if qs.count() == 1:
                obj = qs[0]
                if obj.status == 1:
                    raise serializers.ValidationError('Pipeline already in progress.')

        if tool == 'blastp':
            test_if_pipeline_already_running(Blastp_Result)

        elif tool == 'hhsearch':
            test_if_pipeline_already_running(HHSearch_Result)

        elif tool == 'rpsblast':
            test_if_pipeline_already_running(RPSBlast_Result)

        return data


class ToolsAndDatabasesSerializer(serializers.Serializer):
    hhsearch = serializers.ListField(
        required=False,
        child=serializers.ChoiceField(choices=HHSearch_Result.database_options)
    )
    blastp = serializers.ListField(
        required=False,
        child=serializers.ChoiceField(choices=Blastp_Result.database_options)
    )
    rpsblast = serializers.ListField(
        required=False,
        child=serializers.ChoiceField(choices=RPSBlast_Result.database_options)
    )


class RunAllPhageProteinsAjaxSerializer(serializers.Serializer):
    genome = serializers.CharField(max_length=Genome._meta.get_field('genome_name').max_length)
    rerun = serializers.BooleanField(default=False)
    tools_and_databases = ToolsAndDatabasesSerializer()

    def validate_genome(self, value):
        try:
            Genome.objects.get(genome_name=value)
        except Genome.DoesNotExist:
            raise serializers.ValidationError('Invalid genome name: Does not exist')

        return value


class UploadResultsSerializer(serializers.Serializer):
    tool = serializers.CharField(max_length=20)
    database = serializers.CharField(max_length=20)
    accession = serializers.CharField(max_length=5, min_length=5)
    result = serializers.FileField(required=False)
    status = serializers.IntegerField()

    def validate_accession(self, value):
        # Make sure annotation with accession exists
        if Annotation.objects.filter(pk=int(value, 36)).count() == 0:
            raise serializers.ValidationError('Invalid accession: no annotation with accession found.')

        return value

    def validate_tool(self, value):
        if not value in ['blastp', 'hhsearch', 'rpsblast']:
            raise serializers.ValidationError('Invalid tool: tool not recognized.')

        return value

    def validate_status(self, value):
        if value not in [x[0] for x in search_status_options]:
            raise serializers.ValidationError('Invalid status: not accepted by database')

        return value

    def validate(self, data):
        # Make sure entry in database exists
        database = data['database']
        tool = data['tool']
        annotation = Annotation.objects.get(pk=int(data['accession'], 36))

        def test_if_result_entry_exists(result_class):
            qs = result_class.objects.filter(annotation=annotation, database=database)
            if qs.count() != 1:
                raise serializers.ValidationError(
                    'No {} entry which matches parameters found!'.format(result_class.__name__)
                )

        if tool == 'blastp':
            test_if_result_entry_exists(Blastp_Result)

        elif tool == 'hhsearch':
            test_if_result_entry_exists(HHSearch_Result)

        elif tool == 'rpsblast':
            test_if_result_entry_exists(RPSBlast_Result)

        return data


class DataTablesServerSideSerializer(serializers.Serializer):
    draw = serializers.IntegerField()
    recordsTotal = serializers.IntegerField()
    recordsFiltered = serializers.IntegerField()


class GenomeDataSerializer(serializers.Serializer):
    genome_name = serializers.CharField(max_length=Genome._meta.get_field('genome_name').max_length)
    organism = serializers.CharField(max_length=Genome._meta.get_field('organism').max_length)
    genome_length = serializers.IntegerField()
    num_cds = serializers.IntegerField()
    num_unannotated = serializers.IntegerField()
    num_review = serializers.IntegerField()
    num_green = serializers.IntegerField()
    num_yellow = serializers.IntegerField()
    num_red = serializers.IntegerField()
    num_endo = serializers.IntegerField()
    num_trna = serializers.IntegerField()
    download = serializers.CharField(max_length=100)
    navigator = serializers.CharField(max_length=100)


class GenomeDataListSerializer(DataTablesServerSideSerializer):
    data = GenomeDataSerializer(many=True)


class AnnotationDataSerializer(serializers.Serializer):
    sequence = serializers.CharField(max_length=100)
    annotation = serializers.CharField(max_length=Annotation._meta.get_field('annotation').max_length)
    feature = serializers.CharField(max_length=1000)
    accession = serializers.CharField(max_length=10)
    length = serializers.IntegerField()
    public_notes = serializers.CharField(max_length=Annotation._meta.get_field('public_notes').max_length)
    private_notes = serializers.CharField(max_length=Annotation._meta.get_field('private_notes').max_length)
    flag = serializers.CharField(max_length=100)
    assigned_to = serializers.CharField(max_length=100)
    download = serializers.CharField(max_length=500)
    history = serializers.CharField(max_length=500)
    view_results = serializers.CharField(max_length=500)


class AnnotationDataListSerializer(DataTablesServerSideSerializer):
    data = AnnotationDataSerializer(many=True)
