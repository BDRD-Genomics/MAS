from celery import shared_task
from django.http import Http404, HttpResponseBadRequest
from urllib.error import URLError
import subprocess
import json
from datetime import datetime, timedelta

from MAS.celery import app
from result_viewer.models import *
# from AnnotationToolPipeline.RunSearchForProtein import run_search_for_protein
from AnnotationToolPipeline import RunSearchForProtein
# from AnnotationToolPipeline.RunSearchesForProteins import run_searches_for_proteins
from AnnotationToolPipeline import RunSearchesForProteins


@app.task
def database_maintenance(timeout_hrs=12):
    current_datetime = datetime.now()

    for model in [HHSearch_Result, Blastp_Result, RPSBlast_Result]:
        model.objects.filter(status=1, run_date__lt=current_datetime - timedelta(hours=timeout_hrs)).update(status=2)


@shared_task
def run_single_search(accession, tool, database, site):
    if is_luigi_server_functional() and is_mas_reachable_from_worker(site):
        annotation = Annotation.objects.get(pk=int(accession, 36))

        # Initialize database entry for search
        def create_result_entry(result_class):
            if result_class.objects.filter(annotation=annotation, database=database).exists():
                obj = result_class.objects.get(annotation=annotation, database=database)

                # If already running, do not start a second run
                if obj.status == 1:
                    raise Http404

                obj.status = 1
                obj.run_date = datetime.now()
                obj.save()

            else:
                obj = result_class(
                    annotation=annotation,
                    database=database,
                    status=1,
                    run_date=datetime.now()
                )
                obj.save()

        if tool == 'blastp':
            create_result_entry(Blastp_Result)

        elif tool == 'hhsearch':
            create_result_entry(HHSearch_Result)

        elif tool == 'rpsblast':
            create_result_entry(RPSBlast_Result)

        else:
            raise HttpResponseBadRequest

        # Start pipeline
        # run_search_for_protein(accession, tool, database, site)
        print('Running search pipeline')
        script_path = os.path.abspath(RunSearchForProtein.__file__)
        subprocess.run(['python', script_path, accession, tool, database, site])
        print('Search pipeline finished')


@shared_task
def run_multiple_search(genome_name, rerun, tools_and_databases, site):
    if is_luigi_server_functional() and is_mas_reachable_from_worker(site):
        # Get list of phage's annotations
        genome_obj = Genome.objects.get(genome_name=genome_name)

        # Iterate through phage's annotations, selecting which ones to run searches for
        annotation_objs = Annotation.objects.filter(feature__genome=genome_obj).order_by('feature__start').distinct()

        searches_to_run = []
        for tool in tools_and_databases:
            # Get model for tool
            result_model = None
            if tool == 'hhsearch':
                result_model = HHSearch_Result

            elif tool == 'blastp':
                result_model = Blastp_Result

            elif tool == 'rpsblast':
                result_model = RPSBlast_Result

            else:
                raise HttpResponseBadRequest

            for database in tools_and_databases[tool]:
                codes_to_not_run = [1]
                if not rerun:
                    codes_to_not_run.append(0)

                search_objs = result_model.objects.filter(
                    annotation__in=annotation_objs,
                    database=database,
                    status__in=codes_to_not_run
                )
                annotation_ids_from_filter = {x['annotation_id'] for x in search_objs.values()}

                for annotation in annotation_objs:
                    if annotation.id not in annotation_ids_from_filter:
                        try:
                            obj = result_model.objects.get(annotation=annotation, database=database)
                            obj.status = 1
                            obj.run_date = datetime.now()
                            obj.save()

                        except result_model.DoesNotExist:
                            obj = result_model(
                                annotation=annotation,
                                database=database,
                                status=1,
                                run_date=datetime.now()
                            )
                            obj.save()

                        searches_to_run.append({
                            'tool': tool,
                            'database': database,
                            'accession': annotation.accession
                        })

        if len(searches_to_run) > 0:
            # Due to billiard's open pipe bug, we must start the pipeline through a subprocess command
            # run_searches_for_proteins(searches_to_run, site)
            print('Running search pipeline')
            script_path = os.path.abspath(RunSearchesForProteins.__file__)
            subprocess.run(['python', script_path, site], input=json.dumps(searches_to_run), text=True)
            print('Search pipeline finished')


def is_luigi_server_functional():
    try:
        import luigi
        return luigi.build([])

    except URLError as e:
        print('Celery worker not able to access luigi scheduler!')
        return False


def is_mas_reachable_from_worker(site):
    from AnnotationToolPipeline.AnnotationToolPipeline.global_config import global_config
    from django.urls import reverse
    import requests

    r = requests.get(
        site + reverse('test_connection'),
        auth=(global_config.MAS_USERNAME, global_config.MAS_PASSWORD),
        verify=global_config.MAS_CRT
    )
    if r.status_code != 200:
        message = 'Unable to connect to MAS server from celery worker. ' \
                  'Response status code = {}. Response text = {}.'.format(r.status_code, r.text)
        print(message)
        return False

    else:
        return True
