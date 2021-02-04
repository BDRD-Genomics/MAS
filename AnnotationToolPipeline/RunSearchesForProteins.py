import sys
import json
import os


def run_searches_for_proteins(list_of_searches, mas_server):
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    import luigi
    from AnnotationToolPipeline.AnnotationToolPipeline.pipeline import Run_Pipeline_For_Proteins
    from AnnotationToolPipeline.AnnotationToolPipeline.global_config import global_config

    for search_params in list_of_searches:
        if set(search_params) != {'accession', 'tool', 'database'}:
            raise RuntimeError('Incorrect JSON attributes')

    task = Run_Pipeline_For_Proteins(input_list=list_of_searches, mas_server=mas_server)
    luigi.build([task], workers=global_config.NUM_WORKERS)
    print('PipelineComplete')


if __name__ == '__main__':
    mas_server = sys.argv[-1]
    list_of_searches = json.load(sys.stdin)
    run_searches_for_proteins(list_of_searches, mas_server)