import sys
import os


def run_search_for_protein(annotation_accession, tool, database, mas_server):
    '''
    annotation_accession: Accession of annotation describing protein we are running
    tool: What tool will we use
    database: What database will we search against
    mas_server: Where the site is hosted. Could be ip:port or domain name
    '''
    sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    import luigi
    # Now import pipeline endpoint
    from AnnotationToolPipeline.AnnotationToolPipeline.pipeline import Send_Results_To_MAS

    job_array = [Send_Results_To_MAS(annotation_accession=annotation_accession, database=database, tool=tool, mas_server=mas_server)]
    luigi.build(job_array, workers=1)

    print('PipelineComplete')


if __name__ == '__main__':
    annotation_accession, tool, database, mas_server = sys.argv[1:]
    run_search_for_protein(annotation_accession, tool, database, mas_server)