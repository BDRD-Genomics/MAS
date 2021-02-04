from django.urls import include, path

from result_viewer.api.views import *


urlpatterns = [
    path('run-search/', RunSearchAjaxView.as_view(), name='run_search'),
    path('run-search-for-phage/', RunAllPhageProteinsAjaxView.as_view(), name='run_search_for_phage'),
    path('upload-results/', UploadResultsView.as_view(), name='upload_results'),
    path('test', TestConnectionView.as_view(), name='test_connection'),
    path('get-protein/<str:accession>', GetProtSeqView.as_view(), name='get_protein'),
    path('get-phage-data', GetPhageDataView.as_view(), name='get_phage_data'),
    path('get-annotation-data', GetAnnotationListView.as_view(), name='get_annotation_data'),
    # path('get-phage-data', GetPhageDataView.as_view(), name='get_phage_data'),
    path('api-auth/', include('rest_framework.urls', namespace='rest_framework'))
]