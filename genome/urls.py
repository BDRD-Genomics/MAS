# from django.conf.urls import url, include
# from django.views.decorators.cache import cache_page
from django.urls import path
# from django.views.decorators.csrf import csrf_exempt
from genome import views

app_name = 'genome'

urlpatterns = [
    path('genome/detail/<int:pk>', views.Genome_Detail.as_view(), name='phage_detail'),
    path('genome/list/', views.Genome_List_SS.as_view(), name='phage_list'),
    path('genome/download/fasta/<int:genome_id>', views.genome_download_fasta, name='phage_download_fasta'),
    path('genome/download/deliverables/<int:genome_id>', views.download_deliverables, name='download_deliverables'),
    path('phage-genome/upload/', views.Upload_Phage.as_view(), name='upload_phage'),
    path('bacterial-genome/upload/', views.Upload_Bacterial_Genome.as_view(), name='upload_bacterial_genome'),
    path('custom-genome/upload/', views.Upload_Custom_Genome.as_view(), name='upload_custom_genome'),
    path('genome/delete/', views.Genome_Delete.as_view(), name='phage_delete'),
    path('genome/delete/confirm/', views.Confirm_Genome_Delete.as_view(), name='confirm_phage_delete'),
    path('feature/detail/<int:pk>', views.Feature_Detail.as_view(), name='feature_detail'),
    path('annotation/list/', views.Annotation_List_Serverside.as_view(), name='annotation_list'),
    path('annotation/detail/<int:pk>', views.Annotation_Detail.as_view(), name='annotation_detail'),
    path('annotation/download/<int:annotation_id>', views.annotation_download, name='annotation_download'),
    path('annotation/upload/', views.Upload_Annotation.as_view(), name='upload_annotations'),
    path('annotation/upload/confirm/', views.Confirm_Upload_Annotation.as_view(), name='confirm_upload_annotations'),
    # path('annotation/bulk', views.Annotation_Bulk.as_view(), name='annotation_bulk'),
    # path('annotation/my_annotations', views.My_Annotations.as_view(), name='my_annotations'),
    path('annotation/history/<int:annotation_pk>', views.Annotation_History.as_view(), name='annotation_history'),
    path('annotation/download/excel/', views.download_excel_annotations, name='download_excel_annotations'),
    path('annotation/download/excel_template/', views.download_excel_template, name='download_excel_template'),
    path('annotation/download/excel/unannotated_annotations/', views.download_unannotated_annotations, name='download_unannotated_annotations'),

    # ajax url used for nucleotide sequence drop down
    path('ajax/phage/get', views.Get_Genome.as_view(), name='get_phage'),
    # ajax url used for amino acid sequence drop down
    path('ajax/annotation/get/aa_sequence/', views.Get_AA_Sequence.as_view(), name='get_aa_sequence'),
    path('ajax/feature/get/nucleotide_sequence/', views.Get_Feature_Sequence.as_view(), name='get_feature_sequence'),
]
