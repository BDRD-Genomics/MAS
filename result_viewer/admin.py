from django.contrib import admin

from . import models as result_viewer_models

# Register your models here.
# admin.site.register(result_viewer_models.Genome)
# admin.site.register(result_viewer_models.Feature)
# admin.site.register(result_viewer_models.Annotation)
admin.site.register(result_viewer_models.HHSearch_Result)
admin.site.register(result_viewer_models.Blastp_Result)
admin.site.register(result_viewer_models.RPSBlast_Result)
admin.site.register(result_viewer_models.PDB_Accession_Mapping)