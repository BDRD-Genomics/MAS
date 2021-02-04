from django.contrib import admin

# Register your models here.
from django.contrib import admin

from . import models as genome_models


admin.site.register(genome_models.Genome)
admin.site.register(genome_models.Feature)
admin.site.register(genome_models.Annotation)
# admin.site.register(genome_models.HHSearch_Result)
# admin.site.register(genome_models.Blastp_Result)
# admin.site.register(genome_models.RPSBlast_Result)
