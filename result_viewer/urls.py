from django.urls import path, include

from . import views

urlpatterns = [
    path('result/<str:accession>/<str:navigator>/<str:nav_arg>',
         views.ViewResults.as_view(),
         name='view-results'
         ),
    path('no-result/<str:navigator>/<str:nav_arg>',
         views.ViewNoResults.as_view(),
         name='no-results'
         ),
    path('navigator/flag/<int:flag>', views.FlagNavRedirect.as_view(), name='flag-nav-redirect'),
    path('navigator/phage/<str:phage_name>', views.GenomeNavRedirect.as_view(), name='phage-nav-redirect'),
    path('navigator/assignment/<str:user>', views.AssignmentNavRedirect.as_view(), name='assignment-nav-redirect'),
    path('result/<str:accession>', views.AccessionRedirect.as_view(), name='accession-redirect'),
    path('<str:accession>', views.AccessionRedirect.as_view(), name='accession-redirect'),
    path('preferences/change-theme', views.ChangeTheme.as_view(), name='change-theme'),
    # path('NoProteins-<str:navigator>-<str:nav_arg>', views.ProteinDoesNotExistView.as_view(), name='no-proteins'),
    path('api/', include('result_viewer.api.urls')),
]