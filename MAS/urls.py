"""MAS URL Configuration

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/2.0/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.contrib import admin
from django.conf import settings
from django.urls import path, include
#from django.contrib.auth import urls
from MAS import views
from django.contrib.auth import views as auth_views

urlpatterns = [
    path('', views.login_redirect, name='login_redirect'),
    path('admin/', admin.site.urls),
    # path('blast/', include('blast_viewer.urls')),
    path('viewer/', include('result_viewer.urls')),

    path('accounts/login/', views.MasLoginView.as_view(), name='login'),
    path('accounts/logout/', views.MasLogoutView.as_view(), name='logout'),
    path('accounts/password_change/', auth_views.PasswordChangeView.as_view(), name='password_change'),
    # path('accounts/', include('django.contrib.auth.urls')),

    path('home/', include('home.urls')),
    path('genome/', include('genome.urls'))
    # path('similarity_viewer/', include('similarity_viewer.urls'))
]

if not settings.IN_PRODUCTION:
    import debug_toolbar
    urlpatterns = [path('__debug__/', include(debug_toolbar.urls)),] + urlpatterns
