from django.shortcuts import redirect
from django.contrib.auth import views

from result_viewer.views import MixinForBaseTemplate


def login_redirect(request):
    return redirect('/home/')


#################### Auth View Overrides ####################
# Needed to override auth views to add necessary context
class MasLoginView(MixinForBaseTemplate, views.LoginView):
    pass


class MasLogoutView(MixinForBaseTemplate, views.LogoutView):
    pass
