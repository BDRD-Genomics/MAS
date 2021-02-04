from datetime import datetime, timedelta

from django.views import generic

from result_viewer.views import MixinForBaseTemplate
from result_viewer.models import Annotation


class HomePageView(MixinForBaseTemplate, generic.TemplateView):
    template_name = 'home/index.html'

    def get_context_data(self, **kwargs):
        context = super(HomePageView, self).get_context_data(**kwargs)

        if self.request.user.is_authenticated:
            weekstart = datetime.now() - timedelta(days=datetime.now().weekday(), hours=datetime.now().hour)
            this_weeks_annotations = Annotation.history.filter(
                history_date__gte=weekstart,
                history_user=self.request.user,
                history_type='~'
            )
            context['annotated_this_week'] = this_weeks_annotations.values('id').distinct().count()

            last_weeks_annotations = Annotation.history.filter(
                history_date__lte=weekstart,
                history_date__gte=weekstart - timedelta(days=7),
                history_user=self.request.user,
                history_type='~'
            )
            context['annotated_last_week'] = last_weeks_annotations.values('id').distinct().count()

        return context
