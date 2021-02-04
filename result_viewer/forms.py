from django import forms
from django.db import models
from django.contrib.auth.models import User

# from result_viewer.models import Annotation, Genome
from genome.models import Annotation, Genome

def is_ascii(string):
    '''
    Determine if a given string is ascii
    :param string: string to test
    :return: true/false
    '''
    return len(string) == len(string.encode())


class AnnotationForm(forms.ModelForm):
    go_to_next = forms.BooleanField(label='Go to next protein on submit?', initial=True, required=False)
    assigned_to = forms.ModelChoiceField(queryset=User.objects.filter(groups__name='Bioinformatics'), required=False)

    def clean(self):
        cleaned_data = super().clean()

        for key in cleaned_data:
            value = cleaned_data[key]

            if type(value) == str:
                if not is_ascii(value):
                    raise forms.ValidationError('Non-ascii text detected. All text must be ascii.')

                if value.strip() != value:
                    raise forms.ValidationError('Text may not start or end with whitespace.')

                if key == 'annotation':
                    if any([x in value for x in ['\n', '|', '\t', '\r', '=', ';']]):
                        raise forms.ValidationError('Annotation contains illegal characters (May not contain the '
                                                    'following: \\n, |, \\t, \\r, =, ;)')


    class Meta:
        model = Annotation
        fields = ['annotation', 'flag', 'public_notes', 'private_notes', 'assigned_to']


class GenomeSearchForm(forms.Form):
    search_genome = forms.ModelChoiceField(
        queryset=Genome.objects.order_by('genome_name'),
        empty_label='',
        required=False
    )

    def __init__(self, *args, **kwargs):
        super(GenomeSearchForm, self).__init__(*args, **kwargs)