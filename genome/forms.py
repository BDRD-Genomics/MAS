from django import forms
from django.contrib.auth.models import User
from django.forms import CheckboxInput
from django.db.models.fields.files import FieldFile
from django.core.files.uploadedfile import TemporaryUploadedFile, InMemoryUploadedFile, File
from django.core.exceptions import ValidationError
from django.utils.safestring import mark_safe
from django.urls import reverse_lazy
from django.utils.functional import lazy
from django.utils.text import format_lazy
from django.templatetags.static import static

import re
from copy import deepcopy
import io
from argparse import Namespace

from crispy_forms.helper import FormHelper
from Bio import SeqIO
from Bio.Data.CodonTable import TranslationError
import pandas as pd

from genome import models as genome_models
from genome import views as genome_views
from genome.genomic_loci_conversions import *


def validate_fasta_file(instance):

    try:
        file_text = get_file_handle(instance, mode='r')
        #, alphabet=SingleLetterAlphabet
        #accepts fasta files and tbl files <- not good
        fasta_file = SeqIO.read(file_text, "fasta")
        if len(fasta_file) < 1:
            raise ValidationError("%s has no sequence." % instance.name)

        # Check for non-IUPAC chars
        match = re.search('[^GACTRYSWKMBDHVN]', str(fasta_file.seq).upper())
        if match:
            raise ValidationError('FASTA file sequence contains unaccepted character: {}'.format(match.group()))

    except ValueError as e:
        raise ValidationError('The following error occurred when attempting to validate fasta file: {}'.format(str(e)))


def validate_excel_file(instance):
    try:
        excel_file = get_file_handle(instance)
        file = pd.read_excel(excel_file)
    except:
        raise ValidationError("%s is not an Excel file." % instance.name)


def validate_coordinate_file(instance):
    for line in get_file_handle(instance, mode='r'):
        line_list = line.split()

        if len(line_list) != 3:
            raise ValidationError("%s incorrectly formatted. Each line should contain 3 whitespace separated values: "
                                  "strand (+ or -), start coordinate, and stop coordinate" % instance.name)

        if line_list[0] not in ['+', '-']:
            raise ValidationError("{} incorrectly formatted. First column should be '+' or '-'. Value detected: {}"
                                  "".format(instance.name, line_list[0]))

        try:
            int(line_list[1])
        except ValueError as e:
            raise ValidationError("{} incorrectly formatted. Second column should be an integer. Value detected: {}"
                                  "".format(instance.name, line_list[1]))

        try:
            int(line_list[2])
        except ValueError as e:
            raise ValidationError("{} incorrectly formatted. Third column should be an integer. Value detected: {}"
                                  "".format(instance.name, line_list[2]))


def get_file_handle(instance, mode='rb'):
    re_b = re.compile('b')
    if isinstance(instance, FieldFile):
        instance = instance.file
    if isinstance(instance, TemporaryUploadedFile):
        return open(instance.file.name, mode)
    elif isinstance(instance, InMemoryUploadedFile):
        if re_b.search(mode):
            return deepcopy(instance)
        else:
            return io.TextIOWrapper(deepcopy(instance))
    elif isinstance(instance, File):
        return open(instance.file.name, mode)
    else:
        raise ValidationError('%s may not be a file.' % instance.name)


def parse_prots_from_coords(cds_fh, genome_rec, selected_table):
    for cds_line in cds_fh:
        # Get coordinates
        strand_, start_, stop_ = cds_line.split()
        start, stop, strand = coordinate_file_to_db_standard(int(start_), int(stop_), strand_)

        # produce protein sequence
        prot = get_protein_sequence(start, stop, strand, genome_rec, table=selected_table)

        # Ensure entire sequence translated
        if len(prot) + 1.0 != (1 + int(stop_) - int(start_)) / 3:
            raise TranslationError('"{line}" is not a valid CDS'.format(line=cds_line.strip()))

        yield prot, Namespace(start=start, stop=stop, strand=strand)


def get_set_of_used_organisms():
    return tuple((x, x) for x in genome_models.Genome.objects.all().values_list('organism', flat=True).distinct())

### START classes moved from home.forms in LIMS ###
class CrispyModelForm(forms.ModelForm):
    class Meta:
        abstract = True

    def __init__(self, *args, **kwargs):
        super(CrispyModelForm, self).__init__(*args, **kwargs)
        self.helper = CrispyHorizontalFormHelper()
        self.helper.form_tag = False
        self.helper.disable_csrf = True


class CrispyHorizontalFormHelper(FormHelper):
    def __init__(self, *args, **kwargs):
        super(CrispyHorizontalFormHelper, self).__init__(*args, **kwargs)
        self.label_class = 'col-lg-1'
        self.field_class = 'col-lg-5'
        self.form_class = 'form-horizontal'


class CrispyHorizontalInternalFormHelper(CrispyHorizontalFormHelper):
    def __init__(self, *args, **kwargs):
        super(CrispyHorizontalInternalFormHelper, self).__init__(*args, **kwargs)
        self.form_tag = False
        self.disable_csrf = True
### START classes moved from home.forms in LIMS ###

class BaseAnnotationFormset(forms.BaseInlineFormSet):

    def __int__(self, *args, **kwargs):
        super(BaseAnnotationFormset, self).__init__(*args, **kwargs)
        self.form.helper = CrispyHorizontalInternalFormHelper()
        for form in self.forms:
            form.empty_permitted = False


class Annotation_Form(CrispyModelForm):
    class Meta:
        model = genome_models.Annotation
        fields = (
            'annotation',
            'public_notes',
            'private_notes',
            'flag'
        )

    def save(self, commit=True):
        annotation = super(Annotation_Form, self).save(commit=False)
        if commit:
            annotation.save()
        return annotation


class Bacterial_Genome_Upload_Form(forms.Form):
    user_choices = User.objects.none()

    name = forms.CharField(
        validators=[
            genome_models.validate_genome_name,
            genome_models.validate_duplicate_name
        ],
        max_length=100,
        required=True,
    )

    upload = forms.FileField(
        validators=[validate_fasta_file],
        help_text='Must be a single fasta file!',
    )

    assign_to = forms.ModelChoiceField(
        queryset=user_choices,
        required=False,
    )

    def __init__(self, *args, **kwargs):
        super(Bacterial_Genome_Upload_Form, self).__init__(*args, **kwargs)
        # User choices should be limited to those with permissions to change annotations
        self.fields['assign_to'].queryset = genome_views.get_annotation_editors()


# Does not perform validation to ensure value is one of the selected options
class DynamicChoiceField(forms.ChoiceField):
    def validate(self, value):
        super(forms.ChoiceField, self).validate(value)


class Custom_Genome_Upload_Form(forms.Form):
    user_choices = User.objects.none()

    name = forms.CharField(
        validators=[
            genome_models.validate_genome_name,
            genome_models.validate_duplicate_name
        ],
        max_length=100,
        required=True,
    )

    genome_upload = forms.FileField(
        validators=[validate_fasta_file],
        help_text="Upload a fasta file containing the organism's genome sequence.",
    )

    cds_upload = forms.FileField(
        validators=[validate_coordinate_file],
        help_text=format_lazy(
            'Upload a file containing all CDS positions. This file should contain one line per CDS and each line '
            'should have strain (+ or -), start site (integer), and stop site (integer) separated by whitespace. '
            '<a href="{example_file}" download>Click here</a> to download an example file (Contains CDS for Lambda '
            'phage: <a href={lambda_genbank}">J02459.1</a>)',
            example_file=static('misc/J02459_cds.coords'),
            lambda_genbank='https://www.ncbi.nlm.nih.gov/nuccore/J02459.1/'
        )
    )

    table_choices = (
        (1, '1. The Standard Code'),
        (2, '2. The Vertebrate Mitochondrial Code'),
        (3, '3. The Yeast Mitochondrial Code'),
        (4, '4. The Mold, Protozoan, and Coelenterate Mitochondrial Code and the Mycoplasma / Spiroplasma Code'),
        (5, '5. The Invertebrate Mitochondrial Code'),
        (6, '6. The Ciliate, Dasycladacean and Hexamita Nuclear Code'),
        (9, '9. The Echinoderm and Flatworm Mitochondrial Code'),
        (10, '10. The Euplotid Nuclear Code'),
        (11, '11. The Bacterial, Archaeal and Plant Plastid Code'),
        (12, '12. The Alternative Yeast Nuclear Code'),
        (13, '13. The Ascidian Mitochondrial Code'),
        (14, '14. The Alternative Flatworm Mitochondrial Code'),
        (16, '16. Chlorophycean Mitochondrial Code'),
        (21, '21. Trematode Mitochondrial Code'),
        (22, '22. Scenedesmus obliquus Mitochondrial Code'),
        (23, '23. Thraustochytrium Mitochondrial Code'),
        (24, '24. Rhabdopleuridae Mitochondrial Code'),
        (25, '25. Candidate Division SR1 and Gracilibacteria Code'),
        (26, '26. Pachysolen tannophilus Nuclear Code'),
        (27, '27. Karyorelict Nuclear Code'),
        (28, '28. Condylostoma Nuclear Code'),
        (29, '29. Mesodinium Nuclear Code'),
        (30, '30. Peritrich Nuclear Code'),
        (31, '31. Blastocrithidia Nuclear Code'),
        (33, '33. Cephalodiscidae Mitochondrial UAA-Tyr Code'),
    )
    translation_table = forms.ChoiceField(choices=table_choices)

    organism = DynamicChoiceField(
        choices=lazy(
            get_set_of_used_organisms,
            tuple
        )
    )

    run_trnascan = forms.BooleanField(required=False)

    assign_to = forms.ModelChoiceField(
        queryset=user_choices,
        required=False,
    )

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        # User choices should be limited to those with permissions to change annotations
        self.fields['assign_to'].queryset = genome_views.get_annotation_editors()

    def clean(self):
        cleaned_data = super().clean()
        genome_file = cleaned_data.get("genome_upload")
        cds_file = cleaned_data.get("cds_upload")
        selected_table = cleaned_data.get("translation_table")

        try:
            if genome_file and cds_file:
                genome_fh = get_file_handle(genome_file, mode='r')
                cds_fh = get_file_handle(cds_file, mode='r')

                genome_rec = SeqIO.read(genome_fh, "fasta")

                for prot, cds in parse_prots_from_coords(cds_fh, genome_rec, selected_table):
                    pass

        except (TypeError, TranslationError) as e:
            raise ValidationError(e)

        return cleaned_data


# displayed in phage_upload.html
class Phage_Upload_Form(forms.Form):
    user_choices = User.objects.none()
    error_css_class = 'error'
    required_css_class = 'required'

    name = forms.CharField(
        validators=[
            genome_models.validate_genome_name,
            genome_models.validate_duplicate_name
        ],
        max_length=100,
        required=True,

    )

    upload = forms.FileField(
        validators=[validate_fasta_file],
        help_text='Must be a single fasta file!',
    )

    terminal_repeat = forms.IntegerField(
        required=True,
        initial=0,
    )

    error_choices = [
        ('Yes', 'Yes'),
        ('No', 'No')
    ]

    checkbox = forms.ChoiceField(choices=error_choices, validators=[genome_models.validate_phage_error])

    circularly_permuted = forms.BooleanField(
        # widget=CheckboxInput(attrs={'class': 'confirm_delete'}),
        required=False,
        initial=False
    )

    assign_to = forms.ModelChoiceField(
        queryset=user_choices,
        required=False,
    )

    def __init__(self, *args, **kwargs):
        super(Phage_Upload_Form, self).__init__(*args, **kwargs)
        # User choices should be limited to those with permissions to change annotations
        self.fields['assign_to'].queryset = genome_views.get_annotation_editors()

    def is_valid(self):
        # runs the checks for each field in Phage_Upload_Form
        is_valid = super(Phage_Upload_Form, self).is_valid()

        # checks if the upload information is good
        if self.cleaned_data['terminal_repeat'] != 0 and self.cleaned_data['circularly_permuted']:
            if 'checkbox' in self.errors:
                self.errors.pop('checkbox')
            self.errors.update(Conflict=': Can not have terminal repeat and be cicularly permuted')
            return False

        elif 'name' in self.cleaned_data and 'upload' in self.cleaned_data and self.cleaned_data['terminal_repeat'] > 0:
            # parse genome
            file = get_file_handle(self.cleaned_data['upload'], mode='r')
            genome = SeqIO.read(file, 'fasta').seq.__str__().upper()
            # Ensure both ends of phage genome match at given length
            if genome[:self.cleaned_data['terminal_repeat']] != genome[-self.cleaned_data['terminal_repeat']:]:
                self.errors.update(terminal_repeat=': Incorrect terminal repeat size')
                self.errors.pop('checkbox')
                return False
            else:
                self.errors.pop('checkbox')
                return True

        elif 'name' in self.cleaned_data and 'upload' in self.cleaned_data and self.cleaned_data['terminal_repeat'] == 0:
            self.errors.pop('checkbox')
            return True

        elif 'upload' in self.cleaned_data and 'checkbox' in self.cleaned_data and \
                self.cleaned_data['checkbox'] == 'Yes' and self.cleaned_data['terminal_repeat'] > 0:
            # parse genome
            file = get_file_handle(self.cleaned_data['upload'], mode='r')
            genome = SeqIO.read(file, 'fasta').seq.__str__().upper()
            # Ensure both ends of phage genome match at given length
            if genome[:self.cleaned_data['terminal_repeat']] != genome[-self.cleaned_data['terminal_repeat']:]:
                self.errors.update(terminal_repeat=': Incorrect terminal repeat size')
                self.errors.pop('name')
                return False
            else:
                self.errors.pop('name')
            return True

        elif 'upload' in self.cleaned_data and 'checkbox' in self.cleaned_data and \
                self.cleaned_data['checkbox'] == 'Yes' and self.cleaned_data['terminal_repeat'] == 0:
            self.errors.pop('name')
            return True
        elif 'upload' in self.cleaned_data and 'checkbox' in self.cleaned_data and \
                self.cleaned_data['checkbox'] == 'Yes' and 'circularly_permuted' in self.cleaned_data:
            self.errors.pop('name')
            return True
        else:
            return False


class Genome_Delete(forms.Form):
    genome = forms.ModelMultipleChoiceField(
        queryset=genome_models.Genome.objects.all()
    )

    def __init__(self, *args, **kwargs):
        super(Genome_Delete, self).__init__(*args, **kwargs)
        self.helper = CrispyHorizontalFormHelper()
        self.helper.form_tag = False
        self.helper.disable_csrf = True


class Annotations_Delete(forms.Form):
    annotations = forms.ModelMultipleChoiceField(
        queryset=genome_models.Annotation.objects.all(),
        required=False
    )

    def __init__(self, *args, **kwargs):
        annotations_to_delete = kwargs.pop('annotations_to_delete', None)
        super(Annotations_Delete, self).__init__(*args, **kwargs)
        if annotations_to_delete is not None:
            self.fields['annotations'].queryset = annotations_to_delete


class Confirm_Delete(forms.ModelForm):
    confirm_delete = forms.BooleanField(
        widget=CheckboxInput(attrs={'class': 'confirm_delete'}),
        required=False,
        initial=False
    )

    class Meta:
        model = genome_models.Annotation
        fields = (
            # 'annotation',
        )


class Upload_Annotation(forms.Form):
    upload = forms.FileField(
        validators=[validate_excel_file],
        required=True,
        help_text='Must be a excel file'
    )


class Confirm_Upload_Annotation(forms.Form):

    def __init__(self, *args, **kwargs):
        super(Confirm_Upload_Annotation, self).__init__(*args, **kwargs)
        # self.initial['user_annotation'] = self.kwargs['user_annotation']
        # self.initial['user_public_note'] = self.kwargs['user_public_note']
        # self.initial['user_private_note'] = self.kwargs['user_private_note']
        # self.initial['user_flag'] = self.kwargs['user_flag']

    # select = forms.BooleanField(
    #
    # )
    choices = [('New', 'New'),
               ('Original', 'Original'),
               ('Custom', 'Custom')]

    select_annotation = forms.ChoiceField(
        widget=forms.RadioSelect,
        required=False,
        choices=choices
    )
    select_public_note = forms.ChoiceField(
        widget=forms.RadioSelect,
        required=False,
        choices=choices
    )
    select_private_note = forms.ChoiceField(
        widget=forms.RadioSelect,
        required=False,
        choices=choices
    )
    select_flag = forms.ChoiceField(
        widget=forms.RadioSelect,
        required=False,
        choices=choices
    )

    custom_annotation = forms.CharField(
        widget=forms.TextInput,
        required=False
    )
    custom_public_note = forms.CharField(
        widget=forms.TextInput,
        required=False
    )
    custom_private_note = forms.CharField(
        widget=forms.TextInput,
        required=False
    )
    custom_flag = forms.ChoiceField(
        choices=genome_models.Annotation.flag_options,
        required=False
    )

    user_annotation = forms.CharField(
        # hidden=True
        required=False
        # widget=forms.TextInput
    )
    user_flag = forms.ChoiceField(
        choices=genome_models.Annotation.flag_options,
        # hidden=True
        required=False
    )
    user_public_note = forms.CharField(
        # hidden=True
        required=False
    )
    user_private_note = forms.CharField(
        # hidden=True
        required=False
    )

    db_annotation = forms.CharField(
        # hidden=True
        required=False
        # widget=forms.TextInput
    )
    db_flag = forms.ChoiceField(
        choices=genome_models.Annotation.flag_options,
        # hidden=True
        required=False
    )
    db_public_note = forms.CharField(
        # hidden=True
        required=False
    )
    db_private_note = forms.CharField(
        # hidden=True
        required=False
    )
    db_pk = forms.IntegerField(
        required=True
    )

    class Meta:
        readonly = (
            'db_pk',
            'db_private_note',
            'db_public_note',
            'db_flag',
            'db_annotation',
            'user_private_note',
            'user_public_note',
            'user_flag',
            'user_annotation',
        )
        # model = genome_models.Annotation
        # fields = (
        #     # 'annotation'
        #     # 'private_notes'
        #     # 'public_notes'
        #     # 'flag'
        # )
