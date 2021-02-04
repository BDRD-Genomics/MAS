import subprocess

from django import template
from django.contrib.auth.models import User
from django.conf import settings
from django.template.defaultfilters import stringfilter

from result_viewer.models import Annotation

register = template.Library()

@register.filter
def absolute(value):
    return abs(float(value))


@register.filter
def subtract(value, arg):
    return int(value) - int(arg)


@register.filter(name='accession')
def accession(integer):
    chars = '0123456789ABCDEFGHIJKLMNOPQRSTUVWXYZ'
    integer = abs(integer)
    result = ''
    padding = []

    while integer > 0:
        integer, remainder = divmod(integer, 36)
        padding.insert(0, chars[remainder])
    if len(padding) < 5:
        while len(padding) < 5:
            padding.insert(0, '0')
    for i in padding:
        result = result + i
    return result


@register.filter(name='add_sequence_breaks')
@stringfilter
def add_sequence_breaks(sequence):
    """
    add new line character every 50 characters
    :param string:
    :return:
    """
    # biosequence = Seq.Seq(sequence)
    # splitbiosequence = SeqIO.FastaIO.FastaWriter.
    # word.replace('_', '<wbr>_<wbr>')
    list = []
    newline = '\n'
    count = 0
    for i in sequence:
        if count == 50:
            list.append(i)
            list.append(newline)
            count = 0
        else:
            list.append(i)
            count = count + 1

    newstring = ""
    for i in list:
        newstring = newstring + i

    return newstring


@register.filter(name='has_permission')
def has_group(user, permission_name):
    return user.has_perm(permission_name)


@register.filter(name='get_flag')
def get_flag(flag):
    return Annotation.flag_options[flag][1]


@register.simple_tag()
def get_annotation_flags():
    return Annotation.flag_options


@register.simple_tag()
def get_bioinformatics_users():
    return User.objects.filter(groups__name='Bioinformatics')


@register.simple_tag()
def get_version():
    try:
        result = subprocess.run(
            ['git --git-dir={} describe --tags'.format(settings.GIT_DIR)],
            cwd=settings.BASE_DIR,
            stdout=subprocess.PIPE,
            timeout=2,
            shell=True
        )
        return result.stdout.decode().strip()

    except TimeoutError:
        return ''
