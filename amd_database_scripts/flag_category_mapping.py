'''
Flags are stored in the database as an integer. This file lists the acceptable categories of flags and defines the
integer <-> text mapping for each one.
'''

# If you add to this list, add to the end! indices must be preserved
FLAGS = ['GREEN', 'YELLOW', 'RED', 'REVIEW NAME', 'N/A', 'ORANGE', 'ENDOLYSIN']


def convert_flag_text_to_int(flag_str):
    flag_str = flag_str.upper()
    try:
        return FLAGS.index(flag_str)
    except ValueError:
        msg = '"{}" is not an accepted flag. Accepted flags are: {}'.format(flag_str, FLAGS)
        raise ValueError(msg)


def convert_flag_int_to_text(flag_int):
    try:
        return FLAGS[flag_int]
    except IndexError:
        msg = '{} is not an accepted flag index. Accepted flag indicies are: {}'.format(flag_int, list(range(len(FLAGS))))
        raise ValueError(msg)