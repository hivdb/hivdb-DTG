import re


# TODO, support mutation separator "/"

def parse_mut_str_list(mut_str):

    mut_str_list = [i.strip() for i in mut_str.split(',') if i.strip()]

    return [
        parse_mut_str(i)
        for i in mut_str_list
    ]


def parse_mut_str(mut_str):
    match = re.match(r'^([A-Z])(\d+)([A-Z\_\*\-]*)$', mut_str).groups()
    match = list(match)
    match[1] = int(match[1])

    return dict(zip(['ref', 'pos', 'mut'], match))


def bind_mut_str_list(mutations):

    return ', '.join([
        bind_mut_str(i)
        for i in mutations
    ])


def bind_mut_str(mutation):

    return f"{mutation['ref']}{mutation['pos']}{mutation['mut']}"
