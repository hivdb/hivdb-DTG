from preset.mutation_literal import parse_mut_str_list
from preset.mutation_literal import bind_mut_str_list
from preset.table import group_records_by


def combine_mutation_pattern(mut_pattern_list):
    mut_list = [
        j
        for i in mut_pattern_list
        for j in parse_mut_str_list(i)
    ]

    new_mut_list = []
    for pos, pos_list in group_records_by(mut_list, 'pos').items():
        ref = list(set([
            i['ref']
            for i in pos_list
        ]))

        if len(ref) > 1:
            print(f'warning, ref duplicated, {ref}')

        ref = ref[0]

        muts = list(set(
            j
            for i in pos_list
            for j in i['mut']
        ))

        new_mut_list.append({
            'ref': ref,
            'pos': pos,
            'mut': ''.join(sorted(muts))
        })

    return bind_mut_str_list(new_mut_list)
