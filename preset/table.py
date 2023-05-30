from collections import defaultdict


def group_records_by(records, group_key_list):

    if not isinstance(group_key_list, list):
        group_key_list = [group_key_list]

    group_result = defaultdict(list)
    for r in records:
        primary_key = {}
        for key in group_key_list:
            primary_key[key] = r[key]

        if len(primary_key) == 1:
            primary_key = list(primary_key.values())[0]
        else:
            primary_key = tuple(sorted(primary_key.items()))

        group_result[primary_key].append(r)

    return group_result
