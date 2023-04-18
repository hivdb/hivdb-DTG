#!/usr/bin/env python3.9
import csv
from pathlib import Path
from collections import defaultdict
import re
from operator import itemgetter
from pprint import pprint


WS = Path(__file__).resolve().parent / 'database'


DTG_MAJOR_POSITION = [
    118,
    263,
    148,
    155,
]

DTG_MINOR_POSITION = [
    138,
    232,
    230,
    51,
]


def load_tsv(file_path):
    records = []
    with open(file_path, encoding='utf-8-sig') as fd:
        for record in csv.DictReader(fd, delimiter='\t'):
            records.append(record)

    return records


def load_csv(file_path):
    records = []

    with open(file_path, encoding='utf-8-sig') as fd:
        for record in csv.DictReader(fd):
            records.append(record)
    return records


def dump_csv(file_path, records, headers=None):
    if not records:
        return
    if not headers:
        headers = []
        for rec in records:
            for key in rec.keys():
                if key not in headers:
                    headers.append(key)

    final_records = []
    for rec in records:
        new_rec = {}
        for key in headers:
            new_rec[key] = rec.get(key, '')
        for k, v in rec.items():
            if type(v) == bool:
                new_rec[k] = 'yes' if v else 'no'
        final_records.append(new_rec)

    file_path = Path(file_path)
    file_path.parent.mkdir(exist_ok=True, parents=True)

    with open(file_path, 'w', encoding='utf-8-sig') as fd:
        writer = csv.DictWriter(fd, fieldnames=headers)
        writer.writeheader()
        writer.writerows(final_records)


def dump_csv_raw(file_path, table):

    with open(file_path, 'w', encoding='utf-8-sig') as fd:
        writer = csv.writer(fd)
        writer.writerows(table)


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

    for k, v in group_result.items():
        yield k, v


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


def isolate_has_major_mut(rec):

    if rec['INIMajorDRMs']:
        return True

    minor_muts = parse_mut_str_list(rec['INIMinorDRMs'])
    if any([i['pos'] in DTG_MINOR_POSITION for i in minor_muts]):
        return True

    return False


def is_subpattern_mut_str(mut_str, patterns):
    mut_list = parse_mut_str_list(mut_str)

    for p in patterns:
        if mut_str in p:
            print(f'Found sub or same pattern ({mut_str}) in ({p})')
            return True

        pattern_mut_list = parse_mut_str_list(p)

        check = [
            any([
                is_sub_mutation(mut, p)
                for p in pattern_mut_list
            ])
            for mut in mut_list
        ]

        if all(check):
            print('Sub pattern', mut_str, 'of', p)
            return True

    return False


def is_sub_mutation(mut_a, mut_b):
    if mut_a['pos'] != mut_b['pos']:
        return False
    if not (set(mut_a['mut']).issubset(set(mut_b['mut']))):
        return False
    return True


def prorcess_mixture(record, method):

    major_mut = process_mixture_mut(record, 'INIMajorDRMs', method)
    minor_mut = process_mixture_mut(record, 'INIMinorDRMs', method)

    minor_mut = [
        i
        for i in minor_mut
        if i['pos'] in DTG_MINOR_POSITION
    ]

    drms = major_mut + minor_mut
    drms.sort(key=itemgetter('pos'))

    record['DRMs'] = bind_mut_str_list(drms)
    return record


def process_mixture_mut(record, column, method='if_has_ref_ignore'):
    drm_mut = record[column]

    drm_mut = parse_mut_str_list(drm_mut)

    if method == 'if_has_ref_ignore':
        drm_mut = [
            i
            for i in drm_mut
            if i['ref'] not in i['mut']
        ]
    elif method == 'if_has_ref_ignore_ref':
        [
            i.update({
                'mut': i['mut'].replace(i['ref'], '')
            })
            for i in drm_mut
        ]

    record[column] = bind_mut_str_list(drm_mut)

    return drm_mut


def merge_invivo_muts(pt_iso_list):
    pt_iso_list = [
        prorcess_mixture(i, 'if_has_ref_ignore')
        for i in pt_iso_list
    ]

    pt_iso_list.sort(key=lambda x: len(x['DRMs']))

    merged_by_pattern = {}

    for rec in pt_iso_list[::-1]:
        drms = rec['DRMs']
        if drms in list(merged_by_pattern.keys()):
            print(f'Found sub or same pattern ({drms}) in ('
                  f'{list(merged_by_pattern.keys())})')
            continue
        # if is_subpattern_mut_str(drms, list(merged_by_pattern.keys())):
        #     continue

        merged_by_pattern[drms] = rec

    merged_pt_iso_list = list(merged_by_pattern.values())

    if len(pt_iso_list) > len(merged_pt_iso_list):
        print(
            f"Merge isolates with duplicated"
            f" mutations from patient #{pt_iso_list[0]['PtID']}")

    return merged_pt_iso_list


def adapt_drm_pos_mut(isolates):

    for i in isolates:
        drms = []
        for m in parse_mut_str_list(i['DRMs']):
            # TODO
            if m['pos'] not in DTG_MAJOR_POSITION:
                drms.append(m)
                continue

            m['mut'] = m['mut'].replace(m['ref'], '')

            if m['pos'] == 118:
                m['mut'] = m['mut'].replace('S', '')

            drms.append(m)

        i['DRMs'] = bind_mut_str_list(drms)

    return isolates


def group_by_drm_pattern(isolates):

    drm_pattern = [
        {
            'pattern': k,
            'num_isolate': len(v)
        }
        for k, v in group_records_by(isolates, 'DRMs')
    ]
    return drm_pattern


def group_by_major_pos(drm_pattern):

    for rec in drm_pattern:
        mut_list = parse_mut_str_list(rec['pattern'])

        drm_pos = tuple(set([
            i['pos']
            for i in mut_list
        ]) & set(DTG_MAJOR_POSITION))

        rec['major_pos_list'] = drm_pos
        rec['num_major_pos'] = len(drm_pos)
        rec['pos_list'] = tuple(sorted([
            m['pos']
            for m in mut_list
        ]))

    return group_records_by(drm_pattern, 'major_pos_list')


def combine_drm_by_pos(patterns):

    results = []
    for pos_list, pos_patterns in group_records_by(patterns, 'pos_list'):
        new_pattern = get_combined_pattern(pos_list, pos_patterns)

        if not has_major_pos_combined(new_pattern):
            # pprint(new_pattern)
            results.append(new_pattern)
        else:
            results.extend(pos_patterns)

    return results


def has_major_pos_combined(new_pattern):
    new_mut_list = parse_mut_str_list(new_pattern['pattern'])

    for mut in new_mut_list:
        if ((mut['pos'] in new_pattern['major_pos_list'])
                and len(mut['mut']) > 1):
            return True

    return False


def get_combined_pattern(pos_list, pos_patterns):
    new_pattern = {'pos_list': pos_list}

    new_pattern['num_isolate'] = sum([
        i['num_isolate']
        for i in pos_patterns
    ])
    new_pattern['major_pos_list'] = pos_patterns[0]['major_pos_list']
    new_pattern['num_major_pos'] = pos_patterns[0]['num_major_pos']

    mut_list = [
        parse_mut_str_list(i['pattern'])
        for i in pos_patterns
    ]
    new_mut_list = []
    for pos in pos_list:
        ref = list(set([
            j['ref']
            for i in mut_list
            for j in i
            if j['pos'] == pos
        ]))[0]
        muts = [
            j['mut']
            for i in mut_list
            for j in i
            if j['pos'] == pos
        ]
        new_mut_list.append({
            'ref': ref,
            'pos': pos,
            'mut': ''.join(sorted(set(muts)))
        })
    new_pattern['pattern'] = bind_mut_str_list(new_mut_list)

    return new_pattern


def merge_drm_pattern(drm_pattern):

    drm_pattern_group = defaultdict(list)

    for drm_pos_list, patterns in drm_pattern:
        patterns = combine_drm_by_pos(patterns)

        drm_pattern_group[drm_pos_list] = patterns

    return drm_pattern_group


def sort_pattern_pos(pattern):

    pattern = parse_mut_str_list(pattern)
    pattern.sort(
            key=lambda x:
            (
                DTG_MAJOR_POSITION + DTG_MINOR_POSITION + [x['pos']]
            ).index(x['pos'])
        )

    return bind_mut_str_list(pattern)


def prepare_report(drm_pattern):

    report = []

    for pos in DTG_MAJOR_POSITION:
        for drm_pos_list, patterns in drm_pattern.items():
            if pos not in drm_pos_list:
                continue

            [
                rec.update({'pattern': sort_pattern_pos(rec['pattern'])})
                for rec in patterns
            ]
            patterns.sort(key=lambda x: len(x['pos_list']))

            sub_report = [{
                    'pattern': rec['pattern'],
                    'num_isolate': rec['num_isolate'],
                }
                for rec in patterns
                if not rec.get('processed')
            ]

            [
                rec.update({'processed': True})
                for rec in patterns
            ]

            if sub_report:
                show_sub_report(report, sub_report)

    for drm_pos_list, patterns in drm_pattern.items():
        for rec in patterns:
            if rec.get('processed'):
                continue
            sub_report = [{
                'pattern': rec['pattern'],
                'num_isolate': rec['num_isolate'],
            }]

            show_sub_report(report, sub_report)

    return report


def show_sub_report(report, sub_report):
    report.extend(sub_report)
    report.append({
        'pattern': '----'
    })


def work():
    table_raw = load_tsv(WS / 'geno-rx.dataset.tsv')
    table_major = [
        i
        for i in table_raw
        if isolate_has_major_mut(i)
    ]
    print(f'Including isolates with major DRMs: {len(table_major)}')
    print('*' * 20)

    isolates = []

    for _, pt_iso_list in group_records_by(table_major, 'PtID'):
        if len(pt_iso_list) == 1:
            pt_iso_list = [
                prorcess_mixture(i, 'if_has_ref_ignore_ref')
                for i in pt_iso_list
            ]
            isolates.extend(pt_iso_list)
        else:
            isolates.extend(merge_invivo_muts(pt_iso_list))

    print(f'Isolates after removing duplication: {len(isolates)}')
    print('*' * 20)
    dump_csv(WS / 'unique_isolates.csv', isolates)

    isolates = adapt_drm_pos_mut(isolates)
    drm_pattern = group_by_drm_pattern(isolates)

    drm_pattern = group_by_major_pos(drm_pattern)

    drm_pattern = merge_drm_pattern(drm_pattern)

    dump_csv(WS / 'pathway.csv', prepare_report(drm_pattern))

    print('Finish.')


if __name__ == '__main__':
    work()
