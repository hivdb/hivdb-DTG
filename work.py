#!/usr/bin/env python3.9
from pathlib import Path
from collections import defaultdict
from operator import itemgetter
# from pprint import pprint
from preset.file_formats import load_tsv
from preset.file_formats import load_yaml
from preset.file_formats import dump_csv
from preset.file_formats import load_csv
from preset.mutation_literal import parse_mut_str_list
from preset.mutation_literal import parse_mut_str
from preset.mutation_literal import bind_mut_str_list
from preset.table import group_records_by
from preset.mutation_ops import combine_mutation_pattern
from preset.report import SummaryReport
from collections import OrderedDict
from preset.geno_pheno import get_geno_pheno
from datetime import datetime


WS = Path(__file__).resolve().parent
DB = Path(__file__).resolve().parent / 'database' / 'Archive'

DTG_MAJOR_POSITION = [
    118,
    263,
    148,
    155,
]


def collect_reference_info(table_raw):

    report = SummaryReport()

    report.num_paper = len(set([
        (i['Author'], i['RefYear'], i['Journal'], i['MedlineID'])
        for i in table_raw
    ]))

    report.num_isolate = len(table_raw)
    report.num_patient = len(set([
        i['PtID']
        for i in table_raw
    ]))

    report.num_isolate_complete_mutation = len([
        i
        for i in table_raw
        if i['CompleteMutationListAvailable'] == 'Yes'
    ])
    report.num_patient_complete_mutation = len(set([
        i['PtID']
        for i in table_raw
        if i['CompleteMutationListAvailable'] == 'Yes'
    ]))

    report.num_isolate_partial_mutation = len([
        i
        for i in table_raw
        if i['CompleteMutationListAvailable'] == 'No'
    ])
    report.num_patient_partial_mutation = len(set([
        i['PtID']
        for i in table_raw
        if i['CompleteMutationListAvailable'] == 'No'
    ]))

    report.num_paper_partial_mutation = len(set([
        (i['Author'], i['RefYear'], i['Journal'], i['MedlineID'])
        for i in table_raw
        if i['CompleteMutationListAvailable'] == 'No'
    ]))

    return report


def collect_rx_info(table_raw, table_rx):

    pt_rx = {}

    for ptid, rx_list in group_records_by(table_rx, 'PtID').items():
        [
            i.update({
                'StopDate': datetime.strptime(i['StopDate'], "%Y-%m-%d")})
            for i in rx_list
        ]
        pt_rx[ptid] = rx_list

    for i in table_raw:
        ptid = i['PtID']
        isolateDate = datetime.strptime(i['IsolateDate'], "%Y-%m-%d")
        rx_list = pt_rx[ptid]
        rx_list = [
            i
            for i in rx_list
            if i['StopDate'] <= isolateDate
        ]
        rx_list.sort(key=itemgetter('StopDate'))
        i['Rx'] = '\n'.join([
            i['RegimenName']
            for i in rx_list
        ])


def collect_reference_meta(table_raw, table_meta):

    report = SummaryReport()

    continents = {
        (i['Author'], i['RefYear']): i['Continent']
        for i in table_meta
    }

    categories = {
        (i['Author'], i['RefYear']): i['Category']
        for i in table_meta
    }

    for i in table_raw:
        i['continent'] = continents[(i['Author'], i['RefYear'])]
        i['category'] = categories[(i['Author'], i['RefYear'])]

    for cont, cont_list in group_records_by(table_raw, 'continent').items():
        setattr(report, cont, len(cont_list))

    for cat, cat_list in group_records_by(table_raw, 'category').items():
        setattr(report, cat, len(cat_list))

    return report


def load_main_drm(path):
    main_drm_list = {}
    for category, drm_list in load_yaml(path).items():
        drm_list = [
            parse_mut_str(i)
            for i in drm_list
        ]
        main_drm_list[category] = drm_list

    main_drm_list['nonpolymorphic'] = (
        main_drm_list['SDRM'] + main_drm_list['nonpolymorphic'])

    return main_drm_list


def get_isolate_nonpoly_drm(rec, main_drms):

    complete_muts = rec['CompMutList']
    complete_muts = parse_mut_str_list(complete_muts)

    poly_drms = []
    nonpoly_drms = []
    for mut in complete_muts:
        nonpoly_drms.extend(
            find_drm_from_mutation(mut, main_drms['nonpolymorphic']))

        poly_drms.extend(
            find_drm_from_mutation(mut, main_drms['polymorphic']))

    rec['nonpoly_drms'] = bind_mut_str_list(nonpoly_drms)
    rec['poly_drms'] = bind_mut_str_list(poly_drms)

    drms = nonpoly_drms + poly_drms
    drms.sort(key=itemgetter('pos'))
    rec['drms'] = bind_mut_str_list(drms)

    return rec


def find_drm_from_mutation(mut, drm_list):

    matched_drm = []
    for drm in drm_list:

        if mut['pos'] != drm['pos']:
            continue

        if mut['ref'] != drm['ref']:
            print(f"Warning, ref erro, {mut['pos']}")

        shared_mut = set(mut['mut']) & set(drm['mut'])

        if not (shared_mut):
            continue

        matched_drm.append({
            'ref': mut['ref'],
            'pos': mut['pos'],
            'mut': ''.join(sorted(list(shared_mut)))
        })

    return matched_drm


def get_patient_unique_mutation(pt_iso_list, one_per_person=True):

    selected_pattern = []

    for drm_pattern, rec_list in group_records_by(pt_iso_list, 'drms').items():
        selected_pattern.append(rec_list[0])

    print(
        f"Get unique mutation pattern from PtID {pt_iso_list[0]['PtID']}:"
        f" {len(selected_pattern)}/{len(pt_iso_list)}")

    if one_per_person:
        selected_pattern.sort(key=lambda x: len(parse_mut_str_list(x['drms'])))
        return selected_pattern[-1:]
    else:
        return selected_pattern


def count_drm_pattern_num_isolate(isolates):

    drm_pattern = [
        {
            'drm_pattern': k,
            'num_isolate': len(v)
        }
        for k, v in group_records_by(isolates, 'drms').items()
    ]
    return drm_pattern


def group_by_major_pos(drm_pattern):

    for rec in drm_pattern:
        drm_list = parse_mut_str_list(rec['drm_pattern'])

        major_pos_list = tuple(sorted(list(set([
            i['pos']
            for i in drm_list
        ]) & set(DTG_MAJOR_POSITION))))

        rec['major_pos_list'] = major_pos_list
        rec['num_major_pos'] = len(major_pos_list)

        pos_list = tuple(sorted([
            m['pos']
            for m in drm_list
        ]))
        rec['pos_list'] = pos_list
        rec['num_pos'] = len(pos_list)

    return group_records_by(drm_pattern, 'major_pos_list')


def merge_drm_pattern_same_pos_list(drm_pattern):
    # TODO define the difference between combine and merge

    drm_pattern_group = defaultdict(list)

    for drm_pos_list, patterns in drm_pattern:
        patterns = combine_drm_by_same_pos_list(patterns)

        drm_pattern_group[drm_pos_list] = patterns

    return drm_pattern_group


def combine_drm_by_same_pos_list(patterns):

    results = []
    for pos_list, pos_patterns in group_records_by(patterns, 'pos_list').items():
        new_pattern = get_combined_pattern(pos_list, pos_patterns)

        if not has_major_pos_combined(new_pattern):
            # pprint(new_pattern)
            results.append(new_pattern)
        else:
            results.extend(pos_patterns)

    return results


def get_combined_pattern(pos_list, pos_patterns):
    new_pattern = {
        'pos_list': pos_list,
        'num_pos': len(pos_list),
    }

    new_pattern['num_isolate'] = sum([
        i['num_isolate']
        for i in pos_patterns
    ])
    new_pattern['major_pos_list'] = pos_patterns[0]['major_pos_list']
    new_pattern['num_major_pos'] = pos_patterns[0]['num_major_pos']

    mut_pattern_list = [
        i['drm_pattern']
        for i in pos_patterns
    ]

    new_pattern['drm_pattern'] = combine_mutation_pattern(mut_pattern_list)

    return new_pattern


def has_major_pos_combined(new_pattern):
    new_mut_list = parse_mut_str_list(new_pattern['drm_pattern'])

    for mut in new_mut_list:
        if ((mut['pos'] in new_pattern['major_pos_list'])
                and len(mut['mut']) > 1):
            return True

    return False


def prepare_report(drm_pattern, pos_order):

    report = []

    drm_pattern = [
        k
        for i, j in drm_pattern.items()
        for k in j
    ]

    w_major_pos = [
        i
        for i in drm_pattern
        if i['num_major_pos']
    ]

    prepare_w_major_pos(report, w_major_pos, pos_order)

    wo_major_pos = [
        i
        for i in drm_pattern
        if not i['num_major_pos']
    ]

    prepare_wo_major_pos(report, wo_major_pos, pos_order)

    return report


def prepare_w_major_pos(report, w_major_pos, pos_order):

    w_major_pos = sorted(list(
        group_records_by(w_major_pos, 'num_major_pos').items()),
        key=lambda x: x[0]
        )

    for _, num_major_pos_list in w_major_pos:
        num_major_pos_list = list(
            group_records_by(num_major_pos_list, 'major_pos_list').items())

        num_major_pos_list.sort(
            key=lambda x: min([
                DTG_MAJOR_POSITION.index(i)
                for i in x[0]
            ]))

        for major_pos, major_pos_list in num_major_pos_list:
            major_pos_list.sort(key=lambda x: len(x['pos_list']))

            [
                rec.update({
                    'drm_pattern': sort_pattern_pos(
                        rec['drm_pattern'],
                        rec['major_pos_list'])
                })
                for rec in major_pos_list
            ]

            sub_report = [{
                    'drm_pattern': rec['drm_pattern'],
                    'num_isolate': rec['num_isolate'],
                    'major_pos_list': ', '.join([f'{p}' for p in major_pos]),
                    'num_pos': rec['num_pos'],
                    'pos_list': join_list_str(rec['pos_list'], gap=True),
                    # TODO: a column contains a list of integers, how to sort them in text format?
                }
                for rec in major_pos_list
            ]

            show_sub_report(report, sub_report)


def join_list_str(a_list=[], delimiter=',', gap=False):
    if gap:
        delimiter = f'{delimiter} '

    return delimiter.join([str(i) for i in a_list])


def split_list_str(a_list_str, delimiter=',', vtype=str):

    return [
        vtype(i.strip())
        for i in a_list_str.split(delimiter)]


def prepare_wo_major_pos(report, wo_major_pos, pos_order):

    wo_major_pos = sorted(list(
        group_records_by(wo_major_pos, 'num_pos').items()),
        key=lambda x: x[0]
        )

    for _, num_pos_list in wo_major_pos:

        [
            rec.update({
                'drm_pattern': sort_pattern_pos(
                    rec['drm_pattern'], pos_order)
            })
            for rec in num_pos_list
        ]

        sub_report = [{
                'drm_pattern': rec['drm_pattern'],
                'num_isolate': rec['num_isolate'],
                'major_pos_list': '',
                'num_pos': rec['num_pos'],
                'pos_list': join_list_str(rec['pos_list'], gap=True),
            }
            for rec in num_pos_list
        ]

        show_sub_report(report, sub_report)


def sort_pattern_pos(pattern, pos_order=[]):
    pos_order = list(pos_order)

    pattern = parse_mut_str_list(pattern)
    pattern.sort(key=itemgetter('pos'))
    pattern.sort(
            key=lambda x:
            (pos_order + [x['pos']]).index(x['pos'])
        )

    return bind_mut_str_list(pattern, join_str=' + ')


def show_sub_report(report, sub_report):
    report.extend(sub_report)
    # report.append({
    #     'drm_pattern': '----'
    # })


def get_pos_order(main_drm_list):

    pos_list = [
        i['pos']
        for i in main_drm_list['nonpolymorphic']
        if i['pos'] not in DTG_MAJOR_POSITION
    ]
    pos_list = DTG_MAJOR_POSITION + pos_list + [
        i['pos']
        for i in main_drm_list['polymorphic']
    ]

    return pos_list


def geno_analysis(geno_file, folder, meta, rx):
    table_raw = load_tsv(geno_file)
    main_drm_list = load_main_drm(WS / 'mutations.yml')
    table_meta = load_csv(meta)
    table_rx_history = load_csv(rx)

    table_raw = [
        get_isolate_nonpoly_drm(i, main_drm_list)
        for i in table_raw
    ]

    dump_csv(folder / 'find_drm.csv', table_raw)

    table_nonpoly = [
        i
        for i in table_raw
        if i['nonpoly_drms']
    ]

    print(f'Including isolates with nonpoly DRMs: {len(table_nonpoly)}')
    print('*' * 20)

    report_raw = collect_reference_info(table_raw)
    report_nonpoly = collect_reference_info(table_nonpoly)

    isolates = []
    pt_multiple_isolate = 0
    for _, pt_iso_list in group_records_by(table_nonpoly, 'PtID').items():
        if len(pt_iso_list) == 1:
            isolates.extend(pt_iso_list)
        else:
            pt_multiple_isolate += 1
            isolates.extend(get_patient_unique_mutation(pt_iso_list, True))

    report_nonpoly.patient_with_multiple_isolate = pt_multiple_isolate

    report_meta = collect_reference_meta(isolates, table_meta)
    report_rx = collect_rx_info(isolates, table_rx_history)

    report = report_raw.table() + \
        report_nonpoly.table('W_DRM') + \
        report_meta.table()

    dump_csv(folder / 'reference_info.csv', report)
    print(f'Isolates after removing duplication: {len(isolates)}')
    print('*' * 20)

    dump_csv(folder / 'unique_isolates.csv', isolates)

    drm_pattern = count_drm_pattern_num_isolate(isolates)

    drm_pattern = group_by_major_pos(drm_pattern).items()

    drm_pattern = merge_drm_pattern_same_pos_list(drm_pattern)

    pos_order = get_pos_order(main_drm_list)
    report = prepare_report(drm_pattern, pos_order)

    report.sort(
        key=lambda x:
            [
                x['major_pos_list'],
                x['num_pos'],
                -x['num_isolate']
            ] +
            sorted(split_list_str(x['pos_list'], vtype=int))
        )
    # TIP: fix the order of sorted list, if the num isolate is the same

    dump_csv(folder / 'table 1.csv', report)


def work():
    # get_geno_pheno(
    #     DB / 'May 30, 2023', DB / 'May 30, 2023' / 'geno_rx_pheno.csv')
    geno_analysis(
        geno_file=DB / 'Jun 28, 2023' / 'geno-rx.dataset.tsv',
        folder=DB / 'Jun 28, 2023',
        meta=DB / 'Jun 28, 2023' / 'paper_meta.csv',
        rx=DB / 'Jun 28, 2023' / 'tblRxHistory.csv')


if __name__ == '__main__':
    work()
    print('Finish.')
