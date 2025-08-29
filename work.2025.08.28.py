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
from preset.mutation_literal import bind_mut_str
from preset.table import group_records_by
from preset.mutation_ops import combine_mutation_pattern
from preset.report import SummaryReport
# from preset.geno_pheno import get_geno_pheno
from datetime import datetime
from preset.statistics import calc_contigency_table
from preset.statistics import calc_jaccard
from preset.statistics import calc_spearman
from itertools import combinations
from preset.statistics import calc_holm_Bonferroni_multiply_way
import igraph as ig
from pyvis.network import Network
from dataclasses import dataclass


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
        (i['Author'], i['RefYear']): i['Countries']
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
        main_drm_list.get('SDRM', []) + main_drm_list['nonpolymorphic'])

    return main_drm_list


def get_isolate_nonpoly_drm(
        rec, main_drms,
        no_ref_mixture=True, no_mixture=False):

    complete_muts = rec['CompMutList']
    complete_muts = parse_mut_str_list(complete_muts)

    complete_muts = [
        i
        for i in complete_muts
        if 'X' not in i['mut']
    ]

    if no_ref_mixture:
        [
            i.update({'mut': i['mut'].replace(i['ref'], '')})
            for i in complete_muts
        ]

    if no_mixture:
        complete_muts = [
            i
            for i in complete_muts
            if len(i['mut']) <= 1
        ]
    else:
        mixture_muts = [
            i
            for i in complete_muts
            if len(i['mut']) > 1
        ]
        if len(mixture_muts):
            print('Contain mixture', bind_mut_str_list(mixture_muts))
            print(bind_mut_str_list(complete_muts))

    rec['CompMutList'] = bind_mut_str_list(complete_muts)

    drms = []
    for mut in complete_muts:
        drms.extend(
            find_drm_from_mutation(mut, main_drms['drms']))

    drms.sort(key=itemgetter('pos'))
    rec['drms'] = bind_mut_str_list(drms)

    if not no_mixture:
        rec['has_mixture'] = any(
            len(i['mut']) > 1
            for i in drms
        )

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
            # 'mut': ''.join(sorted(list(shared_mut)))
            'mut': ''.join(sorted([
                i
                for i in mut['mut']
                if i in drm['mut'] or i == mut['ref']
            ]))
        })

    return matched_drm


def get_patient_unique_mutation(table_row):

    isolates = []
    for _, pt_iso_list in group_records_by(table_row, 'PtID').items():
        if len(pt_iso_list) == 1:
            isolates.extend(pt_iso_list)
            continue

        isolate = pt_iso_list[0]
        mutations = [
            bind_mut_str(j)
            for i in pt_iso_list
            for j in parse_mut_str_list(i['CompMutList'])
        ]
        mutations = list(set(mutations))
        isolate['CompMutList'] = ','.join(mutations)
        isolates.append(isolate)

    return isolates

    # for drm_pattern, rec_list in group_records_by(pt_iso_list, ).items():
    #     selected_pattern.append(rec_list[0])


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


def geno_analysis(geno_file, folder, meta=None, rx=None):
    table_raw = load_csv(geno_file)
    main_drm_list = load_main_drm(WS / 'mutations.yml')
    main_drm_list['drms'] = (
        main_drm_list['nonpolymorphic'] + main_drm_list['polymorphic'])

    # table_meta = load_csv(meta)
    # table_rx_history = load_csv(rx)

    print('# Raw isolates: ', len(table_raw))

    isolates = get_patient_unique_mutation(table_raw)

    dump_csv(folder / '1_one_per_pt_isolates.csv', isolates)

    print(f'Isolates after removing duplication: {len(isolates)}')
    print('*' * 20)

    isolates = [
        get_isolate_nonpoly_drm(i, main_drm_list)
        for i in isolates
    ]

    isolates = [
        i
        for i in isolates
        if i['INIMajorDRMs'] or i['INIMinorDRMs']
        # if i['drms']
    ]

    dump_csv(folder / '2_extract_drms.csv', isolates)

    # table_nonpoly = [
    #     i
    #     for i in isolates
    #     if i['nonpoly_drms']
    # ]

    # print(f'Including isolates with nonpoly DRMs: {len(table_nonpoly)}')
    # print('*' * 20)

    return

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


def coevol_analysis(file_path, by_pos=False):
    table = load_csv(file_path)
    drms_seq = [
        i['drms']
        for i in table
        if i['drms']
    ]

    mutations = list(set([
        bind_mut_str({
            'ref': j['ref'],
            'pos': j['pos'],
            'mut': k
        })
        for i in table
        for j in parse_mut_str_list(i['drms'])
        for k in j['mut']
    ]))

    print('#Mutations', len(mutations))
    print(sorted(mutations))

    drms_seq = [
        [i.strip() for i in i.split(',')]
        for i in drms_seq
    ]

    print('# seq', len(drms_seq))
    # main_drm_list = load_main_drm(WS / 'mutations.yml')

    result = []

    for i, j in combinations(mutations, 2):

        posA = parse_mut_str(i)['pos']
        posB = parse_mut_str(j)['pos']
        if posA == posB:
            continue

        if posA > posB:
            i, j = j, i
            posB, posA = posA, posB

        # print(posA, posB)

        # if not ((posA in (263, 155, 118, 148)) or
        #         (posB in (263, 155, 118, 148))):
        #     continue

        listA, listB = get_marked_list(
            drms_seq, parse_mut_str(i), parse_mut_str(j))

        spearman = calc_spearman(listA, listB)
        # jaccard = calc_jaccard(listA, listB)
        contigency = calc_contigency_table(listA, listB)

        resp = {
            'MutA': i,
            'PosA': posA,
            'MutB': j,
            'PosB': posB,
            'Both': contigency['both'],
            'MutA_only': contigency['i_only'],
            'MutB_only': contigency['j_only'],
            'None': contigency['none'],
            'total': contigency['total'],
        }
        resp.update(spearman)
        # resp.update(jaccard)

        result.append(resp)

    # result = [
    #     i
    #     for i in result
    #     if i['Both'] > 1
    # ]

    # assert (len(result) == len(combinations(mutations, 2)))

    dump_csv(file_path.parent / 'mutation_coexist.csv', result)


# TODO, dataclass mutation and mutation mixture
@dataclass
class Mutation:
    ref: str
    pos: int
    mut: list

    @property
    def is_mixture(self):
        return len(self.mut) > 1

    def is_superset(self, mutation):
        if self.ref != mutation.ref:
            print('Ref mutation error.')

        if self.pos != mutation.pos:
            return False

        return set(self.mut).issuperset(set(mutation.mut))

    def contain(self, mutation):
        if self.pos != mutation.pos:
            return False
        if self.ref != mutation.ref:
            print('Ref mutation error.')

        return all([i in self.mut for i in mutation.mut])

    def no_overlap(self, mutation):
        if self.pos != mutation.pos:
            return False
        if self.ref != mutation.ref:
            print('Ref mutation error.')

        return not (set(self.mut) & set(mutation.mut))

    def is_same(self, mutation):
        if self.pos != mutation.pos:
            return False
        if self.ref != mutation.ref:
            print('Ref mutation error.')

        return set(self.mut) == set(mutation.mut)


def get_marked_list(drms_seq, mutA, mutB):
    mutA = Mutation(**mutA)
    mutB = Mutation(**mutB)

    listA = []
    listB = []

    for i in drms_seq:
        check_i = [
            Mutation(**parse_mut_str(j))
            for j in i
        ]

        check_A = [
            j
            for j in check_i
            if j.pos == mutA.pos
        ]

        check_B = [
            j
            for j in check_i
            if j.pos == mutB.pos
        ]

        if (not check_A) and (not check_B):
            listA.append(0)
            listB.append(0)
            continue

        if check_A and (not check_B):
            check_A = check_A[0]
            if check_A.is_same(mutA):
                listA.append(1)
                listB.append(0)
            continue

        if (not check_A) and (check_B):
            check_B = check_B[0]
            if check_B.is_same(mutB):
                listA.append(0)
                listB.append(1)
            # listA.append(0)
            # listB.append(1 if check_B.contain(mutB) else 0)
            continue

        check_A = check_A[0]
        check_B = check_B[0]

        if (
                (check_A.contain(mutA) and not check_A.is_same(mutA))
                or
                (check_B.contain(mutB) and not check_B.is_same(mutB))
        ):
            continue
        if check_A.is_same(mutA) and check_B.is_same(mutB):
            listA.append(1)
            listB.append(1)
        elif check_A.no_overlap(mutA) and check_B.is_same(mutB):
            listA.append(0)
            listB.append(1)
        elif check_A.is_same(mutA) and check_B.no_overlap(mutB):
            listA.append(1)
            listB.append(0)
        elif check_A.no_overlap(mutA) and check_B.no_overlap(mutB):
            listA.append(0)
            listB.append(0)

    return listA, listB


def draw_potential_networks(file_path):

    table = load_csv(file_path)

    # table = [
    #     i
    #     for i in table
    #     if float(i['spearman_rho']) > 0
    # ]

    table = calc_holm_Bonferroni_multiply_way(table, 'spearman_p_value')
    # table = calc_holm_Bonferroni(table, 'jaccard_p_value', 0.05)
    dump_csv(file_path, table)

    draw_by_statistic_value(
        file_path.parent, [
            i
            for i in table
            # if float(i['spearman_p_value']) <= 0.05
            if float(i['HB_spearman_p_value']) <= 0.05
        ],
        'spearman_rho')

    draw_by_statistic_value(
        file_path.parent, [
            i
            for i in table
            # if float(i['jaccard_p_value']) <= 0.05
            if float(i['jaccard_p_value']) <= 0.01
        ],
        'jaccard')

    # draw_by_statistic_p_value(
    #     file_path.parent, table, 'spearman_p_value')

    # draw_by_statistic_p_value(
    #     file_path.parent, table, 'jaccard_p_value')


def draw_by_statistic_value(
        save_folder, table,
        column='spearman_rho', selector=lambda x: float(x) >= 0):

    table = [
        i
        for i in table
        if selector(i[column])
    ]
    print('# edges', len(table))

    save_path = save_folder / f'{column}.svg'

    draw_network(
        save_path, table, column,
        edge_getter=lambda x: int(10 * float(x) + 0.5))


def draw_by_statistic_p_value(
        save_folder, table,
        column='spearman_p_value', selector=lambda x: float(x) < 0.05):

    table = [
        i
        for i in table
        if selector(i[column])
    ]
    print('# edges', len(table))

    save_path = save_folder / f'{column}.svg'

    edge_width_settings = {
        0.05: 1,
        0.005: 2,
        0.0005: 3,
        0.00005: 4,
    }

    draw_network(
        save_path, table, column,
        edge_getter=lambda x: max(
            w
            for c, w in edge_width_settings.items()
            if float(x) <= c
        ))


def draw_network(save_path, table, column, edge_getter):

    edges = get_edge_weight(table, column, edge_getter)

    vertices = sorted(list(set([
        v
        for pair in edges.keys()
        for v in pair
    ])))

    edges = [
        (vertices.index(i), vertices.index(j), w)
        for (i, j), w in edges.items()
    ]

    # graph = ig.Graph()

    # graph.add_vertices(vertices, {'color': 'cyan'})

    # for (vi, vj, w) in edges:
    #     graph.add_edge(vi, vj, width=w, len=14)

    # graph.vs["label"] = graph.vs["name"]

    # G = graph.to_networkx()

    net = Network(
        height='2000px', width='2000px', bgcolor='white', notebook=True)

    # net.from_nx(G)

    for idx, v in enumerate(vertices):
        pos = parse_mut_str(v)['pos']
        if pos in (118, 263, 155, 148):
            net.add_node(idx, label=v, shape='circle', color='cyan')
        else:
            net.add_node(idx, label=v, shape='circle')

    for idx, jdx, w in edges:
        net.add_edge(idx, jdx, width=w*2, color='#9BC0F9')

    print(str(save_path.parent / (save_path.name + '.html')))

    net.show_buttons(filter_=['physics'])

    net.show(str(save_path.parent / (save_path.name + '.html')))

    # graph.write_svg(
    #     str(save_path), width=1000, height=1000,
    #     layout='fr',
    #     font_size=10,
    #     vertex_size=18)


def get_edge_weight(table, column, edge_getter):

    edges = {
        (i['MutA'], i['MutB']): edge_getter(i[column])
        for i in table
    }

    return edges


def work():
    # get_geno_pheno(
    #     DB / 'May 30, 2023', DB / 'May 30, 2023' / 'geno_rx_pheno.csv')
    # geno_analysis(
    #     geno_file=DB / 'Aug 02, 2023' / 'geno-rx.dataset.tsv',
    #     folder=DB / 'Aug 02, 2023',
    #     meta=DB / 'Aug 02, 2023' / 'paper_meta.csv',
    #     rx=DB / 'Aug 02, 2023' / 'tblRxHistory.csv')

    # coevol_analysis(DB / 'Aug 02, 2023' / 'unique_isolates.csv')
    # draw_potential_networks(DB / 'Aug 02, 2023' / 'mutation_coexist.csv')

    geno_analysis(
        geno_file=DB / 'Aug 28, 2025' / 'DTG_InVivo.csv',
        folder=DB / 'Aug 28, 2025',
        meta=None,
        rx=None)

    coevol_analysis(DB / 'Aug 28, 2025' / '2_extract_drms.csv')


if __name__ == '__main__':
    work()
    print('Finish.')
