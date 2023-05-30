from preset.file_formats import load_tsv
from preset.file_formats import dump_csv
from preset.table import group_records_by


def get_geno_pheno(folder, save_path):
    geno_data = load_tsv(folder / 'geno-rx.dataset.tsv')
    pheno_data = load_tsv(folder / 'geno-pheno.dataset.tsv')

    geno_isos = [
        i['IsolateID']
        for i in geno_data
    ]

    pheno_isos = [
        i['IsolateID']
        for i in pheno_data
    ]

    both_isos = set(geno_isos) & set(pheno_isos)

    geno_data = group_records_by(geno_data, 'IsolateID')
    pheno_data = group_records_by(pheno_data, 'IsolateID')

    report = []
    for i in both_isos:
        g_data = geno_data[i]
        p_data = pheno_data[i]
        assert(len(g_data) == 1 and len(p_data) == 1)
        res = g_data[0]
        res.update(p_data[0])
        report.append(res)

    dump_csv(save_path, report)

