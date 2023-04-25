import csv
from pathlib import Path
from ruamel.yaml import YAML


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


yaml = YAML()
yaml.preserve_quotes = True


def load_yaml(file_path):
    return yaml.load(open(file_path))
