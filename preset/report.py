class SummaryReport:

    def __init__(self):
        self.__dict__['summary'] = {}

    def __setattr__(self, name, value):
        self.__dict__['summary'][name] = value

    def table(self, prefix=''):

        table = []
        for k, v in self.__dict__['summary'].items():
            table.append({
                'item': f"{prefix}{'_' if prefix else ''}{k}",
                'value': v
            })

        return table
