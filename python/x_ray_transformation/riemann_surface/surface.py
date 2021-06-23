from .metric import Metric
from .domain import Domain


class Surface:
    def __init__(self, domain, metric):
        if not isinstance(domain, Domain):
            raise Exception(f'Domain is not of object type {Domain}')
        if not issubclass(metric.__class__, Metric):
            raise Exception(f'Metric is not of object type {Metric}')

        self.metric = metric
        self.domain = domain
