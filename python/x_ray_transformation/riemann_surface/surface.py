from .metric.base_metric import BaseMetric
from .domain.base_domain import BaseDomain


class Surface:
    def __init__(self, domain: BaseDomain, metric: BaseMetric):
        if not isinstance(domain, BaseDomain):
            raise Exception(f'Domain is not of object type {BaseDomain}')
        if not issubclass(metric.__class__, BaseMetric):
            raise Exception(f'Metric is not of object type {BaseMetric}')

        self.metric = metric
        self.domain = domain
