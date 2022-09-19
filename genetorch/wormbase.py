from intermine.webservice import Service
import pandas as pd

class Reported_gene:
    def __init__(self,keyword):
        self.keyword = keyword
        self.data = self.query()
        self.gene = []
        for i in self.data:
            self.gene.append(i['symbol'])

    def query(self):
        service = Service("http://im-dev1.wormbase.org/tools/wormmine/service")
        query = service.new_query("Gene")
        query.add_view("primaryIdentifier", "secondaryIdentifier", "symbol")
        query.add_constraint("allele.phenotype.name", "CONTAINS", self.keyword, code="A")
        return query.rows()
