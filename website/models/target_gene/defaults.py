import json

class Defaults:
    def __init__(self):
        self.defaults = {}
        self.defaults["hg19"] = {}
        self.defaults["hg19"]["pos"] = "chr14:35792756-35819812"
        self.defaults["hg19"]["snp"] = "rs11742570"

    def __getitem__(self, item):
        return self.defaults[item]

    def json(self):
        return json.dumps(self.defaults)
