import json

class Defaults:
    def __init__(self):
        self.defaults = {}
        self.defaults["hg19"] = {}
        self.defaults["hg19"]["pos"] = "chr4:83910641-84311040"
        self.defaults["hg19"]["snp"] = "rs149273678"

        self.defaults["mm10"] = {}
        self.defaults["mm10"]["pos"] = "chr1:134054000-134071000"
        self.defaults["mm10"]["snp"] = "rs27106747"

    def __getitem__(self, item):
        return self.defaults[item]

    def json(self):
        return json.dumps(self.defaults)
