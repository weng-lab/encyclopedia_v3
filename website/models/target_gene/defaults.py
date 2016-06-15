import json

class TargetGeneDefaults:
    def __init__(self):
        self.defaults = {}
        self.defaults["hg19"] = {}
        self.defaults["hg19"]["pos"] = "chr9:37102673-37516177"
        self.defaults["hg19"]["snp"] = "rs149273678"

        self.defaults["mm10"] = {}
        self.defaults["mm10"]["pos"] = "chr1:134054000-134071000"
        self.defaults["mm10"]["snp"] = "rs27106747"

    def __getitem__(self, item):
        return self.defaults[item]

    def json(self):
        return json.dumps(self.defaults)
