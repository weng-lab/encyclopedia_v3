#!/usr/bin/env python

import os, json

class AdHoc:
    def __init__(self, toks):
        self.count = toks[0]
	self.term_name = toks[1]
	self.term_id = toks[2]
	self.organism = toks[3]

        self.is_a = ""
        self.part_of = ""
        self.transformed_from = ""
        self.comments = ""

        if len(toks) > 4:
            self.is_a = toks[4]
        if len(toks) > 5:
            self.part_of = toks[5]
        if len(toks) > 6:
            self.transformed_from = toks[6]
        if len(toks) > 7:
            self.comments = toks[7]

    def __repr__(self):
        return "\t".join([self.term_id, self.term_name, self.part_of])

class Ontology:
    def __init__(self):
        self._load_adhoc()
        self._load_ucsc()
        self._load_manual()

    def _load_adhoc(self):
        fnp = os.path.join(os.path.dirname(__file__),
                           "adhoc_slim.tsv.txt")
        with open(fnp) as f:
            lines = [line.rstrip().split('\t') for line in f]
        self.adhoc_lines = [AdHoc(line) for line in lines[1:]]
        self.adhoc = dict((a.term_id, a) for a in self.adhoc_lines if a.term_id)

    def _load_ucsc(self):
        fnp = os.path.join(os.path.dirname(__file__),
                           "ucsc_ontology.txt")
        with open(fnp) as f:
            lines = [line.rstrip().split('\t') for line in f]
        self.ucsc = dict((a[0], a[1]) for a in lines)

    def _load_manual(self):
        fnp = os.path.join(os.path.dirname(__file__),
                           "manual.ontology.txt")
        with open(fnp) as f:
            lines = [line.rstrip().split('\t') for line in f if line]
        self.manual = dict((a[0], a[2]) for a in lines)

    def getTissue(self, epi):
        t = self.getTissueActual(epi).strip()
        if "esc" == t:
            return "ESC"
        if "ipsc" == t:
            return "iPSC"
        if "esc_derived" == t:
            return "ESC derived"
        if "bone element" == t:
            return "bone"
        t = t.replace("gi_", "")
        t = t.replace(" gland", "")
        t = t.replace("large intestine", "intestine")
        t = t.replace(" of body", "")
        t = t.replace("skeletal muscle", "muscle")
        return t

    def getTissueActual(self, epi):
        btn = epi.biosample_term_name
        if "induced pluripotent stem cell" == btn:
            return "ipsc"
        if btn in self.ucsc:
            return self.ucsc[btn]
        if btn in self.manual:
            return self.manual[btn]
        if epi.organ_slims:
            return epi.organ_slims[0]
        bti = epi.biosample_term_id
        if bti in self.adhoc:
            return self.adhoc[bti].part_of
        if "limb" == epi.biosample_term_name:
            return "limb"
        if "neural tube" == epi.biosample_term_name:
            return "brain"
        #print btn, "\t", bti
        return ""

    def getCellType(self, epi):
        bt = epi.biosample_type
        if bt == "PrimaryCulture":
            return "primary cell"
        if bt == "PrimaryCell":
            return "primary cell"
        if bt == "ESCDerived":
            return "ESC derived"
        if bt == "PrimaryTissue":
            return "tissue"
        if bt == "CellLine":
            return "immortalized cell line"
        return bt

def main():
    on = Ontology()
    print on.adhoc["EFO:0002784"]

if __name__ == '__main__':
    main()

