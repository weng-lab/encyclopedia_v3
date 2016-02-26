#!/usr/bin/env python

import os, json
from orangecontrib.bio import ontology

class OntologyAuto:
    def __init__(self):
        d = os.path.realpath(os.path.join(os.path.dirname(__file__), "ontology"))

        uberonFnp = os.path.join(d, "uberon.ext.2016-02Feb-22.obo")
        self.uberon = ontology.OBOOntology(uberonFnp)

        efoFnp = os.path.join(d, "efo.2016-02Feb-22.obo")
        self.efo = ontology.OBOOntology(efoFnp)

    def getTissue(self, epi):
        if epi.organ_slims:
            return epi.organ_slims[0]
        bti = epi.biosample_term_id
        if bti in self.adhoc:
            return self.adhoc[bti].part_of
        return epi.biosample_type

def main():
    on = OntologyAuto()

    if 0:
        print on.uberon.term("UBERON:0001049")
        rts = on.uberon.related_terms("UBERON:0001049")
        for rt in rts:
            if rt[0] == "is_a" or rt[0] == "part_of":
                print rt
                print on.uberon.term(rt[1])
        return

    if 0:
        rts = on.efo.related_terms("EFO:0005233")
        for rt in rts:
            print rt
            print on.efo.term(rt[1])

        print "**************"
        print on.efo.related_terms("EFO:0005292")
        print on.efo.term("EFO:0000322")

    with open("/home/mjp/Dropbox/missing.txt") as f:
        for line in f:
            toks = line.rstrip().split()
            if toks[1].startswith('EFO:'):
                print "*******************************", toks[0], toks[1]
                rts = on.efo.related_terms(toks[1])
                for rt in rts:
                    print rt
                    print on.efo.term(rt[1])

if __name__ == '__main__':
    main()

