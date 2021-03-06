#!/usr/bin/env python

from __future__ import print_function

import os
import sys
import json
import psycopg2
import argparse
import fileinput
import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), "../website/common/"))
from web_epigenomes import WebEpigenomesLoader
from site_info import EnhancersSiteInfo, PromotersSiteInfo

sys.path.append(os.path.join(os.path.dirname(__file__), '../../metadata/utils/'))
from utils import Utils, printWroteNumLines
from dbs import DBS
from metadataws import MetadataWS
from files_and_paths import Datasets


def setupDB(cur):
    for assembly in ["hg19", "mm10"]:
        cur.execute("""
DROP TABLE IF EXISTS bedRanges{assembly};
CREATE TABLE bedRanges{assembly}
(id serial PRIMARY KEY,
chrom text,
startend int4range,
expID text
) """.format(assembly=assembly))


def insertFiles(cur, expID, fnp, assembly):
    peaks = [r for r in fileinput.input(fnp, mode="r", openhook=fileinput.hook_compressed)]
    outF = StringIO.StringIO()
    for peak in peaks:
        toks = peak.rstrip().split('\t')
        outF.write("\t".join([toks[0],
                              "[%s, %s)" % (toks[1], toks[2]),
                              expID]) + "\n")
    outF.seek(0)
    cur.copy_from(outF, "bedRanges" + assembly, '\t',
                  columns=("chrom", "startend", "expID"))
    print("\t", fnp, cur.rowcount)


def test(cur):
    # check chr1:134054000-134071000
    cur.execute("""
SELECT DISTINCT expID
FROM bedRangesmm10
WHERE chrom = 'chr1'
AND startend && int4range(134054000, 134071000)
""")
    print(cur.fetchall())


def build(args, conn, cur):
    setupDB(cur)

    bedFnps = set()

    for siteInfo in [EnhancersSiteInfo, PromotersSiteInfo]:
        epigenomes = WebEpigenomesLoader(args, siteInfo)
        for assembly in ["hg19", "mm10"]:
            for assays in [siteInfo.histMark, "DNase"]:
                epis = epigenomes.GetByAssemblyAndAssays(assembly, assays)
                for epi in epis.epis:
                    for exp in epi.exps():
                        try:
                            bedFnp, bedAssembly = exp.getIDRnarrowPeak(args)
                            if not bedFnp:
                                print(exp.getIDRnarrowPeak(args))
                                print("missing", exp)
                            else:
                                bedFnps.add((exp.encodeID, bedFnp, bedAssembly))
                        except Exception, e:
                            print(str(e))
                            print("bad exp:", exp)

    print("found", len(bedFnps))
    for b in sorted(bedFnps):
        insertFiles(cur, b[0], b[1], b[2])
        conn.commit()

    for assembly in ["mm10", "hg19"]:
        print("indexing", assembly, "chrom")
        cur.execute("""
CREATE INDEX chromIdx{assembly} ON bedRanges{assembly}(chrom);
""".format(assembly=assembly))

        print("indexing", assembly, "startend")
        cur.execute("""
CREATE INDEX rangeIdx{assembly} ON bedRanges{assembly} USING gist (startend);
""".format(assembly=assembly))


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('--process', action="store_true", default=True)
    parser.add_argument('--local', action="store_true", default=False)
    parser.add_argument('--rebuild', action="store_true", default=False)
    args = parser.parse_args()
    return args


def main():
    args = parse_args()

    if args.local:
        dbs = DBS.localAnnotations()
    else:
        dbs = DBS.pgdsn("Annotations")
    dbs["application_name"] = os.path.basename(__file__)

    with psycopg2.connect(**dbs) as conn:
        with conn.cursor() as cur:
            if args.rebuild:
                return build(args, conn, cur)
            # test(cur)


if __name__ == '__main__':
    main()
