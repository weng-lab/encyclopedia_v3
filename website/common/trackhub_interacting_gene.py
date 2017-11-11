import os
import sys
import json
import StringIO

from urls import BIB5

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from common.trackhub import TrackHub
from common.helpers_trackhub import Track, PredictionTrack, BigGenePredTrack, BigWigTrack, officialVistaTrack, bigWigFilters, TempWrap

from common.colors_trackhub import PredictionTrackhubColors, EncodeTrackhubColors, OtherTrackhubColors, GetTrackColorSignal
from common.labs import Labs

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils/'))
from files_and_paths import Datasets
from metadataws import MetadataWS


class TrackHubInteractingGene(TrackHub):
    def __init__(self, args, epigenomes, urlStatus, row):
        super(TrackHubInteractingGene, self).__init__(args, epigenomes, urlStatus, row)

    def makeTrackDb(self):
        trackhubFnp = os.path.join(os.path.dirname(__file__),
                                   "..", "views", "interacting_gene", "trackhub.txt")
        with open(trackhubFnp) as f:
            fileLines = f.read()

        lines = []

        dataset = Datasets.all_human
        m = MetadataWS(dataset)
        exps = m.biosample_term_name("GM12878")

        for exp in sorted(exps, key=lambda x: (x.assay_term_name, x.tf, x.lab)):
            lines += [self.trackhubExp(exp)]

        f = StringIO.StringIO()
        for line in lines:
            if line:
                f.write(line + "\n")

        return fileLines + "\n" + f.getvalue()

    def trackhubExp(self, exp):
        if not exp:
            return ""

        bigWigs = bigWigFilters(self.assembly, exp.files)

        if not bigWigs:
            return ""
        bigWig = bigWigs[0]

        url = bigWig.url
        if not url:
            return ""

        color = GetTrackColorSignal(exp)
        if color:
            color = color.rgb

        desc = " ".join([exp.encodeID, exp.biosample_term_name, Labs.translate(exp.lab),
                         exp.assay_term_name, exp.tf])

        track = BigWigTrack(desc, self.priority, url, color).track()
        self.priority += 1
        return track
