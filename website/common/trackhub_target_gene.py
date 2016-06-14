import os, sys, json
import StringIO

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from common.trackhub import TrackHub
from common.helpers_trackhub import Track, PredictionTrack, BigGenePredTrack, BigWigTrack, officialVistaTrack, bigWigFilters, BIB5, TempWrap

from common.colors_trackhub import PredictionTrackhubColors, EncodeTrackhubColors, OtherTrackhubColors

sys.path.append(os.path.join(os.path.dirname(__file__), '../../../metadata/utils/'))
from files_and_paths import Datasets
from metadataws import MetadataWS

class TrackHubTargetGene(TrackHub):
    def __init__(self, args, epigenomes, urlStatus, row):
        super(TrackHubTargetGene, self).__init__(args, epigenomes, urlStatus, row)

    def makeTrackDb(self):
        trackhubFnp = os.path.join(os.path.dirname(__file__),
                                   "..", "views", "target_gene", "trackhub.txt")
        with open(trackhubFnp) as f:
            fileLines = f.read()

        lines = []

        dataset = Datasets.all_human
        m = MetadataWS(dataset)
        for exp in m.biosample_term_name("GM12878"):
            print(exp)
            lines += [self.trackhubExp(exp)]

        f = StringIO.StringIO()
        map(lambda line: f.write(line + "\n"), lines)

        return fileLines + "\n" + f.getvalue()

    def trackhubExp(self, exp):
        url, name, color = self._getUrln(exp, False)

        if not url:
            return ""
        desc = Track.MakeDesc(name, exp.age, exp.biosample_term_name)
        desc += str(self.priority)

        track = BigWigTrack(desc, self.priority, url, color).track()
        self.priority += 1
        return track

    def _getUrln(self, exp, norm):
        if not exp:
            return None, None, None

        assay = "DNase"
        if exp.isH3K27ac():
            assay = "H3K27ac"
        elif exp.isH3K4me3():
            assay = " H3K4me3"

        bigWigs = bigWigFilters(self.assembly, exp.files)

        if not bigWigs:
            return None, None, None
            raise Exception("missing bigWigs for " + exp.encodeID)
        bigWig = bigWigs[0]

        url = bigWig.url
        if self.urlStatus.find(url) and not self.urlStatus.get(url):
            url = os.path.join(BIB5, "data", bigWig.expID,
                               bigWig.fileID + ".bigWig")

        if norm:
            if "mm10" == self.assembly:
                url = os.path.join(BIB5, "encode_norm", bigWig.expID, bigWig.fileID + ".norm.bigWig")
            else:
                if bigWig.expID.startswith("EN"):
                    url = os.path.join(BIB5, "encode_norm", bigWig.expID, bigWig.fileID + ".norm.bigWig")
                else:
                    url = os.path.join(BIB5, "roadmap_norm/consolidated/",
                                       bigWig.expID,
                                       bigWig.fileID + '-' + assay + ".fc.signal.norm.bigWig")

        if exp.isH3K27ac():
            name = "H3K27ac Signal"
            color = EncodeTrackhubColors.H3K27ac_Signal.rgb
        elif exp.isH3K4me3():
            name = "H3K4me3 Signal"
            color = EncodeTrackhubColors.H3K4me3_Signal.rgb
        elif exp.isDNaseSeq():
            name = "DNase Signal"
            color = EncodeTrackhubColors.DNase_Signal.rgb
        else:
            pass
            #raise Exception("unexpected exp")
        name = "signal"
        color = EncodeTrackhubColors.DNase_Signal.rgb
        return url, name, color
