import os, sys, json

sys.path.append(os.path.join(os.path.dirname(__file__), '../../'))
from common.trackhub import TrackHub

class TrackHubTargetGene(TrackHub):
    def __init__(self, args, epigenomes, urlStatus, row):
        super(TrackHubTargetGene, self).__init__(args, epigenomes, urlStatus, row)

    def makeTrackDb(self):
        trackhubFnp = os.path.join(os.path.dirname(__file__),
                                   "..", "views", "target_gene", "trackhub.txt")
        with open(trackhubFnp) as f:
            return f.read()
