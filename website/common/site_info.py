from common.enums import AssayType

class EnhancersSiteInfo:
    site = "enhancers"
    assayType = AssayType.Enhancer
    histMark = "H3K27ac"

class PromotersSiteInfo:
    site = "promoters"
    assayType = AssayType.Promoter
    histMark = "H3K4me3"
