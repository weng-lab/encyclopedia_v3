from common.enums import AssayType

class EnhancersSiteInfo:
    site = "enhancers"
    assayType = AssayType.Enhancer
    histMark = "H3K27ac"
    name = "enhancer-like"

class PromotersSiteInfo:
    site = "promoters"
    assayType = AssayType.Promoter
    histMark = "H3K4me3"
    name = "promoter-like"

class TargetGeneSiteInfo:
    site = "target_gene"
    assayType = AssayType.TargetGene
    histMark = "TargetGene"

SiteInfos = {}
SiteInfos[EnhancersSiteInfo.assayType] = EnhancersSiteInfo
SiteInfos[PromotersSiteInfo.assayType] = PromotersSiteInfo
SiteInfos[TargetGeneSiteInfo.assayType] = TargetGeneSiteInfo
