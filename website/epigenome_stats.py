class EpigenomeStats:
    def __init__(self, wepigenomes):
        self.wepigenomes = wepigenomes

        mm9b = self.wepigenomes.GetByAssemblyAndAssays("mm9", "Both")
        mm9d = self.wepigenomes.GetByAssemblyAndAssays("mm9", "DNase")
        mm9h = self.wepigenomes.GetByAssemblyAndAssays("mm9", "H3K27ac")
        mm10b = self.wepigenomes.GetByAssemblyAndAssays("mm10", "Both")
        mm10d = self.wepigenomes.GetByAssemblyAndAssays("mm10", "DNase")
        mm10h = self.wepigenomes.GetByAssemblyAndAssays("mm10", "H3K27ac")
        hg19b = self.wepigenomes.GetByAssemblyAndAssays("hg19", "Both")
        hg19d = self.wepigenomes.GetByAssemblyAndAssays("hg19", "DNase")
        hg19h = self.wepigenomes.GetByAssemblyAndAssays("hg19", "H3K27ac")

        self.mouse_both = len(mm9b) + len(mm10b)
        self.human_both = len(hg19b)

        self.mouse_not_both = len(set([x.web_id() for x in mm9d.epis] +
                                      [x.web_id() for x in mm9h.epis] +
                                      [x.web_id() for x in mm10h.epis] +
                                      [x.web_id() for x in mm10h.epis]) -
                                  set([x.web_id() for x in mm9b.epis] +
                                      [x.web_id() for x in mm10b.epis]))

        self.human_not_both = len(set([x.web_id() for x in hg19d.epis] +
                                      [x.web_id() for x in hg19h.epis]) -
                                  set([x.web_id() for x in hg19b.epis]))

