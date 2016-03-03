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
