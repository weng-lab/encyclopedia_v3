class EpigenomeStats:
    def __init__(self, wepigenomes):
        self.wepigenomes = wepigenomes

        hg19 = self.wepigenomes.GetByAssemblyAndAssays("hg19", "TAD")

        self.human = len(hg19)

        self.mouse_not_both = 0

        self.human_not_both = self.human
