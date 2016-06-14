class Labs:
    labs = {"michael-snyder" : "Snyder",
            'barbara-wold' : "Wold",
            'bing-ren' : "Ren",
            'kevin-white' : "White",
            'bradley-bernstein' : "Bernstein",
            'john-stamatoyannopoulos' : "Stam",
            'kevin-struhl' : "Struhl",
            'peggy-farnham' : "Farnham",
            'richard-myers' : "Myers",
            'ross-hardison' : "Hardison",
            'sherman-weissman' : "Weissman",
            'vishwanath-iyer' : "Iyer",
            "gregory-crawford" : "Crawford"}

    @staticmethod
    def translate(lab):
        if lab in Labs.labs:
            return Labs.labs[lab]
        return lab
