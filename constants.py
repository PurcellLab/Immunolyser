from typing import List

class Predictor:
    def __init__(self, short_name, full_name):
        self.short_name = short_name
        self.full_name = full_name

    def __str__(self):
        return self.short_name
    
    def to_dict(self):
        return {"short_name": self.short_name, "full_name": self.full_name}

    def __repr__(self):
        return f"Predictor(short_name='{self.short_name}', full_name='{self.full_name}')"

    @classmethod
    def from_dict(cls, data):
        return cls(data['short_name'], data['full_name'])
    
class Class_One_Predictors:
    MixMHCpred = Predictor("MixMHCpred", "mixMHCpred 3.0")
    MHCflurry = Predictor("MHCflurry", "MHCflurry 2.0")
    NetMHCpan = Predictor("NetMHCpan", "netMHCpan 4.2 b")


class Class_Two_Predictors:
    MixMHC2pred = Predictor("MixMHC2pred", "MixMHC2pred-2.0")
    NetMHCpanII = Predictor("NetMHCpanII", "netMHCIIpan 4.3 e")

def get_all_predictors() -> List[Predictor]:
    """Fetch all predictors dynamically from Class_One_Predictors and Class_Two_Predictors."""
    from constants import Class_One_Predictors, Class_Two_Predictors  # Avoid circular imports
    return [
        getattr(Class_One_Predictors, attr)
        for attr in dir(Class_One_Predictors)
        if isinstance(getattr(Class_One_Predictors, attr), Predictor)
    ] + [
        getattr(Class_Two_Predictors, attr)
        for attr in dir(Class_Two_Predictors)
        if isinstance(getattr(Class_Two_Predictors, attr), Predictor)
    ]

def get_short_name(full_name: str) -> str:
    """Find the short name of a predictor based on its full name."""
    for predictor in get_all_predictors():
        if predictor.full_name == full_name:
            return predictor.short_name
    return None  # Return None if no match is found

class MHC_Class:
    One = "I"
    Two = "II"

APP_VERSION = "2.2.0"
APP_VERSION_DATE = "May 2026"

class ClusteringMethod:
    """Motif-clustering step run on the n-mer peptides of a job.

    GibbsCluster is the long-standing default for both classes. MoDec (Gfeller
    Lab, Racle et al. 2019) is offered for MHC class II only, as an alternative
    that replaces GibbsCluster for that job. The value is persisted per job in
    clustering_method.txt next to mhcclass.txt; a job without that file predates
    MoDec and is, by definition, a GibbsCluster job.
    """
    GIBBS = "gibbscluster"
    MODEC = "modec"
    DEFAULT = GIBBS

    DISPLAY_NAMES = {
        GIBBS: "GibbsCluster 2.0",
        MODEC: "MoDec 1.2",
    }

    # Which methods may be chosen for which MHC class. Single source of truth
    # for the submit-form options, the server-side check in /initialiser and the
    # report page, so a new method is added here and nowhere else.
    ALLOWED_BY_CLASS = {
        MHC_Class.One: (GIBBS,),
        MHC_Class.Two: (GIBBS, MODEC),
    }

    @classmethod
    def normalise(cls, value):
        """Map missing/blank input to the default; return anything else unchanged
        (validation is is_allowed's job, so a bad value is rejected, not coerced)."""
        value = (value or "").strip().lower()
        return value or cls.DEFAULT

    @classmethod
    def is_allowed(cls, method, mhc_class):
        return method in cls.ALLOWED_BY_CLASS.get(mhc_class, ())
