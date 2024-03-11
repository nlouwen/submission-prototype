from forms.structure import StructureMultiple
from forms.min_entry import MinEntryForm
from forms.biological_activity import BioActivityMultiple
from forms.biosynthesis import (
    NRPSForm,
    PKSForm,
    RibosomalForm,
    SaccharideForm,
    TerpeneForm,
    OtherForm,
)
from forms.tailoring import TailoringMultipleForm
from forms.gene_annotation import GeneAnnotationForm


class FormCollection:
    minimal = MinEntryForm
    structure = StructureMultiple
    bioact = BioActivityMultiple

    # Biosynthesis classes
    NRPS = NRPSForm
    PKS = PKSForm
    Ribosomal = RibosomalForm
    Saccharide = SaccharideForm
    Terpene = TerpeneForm
    Other = OtherForm

    tailoring = TailoringMultipleForm
    annotation = GeneAnnotationForm
