from submission.edit.forms.structure import StructureMultiple
from submission.edit.forms.min_entry import MinEntryForm
from submission.edit.forms.biological_activity import BioActivityMultiple
from submission.edit.forms.biosynthesis import (
    NRPSForm,
    PKSForm,
    RibosomalForm,
    SaccharideForm,
    TerpeneForm,
    OtherForm,
)
from submission.edit.forms.biosynthesis_paths import PathMultipleForm
from submission.edit.forms.biosynthesis_modules import ModulesForm
from submission.edit.forms.tailoring import TailoringMultipleForm
from submission.edit.forms.gene_annotation import GeneAnnotationForm


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

    paths = PathMultipleForm
    modules = ModulesForm

    tailoring = TailoringMultipleForm
    annotation = GeneAnnotationForm
