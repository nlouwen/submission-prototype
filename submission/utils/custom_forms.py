""" Collection of custom form classes used throughout the submission system """

from flask import url_for
from wtforms import Form, IntegerField, SelectField, validators

from .custom_fields import TagListField
from .custom_widgets import StructureInput, TextInputWithSuggestions, SelectDefault


class LocationForm(Form):
    """Subform for location entry, use in combination with FormField"""

    start = IntegerField(
        "Start", validators=[validators.Optional(), validators.NumberRange(min=1)]
    )
    end = IntegerField(
        "End", validators=[validators.Optional(), validators.NumberRange(min=2)]
    )


class EvidenceForm(Form):
    method = SelectField(
        "Method",
        choices=[
            "Homology-based prediction",
            "Correlation of genomic and metabolomic data",
            "Gene expression correlated with compound production",
            "Knock-out studies",
            "Enzymatic assays",
            "Heterologous expression",
            "In vitro expression",
        ],
        widget=SelectDefault(),
        validators=[validators.InputRequired()],
    )
    references = TagListField(
        "Citation(s)",
        [validators.InputRequired()],
        description="Comma separated list of references. Accepted formats are: 'doi:10.1016/j.chembiol.2020.11.009', 'PMID:33321099'",
        widget=TextInputWithSuggestions(post_url="/edit/get_references"),
    )


class SubtrateEvidenceForm(Form):
    name = SelectField(
        "Method",
        choices=[
            "Activity assay",
            "ACVS assay",
            "ATP-PPi exchange assay",
            "Enzyme-coupled assay",
            "Feeding study",
            "Heterologous expression",
            "Homology",
            "HPLC",
            "In-vitro experiments",
            "Knock-out studies",
            "Mass spectrometry",
            "NMR",
            "Radio labelling",
            "Sequence-based prediction",
            "Steady-state kinetics",
            "Structure-based inference",
            "X-ray crystallography",
        ],
        widget=SelectDefault(),
    )
    references = TagListField(
        "Citation(s)", widget=TextInputWithSuggestions(post_url="/edit/get_references")
    )


class StructureEvidenceForm(Form):
    method = SelectField(
        "Method",
        choices=[
            "NMR",
            "Mass spectrometry",
            "MS/MS",
            "X-ray crystallography",
            "Chemical derivatisation",
            "Total synthesis",
        ],
        widget=SelectDefault(),
    )
    references = TagListField(
        "Citation(s)", widget=TextInputWithSuggestions(post_url="/edit/get_references")
    )


class FunctionEvidenceForm(Form):
    method = SelectField(
        "Method",
        choices=[
            "Other in vivo study",
            "Heterologous expression",
            "Knock-out",
            "Activity assay",
        ],
        widget=SelectDefault(),
    )
    references = TagListField(
        "Citation(s)", widget=TextInputWithSuggestions(post_url="/edit/get_references")
    )
