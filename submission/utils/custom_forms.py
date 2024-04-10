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
        "Method *",
        choices=[
            "In vitro expression",
            "Heterologous expression",
            "Enzymatic assays",
            "Knock-out studies",
            "Gene expression correlated with compound production",
            "Correlation of genomic and metabolomic data",
            "Homology-based prediction",
        ],
        widget=SelectDefault(),
        validators=[validators.InputRequired()],
    )
    references = TagListField(
        "Citation(s) *",
        [validators.InputRequired()],
        description="Comma separated list of references. Accepted formats are: 'doi:10.1016/j.chembiol.2020.11.009', 'pubmed:33321099'",
        widget=TextInputWithSuggestions(post_url="/edit/get_references"),
    )


class SubtrateEvidenceForm(Form):
    name = SelectField(
        "Method *",
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
        validate_choice=False,
        validators=[validators.InputRequired()],
    )
    references = TagListField(
        "Citation(s) *",
        widget=TextInputWithSuggestions(post_url="/edit/get_references"),
        validators=[validators.InputRequired()],
    )


class StructureEvidenceForm(Form):
    method = SelectField(
        "Method *",
        choices=[
            "NMR",
            "Mass spectrometry",
            "MS/MS",
            "X-ray crystallography",
            "Chemical derivatisation",
            "Total synthesis",
            "Retention time match with authentic standard",
        ],
        description="Technique used to elucidate/verify the structure",
        widget=SelectDefault(),
        validate_choice=False,
        validators=[validators.InputRequired()],
    )
    references = TagListField(
        "Citation(s) *",
        description="Comma separated list of references on this compound using this method",
        widget=TextInputWithSuggestions(post_url="/edit/get_references"),
        validators=[validators.InputRequired()],
    )


class FunctionEvidenceForm(Form):
    method = SelectField(
        "Method *",
        choices=[
            "Other in vivo study",
            "Heterologous expression",
            "Knock-out",
            "Activity assay",
        ],
        widget=SelectDefault(),
        validate_choice=False,
        validators=[validators.InputRequired()],
    )
    references = TagListField(
        "Citation(s) *",
        widget=TextInputWithSuggestions(post_url="/edit/get_references"),
        validators=[validators.InputRequired()],
    )
