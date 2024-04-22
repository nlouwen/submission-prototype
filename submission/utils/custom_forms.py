""" Collection of custom form classes used throughout the submission system """

from markupsafe import Markup
from wtforms import Form, IntegerField, SelectField, validators

from .custom_fields import TagListField
from .custom_widgets import StructureInput, TextInputWithSuggestions, SelectDefault
from .custom_validators import ValidateCitations


def location_form_factory(required: bool = False):
    """Create customized location form

    Args:
        required (bool): flag to add the InputRequired validator

    Returns:
        LocationForm: customized location form
    """

    if required:
        valids = [validators.InputRequired()]
    else:
        valids = [validators.Optional()]

    class LocationForm(Form):
        """Subform for location entry, use in combination with FormField"""

        start = IntegerField(
            "Start", validators=valids + [validators.NumberRange(min=1)]
        )
        end = IntegerField("End", validators=valids + [validators.NumberRange(min=2)])

    return LocationForm


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
            "Synthetic-bioinformatic natural product (syn-BNP)",
            "Homology-based prediction",
        ],
        widget=SelectDefault(),
        validators=[validators.InputRequired()],
    )
    references = TagListField(
        "Citation(s) *",
        [validators.InputRequired(), ValidateCitations()],
        description=Markup(
            "Comma separated list of references. Accepted formats are (in order of preference):<br>"
            "'doi:10.1016/j.chembiol.2020.11.009', 'pubmed:33321099', 'patent:US7070980B2', "
            "'url:https://example.com'.<br>If no publication "
            "is available <u>yet</u>, please use 'doi:pending'."
        ),
        widget=TextInputWithSuggestions(post_url="/edit/get_db_references"),
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
        widget=TextInputWithSuggestions(post_url="/edit/get_db_references"),
        validators=[validators.InputRequired(), ValidateCitations()],
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
            "MicroED",
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
        widget=TextInputWithSuggestions(post_url="/edit/get_db_references"),
        validators=[validators.InputRequired(), ValidateCitations()],
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
        widget=TextInputWithSuggestions(post_url="/edit/get_db_references"),
        validators=[validators.InputRequired(), ValidateCitations()],
    )
