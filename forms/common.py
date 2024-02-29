"""Commonly used forms"""

from typing import Any, Self
from wtforms import (
    Form,
    Field,
    StringField,
    IntegerField,
    validators,
    widgets,
    ValidationError,
    SelectField,
)
import csv
import re


class GeneIdField(StringField):
    """Reusable Gene ID field with validators"""

    def __init__(
        self,
        label: str | None = None,
        validators: list[Any] | None = [
            validators.Optional(),
            validators.Regexp(r"^[^, ]*$", message="Invalid Gene ID"),
        ],
        **kwargs,
    ):
        super(GeneIdField, self).__init__(label, validators, **kwargs)


class MITEEntryId(StringField):
    def __init__(
        self,
        label: str | None = None,
        validators: list[Any] | None = [
            validators.Optional(),
            validators.Regexp(r"^MITE\d{7}$", message="Invalid Gene ID"),
        ],
        **kwargs,
    ):
        super(MITEEntryId, self).__init__(label, validators, **kwargs)


class TagListField(Field):
    """Custom field for comma separated input, use double quotes to enter names containing commas"""

    widget = widgets.TextInput()

    def _value(self):
        if self.data:
            return '"%s"' % '", "'.join(self.data)
        else:
            return ""

    def process_formdata(self, valuelist):
        if valuelist:
            self.data = next(csv.reader(valuelist, skipinitialspace=True))
        else:
            self.data = []


# TODO: convert to factory to customize used validators
class LocationForm(Form):
    """Subform for location entry, use in combination with FormField"""

    start = IntegerField(
        "Start", validators=[validators.InputRequired(), validators.NumberRange(min=1)]
    )
    end = IntegerField(
        "End", validators=[validators.InputRequired(), validators.NumberRange(min=2)]
    )


def validate_geneids(form, field):
    if field.data:
        print(field.data)
        raise ValidationError(str(field.data))


class ValidatTagListRegexp(object):
    def __init__(self, regex, message="Invalid Gene ID(s)"):
        self.regex = regex
        self.message = message

    def __call__(self, form, field):
        for gene_id in field.data:
            print(gene_id)
            if not re.match(self.regex, gene_id):
                print(gene_id)
                raise ValidationError(self.message)


def is_valid_bgc_id(bgc_id: str):
    valid_ids = [f"BGC{num:0>7}" for num in range(1, 2750)]
    if bgc_id not in valid_ids:
        return False
    return True


class SubtrateEvidenceForm(Form):
    name = SelectField(
        "Name",
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
    )
    references = TagListField("Citation(s)")  # TODO: standardize citations


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
    )
    references = TagListField("Citation(s)")  # TODO: standardize citations


class TaxonomyForm(Form):
    name = StringField("Species name")
    ncbitaxid = IntegerField("NCBI TaxId")
