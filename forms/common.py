"""Commonly used forms"""

from typing import Any, Self
from wtforms import (
    Form,
    Field,
    SelectFieldBase,
    StringField,
    IntegerField,
    SelectMultipleField,
    FieldList,
    validators,
    widgets,
    ValidationError,
    SelectField,
)
import csv
import re
from pathlib import Path
from markupsafe import Markup
from flask import url_for


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
        "Start", validators=[validators.Optional(), validators.NumberRange(min=1)]
    )
    end = IntegerField(
        "End", validators=[validators.Optional(), validators.NumberRange(min=2)]
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
    if bgc_id in valid_ids:
        return True
    if Path(f"{bgc_id}_data.json").exists():
        return True
    return False


class EvidenceForm(Form):
    method = SelectField(
        "Method",
        choices=[
            "",
            "Homology-based prediction",
            "Correlation of genomic and metabolomic data",
            "Gene expression correlated with compound production",
            "Knock-out studies",
            "Enzymatic assays",
            "Heterologous expression",
            "In vitro expression",
        ],
        validators=[validators.InputRequired()],
    )
    references = TagListField(
        "Citation(s)",
        [validators.InputRequired()],
        description="Comma separated list of references",
    )


class SubtrateEvidenceForm(Form):
    name = SelectField(
        "Name",
        choices=[
            "",
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
            "",
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


class StringFieldAddBtn(widgets.TextInput):
    def __init__(self, label="add", render_kw={}) -> None:
        super(StringFieldAddBtn, self).__init__()
        self.label = label
        self.render_kw = render_kw

    def __call__(self, field, **kwargs):
        orig = super(StringFieldAddBtn, self).__call__(field, **kwargs)
        add_btn = f"<button class='btn btn-light' {self.html_params(**self.render_kw)}>{self.label}</button>"
        return orig + Markup(add_btn)


class FieldListAddBtn(widgets.SubmitInput):
    def __init__(self, input_type: str | None = None, label="", render_kw={}) -> None:
        super().__init__(input_type)
        self.label = label
        self.render_kw = render_kw

    def __call__(self, field: Field, **kwargs: object) -> Markup:
        kwargs.setdefault("id", field.id)
        kwargs.setdefault("type", self.input_type)
        return Markup(
            f"<button class='btn btn-light' {self.html_params(**kwargs, **self.render_kw)}>{self.label}</button>"
        )


class MultiTextInput(widgets.TextInput):
    def __init__(
        self, input_type: str | None = None, number: int = 1, render_kw: list[dict] = []
    ) -> None:
        super().__init__(input_type)
        self.number = number
        self.render_kw = render_kw

    def __call__(self, field: Field, **kwargs: object) -> Markup:
        kwargs.setdefault("id", field.id)
        kwargs.setdefault("type", self.input_type)

        inputs = Markup("<div>")
        for field_idx in range(self.number):
            render_kw = self.render_kw[field_idx]
            inputs += Markup("<input %s>" % self.html_params(**render_kw, **kwargs))
        inputs += Markup("</div>")
        return inputs


class MultiStringField(StringField):
    def __init__(
        self,
        label=None,
        validators=None,
        description="",
        widget=None,
        render_kw=None,
        field_number: int = 1,
    ) -> None:
        super().__init__(
            label=label,
            validators=validators,
            description=description,
            widget=widget,
            render_kw=render_kw,
        )
        self.number = field_number

    def _value(self):
        if self.data:
            return self.data.split("_")
        else:
            return ""

    def process_formdata(self, valuelist):
        if valuelist:
            self.data = "_".join(valuelist)
        else:
            self.data = []


class MultiCheckboxField(SelectMultipleField):
    """
    A multiple-select, except displays a list of checkboxes.

    Iterating the field will produce subfields, allowing custom rendering of
    the enclosed checkbox fields.
    """

    widget = widgets.ListWidget(prefix_label=False)
    option_widget = widgets.CheckboxInput()


class SelectDefault(widgets.Select):
    def __init__(
        self,
        multiple: bool = False,
        default: tuple[str, str] | None = ("", " --- Select --- "),
        **kwargs,
    ) -> None:
        super().__init__(multiple)
        self.default = default
        self.kwargs = kwargs

    def __call__(self, field: SelectFieldBase, **kwargs: object) -> Markup:
        select = super().__call__(field, **kwargs)
        if self.default:
            val, label = self.default
            insert_idx = select.find(">") + 1
            default_option = Markup(
                f"<option hidden selected disabled {widgets.html_params(value=val, **kwargs, **self.kwargs)}>{label}</option>"
            )
            select = select[:insert_idx] + default_option + select[insert_idx:]
        return select


class TextInputIndicator(widgets.TextInput):
    def __init__(self, input_type: str | None = None) -> None:
        super().__init__(input_type)

    def __call__(self, field: Field, **kwargs: object) -> Markup:
        spinner = Markup(
            f'<img id="spinner" class="htmx-indicator" style="display:inline-block;margin:5px" src="{url_for("static", filename="img/wifi-fade.svg")}" />'
        )
        kwargs.setdefault("style", "width:95%;display:inline-block")
        return super().__call__(field, **kwargs) + spinner
