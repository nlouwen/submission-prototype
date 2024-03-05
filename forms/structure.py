from wtforms import (
    Form,
    StringField,
    FieldList,
    FormField,
    SelectMultipleField,
    validators,
    BooleanField,
    ValidationError,
    SubmitField,
    FloatField,
    SelectField,
)
from forms.common import TagListField
import re


def validate_db_cross(form, field):
    pass


class StructureSingle(Form):
    name = StringField("Compound Name", [validators.DataRequired()])
    synonyms = StringField("Synonyms", [validators.Optional()])
    formula = StringField("Molecular Formula")
    mass = FloatField(
        "Molecular mass", [validators.Optional(), validators.NumberRange(min=0)]
    )
    structure = StringField(
        "SMILES representation",
        [
            validators.Optional(),
            validators.Regexp(regex=re.compile(r"^[\[\]a-zA-Z0-9\@()=\/\\#+.%*-]+$")),
        ],
        render_kw={
            "hx-post": "/render-smiles",
            "hx-trigger": "change, load",
            "hx-swap": "innerHTML",
            "hx-target": "next .struct",
        },
    )
    method = SelectField(  # TODO: add method+citation to 'evidence' formfield, allow multiple
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
    ionType = SelectField(
        "Ion type",
        choices=[
            "",
            "[M+H]+",
            "[M-H]-",
            "[M+Na]+",
            "[2M+Na]+",
            "[2M+H]+",
            "[M+2H]2+",
            "[M+3H]3+",
            "[M+56Fe-2H]+",
            "[M+NH4]+",
            "[M+K]+",
            "[M+H2O+H]+",
            "other",
        ],
        validators=[validators.Optional()],
    )
    classes = SelectMultipleField(
        "Compound class(es)", choices=["alkaloid", "nucleoside", "peptide"]
    )
    cyclic = BooleanField("Cyclic")
    moieties = TagListField("Moieties (csv)")
    references = StringField("Citation(s)")  # TODO: standardize
    db_cross = StringField("Database cross-links")  # TODO: validate input


class StructureMultiple(Form):
    structures = FieldList(FormField(StructureSingle), min_entries=0)
    add = SubmitField("Add another compound")
    submit = SubmitField("Submit")
