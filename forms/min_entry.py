from wtforms import (
    Form,
    StringField,
    IntegerRangeField,
    IntegerField,
    FieldList,
    SelectField,
    SelectMultipleField,
    RadioField,
    FormField,
    validators,
    ValidationError,
    SubmitField,
    Field,
    widgets,
)
from forms.common import TagListField, LocationForm


class MinEntryForm(Form):
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
        )
        references = StringField("Citation", [validators.InputRequired()])

    genome = StringField("Genome identifier", [validators.InputRequired()])
    location = FormField(LocationForm)
    b_class = SelectMultipleField(
        "Biosynthetic class(es)",
        choices=["NRPS", "PKS", "Ribosomal", "Saccharide", "Terpene", "Other"],
    )
    products = TagListField("Product(s)", [validators.InputRequired()])
    evidence = FieldList(FormField(EvidenceForm), min_entries=1)
    add_evidence = SubmitField(
        "Add additional evidence", render_kw={"formnovalidate": True}
    )
    taxonomy = StringField("Species name")
    comments = StringField("Additional comments")

    submit = SubmitField("Submit")
