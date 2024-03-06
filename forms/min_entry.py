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
from forms.common import TagListField, LocationForm, StringFieldAddBtn, FieldListAddBtn


class MinEntryForm(Form):
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
        references = StringField("Citation", [validators.InputRequired()])

    genome = StringField(
        "Genome identifier",
        [validators.InputRequired()],
        description="E.g., AL645882. Only use GenBank accessions, not RefSeq accessions or GI numbers.",
    )
    location = FormField(
        LocationForm,
        description="Start and end coordinates, may be left empty if gene cluster spans entire record.",
    )
    b_class = SelectMultipleField(
        "Biosynthetic class(es)",
        choices=["NRPS", "PKS", "Ribosomal", "Saccharide", "Terpene", "Other"],
        description="Hold ctrl or cmd key to select multiple classes in the case of a hybrid gene cluster. Select all categories that apply: e.g., a polyketide with sugar monomers attached should be both 'Polyketide' and 'Saccharide'.",
    )
    products = TagListField(
        "Product(s)",
        [validators.InputRequired()],
        description='Comma separated list of produced compounds. To enter a compound name containing a comma, encase in double quotes, e.g. "8,9-dihydrolactimidomycin"',
    )
    evidence = FieldList(
        FormField(EvidenceForm),
        min_entries=1,
        description="Type of evidence that shows this gene cluster is responsible for the biosynthesis of the designated molecules.",
        widget=FieldListAddBtn(
            label="Add additional evidence",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_evidence",
                "hx-swap": "beforebegin",
            },
        ),
    )
    # add_evidence = SubmitField(
    #     "Add additional evidence",
    #     render_kw={
    #         "formnovalidate": True,
    #         "hx-post": "/add_evidence",
    #         "hx-swap": "beforeend",
    #         "hx-target": "previous .subform",
    #     },
    # )
    # add_evidence = SubmitField(
    #     "Add additional evidence", render_kw={"formnovalidate": True}
    # )
    taxonomy = StringField(
        "Species name",
        validators=[validators.InputRequired()],
        description="Name of organism including strain identifier, e.g. Streptomyces coelicolor A3(2)",
    )
    comments = StringField("Additional comments (Optional)")

    submit = SubmitField("Submit")
