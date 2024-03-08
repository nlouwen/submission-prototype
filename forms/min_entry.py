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
from forms.common import (
    TagListField,
    LocationForm,
    StringFieldAddBtn,
    FieldListAddBtn,
    MultiTextInput,
    MultiCheckboxField,
    MultiStringField,
    TextInputIndicator,
)


class MinEntryForm(Form):

    class TaxonomyForm(Form):
        ncbi_tax_id = IntegerField("NCBI Taxonomy ID")
        genus = StringField("Genus")
        species = StringField("Species")
        strain = StringField("Strain")

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

    genome = StringField(
        "Genome identifier",
        [validators.InputRequired()],
        widget=TextInputIndicator(),
        description="E.g., AL645882. Only use GenBank accessions, not RefSeq accessions or GI numbers.",
        render_kw={
            "hx-post": "/query_ncbi",
            "hx-trigger": "change",
            "hx-swap": "innerHTML",
            "hx-target": ".subform#taxonomy",
            "hx-indicator": "#spinner",
        },
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
        description="Type of evidence that shows this gene cluster is responsible for the biosynthesis of the designated molecule. Papers highlighting multiple methods can be added under each applicable method using the 'Add additional evidence' button.",
        widget=FieldListAddBtn(
            label="Add additional evidence",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_evidence",
                "hx-swap": "beforebegin",
            },
        ),
    )
    completeness = SelectField(
        "Completeness",
        choices=["", "Complete", "Incomplete", "Unknown"],
        description="Are all genes needed for production of compounds present?",
    )
    taxonomy = FormField(TaxonomyForm)
    comments = StringField("Additional comments (Optional)")

    submit = SubmitField("Submit")
