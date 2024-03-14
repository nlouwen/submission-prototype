from wtforms import (
    Form,
    StringField,
    IntegerField,
    FieldList,
    SelectField,
    SelectMultipleField,
    FormField,
    validators,
    SubmitField,
)
from submission.utils.custom_fields import TagListField, MultiCheckboxField
from submission.utils.custom_forms import LocationForm, EvidenceForm
from submission.utils.custom_widgets import FieldListAddBtn, TextInputIndicator


class MinEntryForm(Form):

    class TaxonomyForm(Form):
        ncbi_tax_id = IntegerField("NCBI Taxonomy ID")
        genus = StringField("Genus")
        species = StringField("Species")
        strain = StringField("Strain")

    class LocusForm(Form):

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
        evidence = FieldList(
            FormField(EvidenceForm),
            min_entries=1,
            description="Type of evidence that shows this gene cluster is responsible for the biosynthesis of the designated molecule. Papers highlighting multiple methods can be added under each applicable method using the 'Add additional evidence' button.",
            widget=FieldListAddBtn(
                label="Add additional evidence",
            ),
        )

    loci = FieldList(
        FormField(LocusForm),
        min_entries=1,
        description="Locus or loci where the gene cluster is located",
        widget=FieldListAddBtn(
            label="Add additional locus",
        ),
    )
    b_class = SelectMultipleField(
        "Biosynthetic class(es)",
        choices=["NRPS", "PKS", "Ribosomal", "Saccharide", "Terpene", "Other"],
        description="Hold ctrl or cmd key to select multiple classes in the case of a hybrid gene cluster. "
        "Select all categories that apply: e.g. a compound resulting from a BGC with both non-ribosomal "
        "and polyketide synthases should be both 'NRPS' and 'PKS', while a BGC with both a polyketide "
        "synthase and a glycosyltransferase would be 'PKS' and 'Saccharide'.",
    )
    products = TagListField(
        "Product(s)",
        [validators.InputRequired()],
        description='Comma separated list of produced compounds. To enter a compound name containing a comma, encase in double quotes, e.g. "8,9-dihydrolactimidomycin"',
    )

    completeness = SelectField(
        "Completeness",
        choices=["", "Complete", "Incomplete", "Unknown"],
        description="Are all genes needed for production of compounds present in the specified locus/loci?",
    )
    taxonomy = FormField(TaxonomyForm)
    comments = StringField("Additional comments (Optional)")

    submit = SubmitField("Submit")
