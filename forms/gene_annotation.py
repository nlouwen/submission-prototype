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
from forms.common import (
    TagListField,
    FieldListAddBtn,
    LocationForm,
    SelectDefault,
    FunctionEvidenceForm,
    SubtrateEvidenceForm,
)


class AddGeneForm(Form):
    gene_id = StringField("Gene identifier")
    location = FieldList(
        FormField(LocationForm),
        description="Locations of exons",
        widget=FieldListAddBtn(
            label="Add exon location",
        ),
    )
    strand = SelectField("Strand", choices=[1, -1], widget=SelectDefault())
    translation = StringField("Translation")


class DeleteGeneForm(Form):
    gene_id = StringField("Gene identifier")
    reason = StringField(
        "Reason",
        description="Rationale why this gene is not a part of this gene cluster",
    )


class AnnotationForm(Form):
    class FunctionForm(Form):
        class MutationPhenotype(Form):
            phenotype = StringField(
                "Phenotype", description="Observed phenotype upon mutation."
            )
            details = StringField("Details")
            references = TagListField("Citation(s)")

        function = SelectField(
            "Function",
            choices=[
                "Activation / processing",
                "Maturation",
                "Precursor",
                "Precursor biosynthesis",
                "Regulation",
                "Resistance/immunity",
                "Scaffold biosynthesis",
                "Tailoring",
                "Transport",
                "Other",
            ],
            widget=SelectDefault(),
        )
        evidence = FieldList(
            FormField(FunctionEvidenceForm),
            min_entries=1,
            widget=FieldListAddBtn(
                label="Add additional evidence",
            ),
        )
        mutation_phenotype = FormField(MutationPhenotype)

    gene_id = StringField("Gene identifier")
    name = StringField("Gene name")  # TODO: duplicated in schema?
    product = StringField("Gene product name")
    functions = FieldList(
        FormField(FunctionForm),
        widget=FieldListAddBtn(
            label="Add additional function",
        ),
    )
    tailoring = None  # get from tailoring form


class DomainForm(Form):
    class SubtrateForm(Form):
        evidence = FieldList(
            FormField(SubtrateEvidenceForm),
            widget=FieldListAddBtn(
                label="Add additional evidence",
            ),
        )
        structure = StringField("Structure SMILES")  # TODO: standardize smiles

    location = FormField(LocationForm)
    name = StringField()
    substrates = FieldList(
        FormField(SubtrateForm),
        widget=FieldListAddBtn(
            label="Add additional substrate",
        ),
    )


class GeneAnnotationForm(Form):
    add = FieldList(
        FormField(AddGeneForm),
        widget=FieldListAddBtn(
            label="Add additional gene",
        ),
    )
    delete = FieldList(
        FormField(DeleteGeneForm),
        widget=FieldListAddBtn(
            label="Add gene to remove",
        ),
    )
    annotations = FieldList(
        FormField(AnnotationForm),
        widget=FieldListAddBtn(
            label="Add additional annotation",
        ),
    )
    domains = FieldList(
        FormField(DomainForm),
        widget=FieldListAddBtn(
            label="Add additional domain",
        ),
    )
    submit = SubmitField("Submit")
