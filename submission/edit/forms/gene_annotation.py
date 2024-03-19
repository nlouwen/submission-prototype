from wtforms import (
    Form,
    StringField,
    FieldList,
    FormField,
    SubmitField,
    SelectField,
    validators,
)
from markupsafe import Markup

from submission.utils.custom_fields import TagListField, GeneIdField
from submission.utils.custom_validators import RequiredIf
from submission.utils.custom_widgets import (
    FieldListAddBtn,
    SelectDefault,
    StructureInput,
    TextInputWithSuggestions,
)
from submission.utils.custom_forms import (
    SubtrateEvidenceForm,
    FunctionEvidenceForm,
    LocationForm,
)


class AddGeneForm(Form):
    gene_id = StringField(
        "Gene identifier", description="The commonly used gene name (e.g. nisA)"
    )
    exons = FieldList(
        FormField(LocationForm),
        label="Exons *",
        description="Location of coding sequences (CDS). Please also include the stop codon in the coordinates",
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add CDS location",
        ),
    )
    strand = SelectField(
        "Strand *",
        choices=[1, -1],
        widget=SelectDefault(),
        description="The directionality of the CDS. '1' indicates forward directionality, '-1' indicates reverse diretionality",
        validate_choice=False,
        validators=[validators.InputRequired()],
    )
    translation = StringField(
        "Translation",
        description="The encoded amino acid sequence in one-letter code. Please omit the stop codon",
    )


class DeleteGeneForm(Form):
    gene_id = GeneIdField("Gene *", validators=[validators.InputRequired()])
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
            references = TagListField(
                "Citation(s)",
                widget=TextInputWithSuggestions(post_url="/edit/get_references"),
                validators=[
                    RequiredIf(
                        "phenotype",
                        message="This field is required when phenotype is filled.",
                    )
                ],
            )

        function = SelectField(
            "Function *",
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
            validate_choice=False,
            validators=[validators.InputRequired()],
            description="The biological function of the gene product. If based on homology, please exclude it from this submission.",
        )
        details = StringField(
            "Details", description="Any additional information on this gene's function"
        )
        evidence = FieldList(
            FormField(FunctionEvidenceForm),
            min_entries=1,
            widget=FieldListAddBtn(
                label="Add additional evidence",
            ),
            description="Evidence describing this gene product function",
        )
        mutation_phenotype = FormField(
            MutationPhenotype, label="Mutation phenotype (Optional)"
        )

    gene_id = GeneIdField("Gene *", validators=[validators.InputRequired()])
    name = StringField("Gene name", description="Commonly used gene name (e.g. scbA)")
    product = StringField(
        "Gene product name",
        description="e.g. acyl-homoserine-lactone acylase or hypothetical protein",
    )
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
        structure = StringField(
            "Substrate structure SMILES *",
            widget=StructureInput(),
            validators=[validators.InputRequired()],
        )  # TODO: standardize smiles

    name = StringField("Domain name *", validators=[validators.InputRequired()])
    location = FormField(
        LocationForm, label="Domain location *"
    )  # TODO: require location
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
        description=Markup(
            "Add extra genes that are <u>missing</u> from the GenBank entry"
        ),
    )
    delete = FieldList(
        FormField(DeleteGeneForm),
        widget=FieldListAddBtn(
            label="Add gene to remove",
        ),
        description=Markup(
            "Specify genes that fall within the locus but are <u>not</u> part of this BGC"
        ),
    )
    annotations = FieldList(
        FormField(AnnotationForm),
        widget=FieldListAddBtn(
            label="Add additional annotation",
        ),
        description=Markup("Add functional gene annotations"),
    )
    domains = FieldList(
        FormField(DomainForm),
        widget=FieldListAddBtn(
            label="Add additional domain",
        ),
        description=Markup(
            "Supply location and substrate specificity for non-module domains. "
            "For NRPS or PKS substrate specificity, refer to the 'Biosynthesis modules' section."
        ),
    )
    submit = SubmitField("Submit")
