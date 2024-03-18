from wtforms import (
    Form,
    StringField,
    HiddenField,
    BooleanField,
    FormField,
    SelectField,
    FieldList,
    validators,
)
from submission.utils.custom_fields import TagListField, GeneIdField
from submission.utils.custom_forms import LocationForm, SubtrateEvidenceForm
from submission.utils.custom_widgets import (
    TextInputWithSuggestions,
    SelectDefault,
    FieldListAddBtn,
    StructureInput,
)


class CondensationDomain(Form):
    # "required": ["type", "gene", "location"]
    _type = HiddenField("condensation")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    subtype = SelectField(
        "Subtype",
        choices=[
            "Dual",
            "Starter",
            "LCL",
            "DCL",
            "Ester bond-forming",
            "Heterocyclization",
        ],
        widget=SelectDefault(),
        validate_choice=False,
    )
    references = TagListField(
        "Citation(s)", widget=TextInputWithSuggestions(post_url="/edit/get_references")
    )


class AdenylationDomain(Form):
    # "required": ["type", "gene", "location"],
    class SubstateForm(Form):
        # "required": ["name", "proteinogenic", "structure"]
        name = StringField("Name")
        proteinogenic = BooleanField("proteinogenic?")
        structure = StringField(
            "structure (SMILES)", widget=StructureInput()
        )  # TODO: standardize smiles input

    _type = HiddenField("adenylation")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    inactive = BooleanField("Inactive?")
    evidence = FieldList(
        FormField(SubtrateEvidenceForm),
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional evidence"),
    )
    precursor_biosynthesis = TagListField("Gene(s) involved in precursor biosynthesis")
    substrates = FieldList(
        FormField(SubstateForm),
        widget=FieldListAddBtn(label="Add additional substrate"),
    )


class CarrierDomain(Form):
    # "required": ["type", "gene", "location"]
    _type = HiddenField("carrier")
    subtype = SelectField(
        "Subtype", choices=["ACP", "PCP"], widget=SelectDefault(), validate_choice=False
    )
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    inactive = BooleanField("Inactive?")
    beta_branching = BooleanField("Beta-branching?")
    evidence = FieldList(
        FormField(SubtrateEvidenceForm),
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional evidence"),
    )


class AminotransferaseDomain(Form):
    _type = HiddenField("aminotransferase")
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location
    inactive = BooleanField("Inactive?")


class CyclaseDomain(Form):
    _type = HiddenField("cyclase")
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location
    references = TagListField(
        "Citation(s)", widget=TextInputWithSuggestions(post_url="/edit/get_references")
    )


class DehydrataseDomain(Form):
    _type = HiddenField("dehydratase")
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location


class EnoylreductaseDomain(Form):
    _type = HiddenField("enoylreductase")
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location


class EpimeraseDomain(Form):
    _type = HiddenField("epimerase")
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location


class HydroxylaseDomain(Form):
    _type = HiddenField("hydroxylase")
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location


class KetoreductaseDomain(Form):
    _type = HiddenField("ketoreductase")
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location
    inactive = BooleanField("Inactive?")
    stereochemistry = SelectField(
        "Stereochemistry",
        choices=["A1", "A2", "B1", "B2", "C1", "C2"],
        widget=SelectDefault(),
        validate_choice=False,
    )
    evidence = FieldList(
        FormField(SubtrateEvidenceForm),
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional evidence"),
    )


class MethyltransferaseDomain(Form):
    _type = HiddenField("methyltransferase")
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location
    subtype = SelectField(
        "Subtype",
        choices=["C", "N", "O", "other"],
        widget=SelectDefault(),
        validate_choice=False,
    )
    details = StringField("Details")


class OtherDomain(Form):
    # "required": ["type", "subtype", "gene", "location"]
    _type = HiddenField("other")
    subtype = StringField("Subtype *", validators=[validators.InputRequired()])
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location


class OxidaseDomain(Form):
    _type = HiddenField("oxidase")
    gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
    location = FormField(LocationForm, label="Location *")  # TODO: require location


class ModificationDomainForm(Form):
    aminotransferase = FieldList(
        FormField(AminotransferaseDomain),
        widget=FieldListAddBtn(label="Add additional domain"),
    )
    cyclase = FieldList(
        FormField(CyclaseDomain), widget=FieldListAddBtn(label="Add additional domain")
    )
    dehydratase = FieldList(
        FormField(DehydrataseDomain),
        widget=FieldListAddBtn(label="Add additional domain"),
    )
    enoylreductase = FieldList(
        FormField(EnoylreductaseDomain),
        widget=FieldListAddBtn(label="Add additional domain"),
    )
    epimerase = FieldList(
        FormField(EpimeraseDomain),
        widget=FieldListAddBtn(label="Add additional domain"),
    )
    hydroxylase = FieldList(
        FormField(HydroxylaseDomain),
        widget=FieldListAddBtn(label="Add additional domain"),
    )
    ketoreductase = FieldList(
        FormField(KetoreductaseDomain),
        widget=FieldListAddBtn(label="Add additional domain"),
    )
    methyltransferase = FieldList(
        FormField(MethyltransferaseDomain),
        widget=FieldListAddBtn(label="Add additional domain"),
    )
    oxidase = FieldList(
        FormField(OxidaseDomain), widget=FieldListAddBtn(label="Add additional domain")
    )
    other = FieldList(
        FormField(OtherDomain), widget=FieldListAddBtn(label="Add additional domain")
    )


class MonomerForm(Form):
    # "required": ["evidence", "name", "structure"]
    name = StringField("Name *", validators=[validators.InputRequired()])
    structure = StringField(
        "Structure (SMILES) *",
        widget=StructureInput(),
        validators=[validators.InputRequired()],
    )  # TODO: standardize smiles
    evidence = FieldList(
        FormField(SubtrateEvidenceForm),
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional evidence"),
    )


class AcyltransferaseForm(Form):
    class SubstrateForm(Form):
        name = SelectField(
            "Name",
            choices=["malonyl-CoA", "methylmalonyl-CoA", "ethylmalonyl-CoA", "other"],
            widget=SelectDefault(),
            validate_choice=False,
        )
        structure = StringField(
            "Structure (SMILES)", widget=StructureInput()
        )  # TODO: standardize smiles
        details = StringField("Details (Optional)")

    _type = HiddenField("acyltransferase")
    gene = GeneIdField()
    location = FormField(LocationForm)
    subtype = SelectField(
        choices=["cis-AT", "trans-AT"], widget=SelectDefault(), validate_choice=False
    )
    inactive = BooleanField("Inactive?")
    substrates = FieldList(
        FormField(SubstrateForm),
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional substrate"),
    )
    evidence = FieldList(
        FormField(SubtrateEvidenceForm),
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional evidence"),
    )
