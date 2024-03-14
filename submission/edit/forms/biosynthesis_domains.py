from wtforms import (
    Form,
    StringField,
    HiddenField,
    BooleanField,
    FormField,
    SelectField,
    FieldList,
)
from submission.utils.custom_fields import TagListField, GeneIdField
from submission.utils.custom_forms import LocationForm, SubtrateEvidenceForm


class CondensationDomain(Form):
    # "required": ["type", "gene", "location"]
    _type = HiddenField("condensation")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    subtype = SelectField(
        "Subtype",
        choices=[
            "",
            "Dual",
            "Starter",
            "LCL",
            "DCL",
            "Ester bond-forming",
            "Heterocyclization",
        ],
    )
    references = StringField("Citation(s)")  # TODO: standardize


class AdenylationDomain(Form):
    # "required": ["type", "gene", "location"],
    class SubstateForm(Form):
        # "required": ["name", "proteinogenic", "structure"]
        name = StringField("Name")
        proteinogenic = BooleanField("proteinogenic?")
        structure = StringField("structure (SMILES)")  # TODO: standardize smiles input

    _type = HiddenField("adenylation")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    inactive = BooleanField("Inactive?")
    evidence = FieldList(FormField(SubtrateEvidenceForm))
    precursor_biosynthesis = TagListField("Gene(s) involved in precursor biosynthesis")
    substrates = FieldList(FormField(SubstateForm))


class CarrierDomain(Form):
    # "required": ["type", "gene", "location"]
    _type = HiddenField("carrier")
    subtype = SelectField("Subtype", choices=["", "ACP", "PCP"])
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    inactive = BooleanField("Inactive?")
    beta_branching = BooleanField("Beta-branching?")
    evidence = FieldList(FormField(SubtrateEvidenceForm))


class AminotransferaseDomain(Form):
    _type = HiddenField("aminotransferase")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    inactive = BooleanField("Inactive?")


class CyclaseDomain(Form):
    _type = HiddenField("cyclase")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    references = StringField("Citation(s)")  # TODO: standardize


class DehydrataseDomain(Form):
    _type = HiddenField("dehydratase")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)


class EnoylreductaseDomain(Form):
    _type = HiddenField("enoylreductase")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)


class EpimeraseDomain(Form):
    _type = HiddenField("epimerase")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)


class HydroxylaseDomain(Form):
    _type = HiddenField("hydroxylase")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)


class KetoreductaseDomain(Form):
    _type = HiddenField("ketoreductase")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    inactive = BooleanField("Inactive?")
    stereochemistry = SelectField(
        "Stereochemistry", choices=["", "A1", "A2", "B1", "B2", "C1", "C2"]
    )
    evidence = FieldList(FormField(SubtrateEvidenceForm))


class MethyltransferaseDomain(Form):
    _type = HiddenField("methyltransferase")
    subtype = SelectField("Subtype", choices=["", "C", "N", "O", "other"])
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)
    details = StringField("Details")


class OtherDomain(Form):
    # "required": ["type", "subtype", "gene", "location"]
    _type = HiddenField("other")
    subtype = StringField("Subtype")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)


class OxidaseDomain(Form):
    _type = HiddenField("oxidase")
    gene = GeneIdField("Gene")
    location = FormField(LocationForm)


class MonomerForm(Form):
    # "required": ["evidence", "name", "structure"]
    evidence = FieldList(FormField(SubtrateEvidenceForm))
    name = StringField("Name")
    structure = StringField("Structure (SMILES)")  # TODO: standardize smiles
