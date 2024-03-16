from wtforms import (
    Form,
    StringField,
    HiddenField,
    BooleanField,
    FormField,
    FieldList,
    SelectField,
)
from submission.utils.custom_fields import TagListField, GeneIdField
from submission.utils.custom_widgets import TextInputWithSuggestions


class ProductForm(Form):
    name = StringField("Name")
    structure = StringField("Structure (SMILES)")  # TODO: standardize smiles
    comment = StringField("Comment")


class StepForm(Form):
    _type = SelectField("Type", choices=["", "enzyme", "module"])
    name = StringField("Name")  # required if type=module
    gene = GeneIdField("Gene")  # required if type=enzyme
    comment = StringField("Comment")


class PathForm(Form):
    products = FieldList(FormField(ProductForm))  # TODO: add btn
    steps = FieldList(FormField(StepForm))  # TODO: add btn
    references = TagListField(
        "Citation(s)", widget=TextInputWithSuggestions(post_url="/edit/get_references")
    )
    isSubcluster = BooleanField("Subcluster?")
    producesPrecursor = BooleanField("Produces precursor?")
