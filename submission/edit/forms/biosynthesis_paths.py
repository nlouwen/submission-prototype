from wtforms import Form, StringField, BooleanField, FormField, FieldList, SubmitField
from submission.utils.custom_fields import TagListField, GeneIdField
from submission.utils.custom_widgets import (
    TextInputWithSuggestions,
    FieldListAddBtn,
    StructureInput,
)

from markupsafe import Markup


class ProductForm(Form):
    name = StringField("Name", description="Name of the product produced by this path")
    structure = StringField("Structure (SMILES)", widget=StructureInput())
    comment = StringField("Comment")


class PathForm(Form):
    steps = StringField(
        "Steps",
        description=Markup(
            "Define the steps in the biosynthetic path using their gene names."
            "<ul class='form-text text-muted'>"
            "<li>Use commas to signify an unordered path, e.g. 'a,b,c', </li>"
            "<li>or a greater-than to indicate ordered paths, e.g. 'a>b>c'. </li>"
            "<li>To indicate a module encase in square brackets, e.g. '[d]'. </li>"
            "<li>To indicate cooperation add a plus sign, e.g. 'a+b'. </li>"
            "<li>Trans-AT can be written as '[d]+c'.</li>"
            "</ul>"
        ),
    )
    products = FieldList(
        FormField(ProductForm),
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional product"),
    )
    references = TagListField(
        "Citation(s)", widget=TextInputWithSuggestions(post_url="/edit/get_references")
    )
    isSubcluster = BooleanField("Is this path carried out by a subcluster?")
    producesPrecursor = BooleanField("Does this path produce a precursor?")


class PathMultipleForm(Form):
    paths = FieldList(
        FormField(PathForm),
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional path"),
    )
    submit = SubmitField("Submit")
