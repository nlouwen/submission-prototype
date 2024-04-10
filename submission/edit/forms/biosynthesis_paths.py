from wtforms import (
    Form,
    StringField,
    BooleanField,
    FormField,
    FieldList,
    SubmitField,
    validators,
)
from submission.utils.custom_fields import (
    TagListField,
    smiles_field_factory,
)
from submission.utils.custom_widgets import (
    TextInputWithSuggestions,
    FieldListAddBtn,
)

from markupsafe import Markup


class ProductForm(Form):
    name = StringField(
        "Name *",
        description="Name of the product produced by this path",
        validators=[validators.InputRequired()],
    )
    structure = smiles_field_factory(label="Structure (SMILES)")
    comment = StringField("Comment")


class PathForm(Form):
    steps = StringField(
        "Steps *",
        description=Markup(
            "Define the steps in the biosynthetic path using their gene names."
            "<ul class='form-text text-muted'>"
            "<li>Use commas to signify an unordered path, e.g. 'exA,exB,exC', </li>"
            "<li>or a greater-than to indicate ordered paths, e.g. 'exA>exB>exC'. </li>"
            "<li>To indicate a module encase in square brackets, e.g. '[exD]'. </li>"
            "<li>To indicate cooperation add a plus sign, e.g. 'exA+exB'. </li>"
            "<li>To indicate two enzymes that can be interchanged add a pipe symbol, e.g. 'exE|exF'"
            "<li>Trans-AT can be written as '[exD]+exC'.</li>"
            "</ul>"
        ),
        validators=[validators.InputRequired()],
    )
    products = FieldList(
        FormField(ProductForm),
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional product"),
    )
    references = TagListField(
        "Citation(s) *",
        widget=TextInputWithSuggestions(post_url="/edit/get_references"),
        validators=[validators.InputRequired()],
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
