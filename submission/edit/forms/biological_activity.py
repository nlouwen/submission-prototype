"""Biological activity for compounds"""

from wtforms import (
    Form,
    StringField,
    BooleanField,
    FloatField,
    SelectField,
    FieldList,
    FormField,
    SubmitField,
    validators,
)
from markupsafe import Markup

from submission.utils.custom_widgets import (
    FieldListAddBtn,
    SelectDefault,
    SubmitIndicator,
    TextInputWithSuggestions,
)
from submission.utils.custom_fields import TagListField
from submission.utils.custom_validators import RequiredIf, ValidateCitations


class AssayForm(Form):
    class ConcentrationForm(Form):
        concentration = FloatField(
            "Concentration",
            description=Markup(
                "Please add the <u>lowest</u> concentration at which the activity was observed. "
                "If it was not observed, add the <u>highest</u> concentration tested."
            ),
            validators=[validators.Optional()],
        )
        # TODO: expand units
        concentration_unit = SelectField(
            "Unit",
            widget=SelectDefault(),
            choices=[
                "mg/mL",
                "µg/mL",
                "ng/mL",
                "pg/mL",
                "milimolar",
                "micromolar",
                "nanomolar",
                "picomolar",
            ],
            validate_choice=False,
            validators=[
                RequiredIf(
                    "concentration",
                    message="This field is required when concentration is filled.",
                )
            ],
        )

    # TODO: show hierarchy, e.g. ionophore has deeper level: chalcophore etc.
    target = SelectField(
        "Property *",
        widget=SelectDefault(),
        choices={
            "chemical properties": [
                "denitrificative",
                "emulsifier",
                "flavor",
                "fluorescent",
                "surfactant",
                "ionophore",
                "chalcophore",
                "lanthanophore",
                "siderophore",
                "zincophore",
                "odorous metabolite",
                "pigment",
                "radical scavenging",
            ],
            "therapeutic properties": [
                "antibacterial",
                "anti-Gram-negative",
                "anti-Gram-positive",
                "anticancer",
                "antineoplastic",
                "antitumor",
                "antifungal",
                "antiinflammatory",
                "antioomycete",
                "antiparasidal",
                "anthelmintic",
                "antiplasmodial",
                "antimalarial",
                "antiprotozoal",
                "anticoccidial",
                "antiproliferative",
                "antitubulin",
                "antiviral",
                "herbicidal",
                "antialgal",
                "immunomodulatory",
                "immunosuppressive",
                "insecticidal",
                "neuroprotective",
                "sodium channel blocking",
                "toxic",
                "cytotoxic",
                "cytostatic",
                "dermatotoxic",
                "DNA-interfering",
                "enterotoxic",
                "hemolytic",
                "hepatotoxic",
                "irritant",
                "neurotoxic",
                "phytotoxic",
                "tumor promoter",
                "vesicant",
            ],
            "cellular processes": [
                "adhesion",
                "biofilm",
                "cell differentiation",
                "cell envelope",
                "cell wall",
                "cell protectant",
                "cyst formation",
                "exopolysaccharide",
                "extracellular capsule",
                "predation",
                "regulatory",
                "inducer",
                "inhibitor",
                "proteasome inhibition",
                "signalling",
                "plant-defense signalling",
                "morphogen",
                "stress response",
                "antioxidant",
                "carbon storage",
                "cold stress",
                "iron reducing",
                "nitrogen reduction",
                "osmolytic",
                "UV protective",
                "swarming motility",
                "virulence factor",
            ],
            "other": ["other"],
        },
        validate_choice=False,
        validators=[validators.InputRequired()],
    )
    references = TagListField(
        "Citation(s) *",
        description="Comma separated list of references highlighted this activity. If references show different concentration, add them separately.",
        widget=TextInputWithSuggestions(post_url="/edit/get_db_references"),
        validators=[validators.InputRequired(), ValidateCitations()],
    )
    obeserved = BooleanField(
        "Observed",
        description=Markup(
            "Untick if this activity was tested for <u>but not observed</u>"
        ),
        render_kw={"checked": True},
    )
    concentration = FormField(ConcentrationForm, label="Concentration (Optional)")


class BioActivityForm(Form):
    compound = StringField("Compound *", validators=[validators.InputRequired()])
    assays = FieldList(
        FormField(AssayForm),
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add assay",
        ),
    )


class BioActivityMultiple(Form):
    activities = FieldList(
        FormField(BioActivityForm),
        widget=FieldListAddBtn(
            label="Add compound",
        ),
    )
    submit = SubmitField("Submit", widget=SubmitIndicator())
