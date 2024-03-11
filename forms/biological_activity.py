"""Biological activity for compounds"""

from wtforms import (
    Form,
    StringField,
    BooleanField,
    FloatField,
    SelectField,
    SelectMultipleField,
    FieldList,
    FormField,
    SubmitField,
    validators,
)
from forms.common import FieldListAddBtn, SelectDefault, TagListField


class AssayForm(Form):
    class ConcentrationForm(Form):
        concentration = FloatField(
            "Concentration",
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
            description="Concentration at which the activity is observed",
        )

    # TODO: show hierarchy, default None
    target = SelectField(
        "Property",
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
                "antineoplastic (duplicate of anticancer)",  # TODO: adress duplicate
                "antitumor (duplicate of anticancer)",  # TODO: adress duplicate
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
                "signalling",
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
    )
    concentration = FormField(ConcentrationForm, label="Concentration (Optional)")
    references = TagListField(
        "Citation(s)",
        description="Comma separated list of references highlighted this activity. If references show different concentration, add them separately.",
    )
    # TODO: remove if no negative data is tracked
    # observed = BooleanField(
    #     "Observed activity",
    #     # description="Leave unticked if the property was not observed at this concentration",
    # )


class BioActivityForm(Form):
    compound = StringField("Compound")
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
    submit = SubmitField("Submit")
