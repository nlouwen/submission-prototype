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


class AssayForm(Form):
    # TODO: show hierarchy, default None
    target = SelectField(
        "Properties",
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
    concentration = FloatField(
        "Concentration",
        validators=[validators.Optional()],
        description="Concentration at which the activity is/isn't observed",
    )
    observed = BooleanField(
        "Observed activity",
        description="Leave unticked if the property was not observed at this concentration",
    )


class BioActivityForm(Form):
    compound = StringField("Compound")
    assays = FieldList(FormField(AssayForm), min_entries=1)
    add = SubmitField(
        "", render_kw={"value": "Add additional activity", "formnovalidate": True}
    )


class BioActivityMultiple(Form):
    activities = FieldList(FormField(BioActivityForm))
    submit = SubmitField("Submit")
