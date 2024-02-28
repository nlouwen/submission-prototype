"""Biological activity for compounds"""

from wtforms import (
    Form,
    StringField,
    BooleanField,
    DecimalField,
    SelectField,
    SelectMultipleField,
    FieldList,
    FormField,
    SubmitField,
)


class AssayForm(Form):
    observed = BooleanField("Observed activity")
    concentration = DecimalField("Concentration")

    # TODO: show hierarchy
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


class BioActivityForm(Form):
    compound = StringField("Compound")  # ideally prefilled
    assay = FieldList(FormField(AssayForm), min_entries=1)
    add = SubmitField(
        "", render_kw={"value": "Add additional activity", "formnovalidate": True}
    )


class BioActivityMultiple(Form):
    activities = FieldList(FormField(BioActivityForm))
    submit = SubmitField("Submit")
