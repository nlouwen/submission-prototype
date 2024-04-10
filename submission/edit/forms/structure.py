from wtforms import (
    Form,
    StringField,
    FieldList,
    FormField,
    SelectMultipleField,
    validators,
    BooleanField,
    SubmitField,
    FloatField,
)

from submission.utils.custom_fields import TagListField, smiles_field_factory
from submission.utils.custom_widgets import (
    FieldListAddBtn,
)
from submission.utils.custom_forms import StructureEvidenceForm


class StructureSingle(Form):
    name = StringField("Compound Name *", [validators.InputRequired()])
    synonyms = TagListField(
        "Synonyms",
        [validators.Optional()],
        description="Synonyms for the compound, separated by commas.",
    )
    structure = smiles_field_factory(  # TODO: crossreference chemical database, only if not in one enter SMILES, mass, formula manually
        label="SMILES representation",
        description="Mandatory for all structurally characterized compounds except for large ones such as most RiPPs and polysaccharides. Chemical structure entered as SMILES string, preferentially isomeric. This can be easily acquired with standard software such as ChemDraw, by, e.g., choosing 'Copy as SMILES'.",
    )
    evidence = FieldList(
        FormField(StructureEvidenceForm),
        description="Evidence linked to elucidation of this compound",
        min_entries=1,
        widget=FieldListAddBtn(label="Add additional evidence"),
    )
    formula = StringField("Molecular Formula")
    mass = FloatField(
        "Molecular mass",
        [validators.Optional(), validators.NumberRange(min=0)],
        description="Monoisotopic mass (Dalton) of the molecule.",
    )
    classes = SelectMultipleField(
        "Compound class(es)",
        description="Hold ctrl or cmd key to select multiple compound classes.",
        choices={
            "Alkaloid": [
                "Amination reaction-derived",
                "Anthranilic acid-derived",
                "Arginine-derived",
                "Guanidine-derived",
                "Histidine-derived",
                "Lysine-derived",
                "Nicotinic acid-derived",
                "Ornithine-derived",
                "Peptide alkaloid",
                "Proline-derived",
                "Purine alkaloid",
                "Serine-derived",
                "Steroidal alkaloid",
                "Tetramate alkaloid",
                "Terpenoid-alkaloid",
                "Tryptophan-derived",
                "Tyrosine-derived",
            ],
            "Shikimic acid-derived": [
                "Aromatic amino acid/simple benzoic acid",
                "Aromatic polyketide",
                "Phenylpropanoid",
                "Terpenoid quinone",
            ],
            "Acetate-derived": [
                "Alkylresorcinol/phloroglucinol polyketide",
                "Chromane polyketide",
                "Cyclic polyketide",
                "Fatty acid",
                "Fatty acid derivate",
                "Linear polyketide",
                "Macrocyclic polyketide",
                "Naphthalene polyketide",
                "Polycyclic polyketide",
                "Polyether polyketide",
                "Xanthone polyketide",
            ],
            "Isoprene-derived": [
                "Atypical terpenoid",
                "Diterpenoid",
                "Hemiterpenoid",
                "Higher terpenoid",
                "Iridoid",
                "Meroterpenoid",
                "Monoterpenoid",
                "Sesquiterpenoid",
                "Steroid",
            ],
            "Peptide": [
                "Beta-lactam",
                "Depsipeptide",
                "Diketopiperazine",
                "Glycopeptide",
                "Glycopeptidolipid",
                "Linear",
                "Lipopeptide",
                "Macrocyclic",
            ],
            "Carbohydrates": [
                "Monosaccharide",
                "Oligosaccharide",
                "Polysaccharide",
                "Nucleoside",
                "Aminoglycoside",
                "Liposaccharide",
                "Glucosinolate",
            ],
            "Glycolysis-derived": ["Butenolides", "γ-Butyrolactones", "Tetronic acids"],
            "Other": ["Lactone", "Ectoine", "Furan", "Phosphonate"],
        },
        render_kw={"style": "height:200px"},
    )
    cyclic = BooleanField("Cyclic Compound?")
    moieties = TagListField(
        "Moieties",
        description="Comma separated list of characteristic and/or noteworthy chemical moieties found in compound.",
        validators=[validators.Optional()],
    )
    db_cross = TagListField(
        "Database cross-reference(s)",
        description="Comma separated list of database cross-references for this compound. Accepted formats: "
        "pubchem:3081434, chebi:29016, chembl:CHEMBL414130, chemspider:6082, npatlas:NPA004746, lotus:Q27102265, gnps:MSV000087858 and cyanometdb:CyanoMetDB_0002",
        validators=[validators.Optional()],
    )  # TODO: validate input


class StructureMultiple(Form):
    structures = FieldList(
        FormField(StructureSingle),
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add additional compound",
        ),
    )
    submit = SubmitField("Submit")
