import re

from wtforms import (
    Form,
    StringField,
    FieldList,
    FormField,
    SelectMultipleField,
    validators,
    BooleanField,
    ValidationError,
    SubmitField,
    FloatField,
    SelectField,
)

from submission.utils.common import TagListField, FieldListAddBtn, StructureInput


def validate_db_cross(form, field):
    pass


class StructureSingle(Form):
    name = StringField("Compound Name", [validators.DataRequired()])
    synonyms = StringField(
        "Synonyms (Optional)",
        [validators.Optional()],
        description="Synonyms for the compound, separated by commas.",
    )
    formula = StringField("Molecular Formula")
    mass = FloatField(
        "Molecular mass",
        [validators.Optional(), validators.NumberRange(min=0)],
        description="Monoisotopic mass (Dalton) of the molecule. Use a dot as a decimal point, not a comma.",
    )
    structure = StringField(  # TODO: crossreference chemical database, only if not in one enter SMILES, mass, formula manually
        "SMILES representation",
        [
            validators.Optional(),
            validators.Regexp(regex=re.compile(r"^[\[\]a-zA-Z0-9\@()=\/\\#+.%*-]+$")),
        ],
        widget=StructureInput(),
        description="Mandatory for all structurally characterized compounds except for large ones such as most RiPPs and polysaccharides. Chemical structure entered as SMILES string, preferentially isomeric. This can be easily acquired with standard software such as ChemDraw, by, e.g., choosing 'Copy as SMILES'.",
    )
    method = SelectField(  # TODO: add method+citation to 'evidence' formfield, allow multiple
        "Method",
        choices=[
            "",
            "NMR",
            "Mass spectrometry",
            "MS/MS",
            "X-ray crystallography",
            "Chemical derivatisation",
            "Total synthesis",
        ],
        description="Technique used to elucidate/verify the structure",
    )
    classes = SelectMultipleField(
        "Compound class(es) (Optional)",
        description="Hold ctrl or cmd key to select multiple compound classes.",
        choices={
            "Alkaloid": [
                "Amination reaction-derived",
                "Anthranilic acid-derived",
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
                "Monoterpenoid",
                "Sesquiterpenoid",
                "Steroid",
            ],
            "Peptide": [
                "Beta-lactam",
                "Depsipeptide",
                "Glycopeptide",
                "Glycopeptidolipid",
                "Linear",
                "Lipopeptide",
                "Macrocyclic",
            ],
        },
    )
    cyclic = BooleanField("Cyclic Compound?")
    moieties = TagListField(
        "Moieties (Optional)",
        description="Characteristic and/or noteworthy chemical moieties found in compound.",
    )
    references = StringField(
        "Citation(s)", description="Comma separated list of references on this compound"
    )  # TODO: standardize
    db_cross = StringField(
        "Database cross-links (Optional)",
        description="Database cross-reference for this compound (pubchem, chebi, chembl, chemspider, npatlas, lotus, gnps, cyanometdb), "
        "e.g. pubchem:3081434 or chebi:29016 or chembl:CHEMBL414130 or chemspider:6082 or npatlas:NPA004746 or lotus:Q27102265 or gnps:MSV000087858 or cyanometdb:CyanoMetDB_0002",
    )  # TODO: validate input


class StructureMultiple(Form):
    structures = FieldList(
        FormField(StructureSingle),
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add additional compound",
        ),
    )
    # add = SubmitField("Add another compound")
    submit = SubmitField("Submit")
