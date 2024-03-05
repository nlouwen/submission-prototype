from wtforms import (
    Form,
    SubmitField,
    StringField,
    SelectField,
    BooleanField,
    SelectMultipleField,
    FieldList,
    FormField,
    IntegerField,
    validators,
)
from forms.common import TagListField, StructureEvidenceForm, LocationForm, TaxonomyForm
from forms.min_entry import MinEntryForm


class AuxEnzymeForm(Form):
    name = StringField("Protein name, e.g. 'NisC'.")
    description = StringField("Brief description of function.")
    databaseIds = TagListField(
        "Database crosslinks, Uniprot or GenBank"
    )  # TODO: db id regexp


class EnzymeForm(Form):
    name = StringField("Protein name, e.g. 'NisB'.")
    description = StringField("Brief description of the enzyme function.")
    databaseIds = TagListField(
        "Database crosslinks, Uniprot or GenBank"
    )  # TODO: db id regexp
    auxiliary_enzymes = FieldList(FormField(AuxEnzymeForm), min_entries=1)
    references = TagListField("Citation(s)")  # TODO: standardize citations


class CompoundForm(Form):
    #  "required": ["name", "evidence"],
    name = StringField(
        "Customarily used compound name.",
        validators=[validators.Regexp(r"^[a-zA-Zα-ωΑ-Ω0-9\[\]'()/&,. +-]+$")],
    )
    synonyms = StringField(
        "Known synonym compound names.",
        validators=[validators.Regexp(r"^[a-zA-Zα-ωΑ-Ω0-9\[\]'()/&,. +-]+$")],
    )
    classes = SelectField("Class", choices=["", "alkaloid", "nucleoside", "peptide"])
    structure = StringField("Structure (SMILES)")  # TODO: standardize smiles
    databaseIds = StringField("Database cross-links")  # TODO: validate input
    evidence = FieldList(FormField(StructureEvidenceForm), min_entries=1)


class GenomicContextForm(Form):
    classes = SelectMultipleField(
        "Biosynthetic class of biosynthetic gene cluster enzyme is associated with, if applicable.",
        choices=["NRPS", "PKS", "Other", "Ribosomal", "Saccharide", "Terpene"],
    )
    accession_genome = StringField(
        "NCBI GenBank genome accession number/ID. RefSeq genomes are prohibited."
    )
    location = FormField(LocationForm)
    taxonomy = FormField(TaxonomyForm)
    # TODO: Auto-fill form based on MIBiG crosslink (if mibig + minimal information is available)
    databaseIds = StringField(
        "MIBiG crosslink", validators=[validators.Regexp(r"^BGC(\d{7,7})$")]
    )
    evidence = FieldList(
        FormField(MinEntryForm.EvidenceForm), min_entries=1
    )  # TODO: move evidenceform to common


class TailoringFunctionForm(Form):
    function = SelectField(
        "Function",
        choices=[
            "",
            "Acetylation",
            "Acylation",
            "Amination",
            "Biaryl bond formation",
            "Carboxylation",
            "Cyclization",
            "Deamination",
            "Decarboxylation",
            "Dehydration",
            "Dehydrogenation",
            "Demethylation",
            "Dioxygenation",
            "Epimerization",
            "FADH2 supply for chlorination",
            "Glycosylation",
            "Halogenation",
            "Heterocyclization",
            "Hydrolysis",
            "Hydroxylation",
            "Macrolactam formation",
            "Methylation",
            "Monooxygenation",
            "Oxidation",
            "Phosphorylation",
            "Prenylation",
            "Reduction",
            "Sulfation",
            "Other",
        ],
    )
    details = StringField("Details")


class HydrogenForm(Form):
    index = IntegerField("The atom index")
    nrHydrogens = IntegerField("The number of explicit hydrogen atoms on this atom")


class ReactionSmartsEvidenceForm(Form):
    evidenceCode = SelectField(
        "Evidence for enzymatic reaction and substrate specificity",
        choices=[
            "",
            "Heterologous expression",
            "In vitro assay",
            "Isothermal titration calorimetry",
            "Knock-out studies",
            "Site-directed mutagenesis",
            "Surface plasmon resonance",
        ],
    )
    references = StringField("Citation(s)")  # TODO: standardize citations


class ReactionSmartsForm(Form):
    reactionSMARTS = StringField(
        "The reaction SMARTS or reaction CXSMARTS string",
        validators=[validators.Regexp(r"^.+>>.+$")],
    )
    isIterative = BooleanField("Is iterative? (modifying all possible substructures).")
    hasLinker = BooleanField(
        "Contains a Link group to represent repeating units of variable length? (e.g. a C-chain of variable length)."
    )
    hasPositionVariationBond = BooleanField(
        "Contains a position variation bond which indicates positional variation of a substituent over multiple atoms (e.g. variable chlorination on an aromatic ring)."
    )
    explicitHydrogen = FieldList(
        FormField(HydrogenForm), min_entries=1
    )  # TODO: add btn
    databaseIds = TagListField("Cross-reference to other databases. (rhea, MITE, EC)")
    evidence = FieldList(
        FormField(ReactionSmartsEvidenceForm), min_entries=1
    )  # TODO: add btn


class ValidatedReactionForm(Form):
    substrate_substructure = StringField(
        "Substrate (sub)structure (SMILES)"
    )  # TODO: standardize smiles
    product_substructure = TagListField(
        "Product (sub)structure(s) (SMILES)"
    )  # TODO: standardize smiles
    isBalanced = BooleanField(
        "Is the validated reaction balanced (i.e. stoichiometric complete)?"
    )
    isAuthentic = BooleanField(
        "Is this substrate-product pair authentic/experimentally verified (i.e. not substructures)"
    )
    isIntermediate = BooleanField(
        "Is this validated reaction an intermediate step (i.e. not the final product)?"
    )
    description = StringField("Additional information about reaction example")
    databaseIds = TagListField("Cross-reference to other databases. (rhea, MITE, EC)")
    evidence = FieldList(
        FormField(ReactionSmartsEvidenceForm), min_entries=1
    )  # TODO: add btn


class ReactionForm(Form):
    tailoring = FieldList(
        FormField(TailoringFunctionForm),
        label="Ontology-derived tailoring/maturation reaction term.",
        min_entries=1,
    )  # TODO: add btn
    description = StringField(
        "Additional information about tailoring/maturation reaction."
    )
    reactionSMARTS = FieldList(
        FormField(ReactionSmartsForm), min_entries=1
    )  # TODO: add btn
    validated_reactions = FieldList(
        FormField(ValidatedReactionForm), min_entries=1
    )  # TODO: add btn


class TailoringForm(Form):
    enzyme = FormField(EnzymeForm)
    # compounds = FieldList(
    #     FormField(CompoundForm),
    #     label="Mature compound(s) (end products) associated with enzyme (were acted on by enzyme).",
    #     min_entries=1,
    # )
    # genomic_context = FormField(GenomicContextForm)
    reactions = FieldList(FormField(ReactionForm), min_entries=1)  # TODO: add btn
    comment = StringField("Any additional information about this entry")
    submit = SubmitField("Submit")
