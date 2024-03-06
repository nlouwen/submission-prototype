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
from forms.common import (
    TagListField,
    StructureEvidenceForm,
    LocationForm,
    TaxonomyForm,
    FieldListAddBtn,
)
from forms.min_entry import MinEntryForm


class AuxEnzymeForm(Form):
    name = StringField("Protein name", description="e.g. 'NisC'.")
    description = StringField(
        "Description", description="Brief description of function."
    )
    databaseIds = TagListField(
        "Database crosslinks", description="Uniprot or GenBank ID of protein"
    )  # TODO: db id regexp


class EnzymeForm(Form):
    name = StringField("Protein name", description="e.g. 'NisB'.")
    description = StringField(
        "Description", description="Brief description of the enzyme function."
    )
    databaseIds = TagListField(
        "Database crosslinks", description="Uniprot or GenBank ID of protein"
    )  # TODO: db id regexp
    auxiliary_enzymes = FieldList(
        FormField(AuxEnzymeForm),
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add auxiliary enzyme",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_aux_enzyme",
                "hx-swap": "beforebegin",
            },
        ),
    )
    references = TagListField(
        "Citation(s)", description="Comma separated references on the protein"
    )  # TODO: standardize citations


class CompoundForm(Form):
    #  "required": ["name", "evidence"],
    name = StringField(
        "Name",
        validators=[validators.Regexp(r"^[a-zA-Zα-ωΑ-Ω0-9\[\]'()/&,. +-]+$")],
        description="Customarily used compound name.",
    )
    synonyms = StringField(
        "Synonyms",
        validators=[validators.Regexp(r"^[a-zA-Zα-ωΑ-Ω0-9\[\]'()/&,. +-]+$")],
        description="Known synonym compound names.",
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
        "SMARTS",
        validators=[validators.Regexp(r"^.+>>.+$")],
        description="The reaction SMARTS or reaction CXSMARTS string",
    )
    isIterative = BooleanField("Iterative? (modifying all possible substructures)")
    hasLinker = BooleanField(
        "Contains a Link group to represent repeating units of variable length? (e.g. a C-chain of variable length)."
    )
    hasPositionVariationBond = BooleanField(
        "Contains a position variation bond which indicates positional variation of a substituent over multiple atoms (e.g. variable chlorination on an aromatic ring)."
    )
    explicitHydrogen = FieldList(
        FormField(HydrogenForm),
        min_entries=0,
        label="Explicit hydrogens",
        widget=FieldListAddBtn(
            label="Add explicit hydrogen",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_hydrogen",
                "hx-swap": "beforebegin",
            },
        ),
    )
    databaseIds = TagListField(
        "Cross-reference",
        description="Cross-reference to other databases. (rhea, MITE, EC), e.g. rhea:16505",
    )
    evidence_sm = FieldList(
        FormField(ReactionSmartsEvidenceForm),
        min_entries=1,
        label="Evidence",
        widget=FieldListAddBtn(
            label="Add additional evidence",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_tail_reaction_smarts_evidence",
                "hx-swap": "beforebegin",
            },
        ),
    )


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
    databaseIds = TagListField(
        "Cross-reference",
        description="Cross-reference to other databases. (rhea, MITE, EC), e.g. rhea:16505",
    )
    evidence_val = FieldList(
        FormField(ReactionSmartsEvidenceForm),
        min_entries=1,
        label="Evidence",
        widget=FieldListAddBtn(
            label="Add additional evidence",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_val_reaction_evidence",
                "hx-swap": "beforebegin",
            },
        ),
    )


class ReactionForm(Form):
    tailoring = FieldList(
        FormField(TailoringFunctionForm),
        label="Ontology-derived tailoring/maturation reaction term.",
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional ontology",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_ontology",
                "hx-swap": "beforebegin",
            },
        ),
    )
    description = StringField(
        "Additional information about tailoring/maturation reaction."
    )
    reaction_smarts = FieldList(
        FormField(ReactionSmartsForm),
        min_entries=1,
        label="Reaction SMARTS",
        widget=FieldListAddBtn(
            label="Add additional reaction SMARTS",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_smarts",
                "hx-swap": "beforebegin",
            },
        ),
    )
    validated_reactions = FieldList(
        FormField(ValidatedReactionForm),
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional validated reaction",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_val_reaction",
                "hx-swap": "beforebegin",
            },
        ),
    )


class TailoringForm(Form):
    enzyme = FormField(EnzymeForm)
    # compounds = FieldList(
    #     FormField(CompoundForm),
    #     label="Mature compound(s) (end products) associated with enzyme (were acted on by enzyme).",
    #     min_entries=1,
    # )
    # genomic_context = FormField(GenomicContextForm)
    reactions = FieldList(
        FormField(ReactionForm),
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add reaction",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_tailoring_reaction",
                "hx-swap": "beforebegin",
            },
        ),
    )
    comment = StringField("Any additional information about this entry")


class TailoringMultipleForm(Form):
    enzymes = FieldList(
        FormField(TailoringForm),
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add enzyme",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_tailoring_enzyme",
                "hx-swap": "beforebegin",
            },
        ),
    )
    submit = SubmitField("Submit")
