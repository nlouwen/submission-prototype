from wtforms import (
    Form,
    SubmitField,
    StringField,
    SelectField,
    BooleanField,
    FieldList,
    FormField,
    IntegerField,
    validators,
)
from submission.utils.custom_fields import TagListField
from submission.utils.custom_widgets import (
    FieldListAddBtn,
    StructureInput,
    SelectDefault,
    TextInputWithSuggestions,
)


class AuxEnzymeForm(Form):
    name = StringField("Protein name", description="e.g. 'NisC'.")
    description = StringField(
        "Description", description="Brief description of function."
    )
    databaseIds = TagListField(
        "Database cross-reference", description="Uniprot or GenBank ID of protein"
    )  # TODO: db id regexp


class EnzymeForm(Form):
    name = StringField("Protein name", description="e.g. 'NisB'.")
    description = StringField(
        "Description", description="Brief description of the enzyme function."
    )
    databaseIds = TagListField(
        "Database cross-reference", description="Uniprot or GenBank ID of protein"
    )  # TODO: db id regexp
    auxiliary_enzymes = FieldList(
        FormField(AuxEnzymeForm),
        min_entries=0,
        label="Auxiliary enzymes (Optional)",
        widget=FieldListAddBtn(
            label="Add auxiliary enzyme",
        ),
    )
    references = TagListField(
        "Citation(s)",
        description="Comma separated references on the protein",
        widget=TextInputWithSuggestions(post_url="/edit/get_references"),
    )


class TailoringFunctionForm(Form):
    function = SelectField(
        "Function",
        choices=[
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
        widget=SelectDefault(),
        validate_choice=False,
    )
    details = StringField("Details")


class HydrogenForm(Form):
    index = IntegerField("The atom index")
    nrHydrogens = IntegerField("The number of explicit hydrogen atoms on this atom")


class ReactionSmartsEvidenceForm(Form):
    evidenceCode = SelectField(
        "Evidence for enzymatic reaction and substrate specificity",
        choices=[
            "Heterologous expression",
            "Inference from genomic data and chemical structure",
            "In vitro assay",
            "Isothermal titration calorimetry",
            "Knock-out studies",
            "Site-directed mutagenesis",
            "Surface plasmon resonance",
        ],
        widget=SelectDefault(),
        validate_choice=False,
    )
    references = TagListField(
        "Citation(s)", widget=TextInputWithSuggestions(post_url="/edit/get_references")
    )


class ReactionSmartsForm(Form):
    reactionSMARTS = StringField(
        "SMARTS",
        validators=[validators.Regexp(r"^.+>>.+$")],
        description="The reaction SMARTS or reaction CXSMARTS string",
    )
    isIterative = BooleanField("Iterative? (modifying all possible substructures)")
    hasFrequencyVariation = BooleanField(
        "Contains a frequency variation to represent repeating units of variable length? (e.g. a C-chain of variable length)."
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
        ),
    )
    databaseIds = TagListField(
        "Cross-reference",
        description='Comma separated cross-references to other databases: Rhea (e.g. "rhea:32647"), MITE (e.g. "MITE0000001"), EC number (e.g. "EC 2.1.1.254")',
    )
    evidence_sm = FieldList(
        FormField(ReactionSmartsEvidenceForm),
        min_entries=1,
        label="Evidence",
        widget=FieldListAddBtn(
            label="Add additional evidence",
        ),
    )


class ValidatedReactionForm(Form):
    substrate_substructure = StringField(
        "Substrate (sub)structure (SMILES)", widget=StructureInput()
    )  # TODO: standardize smiles
    product_substructure = FieldList(
        StringField(label="Product (sub)structure SMILES", widget=StructureInput()),
        label="Product (sub)structure(s)",
        min_entries=1,
        widget=FieldListAddBtn(label="Add product (sub)structure"),
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
    description = StringField(
        "Additional information about reaction example (Optional)"
    )
    databaseIds = TagListField(
        "Cross-reference",
        description='Comma separated cross-references to other databases: Rhea (e.g. "rhea:32647"), MITE (e.g. "MITE0000001"), EC number (e.g. "EC 2.1.1.254")',
    )
    evidence_val = FieldList(
        FormField(ReactionSmartsEvidenceForm),
        min_entries=1,
        label="Evidence",
        widget=FieldListAddBtn(
            label="Add additional evidence",
        ),
    )


class ReactionForm(Form):
    tailoring = FieldList(
        FormField(TailoringFunctionForm),
        label="Ontology-derived tailoring/maturation reaction term.",
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional ontology",
        ),
    )
    description = StringField(
        "Additional information about tailoring/maturation reaction. (Optional)"
    )
    reaction_smarts = FieldList(
        FormField(ReactionSmartsForm),
        min_entries=1,
        max_entries=1,
        label="Reaction SMARTS",
        widget=FieldListAddBtn(render_kw={"style": "display:none"}),
    )
    validated_reactions = FieldList(
        FormField(ValidatedReactionForm),
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional validated reaction",
        ),
    )


class TailoringForm(Form):
    enzyme = FieldList(
        FormField(EnzymeForm),
        label="About the enzyme",
        min_entries=1,
        max_entries=1,
        widget=FieldListAddBtn(render_kw={"style": "display:none"}),
    )
    reactions = FieldList(
        FormField(ReactionForm),
        min_entries=0,
        label="Reactions: One or more substrate (sub)structure -> product (sub)structure pairs that result from application of reaction SMARTS.",
        widget=FieldListAddBtn(
            label="Add reaction",
        ),
    )
    comment = StringField("Any additional information about this entry (Optional)")


class TailoringMultipleForm(Form):
    enzymes = FieldList(
        FormField(TailoringForm),
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add enzyme",
        ),
    )
    submit = SubmitField("Submit")
