import re
from wtforms import (
    Form,
    StringField,
    FieldList,
    FormField,
    SelectField,
    IntegerField,
    SubmitField,
    validators,
)
from submission.utils.custom_fields import (
    GeneIdField,
    ReferenceField,
    TagListField,
    smiles_field_factory,
)
from submission.utils.custom_forms import location_form_factory
from submission.utils.custom_widgets import (
    FieldListAddBtn,
    SubmitIndicator,
    TextInputWithSuggestions,
    SelectDefault,
)
from submission.utils.custom_validators import ValidateTagListRegexp, ValidateCitations


class NRPSForm(Form):
    class ReleaseTypeForm(Form):
        name = SelectField(
            "Release type",
            choices=[
                "Claisen condensation",
                "Hydrolysis",
                "Macrolactamization",
                "Macrolactonization",
                "None",
                "Other",
                "Reductive release",
            ],
            widget=SelectDefault(),
            validate_choice=False,
        )
        details = StringField("Details (Optional)")
        references = ReferenceField(
            "Citation(s)",
            widget=TextInputWithSuggestions(post_url="/edit/get_db_references"),
            validators=[ValidateCitations()],
        )

    class ThioesteraseForm(Form):
        gene = GeneIdField()
        location = FormField(location_form_factory())
        subtype = SelectField(
            "Sub-type",
            choices=["Type I", "Type II"],
            widget=SelectDefault(),
            validate_choice=False,
        )

    subclass = SelectField(
        "Sub-class *",
        choices=["Type I", "Type II", "Type III", "Type IV", "Type V", "Type VI"],
        widget=SelectDefault(),
        validate_choice=False,
        validators=[validators.InputRequired()],
    )
    release_types = FieldList(
        FormField(ReleaseTypeForm),
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional release type",
        ),
    )
    thioesterases = FieldList(
        FormField(ThioesteraseForm),
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional thioesterase",
        ),
    )
    submit = SubmitField("Submit", widget=SubmitIndicator())


class PKSForm(Form):
    subclass = SelectField(
        "Sub-class *",
        choices=[
            "Type I",
            "Type II aromatic",
            "Type II highly reducing",
            "Type II arylpolyene",
            "Type III",
        ],
        widget=SelectDefault(),
        validate_choice=False,
        validators=[validators.InputRequired()],
    )
    cyclases = TagListField(
        "Cyclase(s)",
        [validators.Optional(), ValidateTagListRegexp(r"^[^, ]*$")],
        description="Comma separated list of PKS cyclase gene IDs.",
    )
    starter_unit = None  # TODO: add to schema
    ketide_length = IntegerField("Ketide length", [validators.NumberRange(min=0)])
    iterative = None  # TODO: add to schema
    submit = SubmitField("Submit", widget=SubmitIndicator())


class RibosomalForm(Form):
    class PrecursorForm(Form):
        class CrosslinkForm(Form):
            from_loc = IntegerField(
                "From",
                validators=[validators.Optional(), validators.NumberRange(min=1)],
            )
            to_loc = IntegerField(
                "To", validators=[validators.Optional(), validators.NumberRange(min=2)]
            )
            link_type = SelectField(
                "Type",
                choices=["ether", "thioether", "other"],
                widget=SelectDefault(),
                validate_choice=False,
            )
            details = StringField("Details")

        gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
        core_sequence = StringField(
            "Core sequence *",
            description="Core sequence of precursor in amino acids.",
            validators=[validators.InputRequired()],
        )
        leader_cleavage_location = FormField(
            location_form_factory(), "Leader cleavage location (Optional)"
        )
        follower_cleavage_location = FormField(
            location_form_factory(), "Follower cleavage location (Optional)"
        )
        crosslinks = FieldList(
            FormField(CrosslinkForm),
            "Crosslinks (Optional)",
            min_entries=0,
            widget=FieldListAddBtn(
                label="Add additional crosslink",
            ),
        )
        recognition_motif = StringField("Recognition motif (Optional)")

    subclass = SelectField(
        "Sub-class *",
        description="If unmodified, skip the rest of this form",
        choices=[
            "Unmodified",
            "Atropopeptide",
            "Biarylitide",
            "Bottromycin",
            "Borosin",
            "Crocagin",
            "Cyanobactin",
            "Cyptide",
            "Dikaritin",
            "Epipeptide",
            "Glycocin",
            "Graspetide",
            "Guanidinotide",
            "Head-to-tail cyclized peptide",
            "Lanthipeptide",
            "LAP",
            "Lasso peptide",
            "Linaridin",
            "Methanobactin",
            "Microcin",
            "Microviridin",
            "Mycofactocin",
            "Pearlin",
            "Proteusin",
            "Ranthipeptide",
            "Rotapeptide",
            "Ryptide",
            "Sactipeptide",
            "Spliceotide",
            "Streptide",
            "Sulfatyrotide",
            "Thioamidide",
            "Thiopeptide",
            "Other",
        ],
        widget=SelectDefault(),
        validate_choice=False,
        validators=[validators.InputRequired()],
    )
    # Only if not unmodified
    peptidases = TagListField(
        "Peptidase(s)",
        validators=[ValidateTagListRegexp(r"^[^, ]*$")],
        description="Comma separated list of peptidase gene IDs",
    )
    precursors = FieldList(
        FormField(PrecursorForm),
        "Precursor(s)",
        widget=FieldListAddBtn(
            label="Add additional precursor",
        ),
        description="Note: if the precursor gene is not detected in the genbank entry, please remember to add it in the 'gene annotation' section of the submission system.",
    )
    details = StringField("Details")
    submit = SubmitField("Submit", widget=SubmitIndicator())


class SaccharideForm(Form):
    class GlycosylTranferaseForm(Form):
        gene = GeneIdField("Gene *", validators=[validators.InputRequired()])
        evidence = SelectField(
            "Evidence type *",
            choices=[
                "Sequence-based prediction",
                "Structure-based inference",
                "Knock-out construct",
                "Activity assay",
            ],
            widget=SelectDefault(),
            validate_choice=False,
            validators=[validators.InputRequired()],
        )
        references = ReferenceField(
            "Citation(s) *",
            widget=TextInputWithSuggestions(post_url="/edit/get_db_references"),
            validators=[validators.InputRequired(), ValidateCitations()],
        )
        specificity = smiles_field_factory(
            label="Specificity (SMILES)", show_structure=False
        )

    class SubclusterForm(Form):
        genes = TagListField(
            "Gene(s)",
            [ValidateTagListRegexp(r"^[^, ]*$")],
            description="Comma separated list of subcluster gene IDs",
        )
        specificity = smiles_field_factory(
            label="Specificity (SMILES)", show_structure=False
        )
        references = ReferenceField(
            "Citation(s)",
            widget=TextInputWithSuggestions(post_url="/edit/get_db_references"),
            validators=[ValidateCitations()],
        )

    # TODO: sub-class enum
    # subclass = SelectField("Sub-class", choices=[""])
    glycosyltransferases = FieldList(
        FormField(GlycosylTranferaseForm),
        "Glycosyltransferase(s)",
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional glycosyltransferase",
        ),
    )
    subclusters = FieldList(
        FormField(SubclusterForm),
        "Subcluster(s)",
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add additional subcluster",
        ),
    )
    submit = SubmitField("Submit", widget=SubmitIndicator())


class TerpeneForm(Form):
    subclass = SelectField(
        "Sub-class",
        choices=[
            "Diterpene",
            "Hemiterpene",
            "Monoterpene",
            "Sesquiterpene",
            "Sesterterpene",
            "Triterpene",
        ],
        widget=SelectDefault(),
        validate_choice=False,
    )
    prenyltransferases = TagListField(
        "Prenyltransferase(s)",
        validators=[ValidateTagListRegexp(r"^[^, ]*$")],
        description="Comma separated list of prenyltransferase gene IDs",
    )
    synthases_cyclases = TagListField(
        "Synthase(s)/Cyclase(s)",
        description="Comma separated list of synthase/cyclase gene IDs",
    )
    precursor = SelectField(
        "Precursor",
        choices=["DMAPP", "FPP", "GGPP", "GPP", "IPP"],
        widget=SelectDefault(),
        validate_choice=False,
    )
    submit = SubmitField("Submit", widget=SubmitIndicator())


class OtherForm(Form):
    subclass = SelectField(
        "Sub-class *",
        choices=["aminocoumarin", "cyclitol", "other"],
        widget=SelectDefault(),
        validate_choice=False,
        validators=[validators.InputRequired()],
    )
    details = StringField("Details")
    submit = SubmitField("Submit", widget=SubmitIndicator())


class BioClassesCollection:
    NRPS = NRPSForm
    PKS = PKSForm
    Ribosomal = RibosomalForm
    Saccharide = SaccharideForm
    Terpene = TerpeneForm
    Other = OtherForm


class OperonForm(Form):
    genes = TagListField(
        "Gene(s) *",
        description="Comma separated list of gene IDs forming an operon",
        validators=[validators.InputRequired()],
    )


class OperonMultipleForm(Form):
    operons = FieldList(
        FormField(OperonForm), widget=FieldListAddBtn(label="Add additional operon")
    )
    submit = SubmitField("Submit", widget=SubmitIndicator())
