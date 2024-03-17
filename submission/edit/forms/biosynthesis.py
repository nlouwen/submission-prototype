import re
from wtforms import (
    Form,
    Field,
    StringField,
    FieldList,
    FormField,
    SelectField,
    HiddenField,
    IntegerField,
    SubmitField,
    validators,
)
from submission.utils.custom_fields import (
    GeneIdField,
    TagListField,
)

from submission.utils.custom_forms import LocationForm
from submission.utils.custom_widgets import (
    FieldListAddBtn,
    TextInputWithSuggestions,
    SelectDefault,
)
from submission.utils.custom_validators import ValidateTagListRegexp


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
        )
        details = StringField("Details (Optional)")
        references = TagListField(
            "Citation(s)",
            widget=TextInputWithSuggestions(post_url="/edit/get_references"),
        )

    class ThioesteraseForm(Form):
        gene = GeneIdField()
        location = FormField(LocationForm)
        subtype = SelectField(
            "Sub-type", choices=["Type I", "Type II"], widget=SelectDefault()
        )

    subclass = SelectField(
        "Sub-class",
        choices=["Type I", "Type II", "Type III", "Type IV", "Type V", "Type VI"],
        widget=SelectDefault(),
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
    submit = SubmitField("Submit")


class PKSForm(Form):
    subclass = SelectField(
        "Sub-class",
        choices=[
            "Type I",
            "Type II aromatic",
            "Type II highly reducing",
            "Type II arylpolyene",
            "Type III",
        ],
        widget=SelectDefault(),
    )
    cyclases = TagListField(
        "Cyclase(s)",
        [ValidateTagListRegexp(r"^[^, ]*$")],
        description="Comma separated list of PKS cyclase gene IDs.",
    )
    starter_unit = None  # TODO: add to schema
    ketide_length = IntegerField("Ketide length", [validators.NumberRange(min=0)])
    iterative = None  # TODO: add to schema
    submit = SubmitField("Submit")


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
                "Type", choices=["ether", "thioether", "other"], widget=SelectDefault()
            )
            details = StringField("Details (Optional)")

        gene = GeneIdField()
        core_sequence = StringField(
            "Core sequence", description="Core sequence of precursor in amino acids."
        )
        leader_cleavage_location = FormField(
            LocationForm, "Leader cleavage location (Optional)"
        )
        follower_cleavage_location = FormField(
            LocationForm, "Follower cleavage location (Optional)"
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
        "Sub-class",
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
    )
    # Only if not unmodified
    details = StringField("Details (Optional)")
    peptidases = TagListField(
        "Peptidase(s) (Optional)",
        validators=[ValidateTagListRegexp(r"^[^, ]*$")],
        description="Comma separated list of peptidase gene IDs",
    )
    precursors = FieldList(
        FormField(PrecursorForm),
        "Precursor(s)",
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional precursor",
        ),
    )
    submit = SubmitField("Submit")


class SaccharideForm(Form):
    class GlycosylTranferaseForm(Form):
        gene = GeneIdField()
        evidence = SelectField(
            "Evidence type",
            choices=[
                "Sequence-based prediction",
                "Structure-based inference",
                "Knock-out construct",
                "Activity assay",
            ],
            widget=SelectDefault(),
        )
        references = TagListField(
            "Citation(s)",
            widget=TextInputWithSuggestions(post_url="/edit/get_references"),
        )
        specificity = StringField(
            "Specificity (SMILES)",
            [
                validators.Optional(),
                validators.Regexp(
                    regex=re.compile(r"^[\[\]a-zA-Z0-9\@()=\/\\#+.%*-]+$")
                ),
            ],
        )

    class SubclusterForm(Form):
        genes = TagListField(
            "Gene(s)",
            [ValidateTagListRegexp(r"^[^, ]*$")],
            description="Comma separated list of subcluster gene IDs",
        )
        specificity = StringField(
            "Specificity (SMILES)",
            [
                validators.Optional(),
                validators.Regexp(
                    regex=re.compile(r"^[\[\]a-zA-Z0-9\@()=\/\\#+.%*-]+$")
                ),
            ],
        )
        references = TagListField(
            "Citation(s)",
            widget=TextInputWithSuggestions(post_url="/edit/get_references"),
        )

    # TODO: sub-class enum
    # subclass = SelectField("Sub-class", choices=[""])
    glycosyltransferases = FieldList(
        FormField(GlycosylTranferaseForm),
        "Glycosyltransferase(s)",
        min_entries=0,
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
    submit = SubmitField("Submit")


class TerpeneForm(Form):
    subclass = SelectField(
        "Sub-class",
        choices=[
            "Diterpene",
            "Hemiterpene",
            "Monoterpene",
            "Sesquiterpene",
            "Triterpene",
        ],
        widget=SelectDefault(),
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
    )
    submit = SubmitField("Submit")


class OtherForm(Form):
    subclass = SelectField(
        "Sub-class",
        choices=["aminocoumarin", "cyclitol", "other"],
        widget=SelectDefault(),
    )
    details = StringField("Details")
    submit = SubmitField("Submit")


class BioClassesCollection:
    NRPS = NRPSForm
    PKS = PKSForm
    Ribosomal = RibosomalForm
    Saccharide = SaccharideForm
    Terpene = TerpeneForm
    Other = OtherForm


class OperonForm(Form):
    genes = TagListField("Gene(s) forming operon")


class OperonMultipleForm(Form):
    operons = FieldList(
        FormField(OperonForm), widget=FieldListAddBtn(label="Add additional operon")
    )
