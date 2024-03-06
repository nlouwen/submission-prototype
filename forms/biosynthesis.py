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
from forms.common import (
    GeneIdField,
    LocationForm,
    TagListField,
    ValidatTagListRegexp,
    FieldListAddBtn,
)

# class BiosyntheticClassesForm(Form):
#     b_class = SelectField(
#         "Biosynthetix class",
#         choices=["NRPS", "PKS", "Ribosomal", "Saccharide", "Terpene", "Other"],
#         render_kw={
#             "hx-post": "/get_class",
#             "hx-target": "next div",
#             "hx-trigger": "change",
#             "placeholder": "Select class",
#         },
#     )


class NRPSForm(Form):
    class ReleaseTypeForm(Form):
        name = SelectField(
            "Release type",
            choices=[
                "",
                "Claisen condensation",
                "Hydrolysis",
                "Macrolactamization",
                "Macrolactonization",
                "None",
                "Other",
                "Reductive release",
            ],
        )
        details = StringField("Details")
        references = StringField("Citation")

    class ThioesteraseForm(Form):
        gene = GeneIdField()
        location = FormField(LocationForm)
        subtype = SelectField("Sub-type", choices=["", "Type I", "Type II"])

    subclass = SelectField(
        "Sub-class",
        choices=["", "Type I", "Type II", "Type III", "Type IV", "Type V", "Type VI"],
    )
    release_types = FieldList(
        FormField(ReleaseTypeForm),
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional release type",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_release",
                "hx-swap": "beforebegin",
            },
        ),
    )
    # add_release_type = SubmitField(
    #     "Add release-type", render_kw={"formnovalidate": True}
    # )
    thioesterases = FieldList(
        FormField(ThioesteraseForm),
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional thioesterase",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_thioesterase",
                "hx-swap": "beforebegin",
            },
        ),
    )
    # add_thioesterase = SubmitField(
    #     "Add thioesterase", render_kw={"formnovalidate": True}
    # )
    submit = SubmitField("Submit")


class PKSForm(Form):
    subclass = SelectField(
        "Sub-class",
        choices=[
            "",
            "Type I",
            "Type II aromatic",
            "Type II highly reducing",
            "Type II arylpolyene",
            "Type III",
        ],
    )
    cyclases = TagListField(
        "Cyclase(s)",
        [ValidatTagListRegexp(r"^[^, ]*$")],
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
            link_type = SelectField("Type", choices=["", "ether", "thioether", "other"])
            details = StringField("Details")

        gene = GeneIdField()
        core_sequence = StringField(
            "Core sequence", description="Core sequence of precursor in amino acids."
        )
        leader_cleavage_location = FormField(LocationForm, "Leader cleavage location")
        follower_cleavage_location = FormField(
            LocationForm, "Follower cleavage location"
        )
        crosslinks = FieldList(
            FormField(CrosslinkForm),
            "Crosslinks",
            min_entries=0,
            widget=FieldListAddBtn(
                label="Add additional crosslink",
                render_kw={
                    "formnovalidate": True,
                    "hx-post": "/add_crosslink",
                    "hx-swap": "beforebegin",
                },
            ),
        )
        # add_crosslinks = SubmitField(
        #     "Add crosslink", render_kw={"formnovalidate": True}
        # )
        recognition_motif = StringField("Recognition motif")

    subclass = SelectField(
        "Sub-class",
        choices=[
            "",
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
    )
    # Only if not unmodified
    details = StringField("Details")
    peptidases = TagListField(
        "Peptidase(s)",
        validators=[ValidatTagListRegexp(r"^[^, ]*$")],
        description="Comma separated list of peptidase gene IDs",
    )
    precursors = FieldList(
        FormField(PrecursorForm),
        "Precursor(s)",
        min_entries=1,
        widget=FieldListAddBtn(
            label="Add additional precursor",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_precursor",
                "hx-swap": "beforebegin",
            },
        ),
    )
    # add_precursors = SubmitField("Add precursor", render_kw={"formnovalidate": True})
    submit = SubmitField("Submit")


class SaccharideForm(Form):
    class GlycosylTranferaseForm(Form):
        gene = GeneIdField()
        evidence = SelectField(
            "Evidence type",
            choices=[
                "",
                "Sequence-based prediction",
                "Structure-based inference",
                "Knock-out construct",
                "Activity assay",
            ],
        )
        references = StringField("Citation", [validators.DataRequired()])
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
            [ValidatTagListRegexp(r"^[^, ]*$")],
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
        references = StringField("Citation", [validators.DataRequired()])

    subclass = SelectField("Sub-class", choices=[""])  # TODO: sub-class enum
    glycosyltransferases = FieldList(
        FormField(GlycosylTranferaseForm),
        "Glycosyltransferase(s)",
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add additional glycosyltransferase",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_glycosyltransferase",
                "hx-swap": "beforebegin",
            },
        ),
    )
    # add_glycosyltransferase = SubmitField(
    #     "Add glycosyltransferase", render_kw={"formnovalidate": True}
    # )
    subclusters = FieldList(
        FormField(SubclusterForm),
        "Subcluster(s)",
        min_entries=0,
        widget=FieldListAddBtn(
            label="Add additional subcluster",
            render_kw={
                "formnovalidate": True,
                "hx-post": "/add_subcluster",
                "hx-swap": "beforebegin",
            },
        ),
    )
    # add_subcluster = SubmitField("Add subcluster", render_kw={"formnovalidate": True})
    submit = SubmitField("Submit")


class TerpeneForm(Form):
    subclass = SelectField(
        "Sub-class",
        choices=[
            "",
            "Diterpene",
            "Hemiterpene",
            "Monoterpene",
            "Sesquiterpene",
            "Triterpene",
        ],
    )
    prenyltransferases = TagListField(
        "Prenyltransferase(s)",
        validators=[ValidatTagListRegexp(r"^[^, ]*$")],
        description="Comma separated list of prenyltransferase gene IDs",
    )
    synthases_cyclases = TagListField(
        "Synthase(s)/Cyclase(s)",
        description="Comma separated list of synthase/cyclase gene IDs",
    )
    precursor = SelectField(
        "Precursor", choices=["", "DMAPP", "FPP", "GGPP", "GPP", "IPP"]
    )
    submit = SubmitField("Submit")


class OtherForm(Form):
    subclass = SelectField(
        "Sub-class", choices=["", "aminocoumarin", "cyclitol", "other"]
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


# TODO: separate biosynth components into their own page
class BiosynthesisForm(Form):
    classes = None
    modules = None
    operons = FieldList(FormField(OperonForm))  # add btn
    paths = None
