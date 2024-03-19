from wtforms import (
    Form,
    StringField,
    HiddenField,
    BooleanField,
    FormField,
    FieldList,
    IntegerField,
    SubmitField,
    validators,
)

from submission.utils.custom_fields import TagListField
from submission.utils.custom_widgets import (
    FieldListAddBtn,
)
from submission.edit.forms.biosynthesis_domains import (
    CondensationDomain,
    AdenylationDomain,
    CarrierDomain,
    MonomerForm,
    ModificationDomainForm,
    AcyltransferaseForm,
)


class CalForm(Form):
    _type = HiddenField("cal")
    name = StringField("Name *", validators=[validators.InputRequired()])
    genes = TagListField(
        "Gene(s) *",
        description="Comma separated list of genes in this module",
        validators=[validators.InputRequired()],
    )
    active = BooleanField("Active? *")
    integrated_monomers = FieldList(
        FormField(MonomerForm), widget=FieldListAddBtn(label="Add addional monomer")
    )
    modification_domains = FieldList(
        FormField(ModificationDomainForm),
        min_entries=1,
        max_entries=1,
        widget=FieldListAddBtn(render_kw={"style": "display:none"}),
    )
    comments = StringField("Comments (Optional)")


class NRPS_I_Form(Form):
    # required _type, name, genes, active
    _type = HiddenField("nrps-type1")
    name = StringField("Name *", validators=[validators.InputRequired()])
    genes = TagListField(
        "Gene(s) *",
        description="Comma separated list of genes in this module",
        validators=[validators.InputRequired()],
    )
    active = BooleanField("Active? *")
    c_domain = FieldList(
        FormField(CondensationDomain),
        min_entries=1,
        max_entries=1,
        label="Condensation domain",
        description="If no C domain is present, please skip to 'Adenylation domain'",
        render_kw={"style": "display:none"},
    )
    a_domain = FieldList(
        FormField(AdenylationDomain),
        min_entries=1,
        max_entries=1,
        label="Adenylation domain",
        render_kw={"style": "display:none"},
    )
    carriers = FieldList(
        FormField(CarrierDomain), widget=FieldListAddBtn(label="Add additional carrier")
    )
    integrated_monomers = FieldList(
        FormField(MonomerForm), widget=FieldListAddBtn(label="Add addional monomer")
    )
    modification_domains = FieldList(
        FormField(ModificationDomainForm),
        min_entries=1,
        max_entries=1,
        widget=FieldListAddBtn(render_kw={"style": "display:none"}),
    )
    comments = StringField("Comments (Optional)")


class NRPS_VI_Form(Form):
    _type = HiddenField("nrps-type6")
    name = StringField("Name *", validators=[validators.InputRequired()])
    genes = TagListField(
        "Gene(s) *",
        description="Comma separated list of genes in this module",
        validators=[validators.InputRequired()],
    )
    active = BooleanField("Active? *")
    a_domain = FieldList(
        FormField(AdenylationDomain),
        min_entries=1,
        max_entries=1,
        render_kw={"style": "display:none"},
    )
    carriers = FieldList(
        FormField(CarrierDomain), widget=FieldListAddBtn(label="Add additional carrier")
    )
    integrated_monomers = FieldList(
        FormField(MonomerForm), widget=FieldListAddBtn(label="Add addional monomer")
    )
    modification_domains = FieldList(
        FormField(ModificationDomainForm),
        min_entries=1,
        max_entries=1,
        widget=FieldListAddBtn(render_kw={"style": "display:none"}),
    )
    comments = StringField("Comments (Optional)")


class OtherForm(Form):
    _type = HiddenField("other")
    name = StringField("Name *", validators=[validators.InputRequired()])
    subtype = StringField("Subtype *", validators=[validators.InputRequired()])
    genes = TagListField(
        "Gene(s) *",
        description="Comma separated list of genes in this module",
        validators=[validators.InputRequired()],
    )
    active = BooleanField("Active? *")
    integrated_monomers = FieldList(
        FormField(MonomerForm), widget=FieldListAddBtn(label="Add addional monomer")
    )
    modification_domains = FieldList(
        FormField(ModificationDomainForm),
        min_entries=1,
        max_entries=1,
        widget=FieldListAddBtn(render_kw={"style": "display:none"}),
    )
    comments = StringField("Comments (Optional)")


class PKSIterativeForm(Form):
    _type = HiddenField("pks-iterative")
    name = StringField("Name *", validators=[validators.InputRequired()])
    genes = TagListField(
        "Gene(s) *",
        description="Comma separated list of genes in this module",
        validators=[validators.InputRequired()],
    )
    iterations = IntegerField(
        "Number of iterations *", validators=[validators.InputRequired()]
    )
    active = BooleanField("Active? *")
    ks_domain = None  # TODO: add ketosynthase
    at_domain = FieldList(
        FormField(AcyltransferaseForm),
        min_entries=1,
        max_entries=1,
        render_kw={"style": "display:none"},
    )
    carriers = FieldList(
        FormField(CarrierDomain), widget=FieldListAddBtn(label="Add additional carrier")
    )
    integrated_monomers = FieldList(
        FormField(MonomerForm), widget=FieldListAddBtn(label="Add addional monomer")
    )
    modification_domains = FieldList(
        FormField(ModificationDomainForm),
        min_entries=1,
        max_entries=1,
        widget=FieldListAddBtn(render_kw={"style": "display:none"}),
    )
    comments = StringField("Comments (Optional)")


class PKSModularForm(Form):
    _type = HiddenField("pks-modular")
    name = StringField("Name *", validators=[validators.InputRequired()])
    genes = TagListField(
        "Gene(s) *",
        description="Comma separated list of genes in this module",
        validators=[validators.InputRequired()],
    )
    active = BooleanField("Active? *")
    ks_domain = None  # TODO: add ketosynthase
    at_domain = FieldList(
        FormField(AcyltransferaseForm),
        min_entries=1,
        max_entries=1,
        render_kw={"style": "display:none"},
    )
    carriers = FieldList(
        FormField(CarrierDomain), widget=FieldListAddBtn(label="Add additional carrier")
    )
    integrated_monomers = FieldList(
        FormField(MonomerForm), widget=FieldListAddBtn(label="Add addional monomer")
    )
    modification_domains = FieldList(
        FormField(ModificationDomainForm),
        min_entries=1,
        max_entries=1,
        widget=FieldListAddBtn(render_kw={"style": "display:none"}),
    )
    comments = StringField("Comments (Optional)")


class PKSTransForm(Form):
    _type = HiddenField("pks-trans")
    name = StringField("Name *", validators=[validators.InputRequired()])
    genes = TagListField(
        "Gene(s) *",
        description="Comma separated list of genes in this module",
        validators=[validators.InputRequired()],
    )
    active = BooleanField("Active? *")
    ks_domain = None  # TODO: add ketosynthase
    carriers = FieldList(
        FormField(CarrierDomain), widget=FieldListAddBtn(label="Add additional carrier")
    )
    integrated_monomers = FieldList(
        FormField(MonomerForm), widget=FieldListAddBtn(label="Add addional monomer")
    )
    modification_domains = FieldList(
        FormField(ModificationDomainForm),
        min_entries=1,
        max_entries=1,
        widget=FieldListAddBtn(render_kw={"style": "display:none"}),
    )
    comments = StringField("Comments (Optional)")


class ModulesForm(Form):
    cal = FieldList(
        FormField(CalForm),
        widget=FieldListAddBtn(label="Add additional module"),
        label="Co-enzyme A ligase (CAL)",
    )
    nrps_type1 = FieldList(
        FormField(NRPS_I_Form),
        widget=FieldListAddBtn(label="Add additional module"),
        label="NRPS Type I",
    )
    nrps_type6 = FieldList(
        FormField(NRPS_VI_Form),
        widget=FieldListAddBtn(label="Add additional module"),
        label="NRPS Type VI",
    )
    pks_iterative = FieldList(
        FormField(PKSIterativeForm),
        widget=FieldListAddBtn(label="Add additional module"),
        label="Iterative PKS",
    )
    pks_modular = FieldList(
        FormField(PKSModularForm),
        widget=FieldListAddBtn(label="Add additional module"),
        label="Modular PKS",
    )
    pks_trans_at = FieldList(
        FormField(PKSTransForm),
        widget=FieldListAddBtn(label="Add additional module"),
        label="Trans-AT PKS",
    )
    other = FieldList(
        FormField(OtherForm), widget=FieldListAddBtn(label="Add additional module")
    )
    submit = SubmitField("Submit")
