from wtforms import (
    Form,
    StringField,
    HiddenField,
    BooleanField,
    FormField,
    FieldList,
    IntegerField,
    SubmitField,
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
    name = StringField("Name")
    genes = TagListField("Gene(s)")
    active = BooleanField("Active?")
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
    name = StringField("Name")
    genes = TagListField("Gene(s)")
    active = BooleanField("Active?")
    c_domain = FieldList(
        FormField(CondensationDomain),
        min_entries=1,
        max_entries=1,
        render_kw={"style": "display:none"},
    )
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


class NRPS_VI_Form(Form):
    _type = HiddenField("nrps-type6")
    name = StringField("Name")
    genes = TagListField("Gene(s)")
    active = BooleanField("Active?")
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
    name = StringField("Name")
    subtype = StringField("Subtype")
    genes = TagListField("Gene(s)")
    active = BooleanField("Active?")
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
    name = StringField("Name")
    genes = TagListField("Gene(s)")
    iterations = IntegerField("Number of iterations")
    active = BooleanField("Active?")
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
    name = StringField("Name")
    genes = TagListField("Gene(s)")
    active = BooleanField("Active?")
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
    name = StringField("Name")
    genes = TagListField("Gene(s)")
    active = BooleanField("Active?")
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
        FormField(CalForm), widget=FieldListAddBtn(label="Add additional module")
    )
    nrps_type1 = FieldList(
        FormField(NRPS_I_Form), widget=FieldListAddBtn(label="Add additional module")
    )
    nrps_type6 = FieldList(
        FormField(NRPS_VI_Form), widget=FieldListAddBtn(label="Add additional module")
    )
    pks_iterative = FieldList(
        FormField(PKSIterativeForm),
        widget=FieldListAddBtn(label="Add additional module"),
    )
    pks_modular = FieldList(
        FormField(PKSModularForm), widget=FieldListAddBtn(label="Add additional module")
    )
    pks_trans_at = FieldList(
        FormField(PKSTransForm), widget=FieldListAddBtn(label="Add additional module")
    )
    other = FieldList(
        FormField(OtherForm), widget=FieldListAddBtn(label="Add additional module")
    )
    submit = SubmitField("Submit")
