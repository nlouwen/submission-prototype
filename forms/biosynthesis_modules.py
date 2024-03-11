from wtforms import Form, StringField, HiddenField, BooleanField, FormField, FieldList
from forms.common import TagListField
from forms.biosynthesis_domains import (
    CondensationDomain,
    AdenylationDomain,
    CarrierDomain,
    AminotransferaseDomain,
    CyclaseDomain,
    DehydrataseDomain,
    EnoylreductaseDomain,
    EpimeraseDomain,
    HydroxylaseDomain,
    KetoreductaseDomain,
    MethyltransferaseDomain,
    OtherDomain,
    OxidaseDomain,
    MonomerForm,
)


class CalForm(Form):
    _type = HiddenField("cal")
    # TODO: update?


class NRPS_I_Form(Form):
    # required _type, name, genes, active
    _type = HiddenField("nrps-type1")
    active = BooleanField("Active?")
    c_domain = FormField(CondensationDomain)
    a_domain = FormField(AdenylationDomain)
    carriers = FieldList(FormField(CarrierDomain))  # TODO: include add btn
    comments = StringField("Comments")
    genes = TagListField("Gene(s)")
    integrated_monomers = FieldList(FormField(MonomerForm))  # TODO: add btn
    modification_domains = (
        None  # TODO: array with oneOf 10 domain types, with their own schema
    )
    name = StringField("Name")


class NRPS_VI_Form(Form):
    _type = HiddenField("nrps-type6")
    # TODO: update?


class OtherForm(Form):
    _type = HiddenField("other")
    subtype = StringField("Subtype")
    # TODO: update?


class PKSIterativeForm(Form):
    _type = HiddenField("pks-iterative")
    # TODO: update?


class PKSModularForm(Form):
    _type = HiddenField("pks-modular")
    # TODO: update?


class PKSTransForm(Form):
    _type = HiddenField("pks-trans")
    # TODO: update?
