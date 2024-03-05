from wtforms import Form, StringField, SubmitField, validators, ValidationError
from forms.common import is_valid_bgc_id


def valid_id(form, field):
    # grab existing ids from MIBiG
    if not is_valid_bgc_id(field.data):
        raise ValidationError("Invalid MIBiG ID!")


class SelectExisting(Form):
    accession = StringField(
        "MIBiG ID",
        [validators.DataRequired(), valid_id],
        description="Unique MIBiG identifier starting with 'BGC' or a unique temporary identifier",
    )
    edit = SubmitField("Edit existing entry")
    submit = SubmitField("Submit new entry", render_kw={"formnovalidate": True})
