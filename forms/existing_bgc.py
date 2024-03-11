from wtforms import (
    Form,
    StringField,
    SubmitField,
    validators,
    ValidationError,
    FormField,
    HiddenField,
)
from wtforms.widgets import SubmitInput
from markupsafe import Markup
from forms.common import is_valid_bgc_id, StringFieldAddBtn


def valid_id(form, field):
    # grab existing ids from MIBiG
    if not is_valid_bgc_id(field.data):
        raise ValidationError("Invalid MIBiG ID!")


class SelectExisting(Form):
    accession = StringField(
        "MIBiG ID",
        [validators.DataRequired(), valid_id],
        description="Unique MIBiG identifier starting with 'BGC' or a unique temporary identifier",
        widget=StringFieldAddBtn(
            label="Edit existing entry",
            render_kw={
                "name": "edit",
                "value": "edit",
                "type": "submit",
                "style": "margin-top:5px",
            },
        ),
    )
    # edit = SubmitField("Edit existing entry")
    submit = SubmitField("Submit new entry", render_kw={"formnovalidate": True})
