from wtforms import (
    Form,
    StringField,
    SubmitField,
    validators,
)


class UserDetailsEditForm(Form):
    name = StringField(
        "Name",
        validators=[
            validators.InputRequired(message="Please provide a name")
        ],
    )
    call_name = StringField(
        "What should we call you?",
        validators=[
            validators.InputRequired(message="Please provide a way to address you")
        ],
    )
    orcid = StringField("ORCID")
    organisation = StringField(
        "First affiliation",
        validators=[
            validators.InputRequired(message="Please provide an affiliation")
        ],
    )
    organisation_2 = StringField("Second affiliation (optional)")
    organisation_3 = StringField("Third affiliation (optional)")
    submit = SubmitField("Update user details")
