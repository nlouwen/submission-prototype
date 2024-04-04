from flask_wtf import FlaskForm
from wtforms import (
    BooleanField,
    EmailField,
    SelectMultipleField,
    StringField,
    SubmitField,
    validators,
)

class UserAdd(FlaskForm):
    email = EmailField("Email", validators=[validators.InputRequired(message="Email address required")])
    name = StringField("Name", validators=[validators.InputRequired(message="Name is required")])
    affiliation = StringField("Affiliation", validators=[validators.InputRequired(message="Affiliation is required")])
    active = BooleanField("Active")
    roles = SelectMultipleField("Roles")
    submit = SubmitField("Edit user")


class UserEdit(FlaskForm):
    email = EmailField("Email", validators=[validators.InputRequired(message="Email address required")])
    active = BooleanField("Active")
    roles = SelectMultipleField("Roles")
