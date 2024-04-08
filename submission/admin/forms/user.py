from flask_wtf import FlaskForm
from wtforms import (
    BooleanField,
    EmailField,
    SelectMultipleField,
    StringField,
    SubmitField,
    validators,
)

from submission.models import User


def email_unique(form, field):
    email = field.data
    user = User.query.filter(User.email.like(email)).first()
    if user is not None:
        raise validators.ValidationError(f"{email} already exists in the user database")


class UserAdd(FlaskForm):
    email = EmailField("Email", validators=[validators.InputRequired(message="Email address required"), email_unique])
    name = StringField("Name", validators=[validators.InputRequired(message="Name is required")])
    affiliation = StringField("Affiliation", validators=[validators.InputRequired(message="Affiliation is required")])
    active = BooleanField("Active")
    roles = SelectMultipleField("Roles")
    submit = SubmitField("Edit user")


class UserEdit(FlaskForm):
    email = EmailField("Email", validators=[validators.InputRequired(message="Email address required"), email_unique])
    active = BooleanField("Active")
    roles = SelectMultipleField("Roles")
