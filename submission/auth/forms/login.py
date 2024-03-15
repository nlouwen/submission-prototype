from wtforms import (
    Form,
    PasswordField,
    EmailField,
    BooleanField,
    SubmitField,
    validators,
)


class LoginForm(Form):
    username = EmailField("Username", render_kw={"placeholder": "alice@example.edu"})
    password = PasswordField("Password")
    remember = BooleanField("Remember me")
    submit = SubmitField("Login")


class UserEmailForm(Form):
    username = EmailField(
        "Email",
        validators=[
            validators.InputRequired(message="Please provide an email address")
        ],
        render_kw={"placeholder": "alice@example.edu"},
    )
    submit = SubmitField("Send email")


class PasswordResetForm(Form):
    password = PasswordField("New password")
    confirm = PasswordField(
        "Confirm password",
        validators=[validators.EqualTo("password", "Password mismatch!")],
    )
    submit = SubmitField("Save password")
