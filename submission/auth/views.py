from typing import Union

from flask import render_template, request, redirect, url_for, flash, abort, current_app
from flask_login import login_user, logout_user, login_required
from flask_mail import Message
from werkzeug.wrappers import response

from submission import mail
from submission.auth import bp_auth
from submission.models import User, Token
from .forms.login import LoginForm, UserEmailForm, PasswordResetForm


@bp_auth.route("/login", methods=["GET"])
def login() -> str:
    """Renders login form"""
    form = LoginForm(request.form)
    return render_template("auth/login.html.j2", form=form)


@bp_auth.route("/login", methods=["POST"])
def login_post() -> response.Response:
    """Handles POST requests to the login page to log in users"""
    email = request.form.get("username")
    password = request.form.get("password")
    remember = True if request.form.get("remember") else False

    user: User = User.query.filter(User.email.ilike(email)).first()

    if not user or not user.check_password(password):
        flash("Please check your login details")
        return redirect(url_for("auth.login"))

    login_user(user, remember)

    return redirect(url_for("main.index"))


@bp_auth.route("/logout")
@login_required
def logout() -> response.Response:
    """Logs out current user and redirects to the login page"""
    logout_user()
    return redirect(url_for("auth.login"))


@bp_auth.route("/reset-my-password", methods=["GET", "POST"])
def password_email() -> Union[str, response.Response]:
    """Send an email with password reset link to user provided email"""
    form = UserEmailForm(request.form)
    if request.method == "POST" and form.validate():
        email = form.email.data
        user = User.query.filter_by(email=email).first()

        if not user:
            flash("Unknown email address")
            return redirect(url_for("auth.password_email"))

        token_id = Token.generate_token(user.id, "password_reset")
        # TODO: send email
        mail.send(
            Message(
                subject="Change your MIBiG password",
                recipients=[email],
                body=f"Hello, click this link {current_app.config['BASE_URL']}/auth/reset/{token_id}",
            )
        )
        flash(f"Please check your email")
        return redirect(url_for("auth.login"))

    return render_template("auth/pw_reset_request.html", form=form)


@bp_auth.route("/reset/<token_id>", methods=["GET", "POST"])
def reset_password(token_id: str) -> Union[str, response.Response]:
    """Allow a user to change their password via email provided link

    Arguments:
        token_id (str): uuid token
    """
    token: Token = Token.query.filter_by(token_id=token_id).first()

    if not token or token.purpose != "password_reset":
        abort(403, "Invalid link for password reset")

    if not token.is_created_within(hours=2):
        abort(403, "Token has expired")

    form = PasswordResetForm(request.form)
    if request.method == "POST" and form.validate():
        user = User.query.filter_by(id=token.user_id).first()
        user.password = form.password.data

        token.cleanup_tokens()

        flash("Successfully changed password")
        return redirect(url_for("auth.login"))

    return render_template("auth/pw_reset.html", form=form)
