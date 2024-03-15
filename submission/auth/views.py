from flask import render_template, request, redirect, url_for, flash
from flask_login import login_user, logout_user, login_required

from submission.auth import bp_auth
from submission.models import User
from .forms.login import LoginForm, UserEmailForm, PasswordResetForm


@bp_auth.route("/login", methods=["GET"])
def login():
    form = LoginForm(request.form)
    return render_template("auth/login.html.j2", form=form)


@bp_auth.route("/login", methods=["POST"])
def login_post():
    email = request.form.get("email")
    password = request.form.get("password")
    remember = True if request.form.get("remember") else False

    user = User.query.filter_by(email=email).first()

    if not user or not user.check_password(password):
        flash("Please check your login details")
        return redirect(url_for("auth.login"))

    login_user(user, remember)

    return redirect(url_for("main.profile"))


@bp_auth.route("/logout")
@login_required
def logout():
    logout_user()
    return redirect(url_for("auth.login"))


@bp_auth.route("/reset-my-password", methods=["GET", "POST"])
def password_email():
    """Send an email with pw reset link to user provided email"""

    form = UserEmailForm(request.form)
    if request.method == "POST" and form.validate():
        # TODO: generate token, add email+token to database
        # TODO: send email
        flash("Please check your email")
        return redirect(url_for("auth.login"))

    return render_template("auth/pw_reset_request.html", form=form)


@bp_auth.route("/reset/<token>", methods=["GET", "POST"])
def reset_password(token: str):
    """Allow a user to change their password via email provided link"""

    form = PasswordResetForm(request.form)
    if request.method == "POST" and form.validate():
        # TODO: find user from token, check time past
        # TODO: create new user+pw / change pw of existing user
        return redirect(url_for("auth.login"))

    return render_template("auth/pw_reset.html", form=form)
