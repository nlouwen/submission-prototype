from flask import Blueprint, render_template, request, redirect, url_for, flash
from flask_login import login_user, logout_user, login_required

from models import User

auth = Blueprint("auth", __file__, url_prefix="/auth", template_folder="app/auth/templates")

@auth.route("/login", methods=["GET"])
def login():
    return render_template("auth/login.html.j2")


@auth.route("/login", methods=["POST"])
def login_post():
    email = request.form.get("email")
    password = request.form.get("password")
    remember = True if request.form.get("remember") else False

    user = User.query.filter_by(email=email).first()

    if not user or not user.check_password(password):
        flash("Please check your login details")
        return redirect(url_for("auth.login"))

    login_user(user, remember)

    return redirect(url_for("profile"))


@auth.route("/logout")
@login_required
def logout():
    logout_user()
    return redirect(url_for("auth.login"))
