from pathlib import Path

from flask import render_template, request, redirect, url_for
from flask_login import login_required, current_user

from submission.main import bp_main
from submission.main.forms import SelectExisting
from submission.auth import auth_role
from submission.utils import Storage


@bp_main.route("/", methods=["GET", "POST"])
@login_required
def index():
    """Main page, edit existing entry or start new entry"""
    form = SelectExisting(request.form)
    if request.method == "POST":
        # create new entry
        if form.submit.data:
            bgc_id = Storage.create_new_entry()
            return redirect(url_for("edit.edit_bgc", bgc_id=bgc_id))

        # edit valid existing entry
        if request.form.get("edit") and form.validate():
            bgc_id = form.accession.data
            Storage.create_entry_if_not_exists(bgc_id)

            return redirect(url_for("edit.edit_bgc", bgc_id=bgc_id))

    return render_template("main/index.html", form=form)


@bp_main.route("/delete", methods=["DELETE"])
def delete() -> str:
    """Dummy route to delete any target

    Returns:
        str: Dummy value
    """
    return ""


@bp_main.route("/profile")
@login_required
def profile():
    roles = [role.slug for role in current_user.roles]
    return render_template(
        "main/profile.html.j2", name=current_user.info.name, roles=roles
    )


@bp_main.route("/submitter")
@login_required
@auth_role("submitter")
def submitter():
    return render_template("main/submitter.html.j2", name=current_user.info.name)


@bp_main.route("/reviewer")
@login_required
@auth_role("reviewer")
def reviewer():
    return render_template("main/reviewer.html.j2", name=current_user.info.name)
