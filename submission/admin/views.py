from flask import request, render_template, redirect, url_for
from sqlalchemy import or_


from submission.admin import bp_admin
from submission.auth import auth_role
from submission.extensions import db
from submission.models import Role, User, UserInfo

from submission.admin.forms import UserAdd, UserEdit

@bp_admin.route("/")
def index() -> str:
    return render_template("admin/index.html.j2")


@bp_admin.route("/users", methods=["GET"])
def list_users() -> str:
    users = User.query.join(UserInfo).order_by(UserInfo.name).all()
    return render_template("admin/users.html.j2", users=users)

@bp_admin.route("/users", methods=["POST"])
def search_users() -> str:
    search = request.form["search"]
    users = User.query.join(UserInfo) \
        .filter(or_(UserInfo.name.like(f"{search}%"), User.email.like(f"{search}%"))) \
        .order_by(UserInfo.name).all()
    return render_template("admin/user_search.html.j2", users=users)


@bp_admin.route("/user/<user_id>", methods=["GET"])
def user(user_id: int) -> str:
    user = User.query.get_or_404(user_id)
    return render_template("admin/user_list_line.html.j2", user=user)


@bp_admin.route("/user/<user_id>/edit", methods=["GET", "PUT"])
def user_edit(user_id: int) -> str:
    user = User.query.get_or_404(user_id)
    all_roles = Role.query.order_by(Role.slug).all()
    form = UserEdit(request.form, data={
        "email": user.email,
    })
    form.roles.choices = [(role.slug, role.name) for role in all_roles]

    if form.validate_on_submit():
        user.email = form.email.data
        user.active = form.active.data
        roles = []
        for wanted_role in form.roles.data:
            for role in all_roles:
                if role.slug == wanted_role:
                    roles.append(role)
        user.roles = roles
        db.session.add(user)
        db.session.commit()
        return render_template("admin/user_list_line.html.j2", user=user)

    form.roles.data = [role.slug for role in user.roles]
    form.active.data = user.active

    return render_template("admin/user_edit_form.html.j2", form=form, user=user)


@bp_admin.route("/user/new", methods=["GET", "POST"])
def user_create() -> str:
    form = UserAdd(request.form)
    all_roles = Role.query.order_by(Role.slug).all()
    form.roles.choices = [(role.slug, role.name) for role in all_roles]

    if form.validate_on_submit():
        user = User(email=form.email.data, active=form.active.data, _password="deactivated")
        db.session.add(user)
        db.session.commit()
        name = form.name.data
        info = UserInfo(
            id=user.id,
            alias=UserInfo.generate_alias(),
            name=name,
            call_name=UserInfo.guess_call_name(name),
            organisation=form.affiliation.data
        )
        db.session.add(info)
        for wanted_role in form.roles.data:
            for role in all_roles:
                if role.slug == wanted_role:
                    user.roles.append(role)
        db.session.add(user)
        db.session.commit()
        return redirect(url_for("admin.list_users"))

    return render_template("admin/user_add.html.j2", form=form)
