from flask import Blueprint
from flask_login import current_user, login_required

from submission.auth import auth_role

bp_admin = Blueprint("admin", __file__, url_prefix="/admin")

@bp_admin.before_request
@login_required
@auth_role("admin")
def handle_permissions():
    # just here to enforce admin permissions on all views.
    pass

from submission.admin import views
