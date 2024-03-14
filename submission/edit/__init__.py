from flask import Blueprint

bp_edit = Blueprint("edit", __file__, url_prefix="/edit")

from submission.edit import views
