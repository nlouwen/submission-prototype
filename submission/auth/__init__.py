from flask import Blueprint
from .decorators import auth_role

bp_auth = Blueprint("auth", __file__, url_prefix="/auth")

from submission.auth import views

__all__ = ["auth_role"]
