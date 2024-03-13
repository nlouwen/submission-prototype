from functools import wraps

from flask import redirect, url_for, flash
from flask_login import current_user

def auth_role(role):
    def wrapper(fn):
        @wraps(fn)
        def decorator(*args, **kwargs):
            roles = role if isinstance(role, list) else [role]
            if all(not current_user.has_role(r) for r in roles):
                flash("You don't have permissions to access this content!", "error")
                return redirect(url_for("profile"))
            return fn(*args, **kwargs)

        return decorator

    return wrapper
