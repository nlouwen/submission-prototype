import os
from typing import Optional

from flask import Flask
from dotenv import load_dotenv

from submission.extensions import db, migrate, login_manager, mail
from submission.main import bp_main
from submission.edit import bp_edit
from submission.auth import bp_auth
from submission.models import User


def create_app(test_config: Optional[dict] = None) -> Flask:
    """Flask app factory to create a Flask app instance

    Args:
        test_config (Optional[dict], optional): app config during tests. Defaults to None.

    Returns:
        Flask: Flask app instance
    """

    app = Flask(__name__, instance_relative_config=True)
    app = configure_app(app, test_config)

    app = register_blueprints(app)

    db.init_app(app)
    migrate.init_app(app, db)
    mail.init_app(app)
    login_manager.init_app(app)

    @login_manager.user_loader
    def load_user(user_id):
        return User.query.get(int(user_id))

    return app


def configure_app(app: Flask, test_config: Optional[dict] = None) -> Flask:
    """Add config settings to app

    Args:
        app (Flask): Flask app instance
        test_config (Optional[dict], optional): app config during tests. Defaults to None.

    Returns:
        Flask: Configures Flask app instance
    """
    load_dotenv()
    app.config["SECRET_KEY"] = "IYKYK"
    app.config["SQLALCHEMY_DATABASE_URI"] = "sqlite:///auth.sqlite3"

    app.config["MAIL_SERVER"] = os.getenv("MAIL_SERVER")
    app.config["MAIL_PORT"] = os.getenv("MAIL_PORT")
    app.config["MAIL_USERNAME"] = os.getenv("MAIL_USERNAME")
    app.config["MAIL_PASSWORD"] = os.getenv("MAIL_PASSWORD")
    app.config["MAIL_DEFAULT_SENDER"] = os.getenv("MAIL_DEFAULT_SENDER")
    app.config["MAIL_USE_TLS"] = True
    app.config["MAIL_USE_SSL"] = False

    print(app.config)
    if test_config:
        app.config.from_mapping(test_config)
    else:
        app.config.from_pyfile("config.py", silent=True)

    return app


def register_blueprints(app: Flask) -> Flask:
    """Register all blueprints

    Args:
        app (Flask): Flask app instance

    Returns:
        Flask: Flask app instance with registered blueprints
    """
    app.register_blueprint(bp_main)
    app.register_blueprint(bp_edit)
    app.register_blueprint(bp_auth)
    return app
