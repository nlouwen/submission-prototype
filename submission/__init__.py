import os
from typing import Optional

from flask import Flask


def create_app(test_config: Optional[dict] = None) -> Flask:
    """Flask app factory to create a Flask app instance

    Args:
        test_config (Optional[dict], optional): app config during tests. Defaults to None.

    Returns:
        Flask: Flask app instance
    """

    app = Flask(__name__, instance_relative_config=True)
    app = configure_app(app, test_config)

    return app


def configure_app(app: Flask, test_config: Optional[dict] = None) -> Flask:
    """Add config settigns to app

    Args:
        app (Flask): Flask app instance
        test_config (Optional[dict], optional): app config during tests. Defaults to None.

    Returns:
        Flask: Configures Flask app instance
    """
    app.config["SECRET_KEY"] = "IYKYK"

    if test_config:
        app.config.from_mapping(test_config)
    else:
        app.config.from_pyfile("config.py", silent=True)

    return app
