#!/usr/bin/env python3
# Set up the basic database contents


from submission import create_app, db
from submission.models import Role


def main():
    app = create_app()

    with app.app_context():
        for slug in ("submitter", "reviewer", "admin"):
            role = Role.query.filter_by(slug=slug).first()
            if role is None:
                role = Role(name=f"MIBiG {slug.capitalize()}s", slug=slug)
                db.session.add(role)

        db.session.commit()


if __name__ == "__main__":
    main()
