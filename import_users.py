#!/usr/bin/env python3
# bulk-import users from a mibig3 database dump

from argparse import ArgumentParser, FileType
from csv import DictReader
from dataclasses import dataclass
from typing import Self, TextIO

from sqlalchemy import or_

from submission import create_app, db
from submission.models import User, Role, UserInfo, UserRole

@dataclass
class LegacyUser:
    alias: str
    email: str
    name: str
    call_name: str
    org: str
    public: bool
    active: bool
    orcid: str | None = None
    org_2: str | None = None
    org_3: str | None = None
    reviewer: bool = False

    @classmethod
    def from_legacy(cls, line: str) -> Self:
        parts = line.split("\t")
        if len(parts) != 10:
            raise ValueError(f"Invalid line split: {parts}")

        if "," in parts[2]:
            raise ValueError(f"comma in user name for {parts[2]}")

        return cls(
            alias=parts[0],
            email=parts[1],
            name=parts[2],
            call_name=parts[3],
            org=parts[4],
            # skip password
            public=True if parts[6] == "t" else False,
            # skip gdpr_consent field
            active=True if parts[8] == "t" else False,
        )

    @classmethod
    def from_csv(cls, row: dict[str, str]) -> Self:
        first_and_initials = row["First Name + Initials"].strip()
        call_name = first_and_initials.split(" ")[0]
        last = row["Surname"].strip()
        email = row["E-mail Address"].strip()
        reviewer = row["Reviewer"] == "TRUE"
        public = False
        active = True
        orcid = row["ORCID"].strip() if row["ORCID"] else None
        org1 = row["Affliation/Institution 1 "].strip()
        org2 = row["Affliation/Institution 2"].strip() if row["Affliation/Institution 2"] else None
        org3 = row["Affliation/Institution 3"].strip() if row["Affliation/Institution 3"] else None

        return cls(
            alias="BBBBBBBBBBBBBBBBBBBBBBBB",
            email=email,
            name=f"{first_and_initials} {last}",
            call_name=call_name,
            orcid=orcid,
            org=org1,
            org_2=org2,
            org_3=org3,
            public=public,
            active=active,
            reviewer=reviewer,
        )




def main():
    parser = ArgumentParser()
    parser.add_argument("dumpfile", type=FileType('r', encoding="utf-8"),
                        help="File to read users from")
    parser.add_argument("-m", "--mode", type=str, choices=("legacy", "csv"), default="csv",
                        help="Import the user from legacy database dump or new csv file")
    args = parser.parse_args()

    users = parse_users(args.dumpfile, args.mode)

    load_users(users)


def parse_users(handle: TextIO, mode: str) -> list[LegacyUser]:
    users: list[LegacyUser] = []

    if mode == "legacy":
        for line in handle:
            users.append(LegacyUser.from_legacy(line))
    elif mode == "csv":
        reader = DictReader(handle)
        for row in reader:
            users.append(LegacyUser.from_csv(row))

    return users


def load_users(users: list[LegacyUser]):
    app = create_app()
    with app.app_context():
        reviewers = Role.query.filter_by(slug="reviewer").one()
        for user in users:

            if user.orcid is not None:
                query = User.query.join(UserInfo).filter(or_(
                    User.email == user.email,
                    UserInfo.name == user.name,
                    UserInfo.orcid == user.orcid,
                ))
            else:
                query = User.query.join(UserInfo).filter(or_(
                    User.email == user.email,
                    UserInfo.name == user.name,
                ))

            existing = query.first()
            if existing:
                existing.email = user.email
                existing.info.name = user.name
                existing.info.call_name = user.call_name
                existing.info.organisation = user.org
                existing.info.organisation_2 = user.org_2
                existing.info.organisation_3 = user.org_3
                existing.info.orcid = user.orcid
                existing.active = True
                db.session.add(existing)
                db.session.commit()
                print("updated", existing.info.name)
            else:
                load_user(user)

            if user.reviewer:
                existing = User.query.filter_by(email=user.email).first()
                if not existing:
                    continue
                db.session.add(UserRole(user_id=existing.id, role_id=reviewers.id))
                db.session.commit()


def load_user(legacy: LegacyUser):
    user = User(email=legacy.email, active=legacy.active, _password="deactivated")
    db.session.add(user)
    db.session.commit()

    if legacy.alias == "BBBBBBBBBBBBBBBBBBBBBBBB":
        legacy.alias = UserInfo.generate_alias()

    info = UserInfo(
        id=user.id,
        alias=legacy.alias,
        name=legacy.name,
        call_name=legacy.call_name,
        organisation=legacy.org,
        public=legacy.public,
        orcid=legacy.orcid,
        organisation_2=legacy.org_2,
        organisation_3=legacy.org_3,
    )
    db.session.add(info)
    db.session.commit()


if __name__ == "__main__":
    main()
