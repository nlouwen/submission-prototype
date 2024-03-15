#!/usr/bin/env python3
# bulk-import users from a mibig3 database dump

from argparse import ArgumentParser, FileType
from dataclasses import dataclass
from typing import Self, TextIO

from submission import create_app, db
from submission.models import User, Role, UserInfo

@dataclass
class LegacyUser:
    alias: str
    email: str
    name: str
    call_name: str
    org: str
    public: bool
    active: bool

    @classmethod
    def from_line(cls, line: str) -> Self:
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



def main():
    parser = ArgumentParser()
    parser.add_argument("dumpfile", type=FileType('r', encoding="utf-8"))
    args = parser.parse_args()

    users = parse_users(args.dumpfile)

    load_users(users)


def parse_users(handle: TextIO) -> list[LegacyUser]:
    users: list[LegacyUser] = []
    for line in handle:
        users.append(LegacyUser.from_line(line))

    return users


def load_users(users: list[LegacyUser]):
    app = create_app()
    with app.app_context():
        for user in users:
            load_user(user)


def load_user(legacy: LegacyUser):
    user = User(email=legacy.email, active=legacy.active, _password="deactivated")
    db.session.add(user)
    db.session.commit()
    info = UserInfo(
        id=user.id,
        alias=legacy.alias,
        name=legacy.name,
        call_name=legacy.call_name,
        organisation=legacy.org,
        public=legacy.public,
    )
    db.session.add(info)
    db.session.commit()


if __name__ == "__main__":
    main()
