# MIBiG Submission prototype

- Clone or download repo

`git clone git@github.com:nlouwen/submission-prototype.git`

- Install required packages

`conda env create --name submission -f requirements.txt`

`conda activate submission`

- Set up the database

`flask --app submission db upgrade`

- Populate the database with example users

Open a shell by running `flask --app submission shell`, then past the following:

```
from submission.models import User, Role, UserRole
from submission import db
alice = User(name="Alice", email="alice@example.edu", password="secret")
bob = User(name="Bob", email="bob@example.edu", password="secret")
submitter = Role(name="MIBiG Submitters", slug="submitter")
reviewer = Role(name="MIBiG Reviewers", slug="reviewer")
db.session.add_all([alice, bob, submitter, reviewer])
db.session.commit()
db.session.add_all([UserRole(user_id=alice.id, role_id=submitter.id), UserRole(user_id=alice.id, role_id=reviewer.id), UserRole(user_id=bob.id, role_id=submitter.id)])
db.session.commit()
exit()
```

- Run the flask app

`flask --app submission run --debug`

- Open the submission prototype in your browser at `localhost:5000` or `http://127.0.0.1:5000`
