# MIBiG Submission prototype

- Clone or download repo

`git clone git@github.com:nlouwen/submission-prototype.git`

- Create a virtual python (>v3.11) environment (e.g. with `micromamba` or `conda`:)

`conda create --name submission python==3.11`

`conda activate submission`

- Install the submission app prototype

`pip install -e .`

- Set up the database

`flask --app submission db upgrade`

- Seed the database with the MIBiG roles

`./seed.py`

- Run the flask app

`flask --app submission run --debug`

- Open the submission prototype in your browser at `localhost:5000` or `http://127.0.0.1:5000`
