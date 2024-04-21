from argparse import ArgumentParser
from pathlib import Path

from submission import create_app
from submission.models import NPAtlas


def main():
    parser = ArgumentParser()
    parser.add_argument("dumpfile", type=Path, help="NPAtlas tsv dump file")
    args = parser.parse_args()

    app = create_app()
    with app.app_context():
        NPAtlas.from_tsv_dump(args.dumpfile)


if __name__ == "__main__":
    main()
