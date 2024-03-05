from flask import (
    Flask,
    render_template,
    render_template_string,
    request,
    redirect,
    url_for,
    make_response,
    flash,
)
from forms.existing_bgc import SelectExisting
from forms.structure import StructureMultiple
from forms.min_entry import MinEntryForm
from forms.edit_select import EditSelectForm
from forms.biological_activity import BioActivityMultiple
from forms.biosynthesis import BioClassesCollection
from forms.common import is_valid_bgc_id
from forms.tailoring import TailoringForm
from rdkit import Chem, rdBase
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from werkzeug.datastructures import MultiDict
from pathlib import Path
import json
import csv

app = Flask(__name__)
app.secret_key = "IYKYK"


### Landing page
@app.route("/", methods=["GET", "POST"])
def index():
    """Main page, edit existing entry or start new entry"""
    form = SelectExisting(request.form)
    if request.method == "POST":
        # create new entry
        if form.submit.data:
            bgc_id = create_new_entry()
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))

        # edit valid existing entry
        if form.edit.data and form.validate():
            return redirect(url_for("edit_bgc", bgc_id=form.accession.data))

    return render_template("index.html", form=form)


## TODO: Submit new
@app.route("/WIP")
def work_in_progress():
    return render_template("WIP.html")


## Edit existing
@app.route("/edit/<bgc_id>", methods=["GET", "POST"])
def edit_bgc(bgc_id: str):

    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    form = EditSelectForm(request.form)
    if request.method == "POST":
        for option, chosen in form.data.items():
            if chosen:
                return redirect(url_for(f"edit_{option}", bgc_id=bgc_id))
    return render_template("edit.html", form=form, bgc_id=bgc_id)


# Edit minimal information
@app.route("/edit/<bgc_id>/minimal", methods=["GET", "POST"])
def edit_minimal(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    # try to fill data from existing entry
    if not request.form:
        form = MinEntryForm(MultiDict(read_data(bgc_id).get("Minimal")))
    else:
        form = MinEntryForm(request.form)

    if request.method == "POST" and form.validate():
        # TODO: save to database
        save_data(bgc_id, "Minimal", request.form)
        flash("Submitted minimal entry!")
        return redirect(url_for("edit_bgc", bgc_id=bgc_id))
    return render_template("min_entry.html", form=form, bgc_id=bgc_id)


# Edit structure information
@app.route("/edit/<bgc_id>/structure", methods=["GET", "POST"])
def edit_structure(bgc_id: str):
    """Main page for editing entry structure information"""
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    form = StructureMultiple(request.form)
    if request.method == "POST":
        if form.add.data:
            form.structures.append_entry()
            return render_template("structure.html", form=form, bgc_id=bgc_id)
        elif form.submit.data and form.validate():
            # return form.data  # TODO: save to db
            save_data(bgc_id, "Structure", request.form)
            flash("Submitted structure information!")
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))

    # on GET query db for any products already present
    else:
        # if one (or more) is present prefill
        try:
            products = next(
                csv.reader(
                    [MultiDict(read_data(bgc_id).get("Minimal")).get("products")],
                    skipinitialspace=True,
                )
            )
        except:
            products = [""]

        for product in products:
            form.structures.append_entry(data={"name": product})
    return render_template("structure.html", form=form, bgc_id=bgc_id)


@app.route("/render-smiles", methods=["POST"])
def render_smiles():
    """Render entered smiles string using rdkit"""
    # record which product smiles was entered
    origin = request.headers["Hx-Trigger-Name"]
    smiles = request.form.get(origin).strip()

    if not smiles:
        return ""

    # RDkit implementation
    try:
        mol = Chem.MolFromSmiles(smiles)
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(400, 200)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace("svg:", "")
    except Exception:
        return "Invalid SMILES"

    response = make_response(svg)
    response.content_type = "image/svg+xml"
    return response


## Bio Act
@app.route("/edit/<bgc_id>/bioact", methods=["GET", "POST"])
def edit_activity(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    form = BioActivityMultiple(request.form)
    if request.method == "POST":
        for activity in form.activities:
            if activity.add.data:
                activity.assay.append_entry()
                return render_template(
                    "biological_activity.html", bgc_id=bgc_id, form=form
                )
        if form.submit.data and form.validate():
            # return form.data # TODO: save to db
            save_data(bgc_id, "Bio_activity", request.form)
            flash("Submitted activity information!")
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))
    else:
        # prefill compounds
        try:
            products = products = next(
                csv.reader(
                    [MultiDict(read_data(bgc_id).get("Minimal")).get("products")],
                    skipinitialspace=True,
                )
            )
        except:
            products = [""]

        for product in products:
            form.activities.append_entry(data={"compound": product})
    return render_template("biological_activity.html", bgc_id=bgc_id, form=form)


## Biosynthesis
@app.route("/edit/<bgc_id>/biosynth", methods=["GET", "POST"])
def edit_biosynth(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    # if request.method == "POST":
    #     b_class = request.form.get("b_class")
    #     form = getattr(BioClassesCollection, b_class)(request.form)

    #     return form.data
    # return request.form

    # form = BiosyntheticClassesForm(request.form)
    return render_template("biosynthesis.html", bgc_id=bgc_id)


@app.route("/edit/<bgc_id>/biosynth/<b_class>", methods=["GET", "POST"])
def edit_biosynth_class(bgc_id: str, b_class: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    form = getattr(BioClassesCollection, b_class)(request.form)

    if request.method == "POST":

        # Add extra NRPS fields
        if form.data.get("add_release_type"):
            form.release_types.append_entry()
            return render_template(
                "biosynth_class_specific.html",
                form=form,
                b_class=b_class,
                bgc_id=bgc_id,
            )
        if form.data.get("add_thioesterase"):
            form.thioesterases.append_entry()
            return render_template(
                "biosynth_class_specific.html",
                form=form,
                b_class=b_class,
                bgc_id=bgc_id,
            )

        # Add extra Ribosomal fields
        if form.data.get("precursors"):
            if form.data.get("add_precursors"):
                form.precursors.append_entry()
                return render_template(
                    "biosynth_class_specific.html",
                    form=form,
                    b_class=b_class,
                    bgc_id=bgc_id,
                )
            for precursor in form.precursors:
                if precursor.data.get("add_crosslinks"):
                    precursor.crosslinks.append_entry()
                    return render_template(
                        "biosynth_class_specific.html",
                        form=form,
                        b_class=b_class,
                        bgc_id=bgc_id,
                    )

        # Add extra Saccharide fields
        if form.data.get("add_glycosyltransferase"):
            form.glycosyltransferases.append_entry()
            return render_template(
                "biosynth_class_specific.html",
                form=form,
                b_class=b_class,
                bgc_id=bgc_id,
            )
        if form.data.get("add_subcluster"):
            form.subclusters.append_entry()
            return render_template(
                "biosynth_class_specific.html",
                form=form,
                b_class=b_class,
                bgc_id=bgc_id,
            )

        # if no adds were triggered, validate and process data
        if form.validate():
            save_data(bgc_id, f"BioSynth_{b_class}", request.form)
            flash(f"Submitted {b_class} biosynthesis information!")
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))
            # return form.data

    return render_template(
        "biosynth_class_specific.html", form=form, b_class=b_class, bgc_id=bgc_id
    )


# @app.route("/get_class", methods=["POST"])
# def get_class():
#     if request.headers.get("Hx-Request"):
#         b_class = request.form.get("b_class")
#         form = getattr(BioClassesCollection, b_class)()
#         return render_template("forms_basic.html", form=form)


@app.route("/class_buttons/<bgc_id>", methods=["POST"])
def class_buttons(bgc_id: str):

    # grab classes for current bgc_id
    classes = MultiDict(read_data(bgc_id).get("Minimal")).getlist("b_class")

    if not classes:
        classes = ["NRPS", "PKS", "Ribosomal", "Saccharide", "Terpene", "Other"]
    class_btns = ""
    for cls in classes:
        # class_btns += f"<a href='/edit/{bgc_id}/biosynth/{cls}'><button class='btn btn-light'>{cls}</button><a>"
        class_btns += f"<a class='btn btn-light' style='margin: 5px' role='button' href='/edit/{bgc_id}/biosynth/{cls}'>{cls}<a>"
    return class_btns


## tailoring
@app.route("/edit/<bgc_id>/tailoring", methods=["GET", "POST"])
def edit_tailoring(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    form = TailoringForm(request.form)

    if request.method == "POST" and form.validate():
        save_data(bgc_id, "Tailoring", request.form)
        flash("Submitted tailoring information!")
        return redirect(url_for("edit_bgc", bgc_id=bgc_id))

    return render_template("tailoring.html", bgc_id=bgc_id, form=form)


## utils
# @app.route("/add_entry", methods=["POST"])
# def add_entry():
#     # return dict(request.headers)
#     origin = request.headers["Hx-Trigger"]
#     request.form.get(origin).add_entry()
#     return render_template("forms_basic.html", form=request.form)


@app.route("/add_evidence", methods=["POST"])
def addd():
    form = MinEntryForm(request.form)
    form.evidence.append_entry()

    return render_template_string(
        """{% import 'macros.html' as m %}
        <span>
        {{m.simple_divsubform(field, deletebtn=True)}}</span>""",
        field=form.evidence[-1],
    )


@app.route("/addtst", methods=["POST"])
def addtest():
    r = request
    return ""


@app.route("/delete", methods=["DELETE"])
def delete():
    return ""


@app.route("/submit", methods=["GET", "POST"])
def submit():
    return render_template("base.html")


def save_data(bgc_id: str, section_key: str, req_data: MultiDict):
    """Append data to file containing all answers for one BGC"""
    data = {section_key: [(k, v) for k in req_data for v in req_data.getlist(k)]}
    existing_data = read_data(bgc_id)

    existing_data.update(data)
    with open(f"{bgc_id}_data.json", "w") as outf:
        json.dump(existing_data, outf, sort_keys=True, indent=4)


def read_data(bgc_id: str):
    with open(f"{bgc_id}_data.json", "r") as inf:
        if content := inf.read():
            existing_data = json.loads(content)
        else:
            existing_data = {}
    return existing_data


def create_new_entry():
    max_entry_id = 0
    for filepath in Path(".").glob("new*"):
        if (nr := int(filepath.stem[3:6])) > max_entry_id:
            max_entry_id = nr
    bgc_id = f"new{max_entry_id+1:0>3}"
    Path(f"{bgc_id}_data.json").touch()

    return bgc_id
