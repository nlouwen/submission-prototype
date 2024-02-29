from flask import (
    Flask,
    render_template,
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

app = Flask(__name__)
app.secret_key = "IYKYK"


### Landing page
@app.route("/", methods=["GET", "POST"])
def index():
    form = SelectExisting(request.form)
    if request.method == "POST" and form.validate():
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

    if not request.form:
        form = MinEntryForm(
            MultiDict(
                {
                    "genome": "ABG100001.1",
                    "location-start": 1,
                    "location-end": 4,
                    "products": ["FAKEOMYCIN"],
                }
            )
        )
    else:
        form = MinEntryForm(request.form)
    if request.method == "POST":
        if form.add_evidence.data:
            form.evidence.append_entry()
            return render_template("min_entry.html", form=form, bgc_id=bgc_id)

        if form.validate():
            # return form.data  # TODO: save to db
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
            flash("Submitted structure information!")
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))

    # on GET query db for any products already present
    else:
        # if one (or more) is present prefill
        for c_name, smiles in [
            (
                "EXAMPLEMYCIN",
                "C[C@@H](CCC=C(C)C)[C@@H]1CC[C@]2([C@]1(CC[C@H]3C2=CC[C@@H]4[C@@]3(CC[C@@H](C4(C)C)O)C)C)C",
            )
        ]:
            form.structures.append_entry(data={"name": c_name, "structure": smiles})
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
            flash("Submitted activity information!")
            return redirect(url_for("edit_bgc", bgc_id=bgc_id))
    else:
        # prefill compounds
        for c_name in ["a", "b"]:
            form.activities.append_entry(data={"compound": c_name})
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
    classes = ["NRPS", "PKS", "Ribosomal", "Saccharide", "Terpene", "Other"]
    class_btns = ""
    for cls in classes:
        class_btns += (
            f"<a href='/edit/{bgc_id}/biosynth/{cls}'><button>{cls}</button><a>"
        )
    return class_btns


## tailoring
@app.route("/edit/<bgc_id>/tailoring", methods=["GET", "POST"])
def edit_tailoring(bgc_id: str):
    if not is_valid_bgc_id(bgc_id):
        return "Invalid existing entry!", 404

    form = TailoringForm(request.form)

    if request.method == "POST" and form.validate():
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


@app.route("/delete", methods=["DELETE"])
def delete():
    return ""


@app.route("/submit", methods=["GET", "POST"])
def submit():
    return render_template("base.html")
