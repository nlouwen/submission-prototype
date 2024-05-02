from typing import Union

from flask import make_response, Response
from markupsafe import Markup
from rdkit import Chem
from rdkit.Chem import AllChem, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


def draw_smiles_svg(smiles: str) -> Union[str, Response]:
    """Draw svg image of a SMILES structure string

    Args:
        smiles (str): SMILES/CXSMILES representation of a structure

    Returns:
        Union[str, Response]: error message or structure svg
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        rdDepictor.Compute2DCoords(mol)
        drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace(
            "<svg", "<svg style='max-width:100%;height:100%'"
        )
    except Exception:
        return Markup(
            "<span class='invalid-feedback' style='display:block'>Invalid SMILES</span>"
        )

    response = make_response(svg)
    response.content_type = "image/svg+xml"
    return response


def draw_smarts_svg(smarts: str) -> Union[str, Response]:
    """Draw svg image of SMARTS reaction

    Args:
        smarts (str): SMARTS/CXSMARTS representation of a reaction

    Returns:
        Union[str, Response]: error message or reaction svg
    """
    try:
        rxn = AllChem.ReactionFromSmarts(smarts)
        drawer = rdMolDraw2D.MolDraw2DSVG(-1, -1)
        dopts = drawer.drawOptions()
        dopts.padding = 1e-5
        drawer.DrawReaction(rxn)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText().replace(
            "<svg", "<svg style='max-width:100%;height:100%'"
        )
    except Exception:
        return Markup(
            "<span class='invalid-feedback' style='display:block'>Invalid SMARTS</span>"
        )

    response = make_response(svg)
    response.content_type = "image/svg+xml"
    return response
