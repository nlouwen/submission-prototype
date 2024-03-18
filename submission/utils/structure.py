from typing import Union

from flask import make_response, Response
from rdkit import Chem, rdBase
from rdkit.Chem import Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D


def draw_smiles_svg(smiles: str) -> Union[str, Response]:
    """Draw svg image of a SMILES structure string

    Args:
        smiles (str): SMILES representation of a structure

    Returns:
        str | Response: error message or structure svg
    """
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
