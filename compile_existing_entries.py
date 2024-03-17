import json
from argparse import ArgumentParser
from pathlib import Path


def main():
    """looks for all BGC*.json files in dir and dumps converted jsons in current dir"""
    parser = ArgumentParser()
    parser.add_argument("mibig_json_dir", type=Path)
    args = parser.parse_args()

    output_dir = Path(__file__).parent / "data"
    if not output_dir.exists():
        output_dir.mkdir()

    load_existing_entries(args.mibig_json_dir, output_dir)


def load_existing_entries(data_dir, output_dir):
    for entry in data_dir.glob("BGC*.json"):
        compile_entry(entry, output_dir)


def compile_entry(entry: Path, output_dir: Path):
    full_data: dict = json.load(open(entry))
    data = full_data.get("cluster")

    if data is None:
        return

    entry_data = {
        "Minimal": convert_minimal(data),
        "Structure": convert_structure(data),
        "Bio_activity": convert_bioact(data),
        "Annotation": convert_annotation(data),
    }

    outf = output_dir / f"{entry.stem}_data.json"
    with open(outf, "w") as fp:
        json.dump(entry_data, fp, sort_keys=True, indent=4)


def convert_minimal(data: dict):
    min_data: list[list[str]] = []

    if (locus := data.get("loci")) is not None:
        min_data.append([f"loci-0-genome", locus.get("accession", "")])
        min_data.append([f"loci-0-location-start", str(locus.get("start_coord", ""))])
        min_data.append([f"loci-0-location-end", str(locus.get("end_coord", ""))])

        min_data.append(["completeness", locus.get("completeness", "")])

        for e_idx, evidence in enumerate(locus.get("evidence", "")):
            min_data.append([f"loci-0-evidence-{e_idx}-method", evidence])

    if (compounds := data.get("compounds")) is not None:
        min_data.append(
            [
                "products",
                '"' + '", "'.join([c.get("compound") for c in compounds]) + '"',
            ]
        )

    min_data.append(["taxonomy-ncbi_tax_id", data.get("ncbi_tax_id", "")])

    if (org := data.get("organism_name")) is None:
        genus, species, strain = "", "", ""
    else:
        genus, species, strain = org.split(" ", 2)
    min_data.append(["taxonomy-genus", genus])
    min_data.append(["taxonomy-species", species])
    min_data.append(["taxonomy-strain", strain])

    if (classes := data.get("biosyn_class")) is not None:
        conversion = {
            "Alkaloid": "",  # TODO: alkaloid conversion
            "Polyketide": "PKS",
            "RiPP": "Ribosomal",
            "NRP": "NRPS",
            "Saccharide": "Saccharide",
            "Terpene": "Terpene",
            "Other": "Other",
        }
        for biosynth_class in classes:
            min_data.append(["b_class", conversion.get(biosynth_class, "")])

    # TODO: comments?
    return min_data


def convert_structure(data: dict):
    str_data: list[list[str]] = []

    compounds = data.get("compounds")

    if compounds is None:
        return str_data

    for c_idx, compound in enumerate(compounds):
        str_data.append([f"structures-{c_idx}-name", compound.get("compound", "")])
        str_data.append(
            [f"structures-{c_idx}-synonyms", compound.get("chem_synonyms", "")]
        )
        str_data.append(
            [f"structures-{c_idx}-formula", compound.get("molecular_formula", "")]
        )
        str_data.append([f"structures-{c_idx}-mass", str(compound.get("mol_mass", ""))])
        str_data.append(
            [f"structures-{c_idx}-structure", compound.get("chem_struct", "")]
        )

        str_data.append(
            [
                f"structures-{c_idx}-moieties",
                '"'
                + '", "'.join(
                    [m.get("moiety") for m in compound.get("chem_moieties", "")]
                )
                + '"',
            ]
        )

        str_data.append(
            [
                f"structures-{c_idx}-db_cross",
                '"'
                + '", "'.join([db_id for db_id in compound.get("database_id", "")])
                + '"',
            ]
        )
    return str_data


def convert_bioact(data: dict):
    act_data: list[list[str]] = []

    compounds = data.get("compounds")

    if compounds is None:
        return act_data

    for c_idx, compound in enumerate(compounds):
        act_data.append([f"activities-{c_idx}-compound", compound.get("compound", "")])

        chem_acts = compound.get("chem_acts")

        if chem_acts is None:
            break

        for a_idx, chem_act in enumerate(chem_acts):
            act_data.append(
                [
                    f"activities-{c_idx}-assays-{a_idx}-target",
                    chem_act.get("activity", ""),
                ]
            )
    return act_data


def convert_annotation(data: dict):
    annot_data: list[list[str]] = []

    genes = data.get("genes")

    if genes is None:
        return annot_data

    annotations = genes.get("annotations")

    if annotations is None:
        return annot_data

    for a_idx, annotation in enumerate(annotations):
        annot_data.append([f"annotations-{a_idx}-gene_id", annotation.get("id", "")])
        annot_data.append([f"annotations-{a_idx}-name", annotation.get("name", "")])
        annot_data.append(
            [f"annotations-{a_idx}-product", annotation.get("product", "")]
        )

        functions = annotation.get("functions")

        if functions is None:
            break

        for f_idx, function in enumerate(functions):
            annot_data.append(
                [
                    f"annotations-{a_idx}-functions-{f_idx}-function",
                    function.get("category", ""),
                ]
            )

            evidences = function.get("evidence")

            if evidences is None:
                break

            for e_idx, evidence in enumerate(evidences):
                annot_data.append(
                    [
                        f"annotations-{a_idx}-functions-{f_idx}-evidence-{e_idx}-method",
                        evidence,
                    ]
                )
    return annot_data


if __name__ == "__main__":
    main()
