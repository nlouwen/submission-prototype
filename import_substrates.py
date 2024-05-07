#!/usr/bin/env python3
# populate substrates table

from submission import create_app
from submission.extensions import db
from submission.models import Substrate

PROTEINOGENIC = [
    "ala",
    "arg",
    "asn",
    "asp",
    "cys",
    "gln",
    "glu",
    "gly",
    "his",
    "ile",
    "leu",
    "lys",
    "met",
    "phe",
    "pro",
    "ser",
    "thr",
    "trp",
    "tyr",
    "val",
]


def get_substrate_data():
    """static collection of substrates, based on https://github.com/antismash/antismash/blob/master/antismash/modules/nrps_pks/data/aaSMILES.txt"""
    return [
        ["ala", "NC(C)C(=O)O", "alanine"],
        ["arg", "NC(CCCNC(N)=N)C(=O)O", "arginine"],
        ["asn", "NC(CC(=O)N)C(=O)O", "asparagine"],
        ["asp", "NC(CC(=O)O)C(=O)O", "aspartic acid"],
        ["cys", "NC(CS)C(=O)O", "cysteine"],
        ["gln", "NC(CCC(=O)N)C(=O)O", "glutamine"],
        ["glu", "NC(CCC(=O)O)C(=O)O", "glutamic acid"],
        ["gly", "NCC(=O)O", "glycine"],
        ["his", "NC(CC1=CNC=N1)C(=O)O", "histidine"],
        ["ile", "NC(C(C)CC)C(=O)O", "isoleucine"],
        ["leu", "NC(CC(C)C)C(=O)O", "leucine"],
        ["lys", "NC(CCCCN)C(=O)O", "lysine"],
        ["met", "NC(CCSC)C(=O)O", "methionine"],
        ["phe", "NC(Cc1ccccc1)C(=O)O", "phenylalanine"],
        ["pro", "N1C(CCC1)C(=O)O", "proline"],
        ["ser", "NC(CO)C(=O)O", "serine"],
        ["thr", "NC(C(O)C)C(=O)O", "threonine"],
        ["trp", "NC(CC1=CNc2c1cccc2)C(=O)O", "tryptophan"],
        ["tyr", "NC(Cc1ccc(O)cc1)C(=O)O", "tyrosine"],
        ["val", "NC(C(C)C)C(=O)O", "valine"],
        ["3Me-Glu", "NC(C(C)CC(=O))C(=O)O", "3-methyl-glutamate"],
        ["bAla", "NCCC(=O)O", "beta-alanine"],
        ["aIle", "NC(C(C)CC)C(=O)O", "allo-isoleucine"],
        ["aThr", "CC(O)C(N)C(=O)O", "allo-threonine"],
        ["diOH-Bz", "Oc1c(O)cccc1C(=O)O", "2,3-dihydroxy-benzoic acid"],
        ["bLys", "NCCCC(N)CC(=O)O", "beta-lysine"],
        ["bOH-Tyr", "NC(C(O)c1ccc(O)cc1)C(=O)O", "beta-hydroxy-tyrosine"],
        [
            "Cl2-Hpg",
            "NC(c1cc(Cl)c(O)c(Cl)c1)C(=O)O",
            "3,5-dichloro-4-hydroxy-phenylglycine",
        ],
        ["D-Ala", "NC(C)C(=O)O", "D-alanine"],
        ["D-Hiv", "OC(C(C)C)C(=O)O", "D-2-hydroxyisovalerate"],
        ["D-Hmp", "CCC(C)C(O)C(=O)O", "2-hydroxy-3-methyl-pentanoic acid"],
        ["dhpg", "NC(c1cc(O)cc(O)c1)C(=O)O", "3,5-dihydroxy-phenylglycine"],
        ["Hpr", "N1C(CCCC1)C(=O)O", "pipecolic acid"],
        ["IVal", "NC(CC)(C)C(=O)O", "isovaline"],
        ["OH-Orn", "NC(CCCNO)C(=O)O", "N5-hydroxyornithine"],
        ["Fo-OH-Orn", "NC(CCCN(O)C=O)C(=O)O", "N5-formyl-N5-hydroxyornithine"],
        ["Ac-OH-Orn", "NC(CCCN(O)C(=O)C)C(=O)O", "L-δ-N-acetyl-δ-N-hydroxyornithine"],
        [
            "C10:0-NH2(2)-Ep(9)-oxo(8)",
            "NC(CCCCCC(=O)C1OC1)C(=O)O",
            "2-amino-8-oxo-9,10-decanoate",
        ],
        ["Valol", "NC(C(C)C)CO", "valinol"],
        ["Pgl", "NC(c1ccccc1)C(=O)O", "phenylglycine"],
        ["pPro", "N1CC(CCC)CC1C(=O)O", "4-propyl-proline"],
        ["aad", "NC(CCCC(=O)O)C(=O)O", "2-amino-adipic acid"],
        ["abu", "NC(C(C))C(=O)O", "2-amino-butyric acid"],
        ["bmt", "NC(C(O)C(C)CC=CC)C(=O)O", "4-butenyl-4-methyl threonine"],
        ["cap", "NC(C1CCN=C(N1)N)C(=O)O", "capreomycidine"],
        ["dab", "NC(CCN)C(=O)O", "2,4-diaminobutyric acid"],
        ["dht", "NC(C(=O)C)C(=O)O", "dehydro-threonine"],
        ["hiv", "OC(C(C)C)C(=O)O", "2-hydroxyisovalerate"],
        ["hpg", "NC(c1ccc(O)cc1)C(=O)O", "4-hydroxy-phenylglycine"],
        ["hyv", "NC(C(CO)C)C(=O)O", "4-hydroxy-L-valine"],
        ["hyv-d", "OC(C(C)C)C(=O)O", "2-hydroxy-valeric acid"],
        ["orn", "NC(CCCN)C(=O)O", "ornithine"],
        ["sal", "Oc1ccccc1C(=O)O", "salicylic acid"],
        ["tcl", "NC(CC(C)C(Cl)(Cl)(Cl))C(=O)O", "(4S)-5,5,5-trichloro-leucine"],
        ["LDAP", "NC(CCCC(N)C(=O)O)C(=O)O", "diaminopimelic acid"],
        ["meval", "NC(C(C)(C)C)C(=O)O", "Methyl-Valine"],
        ["alaninol", "NC(C)CO"],
        ["N-(1,1-dimethyl-1-allyl)Trp", "NC(CC1=CN(C(C)(C)C=C)c2c1cccc2)C(=O)O"],
        ["d-lyserg", "CN1CC(C=C2C1CC3=CNC4=CC=CC2=C34)C(=O)O", "D-lysergic acid"],
        ["ser-thr", "NC(C([*])O)C(=O)O", "Serine or Threonine"],
        ["mephe", "NC(C(C)c1ccccc1)C(=O)O", "Cmethyl-phenylalanine"],
        ["hasn", "NC(C(O)C(=O)N)C(=O)O", "hydroxyasparagine"],
        ["s-nmethoxy-trp", "NC(CC1=CN(OC)c2c1cccc2)C(=O)O"],
        ["2S-Hic", "OC(C(O)CC(C)C)C(=O)O", "alpha-hydroxy-isocaproic-acid"],
        ["MeHOval", "O=C(C(C)CC)C(=O)O", "3-Methyl-2-oxovaleric acid"],
        ["2-oxo-isovaleric-acid", "O=C(C(C)C)C(=O)O"],
        ["aoda", "NC(CCCCCC(=O)CC)C(=O)O", "S-2-amino-8-oxodecanoic acid"],
        # ["mal", "CC(=O)", "malonyl-CoA"],
        # ["ohmal", "CC(O)", "malonyl-CoA"],
        # ["ccmal", "C=C", "malonyl-CoA, double-bonded"],
        # ["redmal", "CC", "malonyl-CoA, reduced"],
        # ["me-mal", "C(C)C(=O)", "malonyl-CoA, methylated variant"],
        # ["me-ohmal", "C(C)C(O)", "malonyl-CoA, methylated variant"],
        # ["me-ccmal", "C(C)=C", "malonyl-CoA, methylated variant"],
        # ["me-redmal", "C(C)C", "malonyl-CoA, methylated variant"],
        # ["mxmal", "C(OC)C(=O)", "malonyl-CoA, methoxy variant"],
        # ["ohmxmal", "C(OC)C(O)", "malonyl-CoA, methoxy variant"],
        # ["ccmxmal", "C(OC)=C", "malonyl-CoA, methoxy variant"],
        # ["redmxmal", "C(OC)C", "malonyl-CoA, methoxy variant"],
        # ["emal", "C(CC)C(=O)", "malonyl-CoA, ethyl variant"],
        # ["ohemal", "C(CC)C(O)", "malonyl-CoA, ethyl variant"],
        # ["ccemal", "C(CC)=C", "malonyl-CoA, ethyl variant"],
        # ["redemal", "C(CC)C", "malonyl-CoA, ethyl variant"],
    ]


def create_subtrate_entry(line):
    if len(line) == 2:
        name, structure = line
        fullname = None
    else:
        name, structure, fullname = line

    proteinogenic = name in PROTEINOGENIC

    substrate = Substrate(
        name=name, structure=structure, fullname=fullname, proteinogenic=proteinogenic
    )
    db.session.add(substrate)


def main():
    app = create_app()
    with app.app_context():
        data = get_substrate_data()
        for line in data:
            create_subtrate_entry(line)
        db.session.commit()


if __name__ == "__main__":
    main()
