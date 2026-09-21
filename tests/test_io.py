"""Regression for load_colo1 inert-guide classification (matched-singles)."""

import io as _io

import pandas as pd

from b06.audit import matched_guide_pairs
from b06.io import load_colo1


def _synth_colo1(tmp_path):
    """Mini COLO1-format table.

    - combo guides F1a, F1b, F2a, F2b (two functional genes)
    - inert partners I1, I2 (appear only in single-note rows next to F1*/F2*)
    - a NegativeControls pair of NTC guides
    """
    rows = [
        # doubles (AnchorCombinations): both guides functional
        ("d1", "AnchorCombinations", "F1a", "GENE1", "F2a", "GENE2", 1000, 500, 300, 100, 60, 40),
        ("d2", "AnchorCombinations", "F1b", "GENE1", "F2b", "GENE2", 1000, 500, 300, 100, 60, 40),
        # singles: functional guide + inert partner
        ("s1", "AnchorSingletons", "F1a", "GENE1", "I1", "", 1000, 500, 1000, 500, 250, 125),
        ("s2", "AnchorSingletons", "F2a", "GENE2", "I2", "", 1000, 500, 1000, 500, 250, 125),
        # negative control: both NTC
        ("n1", "NegativeControls", "NTC_A", "", "NTC_B", "", 1000, 500, 1000, 500, 1000, 500),
    ]
    cols = ["ID", "Note", "sgRNA1_ID", "Gene1", "sgRNA2_ID", "Gene2"] + [
        "SIDM00136_CPID2437", "SIDM00136_CPID2440", "SIDM00136_CPID2443",
        "SIDM00136_CPID1020", "SIDM00136_CPID1023", "SIDM00136_CPID1026",
    ]
    df = pd.DataFrame(rows, columns=cols)
    fn = tmp_path / "colo1_mini.txt"
    df.to_csv(fn, sep="\t", index=False)
    return fn


def test_load_colo1_classifies_inert_partners(tmp_path):
    df = load_colo1(str(_synth_colo1(tmp_path)))
    # combo guides stay functional on the a-side
    assert set(df.loc[df.guide_a.isin(["F1a", "F1b", "F2a", "F2b"]), "class_a"]) == {"gene"}
    # inert partners from single rows become non-targeting on the b-side
    assert set(df.loc[df.guide_b.isin(["I1", "I2"]), "class_b"]) == {"non-targeting"}
    # NTC guides are non-targeting on both sides
    assert set(df.loc[df.guide_a.isin(["NTC_A", "NTC_B"]), "class_a"]) == {"non-targeting"}


def test_matched_guide_pairs_uses_inert_singles(tmp_path):
    df = load_colo1(str(_synth_colo1(tmp_path)))
    matched = matched_guide_pairs(df, 0.0)
    # matched guide pairs exist: both sides must have a matched single.
    # F1a and F2a have I1/I2 singles and co-occur in double d1, so (F1a, F2a)
    # is a matched pair; F1b/F2b only appear in d2 and have no single context.
    assert ("F1a", "F2a") in set(map(tuple, matched[["guide_a", "guide_b"]].to_numpy()))
    assert not set(matched["guide_a"]).intersection({"F1b", "F2b"})