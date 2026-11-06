"""Single shared total ordering for DNA->RNA causal statuses (review B3b-13).

Three modules previously kept private copies of this ladder, and evidence_chains
tied PVACSPLICE_SUPPORTED with TARGETED_SPANNING_ONLY while consensus/pipeline
ranked them 4>3, so the "best causal link" could differ between tables for the
same evidence.  One shared mapping is the point of truth; ties are not allowed.
"""

CAUSAL_STATUS_PRIORITY: dict[str, int] = {
    "UNASSESSED": 0,
    "DNA_PREDICTION_ONLY": 1,
    "DNA_RNA_CIS_SUPPORTED": 2,
    "TARGETED_SPANNING_ONLY": 3,
    "PVACSPLICE_SUPPORTED": 4,
    "TARGETED_REQUANT_SUPPORTED": 5,
}
