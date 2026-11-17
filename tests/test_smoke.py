from breadbasket_waves import params


def test_targets_have_sources():
    for t in params.TARGETS:
        assert t.source and t.where


def test_lmf_example_counts_are_consistent():
    # 6/15 divided by 5/20 is 1.6, the worked example in the paper's Methods.
    hits_multi = params.target("lmf_example_hits_multi").value
    n_multi = params.target("lmf_example_n_multi").value
    hits_zero = params.target("lmf_example_hits_zero").value
    n_zero = params.target("lmf_example_n_zero").value
    assert round((hits_multi / n_multi) / (hits_zero / n_zero), 2) == params.target("lmf_example").value
