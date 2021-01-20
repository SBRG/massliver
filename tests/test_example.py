from massliver import module_a


def test_a():
    assert module_a.a() == "a"


def test_b():
    if module_a is not None:
        assert module_a.a() == "a"
    else:
        assert module_a is None
