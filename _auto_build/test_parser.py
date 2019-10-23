from auto_build import run


def test_check_default_is_an_expression():
    assert run.check_if_default_is_expression('a2*w')
    assert run.check_if_default_is_expression('a2/w')
    assert run.check_if_default_is_expression('a2+w')
    assert run.check_if_default_is_expression('a2-w')
    assert run.check_if_default_is_expression('a2^w')


if __name__ == '__main__':
    test_check_default_is_an_expression()