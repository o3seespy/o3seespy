import o3seespy as o3  # for testing only


def test_plain():
    osi = o3.OpenSeesInstance(ndm=2)
    ts = o3.time_series.Linear(osi, factor=1.0)
    o3.pattern.Plain(osi, ts=ts, fact=1.0)


def test_uniform_excitation():
    osi = o3.OpenSeesInstance(ndm=2)
    ts = o3.time_series.Linear(osi, factor=1.0)
    o3.pattern.UniformExcitation(osi, dir=1, accel_series=ts, vel0=1.0, fact=1.0)


def test_multiple_support():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.pattern.MultipleSupport(osi)

