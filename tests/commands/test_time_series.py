import o3seespy as o3  # for testing only


def test_constant():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.time_series.Constant(osi, factor=1.0)


def test_linear():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.time_series.Linear(osi, factor=1.0)


def test_trig():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.time_series.Trig(osi, t_start=1.0, t_end=1.0, period=1.0, factor=1.0, shift=0.0, zero_shift=0.0)


def test_triangle():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.time_series.Triangle(osi, t_start=1.0, t_end=1.0, period=1.0, factor=1.0, shift=0.0, zero_shift=0.0)


def test_rectangular():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.time_series.Rectangular(osi, t_start=1.0, t_end=1.0, factor=1.0)


def test_pulse():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.time_series.Pulse(osi, t_start=1.0, t_end=1.0, period=1.0, width=0.5, shift=0.0, factor=1.0, zero_shift=0.0)

#
# def test_path():
#     osi = o3.OpenSeesInstance(ndm=2)
#     values = [1.0, 1.0]
#     time = [1.0, 1.0]
#     o3.time_series.Path(osi, dt=0.0, values=values, time=time, filepath='', file_time='', factor=1.0, start_time=0.0, use_last="string", prepend_zero="string")

