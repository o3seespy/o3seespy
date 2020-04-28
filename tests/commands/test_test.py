import o3seespy as o3  # for testing only
import pytest


def test_norm_unbalance():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.NormUnbalance(osi, tol=1.0, max_iter=1, p_flag=0, n_type=2, max_incr=None)


def test_norm_disp_incr():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.NormDispIncr(osi, tol=1.0, max_iter=1, p_flag=0, n_type=2)


def test_energy_incr():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.EnergyIncr(osi, tol=1.0, max_iter=1, p_flag=0, n_type=2)


def test_relative_norm_unbalance():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.RelativeNormUnbalance(osi, tol=1.0, max_iter=1, p_flag=0, n_type=2)


def test_relative_norm_disp_incr():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.RelativeNormDispIncr(osi, tol=1.0, max_iter=1, p_flag=0, n_type=2)


def test_relative_total_norm_disp_incr():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.RelativeTotalNormDispIncr(osi, tol=1.0, max_iter=1, p_flag=0, n_type=2)


def test_relative_energy_incr():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.RelativeEnergyIncr(osi, tol=1.0, max_iter=1, p_flag=0, n_type=2)


def test_fixed_num_iter():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.FixedNumIter(osi, max_iter=1, p_flag=0, n_type=2)


def test_norm_disp_and_unbalance():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.NormDispAndUnbalance(osi, tol_incr=1.0, tol_r=1, max_iter=1, p_flag=0, n_type=2, maxincr=-1)


def test_norm_disp_or_unbalance():
    osi = o3.OpenSeesInstance(ndm=2)
    o3.test.NormDispOrUnbalance(osi, tol_incr=1.0, tol_r=1, max_iter=1, p_flag=0, n_type=2, maxincr=-1)

