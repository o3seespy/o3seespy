import numpy as np
import matplotlib.pyplot as plt

import openseespy.opensees as opy
import o3seespy as opw

from tests.conftest import TEST_DATA_DIR

# Load a ground motion
dt = 0.01
rec = np.loadtxt(TEST_DATA_DIR + 'test_motion_dt0p01.txt')

# Define inelastic SDOF
period = 1.0
xi = 0.05
mass = 1.0
f_yield = 1.5  # Reduce this to make it nonlinear
r_post = 0.0

# Initialise OpenSees instance
osi = opw.OpenseesInstance(dimensions=2, state=0)

# Establish nodes
bot_node = opw.node.Node(osi, 0, 0)
top_node = opw.node.Node(osi, 0, 0)

# Fix bottom node
opw.Fix(osi, top_node, opw.cc.FREE, opw.cc.FIXED, opw.cc.FIXED)
opw.Fix(osi, bot_node, opw.cc.FIXED, opw.cc.FIXED, opw.cc.FIXED)
# Set out-of-plane DOFs to be slaved
opw.EqualDOF(osi, top_node, bot_node, [opw.cc.Y, opw.cc.ROTZ])

# nodal mass (weight / g):
opw.Mass(osi, top_node, mass, 0., 0.)

# Define material
k_spring = 4 * np.pi ** 2 * mass / period ** 2
bilinear_mat = opw.uniaxial_material.Steel01(osi, fy=f_yield, e0=k_spring, b=r_post)

# Assign zero length element, # Note: pass actual node and material objects into element
opw.element.ZeroLength(osi, bot_node, top_node, mat_x=bilinear_mat, r_flag=1)

# Define the dynamic analysis
load_tag_dynamic = 1
pattern_tag_dynamic = 1

values = list(-1 * rec)  # should be negative
opy.timeSeries('Path', load_tag_dynamic, '-dt', dt, '-values', *values)
opy.pattern('UniformExcitation', pattern_tag_dynamic, opw.cc.X, '-accel', load_tag_dynamic)

# set damping based on first eigen mode
angular_freq = opy.eigen('-fullGenLapack', 1) ** 0.5
beta_k = 2 * xi / angular_freq
opw.rayleigh.Rayleigh(osi, alpha_m=0.0, beta_k=beta_k, beta_k_init=0.0, beta_k_comm=0.0)

# Run the dynamic analysis
opw.wipe_analysis(osi)

# Run the dynamic analysis
opw.algorithm.Newton(osi)
opw.system.SparseGeneral(osi)
opw.numberer.RCM(osi)
opw.constraints.Transformation(osi)
opw.integrator.Newmark(osi, gamma=0.5, beta=0.25)
opw.analysis.Transient(osi)

opw.test_check.EnergyIncr(osi, tol=1.0e-10, max_iter=10)
analysis_time = (len(values) - 1) * dt
analysis_dt = 0.001
outputs = {
    "time": [],
    "rel_disp": [],
    "rel_accel": [],
    "rel_vel": [],
    "force": []
}

# access underlying openseespy commands to control analysis
while opy.getTime() < analysis_time:

    opy.analyze(1, analysis_dt)
    curr_time = opy.getTime()
    outputs["time"].append(curr_time)
    outputs["rel_disp"].append(opy.nodeDisp(top_node.tag, opw.cc.X))
    outputs["rel_vel"].append(opy.nodeVel(top_node.tag, opw.cc.X))
    outputs["rel_accel"].append(opy.nodeAccel(top_node.tag, opw.cc.X))
    opy.reactions()
    outputs["force"].append(-opy.nodeReaction(bot_node.tag, opw.cc.X))  # Negative since diff node
opy.wipe()
for item in outputs:
    outputs[item] = np.array(outputs[item])


plt.plot(outputs['time'], outputs['rel_disp'], label='o3seespy')
periods = np.array([period])

# Compare closed form elastic solution
from eqsig import sdof
resp_u, resp_v, resp_a = sdof.response_series(motion=rec, dt=dt, periods=periods, xi=xi)
plt.plot(np.arange(len(rec)) * dt, resp_u[0], ls='--', label='Elastic')
plt.legend()
plt.savefig('readme_example.png')
plt.show()


