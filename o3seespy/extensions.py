from inspect import signature
from collections import OrderedDict


def to_commands(op_base_type, parameters):
    para = []
    for i, e in enumerate(parameters):
        if isinstance(e, str):
            e = "'%s'" % e
        elif isinstance(e, float):
            e = '%.6g' % e
            if '.' not in e and 'e' not in e:
                e += '.0'
        para.append(str(e))
        if i > 40:  # avoid verbose print output
            break
    p_str = ', '.join(para)
    return 'opy.%s(%s)' % (op_base_type, p_str)


def to_py_file(osi, ofile='ofile.py'):
    ofile = open(ofile, 'w')
    pstr = 'import openseespy.opensees as opy\n' + '\n'.join(osi.commands)
    pstr += '\nopy.analyze(1, 0.1)\n'
    ofile.write(pstr)
    ofile.close()


def get_o3_kwargs_from_obj(obj, o3_obj, custom=None, overrides=None):
    if custom is None:
        custom = {}
    if overrides is None:
        overrides = {}
    sig = signature(o3_obj)
    kwargs = OrderedDict()
    args = []
    for p in sig.parameters.values():
        if p.name in custom:
            pname = custom[p.name]
        else:
            pname = p.name
        if pname == 'osi':
            continue
        if pname in overrides:
            val = overrides[pname]
        else:
            val = getattr(obj, pname)
        if p.default == p.empty:
            args.append(val)
        else:
            if val is not None:
                kwargs[p.name] = val
    return args, kwargs
