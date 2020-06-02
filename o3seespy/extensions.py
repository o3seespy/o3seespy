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

def _dep_get_fn_name_and_args(line):
    import re
    line = line.replace('.', 'stoppoint')
    fn_match = re.match(r"(?P<function>\w+)\s?\((?P<arg>(?P<args>\w+(,\s?)?)+)\)", line)
    fn_dict = fn_match.groupdict()
    args = [arg.strip() for arg in fn_dict['arg'].split(',')]
    args = [arg.replace('stoppoint', '.') for arg in args]
    return fn_dict['function'], args


def _get_fn_name_and_args(line):
    import re
    fn_name = line.split('(')[0]
    args = re.search(r'\((.*?)\)', line).group(1)
    args = args.replace(' ', '')
    args = args.split(',')
    # line = line.replace('.', 'stoppoint')
    # fn_match = re.match(r"(?P<function>\w+)\s?\((?P<arg>(?P<args>\w+(,\s?)?)+)\)", line)
    # fn_dict = fn_match.groupdict()
    # args = [arg.strip() for arg in fn_dict['arg'].split(',')]
    # args = [arg.replace('stoppoint', '.') for arg in args]
    return fn_name, args


def check_if_opy_lines_consistent(line1, line2, line3=None):
    if line3 is None:
        if line1 == line2:
            return True
        print(line1)
        fn1, args1 = _get_fn_name_and_args(line1)
        print(line2)
        fn2, args2 = _get_fn_name_and_args(line2)
        if fn1 == fn2 and len(args1) == len(args2):
            for i in range(len(args1)):
                if args1[i] == args2[i]:
                    continue
                elif isinstance(args1[i], str) or isinstance(args2, str):  # TODO: currently all str
                    # consistent = 0  # since strings must be identical
                    return False
            return True
        else:
            return False
    else:
        if line1 == line2 and line1 == line3:
            return True
        fn1, args1 = _get_fn_name_and_args(line1)
        fn2, args2 = _get_fn_name_and_args(line2)
        fn3, args3 = _get_fn_name_and_args(line3)
        if fn1 == fn2 and fn1 == fn3 and len(args1) == len(args2) and len(args1) == len(args3):
            for i in range(len(args1)):
                if args1[i] == args2[i] and args1[i] == args3[i]:
                    continue
                elif isinstance(args1[i], str) or isinstance(args2[i], str) or isinstance(args3[i], str):
                    # consistent = 0  # since strings must be identical
                    return False
                else:
                    val1 = float(args1[i])
                    val2 = float(args2[i])
                    val3 = float(args3[i])
                    if val2 - val1 == val3 - val2:
                        continue
                    else:
                        return False
        else:
            return False


def compress_opy_lines(commands):
    slines = 10  # search lines
    latest_rep = 0
    new_commands = []
    for i, com in enumerate(commands[:-slines]):
        if i <= latest_rep:
            continue
        new_commands.append(com)
        print(com)
        dup_detected = 0
        for j in range(1, slines):
            if dup_detected:
                break
            if check_if_opy_lines_consistent(com, commands[i + j]):
                # check all lines in between are repeated
                consistent = 1
                for k in range(1, j):
                    if not check_if_opy_lines_consistent(commands[i + k], commands[i + j + k]):
                        consistent = 0
                if consistent:
                    nr = 1
                    new_rep = i + nr * j
                    while consistent:
                        nr += 1
                        new_rep = i + nr * j
                        for k in range(j):
                            if new_rep + k >= len(commands) - 1:
                                consistent = 0
                            elif commands[i + k] != commands[new_rep + k]:
                                consistent = 0
                    if nr > 3:
                        dup_detected = 1
                        print('FOUND major repetition: ', nr)
                        latest_rep = new_rep
                        new_commands.append(f'for i in range({nr}):')
                        for k in range(j):
                            new_commands.append('    ' + commands[i + k])

    return new_commands


def to_cpy_file(osi, ofile='ofile.py', w_analyze=False):
    commands = compress_opy_lines(osi.commands)
    ofile = open(ofile, 'w')
    pstr = 'import openseespy.opensees as opy\n' + '\n'.join(commands)
    if w_analyze:
        pstr += '\nopy.analyze(1, 0.1)\n'
    ofile.write(pstr)
    ofile.close()


def to_py_file(osi, ofile='ofile.py', w_analyze=False):
    ofile = open(ofile, 'w')
    pstr = 'import openseespy.opensees as opy\n' + '\n'.join(osi.commands)
    if w_analyze:
        pstr += '\nopy.analyze(1, 0.1)\n'
    ofile.write(pstr)
    ofile.close()


def to_tcl_file(osi, ofile='ofile.tcl', w_analyze=False):
    ofile = open(ofile, 'w')
    pstr = '\n'.join(osi.commands)
    if w_analyze:
        pstr += '\nopy.analyze(1, 0.1)\n'
    tcl_str = py2tcl(pstr)
    ofile.write(tcl_str)
    ofile.close()


def get_o3_kwargs_from_obj(obj, o3_obj, custom=None, overrides=None):
    if custom is None:
        custom = {}
    if overrides is None:
        overrides = {}
    sig = signature(o3_obj)
    kwargs = OrderedDict()
    args = []
    sig_vals = sig.parameters.values()
    for p in sig_vals:
        if p.name in custom:
            pname = custom[p.name]
        else:
            pname = p.name
        if pname == 'osi':
            continue
        if pname in overrides:
            val = overrides[pname]
        else:
            try:
                val = getattr(obj, pname)
            except AttributeError as e:
                if p.default == p.empty:
                    raise AttributeError(e)
                else:
                    val = p.default
        if p.default == p.empty:
            args.append(val)
        else:
            if val is not None:
                kwargs[p.name] = val
    return args, kwargs


def has_o3_model_changed(cur_type, prev_type, cur_args, prev_args, cur_kwargs, prev_kwargs):
    import numpy as np
    changed = 0
    if cur_type != prev_type or len(cur_args) != len(prev_args) or len(cur_kwargs) != len(prev_kwargs):
        changed = 1
    else:
        for j, arg in enumerate(cur_args):
            if hasattr(arg, '__len__'):
                if len(arg) != len(prev_args[j]):
                    changed = 1
                    break
                for k, subarg in enumerate(arg):
                    if not np.isclose(subarg, prev_args[j][k]):
                        changed = 1
            elif not np.isclose(arg, prev_args[j]):
                changed = 1
                break
        for pm in cur_kwargs:
            if pm not in prev_kwargs:
                changed = 1
                break
            if hasattr(cur_kwargs[pm], '__len__'):
                if len(cur_kwargs[pm]) != len(prev_kwargs[pm]):
                    changed = 1
                    break
                for k, subarg in enumerate(cur_kwargs[pm]):
                    if not np.isclose(subarg, cur_kwargs[pm][k]):
                        changed = 1
                        break
            elif not np.isclose(cur_kwargs[pm], prev_kwargs[pm]):
                changed = 1
                break
    return changed


def py2tcl(pystr):
    """
    Converts openseespy script to tcl

    Returns
    -------

    """
    # new = '\n'.join(pystr.split()[1:])  # removes the import line
    new = pystr.replace('(', ' ')
    new = new.replace(')', ' ')
    new = new.replace('opy.', '')
    new = new.replace(',', '')
    new = new.replace("'", '')
    new = new.replace('"', '')
    return new
