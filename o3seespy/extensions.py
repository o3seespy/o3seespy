

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
