import os
import re
from collections import OrderedDict
import pandas as pd

ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) + "/../"
w4 = '    '
w8 = '        '
pname_pat = '\``([A-Za-z0-9_\./\\\'-]*)\``'
dtype_pat = '\|([A-Za-z0-9_\./\\-\ ]*)\|'
optype_pat = "\'([A-Za-z0-9_\./\\-]*)\'"

special_words = {
    'lambda': 'lamb',
}

def clean_param_names(params):
    pms = OrderedDict()
    for pm in params:
        new_pm = pm
        dtype_is_obj = False
        if len(pm) == 1 and pm.istitle():
            new_pm = 'big_' + pm.lower()
        else:
            new_pm = convert_camel_to_snake(pm)
        if new_pm in special_words:
            new_pm = special_words[new_pm]

        if len(new_pm) > 4 and new_pm[-4:] == '_tag':
            new_pm = new_pm[:-4]
            dtype_is_obj = True
        pms[pm] = params[pm]
        if dtype_is_obj:
            pms[pm].dtype = 'obj'
        pms[pm].o3_name = new_pm
    # if 'tag' in pms[0]:
    #     pms = pms[1:]
    return pms


def convert_name_to_class_name(name):
    name = name[0].capitalize() + name[1:]
    name = name.replace('_', '')
    return name


def convert_camel_to_snake(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    return re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()


def constructor(base_type, op_type, defaults, op_kwargs):
    kw_pms = []  # TODO: This only works if object should be split into many objects based on kwargs, but
    # some kwargs are just to input a single value.
    for kw in op_kwargs:
        for pm in op_kwargs[kw]:
            kw_pms.append(pm)
    if not len(op_kwargs):
        op_kwargs[''] = []
    para = ['']
    for kw in op_kwargs:
        name_from_kw = kw.replace('-', '')
        name_from_kw = name_from_kw.replace("'", '')
        op_class_name = convert_name_to_class_name(op_type + name_from_kw)
        base_class_name = base_type[0].capitalize() + base_type[1:]
        para.append(f'class {op_class_name}({base_class_name}Base):')
        para.append('')

        pms = clean_param_names(defaults)
        # Determine what is in this class
        cl_pms = []  # TODO: need to add op_kwarg str
        for pm in pms:
            if pms[pm].org_name not in kw_pms or pm in op_kwargs[kw]:
                cl_pms.append(pm)

        # Build definition string
        # check no no default pms are after pms that have defaults
        contains_no_default = False
        for i in range(len(cl_pms)):
            pm = cl_pms[-1-i]
            if pms[pm].default_is_expression and contains_no_default:
                pms[pm].forced_not_optional = True  # TODO: should raise a flag
            if pms[pm].default_value is None:
                contains_no_default = True

        pjoins = []
        for pm in cl_pms:
            o3_name = pms[pm].o3_name
            default = pms[pm].default_value
            if pms[pm].is_flag:
                pjoins.append(f'{pm}=False')
            elif default is not None:
                if pms[pm].forced_not_optional:
                    pjoins.append(f'{o3_name}')
                elif pms[pm].default_is_expression:
                    pjoins.append(f'{o3_name}=None')
                else:
                    pjoins.append(f'{o3_name}={default}')
            else:
                if pms[pm].marker:
                    pjoins.append(f'{o3_name}=None')
                else:
                    pjoins.append(f'{o3_name}')

        pjoined = ', '.join(pjoins)
        para.append(f'    def __init__(self, osi, {pjoined}):')

        # Create init function saving logic
        for i, pm in enumerate(cl_pms):
            o3_name = pms[pm].o3_name
            dtype = pms[pm].dtype
            if dtype == 'float':
                para.append(w8 + f'self.{o3_name} = float({o3_name})')
            elif dtype == 'int':
                para.append(w8 + f'self.{o3_name} = int({o3_name})')
            elif dtype == 'obj':
                para.append(w8 + f'self.{o3_name} = {o3_name}.tag')
            else:
                para.append(w8 + f'self.{o3_name} = {o3_name}')
        para.append(w8 + 'osi.n_mats += 1')
        para.append(w8 + 'self._tag = osi.mats')
        pjoins = []
        need_special_logic = False
        applied_op_warg = False
        for pm in cl_pms:
            o3_name = pms[pm].o3_name
            if pms[pm].marker:
                continue
            if pms[pm].is_flag:
                continue
            if not applied_op_warg and pm in op_kwargs[kw]:
                applied_op_warg = True
                pjoins.append(f"'{kw}'")
            if pms[pm].default_is_expression:
                need_special_logic = True
                break
            if pms[pm].dtype == 'obj':
                pjoins.append(f'self.{o3_name}.tag')
            elif pms[pm].packed:
                pjoins.append('*self.' + o3_name)
            else:
                pjoins.append('self.' + o3_name)
        para.append(w8 + 'self._parameters = [self.op_type, self._tag, %s]' % (', '.join(pjoins)))
        for pm in cl_pms:
            o3_name = pms[pm].o3_name
            if pms[pm].marker:
                para.append(w8 + f"if getattr(self, '{o3_name}') is not None:")
                if pms[pm].packed:
                    para.append(w8 + w4 + f"self._parameters += ['-{pms[pm].marker}', *self.{o3_name}]")
                else:
                    para.append(w8 + w4 + f"self._parameters += ['-{pms[pm].marker}', self.{o3_name}]")
            if pms[pm].is_flag:
                para.append(w8 + f"if getattr(self, '{o3_name}') is not None:")
                para.append(w8 + w4 + f"self._parameters += ['-{o3_name}']")  # TODO: does this always work?
        if need_special_logic:
            sp_logic = False
            sp_pms = []
            for pm in pms:
                if pms[pm].default_is_expression:
                    sp_logic = True
                if sp_logic:
                    sp_pms.append(pm)
            sp_pm_strs = ["'%s'" % pms[pm].o3_name for pm in sp_pms]
            para.append(w8 + f"special_pms = [{', '.join(sp_pm_strs)}]")
            packets = [str(pms[pm].packed) for pm in sp_pms]
            para.append(w8 + f"packets = [{', '.join(packets)}]")
            para.append(w8 + 'for i, pm in enumerate(special_pms):')
            para.append(w8 + w4 + 'if getattr(self, pm) is not None:')
            para.append(w8 + w8 + 'if packets[i]:')
            para.append(w8 + w8 + w4 + 'self._parameters += [*getattr(self, pm)]')
            para.append(w8 + w8 + 'else:')
            para.append(w8 + w8 + w4 + 'self._parameters += [getattr(self, pm)]')
        para.append(w8 + 'self.to_process(osi)')
        para.append('')

        low_op_name = convert_camel_to_snake(op_class_name)
        low_base_name = convert_camel_to_snake(base_class_name)

        # Build test
        tpara = ['', f'def test_{low_op_name}():', w4 + 'osi = o3.OpenseesInstance(dimensions=2)']
        pjoins = []
        for i, pm in enumerate(cl_pms):
            o3_name = pms[pm].o3_name
            default = pms[pm].default_value
            dtype = pms[pm].dtype
            if pms[pm].default_is_expression:
                pjoins.append(f'{o3_name}=None')
            elif default is not None:
                pjoins.append(f'{o3_name}={default}')
            elif dtype == 'float':
                pjoins.append(f'{o3_name}=1.0')
            elif dtype == 'obj':
                pjoins.append(f'{o3_name}=obj')
            else:
                pjoins.append(f'{o3_name}=1')
        pjoint = ', '.join(pjoins)
        tpara.append(w4 + f'o3.{low_base_name}.{op_class_name}(osi, {pjoint})')
        tpara.append('')
        tpara.append('')
    return '\n'.join(para), '\n'.join(tpara)


class Param(object):
    def __init__(self, org_name, default_value, packed=None, dtype=None):
        self.org_name = org_name
        self.o3_name = None
        self.default_value = default_value
        self.packed = packed
        self.dtype = dtype
        self.default_is_expression = False
        self.p_description = ''
        self.marker = None
        self.is_flag = False
        self.forced_not_optional = False  # if default value precedes non default values


def check_if_default_is_expression(defo):
    # if any(re.findall('|'.join(['\*', '\/', '\+', '\-', '\^']), defo)):  # deal with -1, or 1e-1
    #     return True
    try:
        a = float(defo)
        return False
    except ValueError:
        if any(re.findall('|'.join(['\"', "\'"]), defo)):
            return False
        return True


def clean_fn_line(line):
    df = pd.read_csv('overrides.csv')
    defaults = OrderedDict()
    print(line)
    base_type = line.split('.. function:: ')[-1]
    base_type = base_type.split('(')[0]
    optype_res = re.search(optype_pat, line)
    optype = optype_res.group()[1:-1]
    df_op = df[(df['base_type'] == base_type) & (df['op_type'] == optype)]

    print(optype_res)
    inputs_str = line.split(')')[0]
    inputs = inputs_str.split(',')
    inputs = inputs[2:]  # remove class definition and tag
    op_kwargs = OrderedDict()
    markers = OrderedDict()
    flags = []
    cur_kwarg = None
    names_only = []
    for j, inpy in enumerate(inputs):
        name_only = inpy.replace(' ', '')
        name_only = name_only.replace('[', '')
        name_only = name_only.replace(']', '')
        if '*' in name_only:
            name_only = name_only.replace('*', '')
            # name_only = name_only.replace('Args', '')
        name_only = name_only.split('=')[0]
        names_only.append(name_only)
    for j, inpy in enumerate(inputs):
        inpy = inpy.replace(' ', '')
        inpy = inpy.replace('[', '')
        inpy = inpy.replace(']', '')
        if '=' in inpy:
            inp, defo = inpy.split('=')
        else:
            inp = inpy
            defo = None
        if '-' in inp:
            word = inp[2:-1]
            if len(df_op):
                df_p = df_op[df_op['param'] == word]
                if len(df_p):
                    if df_p['marker'].iloc[0] != '':
                        markers[names_only[j + 1]] = word
                        continue
            if len(inputs) > j + 1 and (word == names_only[j + 1] or word + 'Args' == names_only[j + 1] or word + 'Tag' == names_only[j + 1]):
                markers[names_only[j + 1]] = word
                continue
            elif len(inputs) > j + 1 and '-' in names_only[j + 1]:  # TODO: unsure if this is best way to identify flags
                flags.append(word)
                inp = word
            else:
                cur_kwarg = '-' + word
                op_kwargs[cur_kwarg] = []
                continue
        if inp[0] == '*':
            inp = inp[1:]
            packed = True
        else:
            packed = False
        if inp in defaults:
            i = 2
            new_inp = inp + '_%i' % i
            while new_inp in defaults:
                i += 1
                new_inp = inp + '_%i' % i
        else:
            new_inp = inp
        defaults[new_inp] = Param(org_name=inp, default_value=defo, packed=packed)
        if defo is not None and check_if_default_is_expression(defo):
            defaults[inp].default_is_expression = True
        if cur_kwarg is not None:
            op_kwargs[cur_kwarg].append(new_inp)
        for marker in markers:
            defaults[marker].marker = markers[marker]
        for flag in flags:
            defaults[flag].is_flag = True

    return base_type, optype, defaults, op_kwargs


def parse_mat_file(ffp):
    a = open(ffp, encoding="utf8")
    f = a.read()
    lines = f.splitlines()
    doc_str_pms = []
    dtypes = []
    defaults = None
    base_type = None
    optype = None
    op_kwargs = OrderedDict()
    descriptions = []
    for line in lines:
        char_only = line.replace(' ', '')
        char_only = char_only.replace('\t', '')
        if not len(char_only):
            continue
        first_char = char_only[0]
        if first_char == '*':
            continue
        res = re.search(pname_pat, line)
        if res:
            pname = res.group()[2:-2]
            print(pname)
            # if len(res.group()) > 4 and "'-" == res.group()[2:4]:
            #     continue  # op_kwarg
            ei = line.find('|')
            dtype_res = re.search(dtype_pat, line)
            if dtype_res is None:
                continue
            dtype = dtype_res.group()[1:-1]
            des = line[dtype_res.end():]
            des = des.replace('\t', ' ')
            if not len(des):
                continue
            while des[0] == ' ':
                des = des[1:]
                if not len(des):
                    break

            res = re.findall(pname_pat, line[:ei])
            for pm in res:
                if len(pm) > 4 and "'-" == pm[0:2]:
                    pm = pm[2:-1]
                doc_str_pms.append(pm)
                dtypes.append(dtype)
                descriptions.append(des)
        if '.. function:: ' in line:
            if base_type is None:
                base_type, optype, defaults, op_kwargs = clean_fn_line(line)
            else:
                base_type1, optype1, defaults1, op_kwargs1 = clean_fn_line(line)
                assert base_type == base_type1
                assert optype == optype1
                for inp in defaults1:
                    if inp not in defaults:
                        defaults[inp] = defaults1[inp]
                for inp in op_kwargs1:
                    if inp not in op_kwargs:
                        op_kwargs[inp] = op_kwargs1[inp]

    doc_str_pms = doc_str_pms[1:]  # remove mat tag
    dtypes = dtypes[1:]
    print('doc_str: ', doc_str_pms)
    print('fn_inps: ', list(defaults))
    print('op_kwargs: ', list(op_kwargs))
    assert len(doc_str_pms) == len(defaults) + len(op_kwargs), (len(doc_str_pms), (len(defaults), len(op_kwargs)))
    for i, pm in enumerate(doc_str_pms):
        if '-' + pm in op_kwargs:
            continue
        defaults[pm].dtype = dtypes[i]
        print(pm, i)
        defaults[pm].p_description = descriptions[i]

    pstr, tstr = constructor(base_type, optype, defaults, op_kwargs)
    print(pstr, tstr)
    return pstr, tstr
    # if line[3:5] == '``':
    #     para = line[5:]


def parse_all_uniaxial_mat():
    import user_paths as up
    uni_axial_mat_file = open(up.OPY_DOCS_PATH + 'uniaxialMaterial.rst')
    lines = uni_axial_mat_file.read().split('\n')
    collys = {}
    mtype = None
    for line in lines:
        if 'Steel & Reinforcing-Steel Materials' in line:
            mtype = 'steel'
            collys[mtype] = []
            continue
        if 'Concrete Materials' in line:
            mtype = 'concrete'
            collys[mtype] = []
            continue
        if 'Standard Uniaxial Materials' in line:
            mtype = 'standard'
            collys[mtype] = []
            continue
        if 'PyTzQz uniaxial materials' in line:
            mtype = 'pytz'
            collys[mtype] = []
            continue
        if 'Other Uniaxial Materials' in line:
            mtype = 'other'
            collys[mtype] = []
            continue
        if mtype is not None:
            line = line.replace(' ', '')
            line = line.replace('\t', '')
            if ':' in line or '-' in line or '#' in line  or line == '':
                continue
            collys[mtype].append(line)

    floc = ROOT_DIR + 'o3seespy/command/uniaxial_material/'
    for item in collys:
        para = ['from o3seespy.command.uniaxial_material.base_material import UniaxialMaterialBase', '', '']
        tpara = ['import o3seespy as o3  # for testing only', '', '']
        print(item, collys[item])
        for mat in collys[item]:
            if mat == 'steel4':  # consider as a single object, make all values after E0 equal to None, all -ps should be flags.
                continue
            if mat == 'Pinching4':
                continue
            if mat == 'KikuchiAikenHDR' or mat == 'KikuchiAikenLRB':
                continue
            if mat == 'PinchingLimitStateMaterial':
                continue
            open(up.OPY_DOCS_PATH + '%s.rst' % mat)
            ffp = up.OPY_DOCS_PATH + '%s.rst' % mat
            pstr, tstr = parse_mat_file(ffp)
            para.append(pstr)
            tpara.append(tstr)
        with open(floc + f'{item}.py', 'w') as ofile:
            ofile.write('\n'.join(para))
        with open(f'test_{item}.py', 'w') as ofile:
            ofile.write('\n'.join(tpara))




if __name__ == '__main__':
    # parse_mat_file('BoucWen.rst')
    # parse_mat_file('Bond_SP01.rst')
    import user_paths as up
    # parse_mat_file(up.OPY_DOCS_PATH + 'MinMax.rst')
    parse_all_uniaxial_mat()
    # defo = 'a2*k'
    # if any(re.findall('|'.join(['\*', '\/', '\+', '\-', '\^']), defo)):
    #     print('found')