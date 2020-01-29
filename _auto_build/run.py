import os
import re
from collections import OrderedDict
import pandas as pd
import inspect
from inspect import getmembers, isclass
import copy

glob_list = []
ROOT_DIR = os.path.dirname(os.path.abspath(__file__)) + "/../"
w4 = '    '
w8 = '        '
w12 = w8 + w4
w16 = w8 + w8
pname_pat = '\``([A-Za-z0-9_\./\\\'-]*)\``'
dtype_pat = '\|([A-Za-z0-9_\./\\-\ ]*)\|'
optype_pat = "\'([A-Za-z0-9_\./\\-]*)\'"

special_words = {
    'lambda': 'lamb',
    'type': 'otype',
    'as': 'a_s'
}

def clean_param_names(params, base_type):
    pms = OrderedDict()
    for pm in params:
        new_pm = pm
        dtype_is_obj = False
        if len(pm) == 1 and pm.istitle():
            if base_type in ['uniaxialMaterial', 'nDMaterial'] and pm in ['E', 'G', 'K']:
                new_pm = f'{pm.lower()}_mod'
            elif base_type in ['section'] and pm == 'A':
                new_pm = 'area'
            else:
                new_pm = 'big_' + pm.lower()
        else:
            new_pm = convert_camel_to_snake(pm)
        if new_pm in special_words:
            new_pm = special_words[new_pm]

        if len(new_pm) > 4 and new_pm[-4:] == '_tag':
            new_pm = new_pm[:-4]
            dtype_is_obj = True
        if len(new_pm) > 5 and new_pm[-5:] == '_tags':
            new_pm = new_pm[:-5] + 's'
        # if pm == 'eleNodes':

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
    name = name.replace('2d', '2D')
    name = name.replace('3d', '3D')
    if name[0] in '0123456789':
        name = 'N' + name
    return name


def convert_camel_to_snake(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    s1 = re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()
    s1 = re.sub(r'(\d)\_', r'\1', s1)  # removes underscore after number (e.g. 2D)
    return s1


def constructor(base_type, op_type, defaults, op_kwargs, osi_type, cl_name_suf="", obj_blurb=""):
    df_ip = pd.read_csv('force_in_place.csv')
    df_ip = df_ip[(df_ip['base_type'] == base_type) & (df_ip['op_type'] == op_type)]
    if len(op_kwargs) == 1:
        iis = list(op_kwargs)
        item = iis[0]
        pms = op_kwargs[item]
        if len(pms) == 1:
            defaults[pms[0]].marker = item[1:]  # remove '-' since it gets added again
            del op_kwargs[item]

    kw_pms = []  # TODO: This only works if object should be split into many objects based on kwargs, but
    # some kwargs are just to input a single value.
    for kw in op_kwargs:
        for pm in op_kwargs[kw]:
            kw_pms.append(pm)
    if not len(op_kwargs):
        op_kwargs[''] = []
    para = ['']
    tpara = []
    for kw in op_kwargs:
        name_from_kw = kw.replace('-', '')
        name_from_kw = name_from_kw.replace("'", '')
        if op_type is None:
            base_class_name = convert_name_to_class_name(base_type)
            # base_class_name = base_type[0].capitalize() + base_type[1:]
            para.append(f'class {base_class_name}(OpenseesObject):')
            cur_obj_class_name = base_type
            para += build_obj_docstring(base_type, base_class_name, obj_blurb)
            para.append(w4 + f"op_base_type = '{base_type}'")
            para.append('')
        else:
            op_class_name = convert_name_to_class_name(op_type + name_from_kw) + cl_name_suf
            base_class_name = convert_name_to_class_name(base_type)
            # base_class_name = base_type[0].capitalize() + base_type[1:]
            para.append(f'class {op_class_name}({base_class_name}Base):')
            cur_obj_class_name = op_class_name
            para += build_obj_docstring(op_class_name, base_class_name, obj_blurb)
            para.append(w4 + f"op_type = '{op_type}'")
            para.append('')

        pms = clean_param_names(defaults, base_type)
        # Determine what is in this class
        cl_pms = []  # TODO: need to add op_kwarg str
        for pm in pms:
            if pms[pm].org_name not in kw_pms or pm in op_kwargs[kw]:
                cl_pms.append(pm)

        # Build definition string
        # check there are no non-default pms are after pms that have defaults
        contains_no_default = False
        non_defaults_reversed = []
        con_defaults_reversed = []
        for i in range(len(cl_pms)):
            pm = cl_pms[-1-i]
            if pms[pm].default_is_expression and contains_no_default:
                if len(df_ip):
                    pms[pm].forced_not_optional = True  # TODO: should raise a flag
                # raise ValueError
                pass
            if pms[pm].default_value is None and not pms[pm].marker and not pms[pm].is_flag and not pms[pm].depends_on:  # TODO: may need to switch to .has_default in case default=None
                contains_no_default = True
                non_defaults_reversed.append(pm)
            else:
                con_defaults_reversed.append(pm)
        if not len(df_ip):
            inp_pms_order = non_defaults_reversed[::-1] + con_defaults_reversed[::-1]
        else:
            inp_pms_order = list(cl_pms)
        pjoins = ['osi']
        for pm in inp_pms_order:  # TODO: if packed and obj
            o3_name = pms[pm].o3_name
            default = pms[pm].default_value
            if pms[pm].is_flag:
                pjoins.append(f'{o3_name}=False')
            elif default is not None:
                if pms[pm].forced_not_optional:
                    pjoins.append(f'{o3_name}')
                elif pms[pm].default_is_expression:
                    pjoins.append(f'{o3_name}: float=None')
                    pms[pm].o3_default_is_none = True
                else:
                    if pms[pm].marker or pms[pm].depends_on:
                        if pms[pm].dtype in ['str', 'float', 'int']:
                            pjoins.append(f'{o3_name}: {pms[pm].dtype}=None')
                        else:
                            pjoins.append(f'{o3_name}=None')  # cannot have value for marker
                        pms[pm].o3_default_is_none = True
                    else:
                        pjoins.append(f'{o3_name}={default}')
            else:
                if pms[pm].marker or pms[pm].depends_on:
                    if pms[pm].dtype in ['str', 'float', 'int']:
                        pjoins.append(f'{o3_name}: {pms[pm].dtype}=None')
                    else:
                        pjoins.append(f'{o3_name}=None')
                    pms[pm].o3_default_is_none = True
                else:
                    pjoins.append(f'{o3_name}')

        pjoined = ', '.join(pjoins)
        para.append(f'    def __init__(self, {pjoined}):')
        para += build_init_method_docstring(cur_obj_class_name, pms, inp_pms_order)

        # Create init function saving logic
        for i, pm in enumerate(cl_pms):
            o3_name = pms[pm].o3_name
            dtype = pms[pm].dtype
            w_extra = ''
            extra = []
            if pms[pm].o3_default_is_none:
                extra.append(w8 + f'if {o3_name} is None:')
                extra.append(w8 + w4 + f'self.{o3_name} = None')
                extra.append(w8 + f'else:')
                w_extra = w4
            if dtype == 'float':
                para += extra
                para.append(w8 + w_extra + f'self.{o3_name} = float({o3_name})')
            elif dtype == 'int':
                para += extra
                para.append(w8 + w_extra + f'self.{o3_name} = int({o3_name})')
            elif dtype == 'obj':
                para.append(w8 + f'self.{o3_name} = {o3_name}')
            else:
                if pms[pm].list_items_dtype == 'obj':
                    para += extra
                    para.append(w8 + w_extra + f'self.{o3_name} = [x.tag for x in {o3_name}]')
                else:
                    para.append(w8 + f'self.{o3_name} = {o3_name}')
        pjoins = []
        if op_type is not None:
            pjoins.append('self.op_type')
        if osi_type is not None:
            para.append(w8 + f'osi.n_{osi_type} += 1')
            para.append(w8 + f'self._tag = osi.n_{osi_type}')
            pjoins += ['self._tag']
        need_special_logic = False
        applied_op_warg = False
        for pm in cl_pms:
            o3_name = pms[pm].o3_name
            if pms[pm].marker or pms[pm].depends_on:
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
        para.append(w8 + 'self._parameters = [%s]' % (', '.join(pjoins)))
        for pm in cl_pms:
            o3_name = pms[pm].o3_name
            if pms[pm].packed:
                ps = '*'
            else:
                ps = ''

            if pms[pm].marker:
                # para.append(w8 + f"if getattr(self, '{o3_name}') not in [None, '']:")
                para.append(w8 + f"if getattr(self, '{o3_name}') is not None:")
                tt = ''
                if pms[pm].dtype == 'obj':
                    tt = '.tag'
                para.append(w8 + w4 + f"self._parameters += ['-{pms[pm].marker}', {ps}self.{o3_name}{tt}]")

            elif pms[pm].depends_on:
                d_o3 = pms[pms[pm].depends_on].o3_name
                para.append(w8 + f"if getattr(self, '{o3_name}') is not None:")
                para.append(w8 + w4 + f"if getattr(self, '{d_o3}') is None:")
                para.append(w8 + w8 + f"raise ValueError('Cannot set: {o3_name} and not: {d_o3}')")
                para.append(w8 + w4 + f"self._parameters += [{ps}self.{o3_name}]")

            if pms[pm].is_flag:
                para.append(w8 + f"if getattr(self, '{o3_name}'):")
                para.append(w8 + w4 + f"self._parameters += ['-{pms[pm].org_name}']")  # TODO: does this always work?
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
            para.append(w8 + w4 + 'else:')
            para.append(w8 + w8 + 'break')
        para.append(w8 + 'self.to_process(osi)')
        para.append('')

        if op_type is None:
            low_op_name = convert_camel_to_snake(base_class_name)
            op_class_name = base_class_name
        else:
            low_op_name = convert_camel_to_snake(op_class_name)
        low_base_name = convert_camel_to_snake(base_class_name)  # TODO: should actually be file name

        # Build test
        names = {
            'low_op_name': low_op_name,
            'low_base_name': low_base_name,
            'op_class_name': op_class_name
        }
        # if low_base_name == 'element':
        #     tpara1 = build_test_for_element(names, pms, cl_pms)
        # elif low_base_name == 'beam_integration':
        #     tpara1 = build_test_for_beamintegration(names, pms, cl_pms)
        # else:
        tpara1 = build_test_for_generic(names, pms, cl_pms)
        tpara += tpara1
    return '\n'.join(para), '\n'.join(tpara)


def build_test_for_generic(names, pms, cl_pms):
    tpara = [f'def test_{names["low_op_name"]}():', w4 + 'osi = o3.OpenseesInstance(ndm=2)']
    prior_strs = []
    pjoins = ['osi']
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
            if o3_name == 'transf':
                prior_strs.append(w4 + f'{o3_name} = o3.transformation.Linear(osi, [])')
            elif o3_name in ['node', 'i_node']:
                prior_strs.append(w4 + f'{o3_name} = o3.node.Node(osi, 0.0, 0.0)')
            elif o3_name == 'j_node':
                prior_strs.append(w4 + 'j_node = o3.node.Node(osi, 0.0, 1.0)')
            elif o3_name == 'mat':  # TODO: detect ndmaterial
                prior_strs.append(w4 + f'{o3_name} = o3.uniaxial_material.Elastic(osi, 1.0)')
            if o3_name in ['sec', 'sec_i', 'sec_j', 'sec_e', 'section']:
                prior_strs.append(w4 + f'{o3_name} = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)')
            elif o3_name == 'integration':
                prior_strs.append(w4 + 'sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)')
                prior_strs.append(w4 + 'integration = o3.beam_integration.Lobatto(osi, sec, 5)')
            pjoins.append(f'{o3_name}={o3_name}')
        elif dtype == 'listi':
            if o3_name == 'ele_nodes':
                prior_strs.append(w4 + 'coords = [[0, 0], [1, 0], [1, 1], [0, 1]]')
                prior_strs.append(w4 + 'ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]')
            else:
                prior_strs.append(w4 + f'{o3_name} = [1, 1]')
            pjoins.append(f'{o3_name}={o3_name}')
        elif dtype == 'listf':
            prior_strs.append(w4 + f'{o3_name} = [1.0, 1.0]')
            pjoins.append(f'{o3_name}={o3_name}')
        elif dtype == 'str':
            pjoins.append(f'{o3_name}="string"')
        else:
            pjoins.append(f'{o3_name}=1')
    pjoint = ', '.join(pjoins)
    tpara += prior_strs
    tpara.append(w4 + f'o3.{names["low_base_name"]}.{names["op_class_name"]}({pjoint})')
    tpara.append('')
    tpara.append('')

    return tpara


def build_obj_docstring(op_class_name, base_class_name, obj_blurb):
    para = []
    para.append(w4 + '"""')
    para.append(w4 + f'The {op_class_name} {base_class_name} Class')
    para.append(w4 + '')
    para.append(force_line_char_limit(w4 + f'{obj_blurb}', w4))
    para.append(w4 + '"""')
    return para

def force_line_char_limit(line, indent):
    """
    If line is longer than limit then create new line at a space in the text

    :param line:
    :return:
    """
    clim = 120
    if len(line) <= clim:
        return line
    rem = line
    oline = ''
    for i in range(clim):
        j = clim - i
        if rem[j] == ' ':
            oline += rem[:j] + '\n' + indent
            rem = rem[j + 1:]
            if len(rem) < clim:
                oline += rem
                break
    return oline


def build_init_method_docstring(classname, pms, pms_ordered):
    init_blurb = f'Initial method for {classname}'
    dstr = [w8 + '"""', w8 + init_blurb, '', w8 + 'Parameters', w8 + '----------']
    for pm in pms_ordered:
        op_str = ''
        if pms[pm].default_is_expression:
            op_str = f' (default={pms[pm].default_is_expression})'
        dstr.append(w8 + f'{pms[pm].o3_name}: {pms[pm].dtype}{op_str}')
        descr = force_line_char_limit(w12 + f'{pms[pm].p_description.capitalize()}', w12)
        dstr.append(descr)
    dstr.append(w8 + '"""')
    return dstr


class Param(object):
    def __init__(self, org_name, default_value, packed=None, dtype=None):
        self.org_name = org_name
        self.o3_name = None
        self.default_value = default_value
        self.packed = packed
        self.dtype = dtype
        self.list_items_dtype = None
        self.default_is_expression = False
        self.p_description = ''
        self.marker = None
        # self.shared_marker = False  # if marker is used for
        self.is_flag = False
        self.forced_not_optional = False  # if default value precedes non default values
        self.depends_on = None
        self.o3_default_is_none = False


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


suffixes = ['', 'Args', 'Tag', 'Tags', 'MatTag', 'MatTags', 'Flag', 'Vals', 'SeriesTag', 's']


def clean_fn_line(line, has_tag=True):
    df = pd.read_csv('overrides.csv')
    defaults = OrderedDict()
    print(line)
    base_type = line.split('.. function:: ')[-1]
    base_type = base_type.split('(')[0]
    optype_res = re.search(optype_pat, line)
    if optype_res is None:
        optype = None
    else:
        optype = optype_res.group()[1:-1]
    df_op = df[(df['base_type'] == base_type) & (df['op_type'] == optype)]

    print(optype_res)
    inputs_str = line.split(')')[0]
    inputs = inputs_str.split(',')
    if has_tag:
        inputs = inputs[2:]  # remove class definition and tag
    else:
        inputs = inputs[1:]
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
            if len(inputs) > j + 1:
                suf_match = 0
                for suffix in suffixes:
                    if word + suffix == names_only[j + 1]:
                        markers[names_only[j + 1]] = word
                        inp = word
                        suf_match = 1
                        break
                if suf_match:
                    continue
                elif '-' in names_only[j + 1]:  # TODO: unsure if this is best way to identify flags
                    flags.append(word)
                    inp = word
                else:
                    cur_kwarg = '-' + word
                    op_kwargs[cur_kwarg] = []
                    continue
            else:
                flags.append(word)
                inp = word
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
    if "'-doRayleigh'" in defaults and 'rFlag' in defaults:  # e.g. ZeroLength
        del defaults["'-doRayleigh'"]

    return base_type, optype, defaults, op_kwargs


def trim_leading_whitespace(line):
    while line[0] in [' ', '\t']:
        line = line[1:]
        if not len(line):
            break
    return line


def parse_single_file(ffp, osi_type, expected_base_type=None):
    print('process: ', ffp)
    a = open(ffp, encoding="utf8")
    f = a.read()
    lines = f.splitlines()
    doc_str_pms = []
    dtypes = []
    defaults = None
    base_type = None
    optype = None
    defaults1 = None
    base_type1 = None
    op_kwargs1 = {}
    optype1 = None
    op_kwargs = OrderedDict()
    descriptions = {}
    fn_line_counter = 0
    doc_string_open = 0
    title_marker = 0
    two_defs = ''
    pstr = ''
    tstr = ''
    cur_res = []
    obj_blurb = ''
    sub_obj_blurbs = ['']
    ipara = []
    for line in lines:
        if len(line) > 3 and line[:3] == '===':
            title_marker += 1
            continue
        if ' ===' in line or '\t===' in line:
            if not doc_string_open:
                doc_string_open = 1
            else:
                if not fn_line_counter:
                    raise ValueError
                cl_name_suf = ''
                if two_defs == '2Dand3D':
                    cl_name_suf = '2D'
                if len(obj_blurb):
                    cur_obj_blurb = obj_blurb + '\n\n' + w4 + sub_obj_blurbs[0]
                else:
                    cur_obj_blurb = sub_obj_blurbs[0]
                pstr1, tstr1 = refine_and_build(doc_str_pms, dtypes, defaults, op_kwargs, descriptions, optype,
                                                base_type, osi_type, cl_name_suf, cur_obj_blurb)
                pstr += pstr1
                tstr += tstr1
                if two_defs == '2Dand3D':
                    cl_name_suf = '3D'
                    if len(obj_blurb):
                        cur_obj_blurb = obj_blurb + '\n\n' + w4 + sub_obj_blurbs[1]
                    else:
                        cur_obj_blurb = sub_obj_blurbs[1]
                    if defaults1 is None:
                        raise ValueError(defaults, op_kwargs)
                    pstr1, tstr1 = refine_and_build(doc_str_pms, dtypes, defaults1, op_kwargs1, descriptions, optype1,
                                                    base_type1, osi_type, cl_name_suf, cur_obj_blurb)
                    pstr += '\n' + pstr1
                    tstr += tstr1
                # Reset in case two objects in same file
                doc_str_pms = []
                dtypes = []
                doc_string_open = 0
                fn_line_counter = 0
                sub_obj_blurbs = []
                defaults = None
                base_type = None
                optype = None
                op_kwargs = OrderedDict()
                descriptions = {}
            continue
        char_only = line.replace(' ', '')
        char_only = char_only.replace('\t', '')
        if not len(char_only):
            continue

        first_char = char_only[0]
        res = re.search(pname_pat, line)
        if first_char != '*' and res:
            pname = res.group()[2:-2]
            print(pname)
            # if len(res.group()) > 4 and "'-" == res.group()[2:4]:
            #     continue  # op_kwarg
            ei = line.find('|')
            dtype_res = re.search(dtype_pat, line)
            if dtype_res is None:
                continue
            dtype = dtype_res.group()[1:-1]
            if dtype == 'str':
                pass  # TODO: look for options: e.g. (options: ``'beamtop'``, ``'beambot'`` or ``'column'``)
            des = line[dtype_res.end():]
            des = des.replace('\t', ' ')
            if not len(des):
                continue
            # while des[0] == ' ':
            #     des = des[1:]
            #     if not len(des):
            #         break
            des = trim_leading_whitespace(des)

            res = re.findall(pname_pat, line[:ei])
            for pm in res:
                if len(pm) > 4 and "'-" == pm[0:2]:
                    pm = pm[2:-1]
                doc_str_pms.append(pm)
                dtypes.append(dtype)
                descriptions[pm] = des
            cur_res = list(res)
        elif '.. function:: ' in line:
            fn_line_counter += 1
            sub_obj_blurbs.append('')
            if base_type is None:
                base_type, optype, defaults, op_kwargs = clean_fn_line(line, osi_type)
                base_class_name = convert_name_to_class_name(base_type)
                if expected_base_type is not None and expected_base_type != base_class_name:
                    if optype is None:
                        pass
                    else:
                        cam_case = convert_camel_to_snake(base_class_name)
                        new_istr = f'from o3seespy.command.{cam_case} import {base_class_name}Base'
                        if new_istr not in ipara:
                            ipara.append(new_istr)
            else:  # multiple function definitions
                glob_list.append(f'{base_type}-{optype}')

                base_type1, optype1, defaults1, op_kwargs1 = clean_fn_line(line, osi_type)
                df_td = pd.read_csv('two_definitions.csv')
                df_td = df_td[(df_td['base_type'] == base_type) & (df_td['op_type'] == optype)]
                two_defs = df_td['option'].iloc[0]
                if not len(df_td):
                    raise ValueError(f'{base_type}-{optype}')
                assert base_type == base_type1
                assert optype == optype1, (optype, optype1)
                if two_defs == 'combine':
                    for inp in defaults1:
                        if inp not in defaults:
                            defaults[inp] = defaults1[inp]
                    for inp in op_kwargs1:
                        if inp not in op_kwargs:
                            op_kwargs[inp] = op_kwargs1[inp]
        elif doc_string_open and len(cur_res):
            line = trim_leading_whitespace(line)
            # while line[0] == ' ':
            #     line = line[1:]
            #     if not len(line):
            #         break
            if len(line):
                line = ' ' + line
            for pm in cur_res:
                descriptions[pm] += line
        elif len(line) > 7 and '   :noi' in line[:7]:
            # if len(line) > 7 and '   :ref' in line[:7]:
            #     pass
            # else:
            continue
        elif not doc_string_open:
            line = line.replace(':ref:', '')
            if fn_line_counter:
                line = trim_leading_whitespace(line)
                sub_obj_blurbs[fn_line_counter - 1] += line
            elif title_marker == 2:
                line = trim_leading_whitespace(line)
                obj_blurb += line

    if base_type is not None:  # when there are no inputs
        cl_name_suf = ""
        if len(obj_blurb):
            cur_obj_blurb = obj_blurb + '\n\n' + w4 + sub_obj_blurbs[fn_line_counter - 1]
        else:
            cur_obj_blurb = sub_obj_blurbs[fn_line_counter - 1]
        pstr1, tstr1 = refine_and_build(doc_str_pms, dtypes, defaults, op_kwargs, descriptions, optype, base_type,
                                        osi_type, cl_name_suf, cur_obj_blurb)
        pstr += pstr1
        tstr += tstr1
    istr = '\n'.join(ipara)
    return pstr, tstr, istr


def refine_and_build(doc_str_pms, dtypes, defaults, op_kwargs, descriptions, optype, base_type, osi_type, cl_name_suf="", obj_blurb=''):
    if osi_type is not None:
        doc_str_pms = doc_str_pms[1:]  # remove mat tag
        dtypes = dtypes[1:]
    print('doc_str: ', doc_str_pms)
    print('fn_inps: ', list(defaults))
    print('op_kwargs: ', list(op_kwargs))

    for i, pm in enumerate(doc_str_pms):
        if '-' + pm in op_kwargs:
            continue
        if pm not in defaults:
            continue
        if defaults[pm].org_name.endswith('Tag'):
            defaults[pm].dtype = 'obj'
        elif defaults[pm].org_name.endswith('Tags'):
            defaults[pm].list_items_dtype = 'obj'
        else:
            defaults[pm].dtype = dtypes[i]
        defaults[pm].p_description = descriptions[pm]

    if 'eleNodes' in defaults:
        defaults['eleNodes'].list_items_dtype = 'obj'
    hidden_objs = ['iNode', 'jNode', 'sec', 'secI', 'secJ', 'secE']
    for obj in hidden_objs:
        if obj in defaults:
            defaults[obj].dtype = 'obj'

    if ("'-orient'" in defaults or '-orient' in op_kwargs) and 'vecx' in defaults and 'vecyp' in defaults:
        try:
            del defaults["'-orient'"]
        except KeyError:
            del op_kwargs['-orient']
        del defaults['vecx']
        del defaults['vecyp']
        defaults['orient'] = Param(org_name='orient', default_value=None, packed=True)
        defaults['orient'].marker = '-orient'
    if "-orient" in op_kwargs and 'x1' in defaults:  # Element
        del op_kwargs["-orient"]
        del defaults['x1']
        del defaults['x2']
        del defaults['x3']
        if 'yp1' in defaults:
            del defaults['yp1']
            del defaults['yp2']
            del defaults['yp3']
        if 'y1' in defaults:
            del defaults['y1']
            del defaults['y2']
            del defaults['y3']
        defaults['orient'] = Param(org_name='orient', default_value=None, packed=True)
        defaults['orient'].marker = 'orient'
    if "-mass" in op_kwargs and 'm' in defaults:  # Element
        del op_kwargs['-mass']
        defaults['mass'] = copy.deepcopy(defaults['m'])
        defaults['mass'].marker = 'mass'
        del defaults['m']
    if "-mass" in op_kwargs and 'massDens' in defaults:  # Element
        del op_kwargs['-mass']
        defaults['mass'] = copy.deepcopy(defaults['massDens'])
        defaults['mass'].marker = 'mass'
        del defaults['massDens']
    if "-shearDist" in op_kwargs and 'sDratio' in defaults:  # Element
        del op_kwargs['-shearDist']
        defaults['shear_dist'] = copy.deepcopy(defaults['sDratio'])
        defaults['shear_dist'].marker = 'shearDist'
        del defaults['sDratio']
    if "-iter" in op_kwargs and 'maxIter' in defaults:
        del op_kwargs['-iter']
        # defaults['iter'] = copy.deepcopy(defaults['maxIter'])
        defaults['maxIter'].marker = 'iter'
        if 'tol' in defaults:
            defaults['tol'].depends_on = 'maxIter'
    if "-jntOffset" in op_kwargs and 'dI' in defaults:
        del op_kwargs['-jntOffset']
        # defaults['iter'] = copy.deepcopy(defaults['maxIter'])
        defaults['dI'].marker = 'jntOffset'
        if 'dJ' in defaults:
            defaults['dJ'].depends_on = 'dI'
        # del defaults['maxIter']
    #assert len(doc_str_pms) == len(defaults) + len(op_kwargs), (len(doc_str_pms), (len(defaults), len(op_kwargs)))
    # if len(op_kwargs) == 1:
    #     opk = list(op_kwargs)
    #     for item in opk:
    #         defaults[item] = Param(org_name=item, default_value=False, packed=False)
    #         defaults[item].marker = f'-{item}'
    #         del op_kwargs[item]
    from _auto_build import _custom_gen as cust_file

    cust_obj_list = [o[0] for o in getmembers(cust_file) if isclass(o[1])]
    if optype is not None:
        op_class_name = convert_name_to_class_name(optype)
        if op_class_name in cust_obj_list:
            source = inspect.getsource(getattr(cust_file, op_class_name))
            return source + '\n', ""

    pstr, tstr = constructor(base_type, optype, defaults, op_kwargs, osi_type=osi_type, cl_name_suf=cl_name_suf, obj_blurb=obj_blurb)
    print(pstr, tstr)
    return pstr, tstr
    # if line[3:5] == '``':
    #     para = line[5:]

unimats = {
    'Steel & Reinforcing-Steel Materials': 'steel',
    'Concrete Materials': 'concrete',
    'Standard Uniaxial Materials': 'standard',
    'PyTzQz uniaxial materials': 'pytz',
    'Other Uniaxial Materials': 'other'
}


def parse_all_uniaxial_mat():
    import user_paths as up
    from _auto_build import _custom_mat as cust_file

    cust_obj_list = [o[0] for o in getmembers(cust_file) if isclass(o[1])]
    uni_axial_mat_file = open(up.OPY_DOCS_PATH + 'uniaxialMaterial.rst')
    lines = uni_axial_mat_file.read().split('\n')
    collys = {}
    mtype = None
    for line in lines:
        m_found = 0
        for mopt in unimats:
            if mopt in line:
                mtype = unimats[mopt]
                collys[mtype] = []
                m_found = 1
                break
        if m_found:
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
        ipara = []
        print(item, collys[item])
        for mat in collys[item]:
            # if mat == 'PyLiq1':
            #     continue
            if mat == 'steel4':  # consider as a single object, make all values after E0 equal to None, all -ps should be flags.
                continue
            if mat == 'Pinching4':
                continue
            if mat == 'KikuchiAikenHDR' or mat == 'KikuchiAikenLRB':
                continue
            if mat == 'PinchingLimitStateMaterial':
                continue
            mat_name = convert_name_to_class_name(mat)
            if mat_name in cust_obj_list:
                source = inspect.getsource(getattr(cust_file, mat_name))
                print(source)
                para.append('')
                para.append(source)
                continue
            open(up.OPY_DOCS_PATH + '%s.rst' % mat)
            ffp = up.OPY_DOCS_PATH + '%s.rst' % mat
            pstr, tstr, istr = parse_single_file(ffp, osi_type='mat')
            if istr not in ipara and istr != '':
                ipara.append(istr)
            para.append(pstr)
            tpara.append(tstr)
        with open(floc + f'{item}.py', 'w') as ofile:
            ofile.write('\n'.join(ipara))
            if len(ipara):
                ofile.write('\n')
            ofile.write('\n'.join(para))
        with open(f'temp_tests/atest_uniaxial_{item}.py', 'w') as ofile:
            ofile.write('\n'.join(tpara))

ndmats = {
    'Standard Models': 'standard',
    'Materials for Modeling Concrete Walls': 'concrete_walls',
    'Tsinghua Sand Models': 'tsinghua_sand',
    'Contact Materials for 2D and 3D': 'contact',
    'Wrapper material for Initial State Analysis': 'wrapper',
    'UC San Diego soil models': 'uc_san_diego_soil',
    'UC San Diego Saturated Undrained soil': 'uc_san_diego_ud_soil'
}


def parse_all_ndmat():
    import user_paths as up
    from _auto_build import _custom_mat as cust_file

    cust_obj_list = [o[0] for o in getmembers(cust_file) if isclass(o[1])]
    mat_file = open(up.OPY_DOCS_PATH + 'ndMaterial.rst')
    lines = mat_file.read().split('\n')
    collys = {}
    mtype = None
    for line in lines:
        m_found = 0
        for mopt in ndmats:
            if mopt in line:
                mtype = ndmats[mopt]
                collys[mtype] = []
                m_found = 1
                break
        if m_found:
             continue
        if mtype is not None:
            line = line.replace(' ', '')
            line = line.replace('\t', '')
            if ':' in line or '-' in line or '#' in line  or line == '':
                continue
            collys[mtype].append(line)

    floc = ROOT_DIR + 'o3seespy/command/nd_material/'
    for item in collys:
        para = ['from o3seespy.command.nd_material.base_material import NDMaterialBase', '']
        tpara = ['import o3seespy as o3  # for testing only', '', '']
        ipara =[]
        print(item, collys[item])
        for mat in collys[item]:
            if mat == 'PressureDependMultiYield':
                continue
            if mat == 'PressureDependMultiYield02':
                continue
            if mat == 'StressDensityModel':
                continue
            if mat == 'PM4Sand':
                continue
            mat_name = convert_name_to_class_name(mat)
            if mat_name in cust_obj_list:
                source = inspect.getsource(getattr(cust_file, mat))
                print(source)
                para.append('')
                para.append(source)
                continue

            open(up.OPY_DOCS_PATH + '%s.rst' % mat)
            ffp = up.OPY_DOCS_PATH + '%s.rst' % mat
            pstr, tstr, istr = parse_single_file(ffp, osi_type='mat')
            if istr not in ipara and istr != '':
                ipara.append(istr)
            para.append(pstr)
            tpara.append(tstr)
        with open(floc + f'{item}.py', 'w') as ofile:
            ofile.write('\n'.join(ipara))
            if len(ipara):
                ofile.write('\n')
            ofile.write('\n'.join(para))
        with open(f'temp_tests/atest_{item}.py', 'w') as ofile:
            ofile.write('\n'.join(tpara))


eles = {
    'Zero-Length Element': 'zero_length_a',
    'Truss Elements': 'truss',
    'Beam-Column Elements': 'beam_column',
    'Joint Elements': 'joint',
    'Link Elements': 'link',
    'Bearing Elements': 'bearing',
    'Quadrilateral Elements': 'quadrilateral',
    'Triangular Elements': 'triangular',
    'Brick Elements': 'brick',
    'Tetrahedron Elements': 'tetrahedron',
    'UC San Diego u-p element (saturated soil)': 'uc_san_diego_up',
    'Other u-p elements': 'other_up',
    'Contact Elements': 'contact',
    'Cable Elements': 'cable',
    'PFEM Elements': 'pfem',
    'Misc.': 'misc'

}


def parse_all_elements():
    import user_paths as up
    ele_file = open(up.OPY_DOCS_PATH + 'element.rst')
    lines = ele_file.read().split('\n')
    collys = {}
    etype = None
    for line in lines:
        m_found = 0
        for mopt in eles:
            if mopt in line:
                etype = eles[mopt]
                collys[etype] = []
                m_found = 1
                break
        if m_found:
             continue
        if etype is not None:
            line = line.replace(' ', '')
            line = line.replace('\t', '')
            if ':' in line or '-' in line or '#' in line or line == '':
                continue
            collys[etype].append(line)

    floc = ROOT_DIR + 'o3seespy/command/element/'
    for item in collys:
        para = ['from o3seespy.command.element.base_element import ElementBase', '']
        tpara = ['import o3seespy as o3  # for testing only', '', '']
        ipara = []
        print(item, collys[item])
        for ele in collys[item]:
            if ele in ['trussEle', 'corotTruss', 'RJWatsonEqsBearing']:
                continue
            # if ele == 'zeroLengthND':
            #     continue
            # if mat == 'PressureDependMultiYield02':
            #     continue

            # open(up.OPY_DOCS_PATH + '%s.rst' % ele)
            ffp = up.OPY_DOCS_PATH + '%s.rst' % ele
            pstr, tstr, istr = parse_single_file(ffp, osi_type='ele')
            if istr not in ipara and istr != '':
                ipara.append(istr)
            para.append(pstr)
            tpara.append(tstr)
        with open(floc + f'{item}.py', 'w') as ofile:
            ofile.write('\n'.join(ipara))
            if len(ipara):
                ofile.write('\n')
            ofile.write('\n'.join(para))
        with open(f'temp_tests/atest_{item}.py', 'w') as ofile:
            ofile.write('\n'.join(tpara))


def parse_generic_single_file(obj_type, osi_type):
    o3_type = convert_camel_to_snake(obj_type)
    o3_class_type = convert_name_to_class_name(obj_type)
    import user_paths as up
    from _auto_build import _custom_gen as cust_file

    cust_obj_list = [o[0] for o in getmembers(cust_file) if isclass(o[1])]
    uni_axial_mat_file = open(up.OPY_DOCS_PATH + f'{obj_type}.rst')
    lines = uni_axial_mat_file.read().split('\n')
    collys = {o3_type: []}
    mtype = None
    for line in lines:
        if '.. toctree::' in line:
            mtype = o3_type
            continue

        if mtype is not None and len(line) and line[0] in [' ', '\t']:
            line = line.replace(' ', '')
            line = line.replace('\t', '')
            if ':' in line or '-' in line or '#' in line or line == '':
                continue
            collys[mtype].append(line)

    floc = ROOT_DIR + 'o3seespy/command/'
    for item in collys:
        para = ['from o3seespy.base_model import OpenseesObject', '', '']
        para += [f'class {o3_class_type}Base(OpenseesObject):']
        para += [w4 + f'op_base_type = "{obj_type}"', '']
        tpara = ['import o3seespy as o3  # for testing only', '', '']
        ipara = []
        print(item, collys[item])
        for mat in collys[item]:

            mat_name = convert_name_to_class_name(mat)
            if mat_name in cust_obj_list:
                source = inspect.getsource(getattr(cust_file, mat))
                print(source)
                para.append('')
                para.append(source)
                continue

            open(up.OPY_DOCS_PATH + '%s.rst' % mat)
            ffp = up.OPY_DOCS_PATH + '%s.rst' % mat
            pstr, tstr, istr = parse_single_file(ffp, osi_type=osi_type, expected_base_type=o3_class_type)
            if istr not in ipara and istr != '':
                ipara.append(istr)
            para.append(pstr)
            tpara.append(tstr)
        with open(floc + f'{item}.py', 'w') as ofile:
            ofile.write('\n'.join(ipara))
            if len(ipara):
                ofile.write('\n')
            ofile.write('\n'.join(para))
        with open(f'temp_tests/atest_{item}.py', 'w') as ofile:
            ofile.write('\n'.join(tpara))


def test_clean_fn_line():
    ln = ".. function:: element('singleFPBearing', eleTag,*eleNodes,frnMdlTag, Reff, kInit,'-P', PMatTag,'-T', TMatTag,'-My', MyMatTag,'-Mz', MzMatTag,['-orient',[x1, x2, x3], y1, y2, y3],['-shearDist', sDratio],['-doRayleigh'],['-mass', m],['-iter', maxIter, tol])"
    base_type, optype, defaults, op_kwargs = clean_fn_line(ln)
    print(base_type, optype)
    print(defaults)
    print(op_kwargs)
    print(defaults['tol'].default_value)  # TODO: add 'depends_on' based on []


if __name__ == '__main__':
    # parse_single_file('BoucWen.rst')
    # parse_single_file('Bond_SP01.rst')
    import user_paths as up
    #parse_all_ndmat()
    # ps, ts = parse_single_file(up.OPY_DOCS_PATH + 'nonlinearBeamColumn.rst', 'ele')
    all = 1
    # all = 1  # TODO: KikuchiBearing
    # TODO: dettach docstrings - if exists then don't use rst version
    # TODO: add type hinting for default None (w: str = None)
    if not all:
        parse_generic_single_file(obj_type='section', osi_type='sect')
        # print(ts)
        # parse_generic_single_file(obj_type='integrator', osi_type=None)
        # parse_generic_single_file(obj_type='beamIntegration', osi_type='integ')
        # ps, ts, istr = parse_single_file(up.OPY_DOCS_PATH + 'Quad.rst', 'ele')
        # parse_generic_single_file(obj_type='geomTransf', osi_type='transformation')
        # parse_generic_single_file(obj_type='beamIntegration', osi_type='integ')
        # print(ts)
        # parse_single_file(up.OPY_DOCS_PATH + 'UserDefined.rst', 'integ')
        # test_clean_fn_line()

    if all:
        parse_generic_single_file(obj_type='pattern', osi_type='pat')
        parse_generic_single_file(obj_type='timeSeries', osi_type='tseries')
        parse_generic_single_file(obj_type='constraints', osi_type=None)
        parse_generic_single_file(obj_type='integrator', osi_type=None)
        parse_generic_single_file(obj_type='beamIntegration', osi_type='integ')
        parse_generic_single_file(obj_type='section', osi_type='sect')
        parse_generic_single_file(obj_type='geomTransf', osi_type='transformation')
        parse_all_uniaxial_mat()
        parse_all_ndmat()
        parse_all_elements()
    ofile = open('temp.out', 'w')
    ofile.write('\n'.join(glob_list))
    # defo = 'a2*k'
    # if any(re.findall('|'.join(['\*', '\/', '\+', '\-', '\^']), defo)):
    #     print('found')
    # TODO: YamamotoBiaxialHDRorient, YamamotoBiaxialHDRcoRS


