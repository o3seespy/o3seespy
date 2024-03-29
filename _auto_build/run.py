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
    'as': 'a_s',
    'iter': 'max_iter'
}

eles = {
    'Zero-Length Element': 'zero_length',
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


unimats = {
    'Steel & Reinforcing-Steel Materials': 'steel',
    'Concrete Materials': 'concrete',
    'Standard Uniaxial Materials': 'standard',
    'PyTzQz uniaxial materials': 'pytz',
    'Other Uniaxial Materials': 'other'
}


ndmats = {
    'Standard Models': 'standard',
    'Materials for Modeling Concrete Walls': 'concrete_walls',
    'Tsinghua Sand Models': 'tsinghua_sand',
    'Contact Materials for 2D and 3D': 'contact',
    'Wrapper material for Initial State Analysis': 'wrapper',
    'UC San Diego soil models': 'uc_san_diego_soil',
    'UC San Diego Saturated Undrained soil': 'uc_san_diego_ud_soil'
}

extra_heir = {'element': eles,
              'uniaxial_material': unimats,
              'nd_material': ndmats}

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
        elif len(new_pm) > 4 and new_pm[-3:] == 'tag':
            new_pm = new_pm[:-3]
            dtype_is_obj = True
        elif len(new_pm) > 5 and new_pm[-5:] == '_tags':
            new_pm = new_pm[:-5]
            if new_pm[-1] != 's':
                new_pm += 's'
        elif len(new_pm) > 5 and new_pm[-5:-1] == '_tag':  # e.g. Tag1
            suf = new_pm[-1]
            new_pm = new_pm[:-5]
            new_pm += suf
        # if pm == 'eleNodes':
        while '__' in new_pm:
            new_pm = new_pm.replace('__', '_')
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
    s1 = s1.replace('__', '_')
    return s1


def build_additional_methods(obj_type, obj_name):
    if obj_type == 'nDMaterial':
        df = pd.read_csv('markup_lists/nD-set-parameters.csv')
    elif obj_type == 'uniaxialMaterial':
        df = pd.read_csv('markup_lists/uniaxial-set-parameters.csv')
    elif obj_type == 'section':
        df = pd.read_csv('markup_lists/section-set-parameters.csv')
    else:
        return []
    df.fillna('', inplace=True)
    df = df[df['obj'] == obj_name]
    if not len(df):
        return []
    print(obj_type, obj_name)
    o3_names = df['o3_name'].iloc[0].split('-')
    cpp_names = df['cpp_name'].iloc[0].split('-')
    if len(o3_names) == 1 and o3_names[0] == '':
        return []
    para = []
    for i, name in enumerate(o3_names):
        para.append(w4 + f'def set_{name}(self, value, ele=None, eles=None):')
        para.append(w8 + f"self.set_parameter(self.osi, '{cpp_names[i]}', value, ele, eles)")
        para.append('')
    return para


def constructor_of_func(base_name, op_type, defaults, op_kwargs, osi_type, cl_name_suf="", obj_blurb=[]):

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
    for kw in op_kwargs:  # TODO: add a csv of objects that should not split since optional args are not exclusive (e.g. FRPConfinedConcrete02)
        name_from_kw = kw.replace('-', '')
        name_from_kw = name_from_kw.replace("'", '')
        if op_type is None:
            raise ValueError

        pms = clean_param_names(defaults, '')
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
            if pms[pm].default_value is None and not pms[pm].marker and not pms[pm].is_flag and not pms[pm].depends_on:  # TODO: may need to switch to .has_default in case default=None
                contains_no_default = True
                non_defaults_reversed.append(pm)
            else:
                con_defaults_reversed.append(pm)
        inp_pms_order = non_defaults_reversed[::-1] + con_defaults_reversed[::-1]

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
                        elif pms[pm].dtype is not None and 'list' in pms[pm].dtype:
                            pjoins.append(f'{o3_name}: list=None')
                        else:
                            pjoins.append(f'{o3_name}=None')  # cannot have value for marker
                        pms[pm].o3_default_is_none = True
                    else:
                        pjoins.append(f'{o3_name}={default}')
            else:
                if pms[pm].marker or pms[pm].depends_on:
                    if pms[pm].dtype in ['str', 'float', 'int']:
                        pjoins.append(f'{o3_name}: {pms[pm].dtype}=None')
                    elif pms[pm].dtype is not None and 'list' in pms[pm].dtype:
                        pjoins.append(f'{o3_name}: list=None')
                    else:
                        pjoins.append(f'{o3_name}=None')
                    pms[pm].o3_default_is_none = True
                else:
                    pjoins.append(f'{o3_name}')

        pjoined = ', '.join(pjoins)
        low_op_name = convert_camel_to_snake(op_type)
        if 'Returns ' in obj_blurb[0] and 'get_' not in low_op_name:
            low_op_name = f'get_{low_op_name}'
        para.append(f'def {low_op_name}({pjoined}):')
        para += build_fn_docstring(obj_blurb, pms, inp_pms_order)
        import importlib

        # Try to insert the initialisation example from tests
        insert_example = 1
        if insert_example:

            source = None
            try:
                tcommands = importlib.import_module(f"tests.commands.test_{base_name}")
                tname = f'test_{low_op_name}'
                if tname in dir(tcommands):
                    source = inspect.getsource(getattr(tcommands, tname)).splitlines()
            except ModuleNotFoundError:
                pass
            if source is not None:
                para = para[:-1]  # remove triple quote at end of docstring
                para.append('')  # blank line
                para.append(w4 + 'Examples')
                para.append(w4 + '--------')
                para.append(w4 + '>>> import o3seespy as o3')

                for line in source[1:]:  # skip first line which says 'test_'
                    if 'def ' in line:
                        nline = w4 + '>>> ' + '# Example is currently not working'
                        para.append(nline)
                        continue
                    nline = w4 + '>>> ' + line[4:]
                    para.append(nline)
                para.append(w4 + '"""')

        # Create init function saving logic
        for i, pm in enumerate(cl_pms):
            o3_name = pms[pm].o3_name
            dtype = pms[pm].dtype
            w_extra = ''
            extra = []
            if pms[pm].o3_default_is_none and dtype in ['float', 'int']:
                extra.append(w4 + f'if {o3_name} is not None:')
                w_extra = w4
            if dtype == 'float':
                para += extra
                para.append(w4 + w_extra + f'{o3_name} = float({o3_name})')
            elif dtype == 'int':
                para += extra
                para.append(w4 + w_extra + f'{o3_name} = int({o3_name})')
            else:
                if pms[pm].list_items_dtype == 'obj':
                    para += extra
                    if o3_name == 'ele_nodes':
                        pms[pm].o3_pm_name = 'ele_node_tags'
                    para.append(w4 + w_extra + f'{pms[pm].o3_pm_name} = [x.tag for x in {o3_name}]')
        pjoins = []
        # if op_type is not None:
        #     pjoins.append(op_type)
        need_special_logic = False
        applied_op_warg = False
        for pm in cl_pms:
            o3_pm_name = pms[pm].o3_pm_name
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
                pjoins.append(f'{o3_pm_name}.tag')
            elif pms[pm].packed:
                pjoins.append('*' + o3_pm_name)
            else:
                pjoins.append('' + o3_pm_name)
        para.append(w4 + '_parameters = [%s]' % (', '.join(pjoins)))
        for pm in cl_pms:
            o3_pm_name = pms[pm].o3_pm_name
            if pms[pm].packed:
                ps = '*'
            else:
                ps = ''

            if pms[pm].marker:
                # para.append(w8 + f"if getattr(self, '{o3_name}') not in [None, '']:")
                para.append(w4 + f"if {o3_pm_name} is not None:")
                tt = ''
                if pms[pm].dtype == 'obj':
                    tt = '.tag'
                para.append(w4 + w4 + f"_parameters += ['-{pms[pm].marker}', {ps}{o3_pm_name}{tt}]")

            elif pms[pm].depends_on:
                d_o3 = pms[pms[pm].depends_on].o3_name
                para.append(w4 + f"if {o3_pm_name} is not None:")
                para.append(w4 + w4 + f"if {d_o3} is None:")
                para.append(w4 + w8 + f"raise ValueError('Cannot set: {o3_pm_name} and not: {d_o3}')")
                para.append(w4 + w4 + f"_parameters += [{ps}{o3_pm_name}]")

            if pms[pm].is_flag:
                para.append(w4 + f"if {o3_pm_name}:")
                para.append(w4 + w4 + f"_parameters += ['-{pms[pm].org_name}']")  # TODO: does this always work?
        if need_special_logic:
            sp_logic = False
            sp_pms = []
            for pm in pms:
                if pms[pm].default_is_expression:
                    sp_logic = True
                if sp_logic:
                    sp_pms.append(pm)
            sp_pm_strs = ["'%s'" % pms[pm].o3_pm_name for pm in sp_pms]
            para.append(w4 + f"special_pms = [{', '.join(sp_pm_strs)}]")
            packets = [str(pms[pm].packed) for pm in sp_pms]
            para.append(w4 + f"packets = [{', '.join(packets)}]")
            para.append(w4 + 'for i, pm in enumerate(special_pms):')
            para.append(w4 + w4 + 'if pm is not None:')
            para.append(w4 + w8 + 'if packets[i]:')
            para.append(w4 + w8 + w4 + '_parameters += [*pm]')
            para.append(w4 + w8 + 'else:')
            para.append(w4 + w8 + w4 + '_parameters += [pm]')
            para.append(w4 + w4 + 'else:')
            para.append(w4 + w8 + 'break')

        para.append(w4 + f'return osi.to_process("{op_type}", _parameters)')
        para.append('')

        # Build test
        names = {
            'low_op_name': low_op_name,
            'low_base_name': base_name,
            'op_class_name': low_op_name
        }
        # if low_base_name == 'element':
        #     tpara1 = build_test_for_element(names, pms, cl_pms)
        # elif low_base_name == 'beam_integration':
        #     tpara1 = build_test_for_beamintegration(names, pms, cl_pms)
        # else:
        tpara1 = build_test_for_generic(names, pms, cl_pms)
        tpara += tpara1
    return '\n'.join(para), '\n'.join(tpara)


def constructor(base_type, op_type, defaults, op_kwargs, osi_type, cl_name_suf="", obj_blurb=[]):
    df_ip = pd.read_csv('markup_lists/force_in_place.csv')
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
    for kw in op_kwargs:  # TODO: add a csv of objects that should not split since optional args are not exclusive (e.g. FRPConfinedConcrete02)
        name_from_kw = kw.replace('-', '')
        name_from_kw = name_from_kw.replace("'", '')
        if op_type is None:
            base_class_name = convert_name_to_class_name(base_type)
            # base_class_name = base_type[0].capitalize() + base_type[1:]
            para.append(f'class {base_class_name}(OpenSeesObject):')
            cur_obj_class_name = base_type
            para += build_obj_docstring(base_type, base_class_name, obj_blurb)
            para.append(w4 + f"op_base_type = '{base_type}'")
            para.append('')
        else:
            op_class_name = convert_name_to_class_name(op_type + name_from_kw) + cl_name_suf
            base_class_name = convert_name_to_class_name(base_type)
            custom_strs = find_in_custom_defs(base_type, op_class_name)
            if custom_strs is not None:
                return custom_strs
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
                        elif pms[pm].dtype is not None and 'list' in pms[pm].dtype:
                            pjoins.append(f'{o3_name}: list=None')
                        else:
                            pjoins.append(f'{o3_name}=None')  # cannot have value for marker
                        pms[pm].o3_default_is_none = True
                    else:
                        pjoins.append(f'{o3_name}={default}')
            else:
                if pms[pm].marker or pms[pm].depends_on:
                    if pms[pm].dtype in ['str', 'float', 'int']:
                        pjoins.append(f'{o3_name}: {pms[pm].dtype}=None')
                    elif pms[pm].dtype is not None and 'list' in pms[pm].dtype:
                        pjoins.append(f'{o3_name}: list=None')
                    else:
                        pjoins.append(f'{o3_name}=None')
                    pms[pm].o3_default_is_none = True
                else:
                    pjoins.append(f'{o3_name}')

        pjoined = ', '.join(pjoins)
        para.append(f'    def __init__(self, {pjoined}):')
        para += build_init_method_docstring(cur_obj_class_name, pms, inp_pms_order)
        import importlib

        print(base_type)
        if op_type is None:
            low_op_name = convert_camel_to_snake(base_class_name)
            op_class_name = base_class_name
        else:
            low_op_name = convert_camel_to_snake(op_class_name)
        low_base_name = convert_camel_to_snake(base_class_name)  # TODO: should actually be file name

        # Try to insert the initialisation example from tests
        insert_example = 1
        if insert_example:

            source = None
            if low_base_name in extra_heir:
                fdict = extra_heir[low_base_name]
                if low_base_name == 'uniaxial_material':
                    estr = 'uniaxial_'
                else:
                    estr = ''

                for f_old in fdict:
                    f_new = fdict[f_old]
                    tcommands = importlib.import_module(f"tests.commands.{low_base_name}.test_{estr}{f_new}")
                    tname = f'test_{low_op_name}'
                    if tname in dir(tcommands):
                        source = inspect.getsource(getattr(tcommands, tname)).splitlines()
                        break
            else:
                try:
                    tcommands = importlib.import_module(f"tests.commands.test_{low_base_name}")
                    tname = f'test_{low_op_name}'
                    if tname in dir(tcommands):
                        source = inspect.getsource(getattr(tcommands, tname)).splitlines()
                except ModuleNotFoundError:
                    pass
            if source is not None:
                para = para[:-1]  # remove triple quote at end of docstring
                para.append('')  # blank line
                para.append(w8 + 'Examples')
                para.append(w8 + '--------')
                para.append(w8 + '>>> import o3seespy as o3')

                for line in source[1:]:  # skip first line which says 'test_'
                    if 'def ' in line:
                        nline = w8 + '>>> ' + '# Example is currently not working'
                        para.append(nline)
                        continue
                    nline = w8 + '>>> ' + line[4:]
                    para.append(nline)
                para.append(w8 + '"""')

        # Create init function saving logic
        para.append(w8 + 'self.osi = osi')
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
                    if o3_name == 'ele_nodes':
                        pms[pm].o3_pm_name = 'ele_node_tags'
                    para.append(w8 + w_extra + f'self.{pms[pm].o3_pm_name} = [x.tag for x in {o3_name}]')
                    if o3_name == 'ele_nodes':
                        para.append(w8 + w_extra + f'self.{o3_name} = {o3_name}')
                else:
                    para.append(w8 + f'self.{o3_name} = {o3_name}')
        pjoins = []
        if op_type is not None:
            pjoins.append('self.op_type')
        if osi_type is not None:
            if osi_type == 'mat':
                para.append(w8 + 'if osi is not None:')
                para.append(w12 + f'osi.n_{osi_type} += 1')
                para.append(w12 + f'self._tag = osi.n_{osi_type}')
            else:
                para.append(w8 + f'osi.n_{osi_type} += 1')
                para.append(w8 + f'self._tag = osi.n_{osi_type}')
            pjoins += ['self._tag']
        need_special_logic = False
        applied_op_warg = False
        for pm in cl_pms:
            o3_pm_name = pms[pm].o3_pm_name
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
                pjoins.append(f'self.{o3_pm_name}.tag')
            elif pms[pm].packed:
                pjoins.append('*self.' + o3_pm_name)
            else:
                pjoins.append('self.' + o3_pm_name)
        para.append(w8 + 'self._parameters = [%s]' % (', '.join(pjoins)))
        for pm in cl_pms:
            o3_pm_name = pms[pm].o3_pm_name
            if pms[pm].packed:
                ps = '*'
            else:
                ps = ''

            if pms[pm].marker:
                # para.append(w8 + f"if getattr(self, '{o3_name}') not in [None, '']:")
                para.append(w8 + f"if getattr(self, '{o3_pm_name}') is not None:")
                tt = ''
                if pms[pm].dtype == 'obj':
                    tt = '.tag'
                para.append(w8 + w4 + f"self._parameters += ['-{pms[pm].marker}', {ps}self.{o3_pm_name}{tt}]")

            elif pms[pm].depends_on:
                d_o3 = pms[pms[pm].depends_on].o3_name
                para.append(w8 + f"if getattr(self, '{o3_pm_name}') is not None:")
                para.append(w8 + w4 + f"if getattr(self, '{d_o3}') is None:")
                para.append(w8 + w8 + f"raise ValueError('Cannot set: {o3_pm_name} and not: {d_o3}')")
                para.append(w8 + w4 + f"self._parameters += [{ps}self.{o3_pm_name}]")

            if pms[pm].is_flag:
                para.append(w8 + f"if getattr(self, '{o3_pm_name}'):")
                para.append(w8 + w4 + f"self._parameters += ['-{pms[pm].org_name}']")  # TODO: does this always work?
        if need_special_logic:
            sp_logic = False
            sp_pms = []
            for pm in pms:
                if pms[pm].default_is_expression:
                    sp_logic = True
                if sp_logic:
                    sp_pms.append(pm)
            sp_pm_strs = ["'%s'" % pms[pm].o3_pm_name for pm in sp_pms]
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
        if osi_type == 'mat':
            para.append(w8 + 'if osi is None:')
            para.append(w12 + 'self.built = 0')
            para.append(w8 + 'if osi is not None:')
            para.append(w12 + 'self.to_process(osi)')
        else:
            para.append(w8 + 'self.to_process(osi)')
        para.append('')

        para += build_additional_methods(base_type, op_type)

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
    # TODO: detect 2D or 3D
    tpara = [f'def test_{names["low_op_name"]}():', w4 + 'osi = o3.OpenSeesInstance(ndm=2)']
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
            elif o3_name in ['j_node', 'c_node', 'r_node', 'l_node', 'lagr_node', 's_node', 'k_node']:
                prior_strs.append(w4 + f'{o3_name} = o3.node.Node(osi, 0.0, 1.0)')
            elif o3_name == 'mat':  # TODO: detect ndmaterial
                prior_strs.append(w4 + f'{o3_name} = o3.uniaxial_material.Elastic(osi, 1.0)')
            if o3_name in ['sec', 'sec_i', 'sec_j', 'sec_e', 'section']:
                prior_strs.append(w4 + f'{o3_name} = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)')
            elif o3_name == 'integration':
                prior_strs.append(w4 + 'sec = o3.section.Elastic2D(osi, 10.0, 1.0, 1.0)')
                prior_strs.append(w4 + 'integration = o3.beam_integration.Lobatto(osi, sec, 5)')
            pjoins.append(f'{o3_name}={o3_name}')
        elif dtype == 'list' and pms[pm].list_items_dtype == 'obj':
            if o3_name == 'ele_nodes':
                prior_strs.append(w4 + 'coords = [[0, 0], [1, 0], [1, 1], [0, 1]]')
                prior_strs.append(w4 + 'ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(len(coords))]')
            else:
                prior_strs.append(w4 + f'{o3_name} = [1, 1]')
            pjoins.append(f'{o3_name}={o3_name}')
        elif dtype == 'list' and pms[pm].list_items_dtype == 'int':
            if o3_name == 'ele_nodes':
                prior_strs.append(w4 + 'coords = [[0, 0], [1, 0], [1, 1], [0, 1]]')
                prior_strs.append(w4 + 'ele_nodes = [o3.node.Node(osi, *coords[x]) for x in range(4)]')
            else:
                prior_strs.append(w4 + f'{o3_name} = [1, 1]')
            pjoins.append(f'{o3_name}={o3_name}')
        elif dtype == 'list' and pms[pm].list_items_dtype == 'float':
            prior_strs.append(w4 + f'{o3_name} = [1.0, 1.0]')
            pjoins.append(f'{o3_name}={o3_name}')
        elif dtype == 'str':
            pjoins.append(f'{o3_name}="string"')
        elif dtype == 'bool':
            pjoins.append(f'{o3_name}=False')
        else:
            pjoins.append(f'{o3_name}=1')
    pjoint = ', '.join(pjoins)
    tpara += prior_strs
    tpara.append(w4 + f'o3.{names["low_base_name"]}.{names["op_class_name"]}({pjoint})')
    tpara.append('')
    tpara.append('')

    return tpara


def build_obj_docstring(op_class_name, base_class_name, obj_blurb):
    use_raw = 0
    joined_obj_blurb = ''.join(obj_blurb)
    if 'math:' in joined_obj_blurb:
        use_raw = 1

    para = []
    if use_raw:
        para.append(w4 + 'r"""')
    else:
        para.append(w4 + '"""')
    para.append(w4 + f'The {op_class_name} {base_class_name} Class')
    para.append(w4 + '')
    if '.. math::' in joined_obj_blurb:
        lines = []
        mline = False
        for i, line in enumerate(obj_blurb):
            line = line.replace('\t', '    ')
            if '.. math::' in line:
                para.append(force_line_char_limit(w4 + f'{"".join(lines)}', w4))  # parse all text above math.
                para.append('\n    .. math::\n')
                mline = True  # define a new line for math
                lines = []
            elif mline and w4 + ' ' == line[:5]:  # copy verbose if line is indented
                para.append(line)
            else:
                if mline:
                    para.append('\n')
                    mline = False
                lines.append(trim_leading_whitespace(line))
        if len(lines):
            para.append(force_line_char_limit(w4 + f'{"".join(lines)}', w4))  # parse all remaining text.
    else:
        lines = []
        for line in obj_blurb:
            lines.append(trim_leading_whitespace(line))

        para.append(force_line_char_limit(w4 + f'{"".join(lines)}', w4))
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


def clean_docstring_content(doc_str):
    # doc_str = doc_str.replace('\\l', '\\\\l')  # Now using raw strings
    # doc_str = doc_str.replace('\\s', '\\\\s')
    # doc_str = doc_str.replace('\\e', '\\\\e')
    # doc_str = doc_str.replace('\\d', '\\\\d')
    # doc_str = doc_str.replace('\\D', '\\\\D')
    doc_str = doc_str.replace('tag', 'object')
    doc_str = doc_str.replace('Tag', 'Object')
    doc_str = doc_str.replace('(optional)', '')  # Note this is now dealt with in the type definition
    doc_str = doc_str.replace('uniaxialmaterial', 'uniaxial_material')
    return doc_str


def build_init_method_docstring(classname, pms, pms_ordered):
    init_blurb = f'Initial method for {classname}'
    dstr = [w8 + '"""', w8 + init_blurb, '', w8 + 'Parameters', w8 + '----------',
            w8 + f'osi: o3seespy.OpenSeesInstance']
    for pm in pms_ordered:
        op_str = ''
        pdes = clean_docstring_content(pms[pm].p_description.capitalize())
        if pms[pm].default_is_expression:
            op_str = f' (default={pms[pm].default_is_expression})'
        if pms[pm].o3_default_is_none or pms[pm].default_value:
            ostr = ', optional'
        else:
            ostr = ''
        dstr.append(w8 + f'{pms[pm].o3_name}: {pms[pm].dtype}{op_str}{ostr}')
        descr = force_line_char_limit(w12 + f'{pdes}', w12)
        dstr.append(descr)
    dstr.append(w8 + '"""')
    if 'math:' in ''.join(dstr):
        dstr[0] = w8 + 'r"""'
    return dstr


def build_fn_docstring(fn_blurb, pms, pms_ordered):
    use_raw = 0
    joined_fn_blurb = ''.join(fn_blurb)
    if 'math:' in joined_fn_blurb:
        use_raw = 1

    para = []
    if use_raw:
        para.append(w4 + 'r"""')
    else:
        para.append(w4 + '"""')
    if '.. math::' in joined_fn_blurb:
        lines = []
        mline = False
        for i, line in enumerate(fn_blurb):
            line = line.replace('\t', '    ')
            if '.. math::' in line:
                para.append(force_line_char_limit(w4 + f'{"".join(lines)}', w4))  # parse all text above math.
                para.append('\n    .. math::\n')
                mline = True  # define a new line for math
                lines = []
            elif mline and w4 + ' ' == line[:5]:  # copy verbose if line is indented
                para.append(line)
            else:
                if mline:
                    para.append('\n')
                    mline = False
                lines.append(trim_leading_whitespace(line))
        if len(lines):
            para.append(force_line_char_limit(w4 + f'{"".join(lines)}', w4))  # parse all remaining text.
    else:
        lines = []
        for line in fn_blurb:
            lines.append(trim_leading_whitespace(line))

        para.append(force_line_char_limit(w4 + f'{"".join(lines)}', w4))

    para += ['', w4 + 'Parameters', w4 + '----------',
            w4 + f'osi: o3seespy.OpenSeesInstance']
    for pm in pms_ordered:
        op_str = ''
        pdes = clean_docstring_content(pms[pm].p_description.capitalize())
        if pms[pm].default_is_expression:
            op_str = f' (default={pms[pm].default_is_expression})'
        if pms[pm].o3_default_is_none or pms[pm].default_value:
            ostr = ', optional'
        else:
            ostr = ''
        para.append(w4 + f'{pms[pm].o3_name}: {pms[pm].dtype}{op_str}{ostr}')
        descr = force_line_char_limit(w12 + f'{pdes}', w12)
        para.append(descr)
    para.append(w4 + '"""')
    if use_raw:
        para[0] = w4 + 'r"""'
    return para


class Param(object):
    def __init__(self, org_name, default_value, packed=None, dtype=None):
        self.org_name = org_name
        self.o3_name = None  # input name for object (e.g. ele_nodes)
        self._o3_pm_name = None  # name passed to parameters (e.g. ele_node_tags)
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

    @property
    def o3_pm_name(self):
        if self._o3_pm_name is None:
            return self.o3_name
        return self._o3_pm_name

    @o3_pm_name.setter
    def o3_pm_name(self, name):
        self._o3_pm_name = name

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


suffixes = ['', 'Args', 'Tag', 'Tags', 'Tags', 'Tag1', 'Tag2', 'Tag3', 'Tag4', 'MatTag', 'MatTags',
            'Flag', 'Vals', 'SeriesTag', 's', 'Points', 'Code', 'PerLength']


def clean_fn_line(line, has_tag=True):
    df = pd.read_csv('markup_lists/overrides.csv')
    defaults = OrderedDict()
    base_type = line.split('.. function:: ')[-1]
    base_type = base_type.split('(')[0]
    optype_res = re.search(optype_pat, line)
    if optype_res is None:
        optype = None
    else:
        optype = optype_res.group()[1:-1]
    df_op = df[(df['base_type'] == base_type) & (df['op_type'] == optype)]

    inputs_str = line.split(')')[0]
    inputs1 = inputs_str.split(',')
    mstr = r'(?![^)(]*\([^)(]*?\)\)),(?![^\[]*\])'
    inputs = re.split(mstr, inputs_str)
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
        name_only = name_only.replace('<', '')  # TODO: detect and group
        name_only = name_only.replace('>', '')
        if '*' in name_only:
            name_only = name_only.replace('*', '')
            # name_only = name_only.replace('Args', '')
        name_only = name_only.split('=')[0]
        names_only.append(name_only)
    for j, inpy in enumerate(inputs):
        inpy = inpy.replace(' ', '')
        inpy = inpy.replace('[', '')
        inpy = inpy.replace(']', '')
        inpy = inpy.replace('<', '')  # TODO: detect and group
        inpy = inpy.replace('>', '')
        if '=' in inpy:
            print(line)
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
                if word == names_only[j + 1].upper():
                    markers[names_only[j + 1]] = word
                    continue
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
                # elif 'doRayleigh' in word:
                #     flags.append(word)
                #     inp = word
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
    # if "doRayleigh" in defaults and 'rFlag' in defaults:  # e.g. ZeroLength
    #     del defaults["doRayleigh"]
    # if "doRayleigh" in defaults:  # TODO: not working
    #     defaults["doRayleigh"].is_flag = True


    return base_type, optype, defaults, op_kwargs


def trim_leading_whitespace(line):
    while line[0] in [' ', '\t']:
        line = line[1:]
        if not len(line):
            break
    return line


def parse_single_file(ffp, osi_type, expected_base_type=None, multi_def=False, is_fn=False, fname=''):
    a = open(ffp)
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
    obj_blurb = []
    sub_obj_blurbs = [[]]
    ipara = []
    if osi_type is None:
        has_tag = False
    else:
        has_tag = True
    for ss, line in enumerate(lines):
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
                df_td = pd.read_csv('markup_lists/two_definitions.csv')
                df_td = df_td[(df_td['base_type'] == base_type) & (df_td['op_type'] == optype)]
                if len(df_td):
                    two_defs = df_td['option'].iloc[0]
                if two_defs == '2Dand3D':
                    cl_name_suf = '2D'
                if len(obj_blurb):
                    cur_obj_blurb = obj_blurb + ['\n\n' + w4] + sub_obj_blurbs[0]
                else:
                    cur_obj_blurb = sub_obj_blurbs[0]
                if two_defs != 'doubleUp' or optype1 is not None:
                    pstr1, tstr1 = refine_and_build(doc_str_pms, dtypes, defaults, op_kwargs, descriptions, optype,
                                            base_type, osi_type, cl_name_suf, cur_obj_blurb, is_fn=is_fn, fname=fname)
                    pstr += pstr1
                    tstr += tstr1

                if multi_def:
                    base_type = None
                    doc_string_open = 0
                if two_defs == '2Dand3D':
                    cl_name_suf = '3D'
                    if len(obj_blurb):
                        cur_obj_blurb = obj_blurb + ['\n\n' + w4] + sub_obj_blurbs[1]
                    else:
                        cur_obj_blurb = sub_obj_blurbs[1]
                    if defaults1 is None:
                        raise ValueError(defaults, op_kwargs)
                    pstr1, tstr1 = refine_and_build(doc_str_pms, dtypes, defaults1, op_kwargs1, descriptions, optype1,
                                                    base_type1, osi_type, cl_name_suf, cur_obj_blurb, is_fn=is_fn,
                                                    fname=fname)
                    pstr += '\n' + pstr1
                    tstr += tstr1
                # Reset in case two objects in same file
                if two_defs != 'doubleUp':
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
        if '``<specific parameter args>``' in line:
            if 'specificparameterargs' in defaults:
                defaults['p_args'] = defaults['specificparameterargs']
                del defaults['specificparameterargs']
            doc_str_pms.append('p_args')
            dtypes.append('unk')
            descriptions['p_args'] = trim_leading_whitespace(line.replace('``<specific parameter args>``', ''))
            continue
        # rline = line.replace('<', '')
        # rline = rline.replace('>', '')
        rline = line
        res = re.search(pname_pat, rline)
        if first_char != '*' and res:
            pname = res.group()[2:-2]
            # if len(res.group()) > 4 and "'-" == res.group()[2:4]:
            #     continue  # op_kwarg
            ei = rline.find('|')
            dtype_res = re.search(dtype_pat, rline)
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

            res = re.findall(pname_pat, rline[:ei])
            for pm in res:
                if len(pm) > 4 and "'-" == pm[0:2]:
                    pm = pm[2:-1]
                doc_str_pms.append(pm)
                dtypes.append(dtype)
                descriptions[pm] = des
            cur_res = list(res)
        elif '.. function:: ' in line:
            fn_line_counter += 1
            sub_obj_blurbs.append([])
            if base_type is None:
                base_type, optype, defaults, op_kwargs = clean_fn_line(line, has_tag)
                base_class_name = convert_name_to_class_name(base_type)
                if expected_base_type is not None and expected_base_type != base_class_name:
                    if optype is None:
                        pass
                    else:
                        cam_case = convert_camel_to_snake(base_class_name)
                        if is_fn is False:
                            new_istr = f'from o3seespy.command.{cam_case} import {base_class_name}Base'
                            if new_istr not in ipara:
                                ipara.append(new_istr)
                if multi_def and optype is None:
                    base_type = None
                    doc_string_open = 0
            else:  # multiple function definitions
                glob_list.append(f'{base_type}-{optype}')

                base_type1, optype1, defaults1, op_kwargs1 = clean_fn_line(line, has_tag)
                df_td = pd.read_csv('markup_lists/two_definitions.csv')
                df_td = df_td[(df_td['base_type'] == base_type) & (df_td['op_type'] == optype)]
                two_defs = df_td['option'].iloc[0]
                if not len(df_td):
                    raise ValueError(f'{base_type}-{optype}')
                assert base_type == base_type1
                #assert optype == optype1, (optype, optype1)
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
                if len(pm) > 4 and "'-" == pm[0:2]:
                    pm = pm[2:-1]
                descriptions[pm] += line
        elif len(line) > 7 and '   :noi' in line[:7]:
            # if len(line) > 7 and '   :ref' in line[:7]:
            #     pass
            # else:
            continue
        elif not doc_string_open:
            line = line.replace(':ref:', '')
            if fn_line_counter:
                sub_obj_blurbs[fn_line_counter - 1].append(line)
            elif title_marker == 2:
                obj_blurb.append(line)

    if base_type is not None:  # when there are no inputs
        cl_name_suf = ""
        if len(obj_blurb):
            cur_obj_blurb = obj_blurb + ['\n\n' + w4] + sub_obj_blurbs[fn_line_counter - 1]
        else:
            cur_obj_blurb = sub_obj_blurbs[fn_line_counter - 1]
        pstr1, tstr1 = refine_and_build(doc_str_pms, dtypes, defaults, op_kwargs, descriptions, optype, base_type,
                                        osi_type, cl_name_suf, cur_obj_blurb, is_fn=is_fn, fname=fname)
        pstr += pstr1
        tstr += tstr1
    istr = '\n'.join(ipara)
    return pstr, tstr, istr


def find_in_custom_defs(base_type, optype):
    from _auto_build import _custom_gen as cust_file

    cust_obj_list = [o[0] for o in getmembers(cust_file) if isclass(o[1])]
    cust_parent_obj_list = [str(o[1].__bases__[0].__name__) for o in getmembers(cust_file) if isclass(o[1])]
    if optype is not None:
        op_class_name = convert_name_to_class_name(optype)
        if op_class_name in cust_obj_list:
            indy = cust_obj_list.index(op_class_name)
            base_class_name = convert_name_to_class_name(base_type)
            parent_class = f'{base_class_name}Base'
            if parent_class == cust_parent_obj_list[indy]:
                source = inspect.getsource(getattr(cust_file, op_class_name))
                return '\n' + source + '\n', ""
    return None


def refine_and_build(doc_str_pms, dtypes, defaults, op_kwargs, descriptions, optype, base_type, osi_type, cl_name_suf="", obj_blurb=[],
                     is_fn=False, fname=''):
    custom_str = find_in_custom_defs(base_type, optype)
    if custom_str is not None:
        return custom_str
    if osi_type is not None:
        doc_str_pms = doc_str_pms[1:]  # remove mat tag
        dtypes = dtypes[1:]
    verbose = 0
    if verbose:
        print('doc_str: ', doc_str_pms)
        print('fn_inps: ', list(defaults))
        print('op_kwargs: ', list(op_kwargs))

    for i, pm in enumerate(doc_str_pms):
        if '-' + pm in op_kwargs:
            continue
        if pm not in defaults:
            continue
        if defaults[pm].is_flag:
            defaults[pm].dtype = 'bool'
        elif defaults[pm].org_name.endswith('Tag'):
            defaults[pm].dtype = 'obj'
        elif defaults[pm].org_name.endswith('Tags'):
            defaults[pm].list_items_dtype = 'obj'
            defaults[pm].dtype = 'list'
        elif defaults[pm].org_name[:-1].endswith('Tag'):  # e.g. Tag1, Tag2
            defaults[pm].dtype = 'obj'
        elif 'listi' in dtypes[i]:
            defaults[pm].dtype = 'list'
            defaults[pm].list_items_dtype = 'int'
        elif 'listf' in dtypes[i]:
            defaults[pm].dtype = 'list'
            defaults[pm].list_items_dtype = 'float'
        else:
            defaults[pm].dtype = dtypes[i]
        defaults[pm].p_description = descriptions[pm]
        if descriptions[pm].startswith('Optional: '):
            defaults[pm].default_value = None

    if 'eleNodes' in defaults:
        defaults['eleNodes'].list_items_dtype = 'obj'
    hidden_objs = ['iNode', 'jNode', 'kNode', 'cNode', 'lNode', 'rNode', 'sNode', 'lagr_node',
                   'crdTransf', 'sec', 'secI', 'secJ', 'secE']
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
        defaults['orient'].dtype = 'list'
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
    if "-maxIter" in op_kwargs and 'iter' in defaults:
        del op_kwargs['-maxIter']
        # defaults['iter'] = copy.deepcopy(defaults['maxIter'])
        defaults['iter'].marker = 'iter'
        if 'tol' in defaults:
            defaults['tol'].depends_on = 'iter'
    if "-jntOffset" in op_kwargs and 'dI' in defaults:
        del op_kwargs['-jntOffset']
        # defaults['iter'] = copy.deepcopy(defaults['maxIter'])
        defaults['dI'].marker = 'jntOffset'
        if 'dJ' in defaults:
            defaults['dJ'].depends_on = 'dI'
    if "-doRayleigh" in op_kwargs and 'rFlag' in defaults:
        del op_kwargs['-doRayleigh']
        # defaults['iter'] = copy.deepcopy(defaults['maxIter'])
        defaults['rFlag'].marker = 'doRayleigh'
        # del defaults['maxIter']
    #assert len(doc_str_pms) == len(defaults) + len(op_kwargs), (len(doc_str_pms), (len(defaults), len(op_kwargs)))
    # if len(op_kwargs) == 1:
    #     opk = list(op_kwargs)
    #     for item in opk:
    #         defaults[item] = Param(org_name=item, default_value=False, packed=False)
    #         defaults[item].marker = f'-{item}'
    #         del op_kwargs[item]
    if is_fn:
        pstr, tstr = constructor_of_func(fname, base_type, defaults, op_kwargs, osi_type=osi_type, cl_name_suf=cl_name_suf,
                                 obj_blurb=obj_blurb)
    else:
        pstr, tstr = constructor(base_type, optype, defaults, op_kwargs, osi_type=osi_type, cl_name_suf=cl_name_suf, obj_blurb=obj_blurb)
    return pstr, tstr
    # if line[3:5] == '``':
    #     para = line[5:]


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
            if ':' in line or '-' in line or '#' in line or line == '':
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
            print(mat)
            pstr, tstr, istr = parse_single_file(ffp, osi_type='mat')
            if istr not in ipara and istr != '':
                ipara.append(istr)
            para.append(pstr)
            tpara.append(tstr)
        print(item, collys[item])
        with open(floc + f'{item}.py', 'w') as ofile:
            ofile.write('\n'.join(ipara))
            if len(ipara):
                ofile.write('\n')
            ofile.write('\n'.join(para))
        diry = 'temp_tests/uniaxial/'
        if not os.path.exists(diry):
            os.makedirs(diry)
        with open(f'temp_tests/uniaxial/atest_{item}.py', 'w') as ofile:
            ofile.write('\n'.join(tpara))


def parse_all_ndmat():
    import user_paths as up

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
            # if mat == 'PressureDependMultiYield':
            #     continue
            # if mat == 'PressureDependMultiYield02':
            #     continue
            # if mat == 'StressDensityModel':
            #     continue
            if mat == 'PM4Sand':
                continue
            mat_name = convert_name_to_class_name(mat)

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
        diry = 'temp_tests/nd_material/'
        if not os.path.exists(diry):
            os.makedirs(diry)
        with open(f'temp_tests/nd_material/atest_{item}.py', 'w') as ofile:
            ofile.write('\n'.join(tpara))


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
            # if ele in ['trussEle', 'corotTruss', 'RJWatsonEqsBearing']:
            #     continue
            # if ele == 'zeroLengthND':
            #     continue
            # if mat == 'PressureDependMultiYield02':
            #     continue

            # open(up.OPY_DOCS_PATH + '%s.rst' % ele)
            ffp = up.OPY_DOCS_PATH + '%s.rst' % ele
            print(ele)
            pstr, tstr, istr = parse_single_file(ffp, osi_type='ele')
            if istr not in ipara and istr != '':
                ipara.append(istr)
            para.append(pstr)
            tpara.append(tstr)
        print(item, collys[item])
        with open(floc + f'{item}.py', 'w') as ofile:
            ofile.write('\n'.join(ipara))
            if len(ipara):
                ofile.write('\n')
            ofile.write('\n'.join(para))
        diry = 'temp_tests/elements/'
        if not os.path.exists(diry):
            os.makedirs(diry)
        with open(f'temp_tests/elements/atest_{item}.py', 'w') as ofile:
            ofile.write('\n'.join(tpara))


def str_of_set_parameter_base_method():
    return [w4 + 'def set_parameter(self, osi, pstr, value, ele, eles):',
        w8 + 'from o3seespy import set_parameter',
        w8 + 'if ele is not None:',
        w12 + 'set_parameter(osi, value=value, eles=[ele], args=[pstr, 1])',
        w8 + 'if eles is not None:',
        w12 + 'set_parameter(osi, value=value, eles=eles, args=[pstr, 1])', '']


def parse_generic_single_file(obj_type, osi_type, extras=None, multi_def=False, is_fn=False):
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
    if extras is not None:
        collys[obj_type] = extras

    floc = ROOT_DIR + 'o3seespy/command/'
    for item in collys:
        if is_fn is False:
            para = ['from o3seespy.base_model import OpenSeesObject', '', '']
            para += [f'class {o3_class_type}Base(OpenSeesObject):']
            para += [w4 + f'op_base_type = "{obj_type}"', '']
            if obj_type == 'section':
                para += str_of_set_parameter_base_method()
        else:
            para = []
        tpara = ['import o3seespy as o3  # for testing only', '', '']
        ipara = []
        for mat in collys[item]:

            mat_name = convert_name_to_class_name(mat)
            if mat_name in cust_obj_list:
                source = inspect.getsource(getattr(cust_file, mat))
                para.append('')
                para.append(source)
                continue

            open(up.OPY_DOCS_PATH + '%s.rst' % mat)
            ffp = up.OPY_DOCS_PATH + '%s.rst' % mat
            pstr, tstr, istr = parse_single_file(ffp, osi_type=osi_type, expected_base_type=o3_class_type, multi_def=multi_def, is_fn=is_fn, fname=obj_type)
            if istr not in ipara and istr != '':
                ipara.append(istr)
            para.append(pstr)
            tpara.append(tstr)
        name = item
        # if name == 'test':
        #     name = 'check_test'
        with open(floc + f'{name}.py', 'w') as ofile:
            ofile.write('\n'.join(ipara))
            if len(ipara):
                ofile.write('\n')
            ofile.write('\n'.join(para))
        with open(f'temp_tests/atest_{name}.py', 'w') as ofile:
            ofile.write('\n'.join(tpara))


def a_test_clean_fn_line():
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

    all = 0
    # all = 1  # TODO: KikuchiBearing
    # TODO: dettach docstrings - if exists then don't use rst version
    # TODO: add type hinting for default None (w: str = None)
    # TODO: Add SteelZ01 and others that are not in the openseespy docs
    if not all:
        # parse_generic_single_file(obj_type='section', osi_type='sect')
        # print(ts)
        # parse_generic_single_file(obj_type='integrator', osi_type=None)
        # parse_generic_single_file(obj_type='beamIntegration', osi_type='integ')
        # ps, ts, istr = parse_single_file(up.OPY_DOCS_PATH + 'Quad.rst', 'ele')
        # parse_generic_single_file(obj_type='geomTransf', osi_type='transformation')
        # parse_generic_single_file(obj_type='beamIntegration', osi_type='integ')
        # print(ts)
        pstr, tstr, istr = parse_single_file(up.OPY_DOCS_PATH + 'Joint2D.rst', osi_type='ele')
        print(pstr)
        # parse_all_ndmat()
        # parse_generic_single_file(obj_type='senscmds', osi_type='senscmds')
        # p = parse_single_file(up.OPY_DOCS_PATH + 'SFI_MVLEM.rst', osi_type='ele')
        # print(p[0])
        #
        # p = parse_generic_single_file(obj_type='senscmds', osi_type=None, is_fn=True)
        # p = parse_generic_single_file(obj_type='parallelcmds', osi_type=None, is_fn=True)
        # p = parse_generic_single_file(obj_type='reliabilitycmds', osi_type=None, is_fn=True)
        # p = parse_generic_single_file(obj_type='fsicmds', osi_type=None, is_fn=True)


        # p = parse_single_file(up.OPY_DOCS_PATH + 'elastomericBearingBoucWen.rst', osi_type=None)
        # print(p[0])
        # parse_all_elements()
        # pstr, tstr, istr = parse_single_file(up.OPY_DOCS_PATH + 'PathTs.rst', 'tseries')
        # print(pstr)
        # a_test_clean_fn_line()
        # parse_generic_single_file(obj_type='section', osi_type='sect')
        # parse_generic_single_file(obj_type='patch', osi_type=None, extras=['patch'], multi_def=True)
        # ffp = up.OPY_DOCS_PATH + '%s.rst' % 'patch'
        # pstr, tstr, istr = parse_single_file(ffp, osi_type='mat', expected_base_type='patch', multi_def=True)
        # print(pstr, tstr, istr)
    if all:
        parse_generic_single_file(obj_type='frictionModel', osi_type='ele')
        parse_generic_single_file(obj_type='pattern', osi_type='pat')
        parse_generic_single_file(obj_type='timeSeries', osi_type='tseries')
        parse_generic_single_file(obj_type='constraints', osi_type=None)
        parse_generic_single_file(obj_type='integrator', osi_type=None)
        parse_generic_single_file(obj_type='test', osi_type=None)
        parse_generic_single_file(obj_type='beamIntegration', osi_type='integ')
        parse_generic_single_file(obj_type='section', osi_type='sect')  # TODO: deal with duplicate in Fiber
        parse_generic_single_file(obj_type='geomTransf', osi_type='transformation')
        parse_generic_single_file(obj_type='patch', osi_type=None, extras=['patch'], multi_def=True)
        parse_generic_single_file(obj_type='layer', osi_type=None, extras=['layer'], multi_def=True)
        parse_generic_single_file(obj_type='system', osi_type=None)
        parse_generic_single_file(obj_type='mesh', osi_type='mesh')
        parse_generic_single_file(obj_type='senscmds', osi_type='senscmds', is_fn=True)
        # parse_generic_single_file(obj_type='system', osi_type=None)
        # TODO: deal with senscmds.rst
        parse_all_uniaxial_mat()
        parse_all_ndmat()
        parse_all_elements()

        # load old versions
        a = open('req_objs_that_are_deprecated/required_system.txt', 'r')
        b = a.read().splitlines()
        a.close()
        ffp = ROOT_DIR + 'o3seespy/command/system.py'
        a = open(ffp, 'a')
        a.write('\n'.join(b))
        a.close()

    ofile = open('temp.out', 'w')
    ofile.write('\n'.join(glob_list))
    # ps, ts, istr = parse_single_file(up.OPY_DOCS_PATH + 'block2D.rst', 'block')
    #ps, ts, istr = parse_single_file(up.OPY_DOCS_PATH + 'region.rst', 'block')
    ##print(ts)
    # defo = 'a2*k'
    # if any(re.findall('|'.join(['\*', '\/', '\+', '\-', '\^']), defo)):
    #     print('found')
    # TODO: YamamotoBiaxialHDRorient, YamamotoBiaxialHDRcoRS


