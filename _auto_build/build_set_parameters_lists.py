# not needed - Read the file: OpenSeesNDMaterialCommands.cpp
# Search in files, material/nd/*.cpp
# First search if file name matches else look through all files
# looking for "<material-name>::getCopy
# build map file of object to .cpp file

# Then find: 'CycLiqCPSP::setParameter' and inside search for:
#   if (strcmp(argv[0],"updateMaterialStage") == 0) {
import pandas as pd
import user_paths as up
import re


def convert_camel_to_snake(name):
    s1 = re.sub('(.)([A-Z][a-z]+)', r'\1_\2', name)
    s1 = re.sub('([a-z0-9])([A-Z])', r'\1_\2', s1).lower()
    s1 = re.sub(r'(\d)\_', r'\1', s1)  # removes underscore after number (e.g. 2D)
    return s1

md = {
    'e': 'e_mod',
    'g': 'g_mod',
    'shear_modulus': 'g_mod',
    'bulk_modulus': 'bulk_mod',
    'poisson_ratio': 'nu',

}


def create():
    otypes = [
        ('nD', 'nDMaterials', 'OpenSeesNDMaterialCommands', 'Material',),
        ('uniaxial', 'uniaxialMaterials', 'OpenSeesUniaxialMaterialCommands', 'Material'),
        ('section', 'function', 'OpenSeesSectionCommands', 'Section'),
        # ('yieldSurface', 'OpenSeesYSCommands')  # Not yet in Python
    ]
    for j in range(len(otypes)):
        folder = otypes[j][0]
        df = pd.read_csv(f'markup_lists/{folder}-ffps.csv')
        od = {'obj': [],
              'o3_name': [],
              'cpp_name': []}
        for i in range(len(df)):
            ffp = df['ffp'].iloc[i]
            if ffp == '-':
                continue
            obj = df['obj'].iloc[i]
            print(up.OPENSEES_SOURCE_PATH + ffp)
            lines = open(up.OPENSEES_SOURCE_PATH + ffp, 'r', encoding="utf8", errors='ignore').read().splitlines()
            on = 0
            if df['ext'].iloc[i] == '-':
                objname = obj
            else:
                objname = f'{obj}{df["ext"].iloc[i]}'
            params = []
            for line in lines:

                if on:
                    if f'int {objname}' in line:
                        on = 0
                    elif 'strcmp(argv[0],' in line:
                        param = line.split('strcmp(argv[0],')[-1]
                        param = param.split('"')[1]
                        params.append(param)
                elif f'{objname}::setParameter' in line:
                    on = 1
                    continue
            snames = []
            for pm in params:
                print(obj, pm)
                sname = convert_camel_to_snake(pm)
                if sname in md:
                    sname = md[sname]
                else:
                    if pm[0] == 'G' and pm[1] in ['xyz']:
                        sname = 'g_mod_' + pm[1:]
                    if pm[0] == 'E' and pm[1] in ['xyz']:
                        sname = 'e_mod_' + pm[1:]
                snames.append(sname)  # TODO: deal with duplicate name
            od['obj'].append(obj)
            od['o3_name'].append('-'.join(snames))
            od['cpp_name'].append('-'.join(params))
            if not len(params):
                print('None found: ', obj)
        df = pd.DataFrame.from_dict(od)
        df.to_csv(f'markup_lists/{folder}-set-parameters.csv', index=False)



if __name__ == '__main__':
    create()
