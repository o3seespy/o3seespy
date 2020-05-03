import glob
import user_paths as up
import os
import pandas as pd


def build():
    # folder name, map name, interpreter file name
    otypes = [
        ('nD', 'nDMaterials', 'OpenSeesNDMaterialCommands', 'Material',),
        ('uniaxial', 'uniaxialMaterials', 'OpenSeesUniaxialMaterialCommands', 'Material'),
        ('section', 'function', 'OpenSeesSectionCommands', 'Section'),
        # ('yieldSurface', 'OpenSeesYSCommands')  # Not yet in Python
              ]
    for j in range(len(otypes)):
        folder = otypes[j][0]
        map_name = otypes[j][1]
        obj_fname = otypes[j][2]
        ext_name = otypes[j][3]

        ifile = open(up.OPENSEES_SOURCE_PATH + f'interpreter/{obj_fname}.cpp')
        lines = ifile.read().splitlines()
        mstr = f'{map_name}Map.insert(std::make_pair("'
        objs = []
        for line in lines:
            if mstr in line:
                obj = line.split(mstr)[-1]
                obj = obj.split('"')[0]
                objs.append(obj)
        print(objs)
        nd_dir = up.OPENSEES_SOURCE_PATH + f'material/{folder}/'
        subdirs = [x[0] for x in os.walk(nd_dir)]
        ffps = []
        for sdir in subdirs:
            ffps += glob.glob(sdir + '/*.cpp')
        fnames = []
        for ffp in ffps:
            print(ffp)
            fnames.append(ffp.split('/')[-1])
        matched_ffps = []
        exts = []
        two_d = []
        three_d = []
        for obj in objs:
            if f'{obj}.cpp' in fnames:
                print('found: ', obj)
                ind = fnames.index(f'{obj}.cpp')
                matched_ffps.append(ffps[ind])
                exts.append('-')
            elif f'{obj}{ext_name}.cpp' in fnames:
                print('found: ', obj)
                ind = fnames.index(f'{obj}{ext_name}.cpp')
                matched_ffps.append(ffps[ind])
                exts.append(ext_name)
            elif f'{obj}{ext_name}2d.cpp' in fnames:
                print('found: ', obj)
                ind = fnames.index(f'{obj}{ext_name}2d.cpp')
                matched_ffps.append(ffps[ind])
                exts.append(f'{ext_name}2d')
            elif f'{obj}{ext_name}3d.cpp' in fnames:
                print('found: ', obj)
                ind = fnames.index(f'{obj}{ext_name}3d.cpp')
                matched_ffps.append(ffps[ind])
                exts.append(f'{ext_name}3d')
            else:
                print('NOT FOUND: ', obj)
                matched_ffps.append('-')
                exts.append('')
        for i, ffp in enumerate(matched_ffps):
            if ffp != '-':
                matched_ffps[i] = ffp.replace(up.OPENSEES_SOURCE_PATH, '')

        od = {'obj': objs,
              'ffp': matched_ffps,
              'ext': exts}
        df = pd.DataFrame.from_dict(od)
        df.to_csv(f'markup_lists/{folder}-ffps.csv', index=False)


if __name__ == '__main__':
    build()
