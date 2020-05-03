import glob
import user_paths as up
import os
import pandas as pd


def build():
    ifile = open(up.OPENSEES_SOURCE_PATH + 'interpreter/OpenSeesNDMaterialCommands.cpp')
    lines = ifile.read().splitlines()
    mstr = 'nDMaterialsMap.insert(std::make_pair("'
    objs = []
    for line in lines:
        if mstr in line:
            obj = line.split(mstr)[-1]
            obj = obj.split('"')[0]
            objs.append(obj)
    print(objs)
    nd_dir = up.OPENSEES_SOURCE_PATH + 'material/nD/'
    subdirs = [x[0] for x in os.walk(nd_dir)]
    ffps = []
    for sdir in subdirs:
        ffps += glob.glob(sdir + '/*.cpp')
    fnames = []
    for ffp in ffps:
        print(ffp)
        fnames.append(ffp.split('/')[-1])
    matched_ffps = []
    mats = []
    for obj in objs:
        if f'{obj}.cpp' in fnames:
            print('found: ', obj)
            ind = fnames.index(f'{obj}.cpp')
            matched_ffps.append(ffps[ind])
            mats.append(0)
        elif f'{obj}Material.cpp' in fnames:
            print('found: ', obj)
            ind = fnames.index(f'{obj}Material.cpp')
            matched_ffps.append(ffps[ind])
            mats.append(1)
        else:
            print('NOT FOUND: ', obj)
            matched_ffps.append('-')
            mats.append(0)
    for i, ffp in enumerate(matched_ffps):
        if ffp != '-':
            matched_ffps[i] = ffp.replace(up.OPENSEES_SOURCE_PATH, '')

    od = {'obj': objs,
          'ffp': matched_ffps,
          'mat': mats}
    df = pd.DataFrame.from_dict(od)
    df.to_csv('nd-material-ffps.csv', index=False)


if __name__ == '__main__':
    build()
