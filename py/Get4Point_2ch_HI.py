import h5py
import argparse
import os

parser = argparse.ArgumentParser()
parser.add_argument('--init', type=str, default='', help='Reference mu=0 (Q and A) data folder')
args = parser.parse_args()

if args.init == '':
    init = None
else:
    init = args.init

mu_scales = [0.0, 0.2, 0.4, 0.6]
work_dir = os.getcwd()
for mu in mu_scales:
    folder = os.path.join(work_dir, '%.2f'%mu)
    if not os.path.exists(folder):
        os.mkdir(folder)
    os.chdir(folder)
    cmd = f'python3 -m GetMuFromFit_2ch {mu}'
    print('Running ' + cmd)
    os.system(cmd)
    # print(os.getcwd())
    os.chdir(work_dir)
    