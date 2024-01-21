import os
import h5py


def iterate(fname, T, mQ, mG, G, G1, L, screen, suppress, init='', mode='LO', force=False):
    exists = 0

    if os.path.exists(fname):
        f = h5py.File(fname, 'r')
        if f.attrs['status'] == 'DONE':
            f.close()
            exists = 1
            # init_arg = '--init data_single_%i.hdf5'%(int(1e3*T))
            ret_code = 0
        else:
            f.close()
    
    if not exists or force:
        cmd = f'python3 -m tmqgp_iterate_single {T} {mQ} {mG} {G} {G1} {L} {screen} {suppress} --showtime --mode {mode} '
        if init == '':
            init_arg = ''
        else:
            init_arg = '--init ' + init
        cmd += init_arg

        cmd += ' --file ' + fname

        cmd += ' --save_iter'

        print('Running ' + cmd)
        ret_code = os.system(cmd)

    return ret_code

def thermo(fname):
    print(fname)
    # TODO: very crude, redo
    if os.path.exists('th_' + fname):
        return 0
    cmd_th = f'python3 -m tmqgp_thermo_single {fname}'
    print('Running ' + cmd_th)
    return os.system(cmd_th)



def iterate_mu(fname, mu, T, mQ, mG, G, G1, L, screen, suppress, init='', force=False):
    exists = 0

    if os.path.exists(fname):
        f = h5py.File(fname, 'r')
        if f.attrs['status'] == 'DONE':
            f.close()
            exists = 1
            # init_arg = '--init data_single_%i.hdf5'%(int(1e3*T))
            ret_code = 0
        else:
            f.close()
    
    if not exists or force:
        cmd = f'python3 -m tmqgp_iterate_single_mu {mu} {T} {mQ} {mG} {G} {G1} {L} {screen} {suppress} '
        if init == '':
            init_arg = ''
        else:
            init_arg = '--init ' + init
        cmd += init_arg

        cmd += ' --file ' + fname

        cmd += ' --save_iter'

        print('Running ' + cmd)
        ret_code = os.system(cmd)

    return ret_code

def thermo_mu(fname):
    print(fname)

    if os.path.exists('th_' + fname):
        return 0
    cmd_th = f'python3 -m tmqgp_thermo_single_mu {fname}'
    print('Running ' + cmd_th)
    return os.system(cmd_th)