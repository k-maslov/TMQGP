from mpi4py import MPI
import numpy as np
import QuarkTM
import TMQGP as tm

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

sendbuf = None
if rank == 0:
    sendbuf = np.array([([i for i in range(j*10, (j+1)*10)]) for j in range(size)], dtype='i')

recvbuf = np.empty(10, dtype='i')
comm.Scatter(sendbuf, recvbuf, root=0)

print('Rank %i, data:'%rank, recvbuf)

out = np.array(recvbuf) * 2

print('Rank %i, out:'%rank, out)

out = np.array(recvbuf) * 2
res_buf = None
if rank == 0:
    res_buf = np.empty([size, 10], dtype='i')

comm.Gather(out, res_buf, root=0)

if rank == 0:
    print(res_buf)