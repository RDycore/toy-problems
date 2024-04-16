import sys, getopt, os, subprocess
my_env = os.environ.copy()
PETSC_DIR = my_env['PETSC_DIR']
sys.path.append(PETSC_DIR + '/lib/petsc/bin')
import PetscBinaryIO
import numpy as np
from matplotlib import pyplot as plt

def main(argv):

	ex = 'ex1'
	Nx = 200
	Ny = 200

	opts, args = getopt.getopt(argv,"h:e:x:y:",["ifile=","ofile="])
	for opt, arg in opts:
		if opt == '-h':
			print ('show_results.py -ex <example number> -x <Number of grid in x> -y <Number of grid in y>')
			sys.exit()
		elif opt in ("-e","--exmaple"):
			ex = arg
		elif opt in ("-x","--Nx"):
			Nx = int(arg)
		elif opt in ("-y","--Ny"):
			Ny = int(arg)

	print('Showing results from ' + ex)

	if os.path.exists(ex):
		os.remove(ex)

	bashCommand = "make " + ex
	process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, env=my_env)
	output, error = process.communicate()
	if not error:
		print(ex + ' is compiled!')
	else:
		print(output)

	if ex == 'ex2e':
		bashCommand = "mpiexec -n 2 ./ex2b -dt 0.04 -Nt 180 -Nx " + str(Nx) + " -Ny " + str(Ny)  + " -savef"
		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, env=my_env)
		output, error = process.communicate()

		io = PetscBinaryIO.PetscBinaryIO()

		fh = './outputs/ex2b_output_dt_0.040000_final_solution.dat'
		ob = io.readObjectType(fh)
		if ob == 'Vec':
			data = io.readVec(fh)
			data = data[1:]
		data = data.reshape((int(len(data)/3),3))
		h2b    = data[:,0].reshape((Nx,Ny))

		bashCommand = "mpiexec -n 2 ./ex2e -dt 0.4 -Nt 18 -Nx " + str(Nx) + " -Ny " + str(Ny)  + " -savef"
		process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, env=my_env)
		output, error = process.communicate()

		io = PetscBinaryIO.PetscBinaryIO()

		fh = './outputs/ex2e_output_dt_0.400000_final_solution.dat'
		ob = io.readObjectType(fh)
		if ob == 'Vec':
			data = io.readVec(fh)
			data = data[1:]
		data = data.reshape((int(len(data)/3),3))
		h2e    = data[:,0].reshape((Nx,Ny))
		
		fig, (ax1,ax2,ax3) = plt.subplots(figsize=(9, 3), ncols=3)
		pos1 = ax1.imshow(h2b.transpose(),vmin=5, vmax=10, cmap='cool', aspect='auto')
		fig.colorbar(pos1, ax=ax1)
		ax1.set_title('ex2b dt = 0.04 [s]')
		pos2 = ax2.imshow(h2e.transpose(),vmin=5, vmax=10, cmap='cool', aspect='auto')
		fig.colorbar(pos2, ax=ax2)
		ax2.set_title('ex2e dt = 0.4 [s]')
		pos3 = ax3.imshow(h2e.transpose()-h2b.transpose(),vmin=-1, vmax=1, cmap='bwr', aspect='auto')
		fig.colorbar(pos3, ax=ax3)
		ax3.set_title('ex2e - ex2b')
		plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])