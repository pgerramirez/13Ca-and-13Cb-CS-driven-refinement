import optparse
import datetime
import sys
from rosetta import *
init()

##### OPTIONS ######
parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',
    default = '1HKO',    # default example PDB
    help = 'the PDB file containing the protein to refine')
parser.add_option('--rmsd_pdb', dest = 'rmsd_pdb',
    default = '1CYO.clean',
    help = 'the pdb file containing the structure to compare rmsd with')
parser.add_option("--iterations", dest="iterations",
	default = "6",
	help = "the number of final structures that will be obtained")
(options,args) = parser.parse_args()

pdb_filename = str(options.pdb_filename)
rmsd_pdb = str(options.rmsd_pdb)
iterations = int(options.iterations)

###### Functions ######

def start(pdb_filename):
	pose = Pose()
	pdb_file = "%s.pdb" % pdb_filename # May have to change this for multi-conformation .pdb files
	pose_from_pdb(pose, pdb_file)
	pymover = PyMOL_Mover()
	
	scorefxn= create_score_function('talaris2013')
	#scorefxn.set_weight(rg , 1)
	scorefxn(pose)
	pymover.update_energy(True)
	pymover.send_energy(pose)
	pymover.keep_history(True)

	aarmsd = all_atom_rmsd(pose, rmsd_pose)
	carmsd = native_CA_rmsd(pose, rmsd_pose)
	gdt = CA_gdtmm(pose, rmsd_pose)

	fe.write("%*s%*s scorefxn:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "START", 4, "0", 14, str(scorefxn(pose)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))

	pose.pdb_info().name('refine_start')
	pymover.apply(pose)
	pymover.send_energy(pose)

	pose.dump_pdb("refine_start.pdb")
	pdb_file = "refine_start.pdb"
	return pdb_file


def minimize(pdb_file, counter):
	pose = Pose()
	pose_from_pdb(pose, pdb_file)
	pymover = PyMOL_Mover()
	
	scorefxn= create_score_function('talaris2013')
	#scorefxn.set_weight(rg , 1)
	pymover.apply(pose)
	scorefxn(pose)
	pymover.update_energy(True)
	pymover.send_energy(pose)
	pymover.keep_history(True)



	minmover = MinMover()
	movemap = MoveMap()
	minmover.movemap(movemap)
	minmover.score_function(scorefxn)
	minmover.apply(pose)

	aarmsd = all_atom_rmsd(pose, rmsd_pose)
	carmsd = native_CA_rmsd(pose, rmsd_pose)
	gdt = CA_gdtmm(pose, rmsd_pose)

	fe.write("%*s%*s scorefxn:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "minimize", 4, counter, 14, str(scorefxn(pose)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))

	pose.dump_pdb("refine_minimize%s.pdb" % counter)
	pdb_file = "refine_minimize%s.pdb" % counter
	return pdb_file

def fast_relaxation(pdb_file, counter):

	pose = Pose()
	pose_from_pdb(pose, pdb_file)
	pymover = PyMOL_Mover()
	
	scorefxn= create_score_function('talaris2013')
	scorefxn.set_weight(rg , 1)
	pymover.apply(pose)
	scorefxn(pose)
	pymover.update_energy(True)
	pymover.send_energy(pose)
	pymover.keep_history(True)


	movemap = MoveMap()
	movemap.set_bb(True)
	movemap.set_chi(True)
	

	relax = FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.set_movemap(movemap)

	relax.apply(pose)

	aarmsd = all_atom_rmsd(pose, rmsd_pose)
	carmsd = native_CA_rmsd(pose, rmsd_pose)
	gdt = CA_gdtmm(pose, rmsd_pose)
	fe.write("%*s%*s scorefxn:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "fast_relaxation", 4, counter, 14, str(scorefxn(pose)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))
	pose.dump_pdb("refine_relax%s.pdb" % counter)
	pdb_file = "refine_relax%s.pdb" % counter


	return pdb_file

###### START #####


rmsd_file = "%s.pdb" % rmsd_pdb
rmsd_pose = Pose()
pose_from_pdb(rmsd_pose, rmsd_file)
rmsd_pose.dump_pdb("%s_TRANSFORMED.pdb" % rmsd_pdb)


fe = open('rosetta.log', 'w')
now = datetime.datetime.now()
fe.write("%s --pdb_filename %s --rmsd_pdb %s --iterations %s\n" % (str(now), pdb_filename, rmsd_pdb, iterations))

pdb_file = start(pdb_filename)

counter = 0

for i in range(0, iterations):

	pdb_file = minimize(pdb_file, counter)

	pdb_file = fast_relaxation(pdb_file, counter)

	pdb_file = "refine_start.pdb" #To start all over again

	counter += 1

now = datetime.datetime.now()
fe.write(str(now))
fe.close()
sys.exit("**Number of iterations completed, protocol finished.**")


