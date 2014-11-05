import optparse
import datetime
import random
import sys
#import numpy
from itertools import groupby
from operator import itemgetter
from rosetta import *
from rosetta.protocols.loops.loop_mover.refine import *
from rosetta.protocols.loops.loop_closure.ccd import *
import cheshift as cs
init()
#init(extra_options = "-constant_seed")
#random.seed(100)


path_db = 'CS_DB'
path_triple = 'TRIPLE_5'
tolerance = 1

def get_groups(walkers, max_lenght_loop, min_length_loop):
    group_walkers = []
    grouped = []
    for k, g in groupby(enumerate(walkers), lambda (i,x):i-x):
        grouped.append(map(itemgetter(1), g))
    
    #print grouped
    for i in range(0, len(grouped)):
        for j in range(0, len(grouped[i]), max_lenght_loop):
        	if len(grouped[i][j:j+max_lenght_loop]) >= min_length_loop:
        		group_walkers.append(grouped[i][j:j+max_lenght_loop])
    return group_walkers

# def get_groups(walkers, max_lenght_loop):
#     group_walkers = []
#     grouped = []
#     for k, g in groupby(enumerate(walkers), lambda (i,x):i-x):
#         grouped.append(map(itemgetter(1), g))
    
#     for i in range(0, len(grouped)):
#         for j in range(0, len(grouped[i]), max_lenght_loop):
#             group_walkers.append(grouped[i][j:j+max_lenght_loop])
#     return group_walkers

def GiveWalkers(pdb_file):
    reference = 0.
    walkers = cs.CheShift(pdb_file, cs_exp, reference, Db, path_triple, tolerance)
    return walkers

def minimize(pdb_file, counter):
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

	aarmsd = all_atom_rmsd(pose, rmsd_pose)
	carmsd = native_CA_rmsd(pose, rmsd_pose)
	gdt = CA_gdtmm(pose, rmsd_pose)

	movemap = MoveMap()
	walkers = GiveWalkers(pdb_file)
	fd.write("**WALKERS MIN %s : %s\n" % (counter, walkers))
	walkerlist.append(len(walkers))
	fe.write("%*s%*s scorefxn:%*s flaws:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "minimize", 4, counter, 14, str(scorefxn(pose)), 3, str(len(walkers)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))

	minmover = MinMover()
	minmover.score_function(scorefxn)
	minmover.movemap(movemap)
	minmover.apply(pose)

	pose.pdb_info().name('refine_minimize_%s' % counter)
	#pymover.apply(pose)
	pymover.send_energy(pose)

	pose.dump_pdb("refine_minimize%s.pdb" % counter)
	pdb_file = "refine_minimize%s.pdb" % counter
	return pdb_file

def refine_chi(pdb_file, counter,
        kT = 1.0, cycles = 9,
        jobs = 1, job_output = 'refine_chi'):
	pose = Pose()
	pose_from_pdb(pose, pdb_file)

	starting_pose = Pose()
	starting_pose.assign(pose)

	scorefxn = create_score_function('talaris2013') # scorefxn = get_fa_scorefxn() #
	scorefxn.set_weight(rg , 1)

	aarmsd = all_atom_rmsd(pose, rmsd_pose)
	carmsd = native_CA_rmsd(pose, rmsd_pose)
	gdt = CA_gdtmm(pose, rmsd_pose)

	movemap = MoveMap()
	movemap.set_bb(False)
	movemap.set_chi(False)
	walkers = GiveWalkers(pdb_file)
	fd.write("**WALKERS CHI %s : %s\n" % (counter, walkers))
	walkerlist.append(len(walkers))
	fe.write("%*s%*s scorefxn:%*s flaws:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "refine_chi", 4, counter, 14, str(scorefxn(pose)), 3, str(len(walkers)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))
	

	for residue in walkers:
		movemap.set_chi(residue, True)

	minmover = MinMover()
	minmover.movemap(movemap)
	minmover.score_function(scorefxn)

	to_pack = standard_packer_task(starting_pose)
	to_pack.restrict_to_repacking()    # prevents design, packing only
	length = pose.total_residue()
	vector = rosetta.utility.vector1_bool()
	for i in range (1, length+1):
		if i in walkers:
			vector.append(True)
		else:
			vector.append(False)
	to_pack.restrict_to_residues(vector)

	# to_pack.temporarily_fix_everything()
	# for residue in walkers:
	# 	to_pack.temporarily_set_pack_residue(residue, True)
	
	to_pack.or_include_current(True)    # considers the original sidechains
	packmover = PackRotamersMover(scorefxn, to_pack)


	scorefxn(pose)
	pymover = PyMOL_Mover()
	pymover.update_energy(True)
	pymover.apply(pose)
	pymover.send_energy(pose)
	pymover.keep_history(True)

	combined_mover = SequenceMover()
	combined_mover.add_mover(packmover)
	combined_mover.add_mover(minmover)
	combined_mover.add_mover(pymover)

	mc = MonteCarlo(pose, scorefxn, kT)    # must reset for each trajectory!

	trial = TrialMover(combined_mover, mc)

	chi_refinement = RepeatMover(trial, cycles)

	jd = PyJobDistributor(job_output, jobs, scorefxn)
	jd.native_pose = starting_pose

	scores = [0]*(jobs + 1)
	scores[0] = scorefxn(starting_pose)

	cycle_counter = 0    # for exporting to PyMOL
	while not jd.job_complete:
		# -reload the starting pose
		pose.assign(starting_pose)
		cycle_counter += 1
		pose.pdb_info().name(job_output + '_' + str(cycle_counter))
		# -reset the MonteCarlo object (sets lowest_score to that of p)
		mc.reset(pose)

		chi_refinement.apply(pose)

		mc.recover_low(pose)
		jd.output_decoy(pose)

		pose.pdb_info().name( job_output + '_' + str(cycle_counter) + '_final')
		pymover.apply(pose)
		pymover.send_energy(pose)    # see the total score in color

		# -store the final score for this trajectory
		scores[cycle_counter] = scorefxn(pose)
		#print scorefxn.show(pose)
	# Final print
	print 'Original Score\t:\t' , scores[0]
	for i in range(1, len(scores)):    # print out the job scores
		print job_output + '_' + str(i) + '\t:\t', scores[i]
	pdb_file = "%s_1.pdb" % (job_output)
	return pdb_file


def refine_bb(pdb_file, counter,
        kT = 1.0, smallmoves = 3, shearmoves = 5,
        backbone_angle_max = 7, cycles = 9,
        jobs = 1, job_output = 'refine_bb'):

	pose = Pose()
	pose_from_pdb(pose, pdb_file)

	starting_pose = Pose()
	starting_pose.assign(pose)

	#scorefxn_low = create_score_function('cen_std')
	scorefxn= create_score_function('talaris2013') # scorefxn = get_fa_scorefxn() #
	scorefxn.set_weight(rg , 1)

	aarmsd = all_atom_rmsd(pose, rmsd_pose)
	carmsd = native_CA_rmsd(pose, rmsd_pose)
	gdt = CA_gdtmm(pose, rmsd_pose)

	movemap = MoveMap()
	movemap.set_bb(False)
	movemap.set_chi(False)
	walkers = GiveWalkers(pdb_file)
	fd.write("**WALKERS -BB %s : %s\n" % (counter, walkers))
	walkerlist.append(len(walkers))
	fe.write("%*s%*s scorefxn:%*s flaws:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "refine_bb", 4, counter, 14, str(scorefxn(pose)), 3, str(len(walkers)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))

	for residue in walkers:
		movemap.set_bb(residue, True)
	for residue in walkers:
		movemap.set_chi(residue, True)

	smallmover = SmallMover(movemap, kT, smallmoves)
	smallmover.angle_max(backbone_angle_max)
	shearmover = ShearMover(movemap, kT, shearmoves)
	shearmover.angle_max(backbone_angle_max)

	minmover = MinMover()
	minmover.movemap(movemap)
	minmover.score_function(scorefxn)

	to_pack = standard_packer_task(starting_pose)
	to_pack.restrict_to_repacking()    # prevents design, packing only

	length = pose.total_residue()
	vector = rosetta.utility.vector1_bool()
	for i in range (1, length+1):
		if i in walkers:
			vector.append(True)
		else:
			vector.append(False)
	to_pack.restrict_to_residues(vector)

	# to_pack.temporarily_fix_everything()
	# for residue in walkers:
	#  	to_pack.temporarily_set_pack_residue(residue, True)

	to_pack.or_include_current(True)    # considers the original sidechains
	packmover = PackRotamersMover(scorefxn, to_pack)

	scorefxn(pose)
	pymover = PyMOL_Mover()
	pymover.update_energy(True)
	pymover.apply(pose)
	pymover.send_energy(pose)
	pymover.keep_history(True)

	combined_mover = SequenceMover()
	combined_mover.add_mover(smallmover)
	combined_mover.add_mover(shearmover)
	combined_mover.add_mover(minmover)
	combined_mover.add_mover(packmover)
	combined_mover.add_mover(pymover)

	mc = MonteCarlo(pose, scorefxn, kT)    # must reset for each trajectory!
	trial = TrialMover(combined_mover, mc)

	bb_refinement = RepeatMover(trial, cycles)

	jd = PyJobDistributor(job_output, jobs, scorefxn)
	jd.native_pose = starting_pose

	scores = [0]*(jobs + 1)
	scores[0] = scorefxn(starting_pose)

	cycle_counter = 0    # for exporting to PyMOL
	while not jd.job_complete:
		# -reload the starting pose
		pose.assign(starting_pose)
		cycle_counter += 1
		pose.pdb_info().name(job_output + '_' + str(cycle_counter))
		# -reset the MonteCarlo object (sets lowest_score to that of p)
		mc.reset(pose)

		bb_refinement.apply(pose)
		
		mc.recover_low(pose)
		jd.output_decoy(pose)

		pose.pdb_info().name( job_output + '_' + str(cycle_counter) + '_final')
		pymover.apply(pose)
		pymover.send_energy(pose)    # see the total score in color

		# -store the final score for this trajectory
		scores[cycle_counter] = scorefxn(pose)
		#print scorefxn.show(pose)
	# Final print
	print 'Original Score\t:\t' , scores[0]
	for i in range(1, len(scores)):    # print out the job scores
		print job_output + '_' + str(i) + '\t:\t', scores[i]
	pdb_file = "%s_1.pdb" % job_output
	return pdb_file

def low_res(pdb_file, counter,
        kT = 1.0, smallmoves = 3, shearmoves = 5,
        backbone_angle_max = 7, cycles = 9):

	pose = Pose()
	pose_from_pdb(pose, pdb_file)
	starting_pose = Pose()
	starting_pose.assign(pose)
	scorefxn_high = create_score_function('talaris2013')
	scorefxn_low = create_score_function('cen_std')
	scorefxn_low.apply_patch_from_file('score4L')
	#scorefxn_low.set_weight(rg , 50)
	to_centroid = SwitchResidueTypeSetMover('centroid')
	to_fullatom = SwitchResidueTypeSetMover('fa_standard')
	walkers = GiveWalkers(pdb_file)
	fd.write("**WALKERS LOW %s : %s\n" % (counter, walkers))
	walkerlist.append(len(walkers))

	aarmsd = all_atom_rmsd(pose, rmsd_pose)
	carmsd = native_CA_rmsd(pose, rmsd_pose)
	gdt = CA_gdtmm(pose, rmsd_pose)

	fe.write("%*s%*s scorefxn:%*s flaws:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "low_res", 4, counter, 14, str(scorefxn(pose)), 3, str(len(walkers)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))

	movemap = MoveMap()
	movemap.set_bb(False)
	movemap.set_chi(False)
	for residue in walkers:
		movemap.set_bb(residue, True)
	for residue in walkers:
		movemap.set_chi(residue, True)

	to_centroid.apply(pose)

	pymover = PyMOL_Mover()
	scorefxn_low(pose)
	#pymover.update_energy(True)
	pymover.apply(pose)
	#pymover.send_energy(pose)
	#pymover.keep_history(True)

	smallmover = SmallMover(movemap, kT, smallmoves)
	smallmover.angle_max(backbone_angle_max)
	shearmover = ShearMover(movemap, kT, shearmoves)
	shearmover.angle_max(backbone_angle_max)

	minmover = MinMover()
	minmover.movemap(movemap)
	minmover.score_function(scorefxn_low)

	combined_mover = SequenceMover()
	combined_mover.add_mover(smallmover)
	combined_mover.add_mover(shearmover)
	combined_mover.add_mover(minmover)

	mc = MonteCarlo(pose, scorefxn_low, kT)    # must reset for each trajectory!
	trial = TrialMover(combined_mover, mc)

	bb_refinement = RepeatMover(trial, cycles)

	bb_refinement.apply(pose)
	mc.recover_low(pose)

	to_fullatom.apply(pose)
	recover_sidechains = ReturnSidechainMover(starting_pose)
	recover_sidechains.apply(pose)


	pose.pdb_info().name( 'lowres_' + str(counter))
	pymover.apply(pose)

	pose.dump_pdb("refine_lowres%s.pdb" % counter)
	pdb_file = "refine_lowres%s.pdb" % counter
	return pdb_file


def loop_modeler(pdb_file, counter):
	pose = Pose()
	pose_from_pdb(pose, pdb_file)
	starting_pose = Pose()
	starting_pose.assign(pose)
	scorefxn_high = create_score_function('talaris2013')
	scorefxn_low = create_score_function('cen_std') # scorefxn = get_fa_scorefxn() #
	scorefxn_low.apply_patch_from_file('score4L')
	#scorefxn_low.set_weight(chainbreak, 1)
	#scorefxn_low.set_weight(rg , 50)
	to_centroid = SwitchResidueTypeSetMover('centroid')
	to_fullatom = SwitchResidueTypeSetMover('fa_standard')

	pymov = PyMOL_Mover()
	scorefxn_high(starting_pose)    # for exporting the scores
	pymov.apply(starting_pose)
	pymov.send_energy(starting_pose)

	aarmsd = all_atom_rmsd(pose, rmsd_pose)
	carmsd = native_CA_rmsd(pose, rmsd_pose)
	gdt = CA_gdtmm(pose, rmsd_pose)

	walkers = GiveWalkers(pdb_file)
	fd.write("**WALKERS LOO %s : %s\n" % (counter, walkers))
	walkerlist.append(len(walkers))
	fe.write("%*s%*s scorefxn:%*s flaws:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "loop_modeler", 4, counter, 14, str(scorefxn(pose)), 3, str(len(walkers)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))

	group_walkers = get_groups(walkers, 5, 3)
	print group_walkers


	cycle_counter = 0
	for i in range(0, len(group_walkers)):
		looplist = group_walkers.pop(random.randrange(0, len(group_walkers)))
		print looplist

		loop_begin = looplist[0]
		loop_end = looplist[-1]
		loop_cutpoint = (loop_begin + loop_end) / 2
		my_loop = Loop(loop_begin, loop_end, loop_cutpoint)
		set_single_loop_fold_tree(pose, my_loop)
		add_cutpoint_variants(pose)

		loops = Loops()
		loops.add_loop(my_loop)

		movemap = MoveMap()
		movemap.set_bb(False)
		movemap.set_chi(False)
		movemap.set_bb_true_range(loop_begin, loop_end)
		movemap.set_chi_true_range(loop_begin, loop_end)

		smallmover = SmallMover(movemap, kT, smallmoves)
		smallmover.angle_max(backbone_angle_max)
		shearmover = ShearMover(movemap, kT, shearmoves)
		shearmover.angle_max(backbone_angle_max)

		minmover = MinMover()
		minmover.movemap(movemap)
		minmover.score_function(scorefxn_low)

		to_pack = standard_packer_task(starting_pose)
		to_pack.restrict_to_repacking()    # prevents design, packing only

		length = pose.total_residue()
		vector = rosetta.utility.vector1_bool()
		for i in range (1, length+1):
			if i in looplist:
				vector.append(True)
			else:
				vector.append(False)
		to_pack.restrict_to_residues(vector)

		to_pack.or_include_current(True)    # considers the original sidechains
		packmover = PackRotamersMover(scorefxn_low, to_pack)

		combined_mover = SequenceMover()
		combined_mover.add_mover(smallmover)
		combined_mover.add_mover(shearmover)
		combined_mover.add_mover(minmover)
		#combined_mover.add_mover(packmover)

		to_centroid.apply(pose)
		mc = MonteCarlo(pose, scorefxn_low, kT)    # must reset for each trajectory!
		trial = TrialMover(combined_mover, mc)
		bb_refinement = RepeatMover(combined_mover, 9)
		print "STARTING"
		print "BB_REFINEMENT"
		bb_refinement.apply(pose)
		print "DONE"
		mc.recover_low(pose)
		loop_refine = LoopMover_Refine_CCD(loops)
		#loop_refine.apply(pose)
		#print "loop refine done"
		ccd_closure = CcdLoopClosureMover(my_loop, movemap)
		ccd_closure.apply(pose)
		print "CCD CLOSURE DONE"
		to_fullatom.apply(pose)
		recover_sidechains = ReturnSidechainMover(starting_pose)
		recover_sidechains.apply(pose)
		pose.pdb_info().name( 'looprefine_' + str(counter) + '_' + str(cycle_counter))
		pymov.apply(pose)


		cycle_counter += 1

	pose.pdb_info().name( 'looprefine_' + str(counter) + '_final')
	pymov.apply(pose)
	scorefxn_high(pose)
	pymov.send_energy(pose)
	pose.dump_pdb("refine_loop%s.pdb" % counter)
	pdb_file = "refine_loop%s.pdb" % counter
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

	aarmsd = all_atom_rmsd(pose, rmsd_pose)
	carmsd = native_CA_rmsd(pose, rmsd_pose)
	gdt = CA_gdtmm(pose, rmsd_pose)

	movemap = MoveMap()
	movemap.set_bb(False)
	movemap.set_chi(False)
	walkers = GiveWalkers(pdb_file)
	fd.write("**WALKERS REL %s : %s\n" % (counter, walkers))
	walkerlist.append(len(walkers))
	fe.write("%*s%*s scorefxn:%*s flaws:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "fast_relaxation", 4, counter, 14, str(scorefxn(pose)), 3, str(len(walkers)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))
	for residue in walkers:
		movemap.set_bb(residue, True)
	for residue in walkers:
		movemap.set_chi(residue, True)

	######## Esto es para permitir el relajamiento de la estructura
	#if counter % 5 == 0:
		#movemap.set_bb(True)
		#movemap.set_chi(True)
	######## Probablemente es mala idea

	relax = FastRelax()
	relax.set_scorefxn(scorefxn)
	relax.set_movemap(movemap)

	relax.apply(pose)

	# pose.pdb_info().name( 'relax_' + str(counter) + '_final')
	# pymover.apply(pose)
	# pymover.send_energy(pose)

	pose.dump_pdb("refine_relax%s.pdb" % counter)
	pdb_file = "refine_relax%s.pdb" % counter
	return pdb_file

##### OPTIONS ######
parser = optparse.OptionParser()
parser.add_option('--pdb_filename', dest = 'pdb_filename',
    default = '1HKO',    # default example PDB
    help = 'the PDB file containing the protein to refine')
parser.add_option('--cs_exp', dest = 'cs_exp',
    default = '4803',
    help = 'the experimental chemical shifts files for the structure to refine')
parser.add_option('--rmsd_pdb', dest = 'rmsd_pdb',
    default = '1CYO.clean2',
    help = 'the pdb file containing the structure to compare rmsd with')
# custom refinement options
parser.add_option('--kT', dest='kT',
    default = '1.0',
    help = 'the \"temperature\" of the sample refinement protocol')
parser.add_option( '--smallmoves', dest='smallmoves',
    default = '3',
    help = 'the number of times SmallMover is applies in\
        the custom refinement protocol' )
parser.add_option('--shearmoves', dest='shearmoves',
    default = '5',
    help = 'the number of times ShearMover is applied in\
        the custom refinement protocol' )
parser.add_option( '--backbone_angle_max', dest='backbone_angle_max',
    default = '7',
    help = 'the maximum angle perturbation of SmallMover and ShearMover in\
        the custom refinement protocol')
parser.add_option("--iterations", dest="iterations",
	default = "6",
	help = "the number of times the structure will cycle between cheshift\
		and pyrosetta")
parser.add_option("--finals", dest="finals",
	default = "1",
	help = "the number of times the final structures that will be\
		obtained")
parser.add_option('--cycles', dest='cycles',
    default = '9',
    help = 'the number of refinement rounds (small, shear, min, pack) in\
        the sample refinement protocol')
# PyJobDistributor options
parser.add_option('--jobs', dest='jobs',
    default = '1',    # default to single trajectory for speed
    help = 'the number of jobs (trajectories) to perform')
# parser.add_option('--job_output', dest = 'job_output',
#     default = 'refine',    # if a specific output name is desired
#     help = 'the name preceding all output, output PDB files and .fasc')

(options,args) = parser.parse_args()

# PDB file option
cs_exp = str(options.cs_exp)
pdb_filename = str(options.pdb_filename)
rmsd_pdb = str(options.rmsd_pdb)
# custom refinement options
kT = float(options.kT)
smallmoves = int(options.smallmoves)
shearmoves = int(options.shearmoves)
backbone_angle_max = int(options.backbone_angle_max)
cycles = int(options.cycles)
iterations = int(options.iterations)
finals = int(options.finals)
# JobDistributor options
jobs = int(options.jobs)
#job_output = options.job_output

reference, outliers, Db, ok = cs.setup_CheShift(cs_exp, pdb_filename, path_db)

pdb_file = "%s_00000.pdb" % pdb_filename
rmsd_file = "%s.pdb" % rmsd_pdb
rmsd_pose = Pose()
pose_from_pdb(rmsd_pose, rmsd_file)
rmsd_pose.dump_pdb("%s_TRANSFORMED.pdb" % rmsd_pdb)

counter = 0
final_structures = 0

fd = open('rosettaflaws.log', 'w')
fe = open('rosetta.log', 'w')
now = datetime.datetime.now()
fe.write("%s --pdb_filename %s --cs_exp %s --rmsd_pdb %s --iterations %s\n" % (str(now), pdb_filename, cs_exp, rmsd_pdb,iterations))

walkerlist = []

pdb_file = minimize(pdb_file, counter)

for i in range(0, iterations):
	print walkerlist

	if i >= 4:
		if walkerlist[-1] == walkerlist[-2] == walkerlist[-3] == walkerlist[-4] == walkerlist[-5] == walkerlist[-6] == walkerlist[-7] == walkerlist[-8] == walkerlist[-9] == walkerlist[-10] == walkerlist[-11]:
			pose = Pose()
			pose_from_pdb(pose, pdb_file)
			scorefxn = create_score_function('talaris2013')
			scorefxn.set_weight(rg , 1)

			aarmsd = all_atom_rmsd(pose, rmsd_pose)
			carmsd = native_CA_rmsd(pose, rmsd_pose)
			gdt = CA_gdtmm(pose, rmsd_pose)

			walkers = GiveWalkers(pdb_file)

			pose.dump_pdb("reFINAL_%s_%s.pdb" % (pdb_filename, counter))

			pdb_file = "%s_00000.pdb" % pdb_filename

			fd.write("**WALKERS FIN %s : %s\n" % (counter, walkers))
			walkerlist.append(len(walkers))
			fe.write("%*s%*s scorefxn:%*s flaws:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "FINAL", 4, counter, 14, str(scorefxn(pose)), 3, str(len(walkers)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))
			fe.write("**Number of flaws unchanged in last 10 movements, protocol restarting.**\n")

			pdb_file = minimize(pdb_file, counter)

	#print(pdb_file)

	#pdb_file = fast_relaxation(pdb_file, counter)
	#pdb_file = loop_modeler(pdb_file, counter)

	job_output = "refine_chi%s" % counter
	pdb_file = refine_chi(pdb_file, counter,
	kT, cycles,	jobs, job_output)

	
	#print(pdb_file)
	#pdb_file = low_res(pdb_file, counter,
	#kT, smallmoves, shearmoves, backbone_angle_max, cycles)

	job_output = "refine_bb%s" % counter
	pdb_file = refine_bb(pdb_file, counter,
	kT, smallmoves, shearmoves, backbone_angle_max, cycles,	jobs, job_output)

	#if counter % 2 == 0 and counter != 0:
		#pdb_file = fast_relaxation(pdb_file, counter)
	pdb_file = fast_relaxation(pdb_file, counter)

	counter += 1

# for i in range(0, iterations):
# 	if i >= 4:
# 		if walkerlist[-1] == walkerlist[-2] == walkerlist[-3] == walkerlist[-4] == walkerlist[-5] == walkerlist[-6] == walkerlist[-7] == walkerlist[-8] == walkerlist[-9] == walkerlist[-10] == walkerlist[-11]:
# 			fd.close()
# 			sys.exit("**Sin cambios en las ultimas 11 lineas, protocolo terminado.**")

# 	pdb_file = fast_relaxation(pdb_file, counter)

# 	pdb_file = low_res(pdb_file, counter,
#  	kT, smallmoves, shearmoves, backbone_angle_max, cycles)

#  	counter += 1

pose = Pose()
pose_from_pdb(pose, pdb_file)
scorefxn= create_score_function('talaris2013')
scorefxn.set_weight(rg , 1)

aarmsd = all_atom_rmsd(pose, rmsd_pose)
carmsd = native_CA_rmsd(pose, rmsd_pose)
gdt = CA_gdtmm(pose, rmsd_pose)

pymover = PyMOL_Mover()
pymover.apply(rmsd_pose)

#pose.dump_pdb("reFINAL_%s.pdb" % pdb_filename)

walkers = GiveWalkers(pdb_file)
fd.write("**WALKERS FIN %s : %s\n" % (counter, walkers))
walkerlist.append(len(walkers))
fe.write("%*s%*s scorefxn:%*s flaws:%*s carmsd:%*s aarmsd:%*s gdt:%*s\n" % (15, "FINAL", 4, counter, 14, str(scorefxn(pose)), 3, str(len(walkers)), 14, str(carmsd), 14, str(aarmsd), 14, str(gdt)))
fd.close()
fe.write("**Number of iterations completed, protocol finished.**\n")
now = datetime.datetime.now()
fe.write(str(now))
fe.close()
sys.exit("**Number of iterations completed, protocol finished.**")