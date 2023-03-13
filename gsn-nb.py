from openmm.app import *
from openmm import *
from openmm.unit import *
from pdbfixer import *
from sys import stdout

# Load in the PDB strucure
# pdb = PDBFile('6H1F.pdb')

# Retrieve it from the RCSB and fix missing atoms
fixer = PDBFixer(pdbid="6H1F")

# Add missing (unresolved) residues. We don't want to model anything.
fixer.findMissingResidues()
fixer.missingResidues = {}
# fixer.addMissingResidues()


# Add missing (unresolved) atoms
fixer.findMissingAtoms()
fixer.addMissingAtoms()

# Protonate (roughly) at chosen pH
fixer.addMissingHydrogens(pH=7.0)

fixer.addSolvent(boxSize=10 * Vec3(1, 1, 1))

PDBFile.writeFile(fixer.topology, fixer.positions, open("6H1F-fixed.pdb", "w"))


# There is an "SCN" residue to remove
modeller = Modeller(fixer.topology, fixer.positions)

res_SCN = [r for r in modeller.topology.residues() if r.name == "SCN"]
modeller.delete(res_SCN)

PDBFile.writeFile(
    modeller.topology, modeller.positions, open("6H1F-modelled.pdb", "w"), keepIds=True
)


forcefield = ForceField("amber14-all.xml", "amber14/tip3pfb.xml")


system = forcefield.createSystem(
    modeller.topology,
    nonbondedMethod=PME,
    nonbondedCutoff=1 * nanometer,
    constraints=HBonds,
)

integrator = LangevinMiddleIntegrator(300 * kelvin, 1 / picosecond, 0.004 * picoseconds)


# Pressure, Temperature (only used for calculation),
# Frequency (how frequently the system should update the box size)
barostat = MonteCarloBarostat(1.0 * atmosphere, 300.0 * kelvin, 25)

# The barostat is added directly to the system rather than the larger
# simulation.
system.addForce(barostat)


# Combines the molecular topology, system, and integrator
# to begin a new simulation.
simulation = Simulation(modeller.topology, system, integrator)
simulation.context.setPositions(modeller.positions)

# Perform local energy minimization
print("Minimizing energy...")
simulation.minimizeEnergy(maxIterations=200)

# Write the minimized coordinates

PDBFile.writeFile(
    simulation.topology,
    simulation.context.getState(getPositions=True).getPositions(),
    open("6H1F-minimized.pdb", "w"),
    keepIds=True,
)


# Write the trajectory to a file called "output.pdb"
simulation.reporters.append(
    DCDReporter("output.dcd", reportInterval=1000, enforcePeriodicBox=True)
)


Nsteps = 5000

# Report infomation to the screen as the simulation runs
simulation.reporters.append(
    StateDataReporter(
        stdout,
        100,
        step=True,
        time=True,
        potentialEnergy=True,
        kineticEnergy=True,
        totalEnergy=True,
        temperature=True,
        volume=True,
        density=True,
        progress=True,
        remainingTime=True,
        speed=True,
        elapsedTime=True,
        separator=" ",
        totalSteps=Nsteps,
    )
)


# Run the simulation for 1000 timesteps
print("Running simulation...")
simulation.step(Nsteps)
