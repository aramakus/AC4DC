import sys
import numpy as np


def main():
  if len(sys.argv) < 2:
    print('No input .pdb file provided.')
    return 1
  
  if sys.argv[1].partition('.')[-1] != 'pdb':
    print('No input .pdb file provided.')
    return 1

  # Read Harry Format pdb file.
  atoms = read_pdb(sys.argv[1])

  atom_cnt = 0
  # Output atomic composition.
  for k, v in atoms.items():
    atom_cnt += len(v)
    print("{} : {} atoms".format(k, len(v)))


  # Output unit cell dimensions.
  uc_bounds = []
  for i in range(3):
    tmp = [0, 0]
    for k in atoms:
      atoms[k] = np.array(atoms[k])
      tmp[0] = min(tmp[0], np.min(atoms[k][:, i]))
      tmp[1] = max(tmp[1], np.max(atoms[k][:, i]))

    uc_bounds.append(tmp)
    print("{}:  min {} ; max {}".format(i, uc_bounds[i][0], uc_bounds[i][1]))

  feedin_dict = {
    'boundary' : uc_bounds,
    'atoms' : atoms
  }

  make_plasma_input(**feedin_dict)
  make_MD_input(**feedin_dict)

  return 0


def read_pdb(name):
  atoms = {}
  uc_bound = [[-0.5*79.0, 0.5*79.0], [-0.5*79.0, 0.5*79.0], [-19.0, 19.0]]
  traslations = [v[1] - v[0] for v in uc_bound]
  atom_names = ["H", "C", "N", "O", "P", "S"]
  file = open(sys.argv[1])

  for i, line in enumerate(file):
    tmp = line.split()
    """
    if tmp[0] == "TER":
      # Atoms already recorded are inside a unit cell.
      # Determine it's boundaries and drop atoms beyond it.
      for i in range(3):
        boundary = [0, 0]
        for k in atoms:
          axis = np.array(atoms[k])
          boundary[0] = min(boundary[0], np.min(axis[:, i]))
          boundary[1] = max(boundary[1], np.max(axis[:, i]))
        uc_bound.append(boundary)
    """

    if tmp[-1] not in atom_names: continue
    name = tmp[-1]
    x = float(tmp[-6])
    y = float(tmp[-5])
    z = float(tmp[-4])
    if len(uc_bound) > 0:
      if x < uc_bound[0][0] or x > uc_bound[0][1]: continue
      if y < uc_bound[1][0] or y > uc_bound[1][1]: continue
      if z < uc_bound[2][0] or z > uc_bound[2][1]: continue

    if name in atoms:
      atoms[name].append([x, y, z])
    else:
      atoms[name] = [[x, y, z]]

  file.close()
  return atoms


def make_plasma_input(**kwargs):
  uc_bounds = kwargs.get('boundary')
  atoms = kwargs.get('atoms')

  # Creare plasma input file template from Harry's pdb.
  plasma_file = open('Plasma.txt', 'w')
  plasma_file.write("""// Template AC4DC input file created from file.pdb\n\n""")
  plasma_file.write("""#ATOMS\n""")
  for k, v in atoms.items():
    plasma_file.write("""{} {}\n""".format(k, len(v)))
  
  # Create a volume section.
  plasma_file.write("""\n#VOLUME\n""")
  tmp = 1
  for coord in uc_bounds: tmp *= coord[1] - coord[0]
  plasma_file.write("""%.2f      // Volume per molecule in Angstrom^3.\n""" % tmp)
  plasma_file.write("""250          // Radius of a sample in Angstrom. Used for effective escape rate of photo-electrons.\n\n""")

  # Fill in th rest with default values.
  plasma_file.write("""#PULSE\n""")
  plasma_file.write("""8000         // Photon energy in eV.\n""")
  plasma_file.write("""15           // Pulse width in femtoseconds (defined as FWHM for Gaussian pulse).\n""")
  plasma_file.write("""10            // Pulse fluence in 10^4 * J/cm^2.\n\n""")

  plasma_file.write("""#NUMERICAL\n""")
  plasma_file.write("""4000         // Initial guess for number of time step points for rate equation.\n""")
  plasma_file.write("""If the value is 0, program skips rate equation soluving step.\n""")
  plasma_file.write("""12           // Number of threads in OpenMP.\n\n""")


  plasma_file.write("""#OUTPUT\n""")
  plasma_file.write("""4000         // Number of time steps in the output files.\n""")
  plasma_file.write("""N            // Write atomic charges in a separate file (Y/N)?\n""")
  plasma_file.write("""Y            // Write intensity in a separate file (Y/N)?\n""")
  plasma_file.write("""Y            // Write data for molecular dynamics (MD) in a separate file (Y/N)?""")

  plasma_file.close()


def make_MD_input(**kwargs):
  uc_bounds = kwargs.get('boundary')
  atoms = kwargs.get('atoms')

 # Creare plasma input file template from Harry's pdb.
  md_file = open('MD.xyz', 'w')
  num_atoms = 0
  for k, v in atoms.items():
    num_atoms += len(v)
  md_file.write("""{}\n""".format(num_atoms))
  for elem in uc_bounds:
    md_file.write("""%.3f   """ % (elem[1] - elem[0]))

  # Atoms.
  x_shift = uc_bounds[0][0]
  y_shift = uc_bounds[1][0]
  z_shift = uc_bounds[2][0]

  for name in atoms:
    for atom in atoms[name]:
      md_file.write("""\n{}  {:.3f}  {:.3f}  {:.3f}""".format(name, atom[0] - x_shift, atom[1] - y_shift, atom[2] - z_shift))



if __name__ == "__main__":
  main()