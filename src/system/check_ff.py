import openmm.app as app

def check_residue(resname, xml_files):
    print(f"Checking {resname} in {xml_files}...")
    try:
        ff = app.ForceField(*xml_files)
        # Try to get a template
        # ForceField._templates is a dict: {name: template}
        # But it might be private.
        # We can iterate over templates.
        
        found = False
        for name, template in ff._templates.items():
            if name == resname:
                print(f"FOUND {resname}!")
                found = True
                # Print atoms
                print("Atoms:", [atom.name for atom in template.atoms])
                break
        
        if not found:
            print(f"{resname} not found.")
            
    except Exception as e:
        print(f"Error loading {xml_files}: {e}")

# Check common lipid files
xmls = ['amber14-all.xml', 'amber14/lipid17.xml', 'charmm36.xml']
# Note: charmm36.xml might need to be downloaded or is in openmm?
# OpenMM usually has amber14, charmm36.

# Stearate might be STE, SA, STEAR.
# Also check for generic fatty acids.
check_residue('STE', ['amber14-all.xml'])
check_residue('SA', ['amber14-all.xml'])
check_residue('STR', ['amber14-all.xml'])

# Check lipid17 explicitly
check_residue('STE', ['amber14/lipid17.xml'])
check_residue('SA', ['amber14/lipid17.xml'])

# Check CHARMM36 if available
try:
    check_residue('STE', ['charmm36.xml'])
except:
    pass
