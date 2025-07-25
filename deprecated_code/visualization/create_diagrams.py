from ccdc import io
from ccdc.diagram import DiagramGenerator
from ccdc.io import EntryReader
from ccdc.molecule import Molecule

diagram_generator = DiagramGenerator()

diagram_generator.settings.font_size = 12
diagram_generator.settings.line_width = 1.6
diagram_generator.settings.image_width = 500
diagram_generator.settings.image_height = 500

# csd_reader = EntryReader('CSD')

# mol = csd_reader.molecule('TITTUO')

smiles = "O=[U]123(=O)(OSc4ccc(O)cc4C=N1C(-Ns2=C(P)c1ccbcc1P3)SCCCCC)[NH]CCC"  
mol = Molecule.from_string(smiles, "smiles") 

img = diagram_generator.image(mol)

img.save("aspirin.png")
 
# abebuf = csd_reader.entry('ABEBUF')
# img = diagram_generator.image(abebuf)


# Highlighting
# from ccdc.search import SubstructureSearch, SMARTSSubstructure

# searcher = SubstructureSearch()
# sub_id = searcher.add_substructure( SMARTSSubstructure('c1ncccc1') )
# hits = searcher.search(abebuf.molecule)
# selection = hits[0].match_atoms()
# diagram_generator.settings.element_coloring = False
# img = diagram_generator.image(abebuf, highlight_atoms=selection)
# diagram_generator.settings.element_coloring = True
# abahui = csd_reader.molecule('ABAHUI')
# diagram_generator.settings.detect_intra_hbonds = True
# diagram_generator.settings.shrink_symbols = False
# diagram_generator.settings.return_type = 'SVG'
# image = diagram_generator.image(abahui)