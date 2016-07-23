
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
def get_unipro( ps ):
    uniq = set( [ Chem.MolToSmiles( mol[0], isomericSmiles =True ) for mol in ps ] )
    mols = [ Chem.MolFromSmiles( smi ) for smi in uniq ]
    return mols

def get_unimol( mols ):
    uniq = set( [ Chem.MolToSmiles( mol, isomericSmiles =True ) for mol in mols ] )
    mols = [ Chem.MolFromSmiles( smi ) for smi in uniq ]
    return mols
c5 = Chem.MolFromSmiles( 'C1CCCC1' )
c6 = Chem.MolFromSmiles( 'C1CCCCC1' )
c7 = Chem.MolFromSmiles( 'C1CCCCCC1' )
c45 = Chem.MolFromSmiles( 'C2CC1CCC1C2' )
c46 = Chem.MolFromSmiles( 'C2CCC1CCC1C2' )
c55 = Chem.MolFromSmiles( 'C1CC2CCCC2C1' )
c56 = Chem.MolFromSmiles( 'C2CCC1CCCC1C2' )
c66 = Chem.MolFromSmiles( 'C1CCC2CCCCC2C1' )
sp35 = Chem.MolFromSmiles( 'C1CC11CCCC1' )
sp36 =  Chem.MolFromSmiles( 'C1CC11CCCCC1' )
sp45 = Chem.MolFromSmiles( 'C1CCC11CCCC1' )
sp46 =  Chem.MolFromSmiles( 'C1CCC11CCCCC1' )
sp55 = Chem.MolFromSmiles( 'C1CCC11CCCC1' )
sp55 = Chem.MolFromSmiles( 'C1CCCC11CCCC1' )
sp65 = Chem.MolFromSmiles( 'C1CCCC11CCCCC1' )
sp66 =  Chem.MolFromSmiles( 'C1CCCCC11CCCCC1' )

ringsets = [
    c5, c6, c7,
    c45, c46, c55, c56, c66,
    sp35, sp36, sp45, sp46, sp55, sp65, sp66]
rxn = AllChem.ReactionFromSmarts( '[CH2&$(C(CC)(CC)):2]>>[NH:2]' )
monoamine = []
for ring in ringsets:
    ps = rxn.RunReactants( (ring,) )
    res = get_unipro( ps )
    monoamine.extend( res )
diamine = []
for ring in monoamine:
    ps = rxn.RunReactants( (ring,) )
    res = get_unipro( ps )
    diamine.extend( res )
diamine = get_unimol( diamine )
len(diamine)

Draw.MolsToGridImage(diamine, molsPerRow= 5 )
Draw.MolsToGridImage(diamine, molsPerRow= 5 )
