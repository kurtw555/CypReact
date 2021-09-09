/**
 * Authors: Yannick Djoumbou Feunang
 * Class Description:
 */


package reactantpredictor.utils;

import org.openscience.cdk.CDKConstants;
import org.openscience.cdk.aromaticity.Aromaticity;
import org.openscience.cdk.aromaticity.ElectronDonation;
import org.openscience.cdk.aromaticity.Kekulization;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.ConnectivityChecker;
import org.openscience.cdk.graph.Cycles;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.layout.StructureDiagramGenerator;
import org.openscience.cdk.qsar.descriptors.molecular.WeightDescriptor;
import org.openscience.cdk.qsar.result.IDescriptorResult;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import ambit2.smarts.query.SMARTSException;
import ambit2.smarts.query.SmartsPatternCDK;

public class MoleculeExplorer {

	
	/**
	 * This function applies some preprocessing operations, such as setting the
	 * flag of atoms from aromatic rings to "ISAROMATIC", and kelulizing
	 * molecules.
	 * 
	 * @param molecule
	 *            : A molecule of interest
	 * @return : A processed molecule (AtomContainer)
	 * @throws CDKException
	 */
	public static IAtomContainer preprocessContainer(IAtomContainer molecule)
			throws CDKException {
		
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		CDKHydrogenAdder.getInstance(molecule.getBuilder()).addImplicitHydrogens(molecule);
		 
	    Aromaticity aromaticity = new Aromaticity( ElectronDonation.cdk(), Cycles.cdkAromaticSet());
//		Aromaticity aromaticity = new Aromaticity(ElectronDonation.daylight(), Cycles.all());
		
		for (IBond bond : molecule.bonds()) {
			if (bond.getFlag(CDKConstants.SINGLE_OR_DOUBLE)) {
				bond.setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(0).setFlag(CDKConstants.ISAROMATIC, true);
				bond.getAtom(1).setFlag(CDKConstants.ISAROMATIC, true);

			} 
		}
//		aromaticity.apply(molecule);

		Kekulization.kekulize(molecule);
		
		
		StructureDiagramGenerator sdg = new StructureDiagramGenerator();
		sdg.setMolecule(molecule);
		sdg.generateCoordinates();
		
		IAtomContainer layedOutMol = sdg.getMolecule();

//		StringWriter w2 = new StringWriter();
//		MDLWriter mw2 = new MDLWriter(w2);
//		mw2.write(layedOutMol);	
//		System.out.println("After preprocessing\n" + w2.toString() + "\n\n");

		return layedOutMol;
	}
	public static boolean isMixture(IAtomContainer molecule) throws CDKException{
		// compound is not a mixture (checkConnectivity returns 2 or more atomContainers)
		boolean mixture = ConnectivityChecker.partitionIntoMolecules(molecule).getAtomContainerCount()>1;
		return mixture;	
	}
	
	public static boolean containsCarbon(IAtomContainer molecule) {
		boolean carbon = false;
		for(IAtom at : molecule.atoms()){
			if(at.getAtomicNumber() == 6){
				carbon = true;
				break;
			}
		}	
		return carbon;
	}	
	public static boolean isEtherLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		String constraints = "[$([#8;X2][#6;A;H2X4]!@-[#6;A;X4](!@-[!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([#8]!@-[#6;A;X4](!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])]),"
				+ "$([!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])][#6;A;H2X4]!@-[#6;A;X3](!@=[O;X1])!@-[#6;A;H2X4][!#1!#6;$([OX2H1]),$([OX2]-[CX4H2]-[#6;A])])"
				+ "]";		
		SmartsPatternCDK pattern = new SmartsPatternCDK(constraints);
		b = pattern.hasSMARTSPattern(molecule)>0;
			
		return b;
		
	}
	public static boolean isGlyceroLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		String constraints = "[$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#6;A;H2X4R0][#8]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])]),$([OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])][#6;A;H2X4R0][#6;A;H1X4R0]([#6;A;H2X4R0][OX1-,OX2H1,$([OX2](-[CX4])[CX4;R0;H2]),$([OX2]-[CX3]=[CX3]),$([OX2]-[CX3](=O)-[#1,#6])])[#8;X2]-[CX4,$([CX3]=[CX3]),$([CX3](=O)-[#1,#6])])]"	;
		SmartsPatternCDK pattern = new SmartsPatternCDK(constraints);
		b = pattern.hasSMARTSPattern(molecule)>0;
			
		return b;
		
	}	
	public static boolean isGlycerophosphoLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		String constraints =  				"[#8]P([!#1!#6;OX1-,OX2H1,$([O]-[#6])])(=O)[#8;R0][#6;A;H2X4R0][#6;A;H1X4R0]([R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])])[#6;A;H2X4R0][R0;OX1-,OX2H1,$([OX2]-[#1,CX4,$(C(=O)-[#1,#6])])]";
		SmartsPatternCDK pattern = new SmartsPatternCDK(constraints);
		b = pattern.hasSMARTSPattern(molecule)>0;
			
		return b;
		
	}

	public static boolean isSphingoLipid(IAtomContainer molecule) throws SMARTSException{	
		
		boolean b = true;
		String constraints = "[$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8]),$([#6;A;H3X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4][#6;A;H2X4]-[#6;A;H1X4]=[#6;A;H1X4]-[#6;A;H1X4]([#8;A;X2H1,X1-])[#6;A;H1X4]([!#1!#6;NX4H3+,NX3H2,$([#7](-[#1])-[#6]([#8])-[#6,#1])])[#6;A;H2X4][#8])]";
		SmartsPatternCDK pattern = new SmartsPatternCDK(constraints);
		b = pattern.hasSMARTSPattern(molecule)>0;
			
		return b;
		
	}
	/**
	 * This function will check if the input molecule is a possible reactant or not.
	 * It returns false if the molecule matches any of the listed conditions in the if clause	
	 * @param molecule
	 * @return
	 * @throws CDKException
	 * @throws SMARTSException
	 */
	public boolean isInvalidCandidate(IAtomContainer oneMole) throws CDKException, SMARTSException{
		boolean invalid = false;		
		IAtomContainer molecule = preprocessContainer(oneMole);
        IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);		
		Double weight =  MolecularFormulaManipulator.getMajorIsotopeMass(formula);
		if((weight > 1500.0 || isMixture(molecule) || 
			isGlycosylatedCompound(molecule) ||
			isGlutathioneConjugate(molecule)||
			isSulfatedCompound(molecule) ||
			isGlycinatedCompound(molecule) ||
			isTaurinatedCompound(molecule) ||
			isGlutamateConjugate(molecule) ||
			isAcylcarnitineConjugate(molecule) ||
			isCysteinylGlycineConjugate(molecule) ||					
			isAcylCoAConjugate(molecule) ||
			isTetrapyrrole(molecule) ||
			isSaccharide(molecule) ||
			isEtherLipid(molecule) ||
			isGlyceroLipid(molecule) ||
			isGlycerophosphoLipid(molecule) ||
			isGlycerol_3_Phosphate_Inositol(molecule) ||
			isSphingoLipid(molecule) )){
					
				invalid = true;
				return invalid;
		}
	
//		IAtomContainer pmol = preprocessContainer(molecule);
//		
//		if(isMixture(pmol)){ //|| containsCarbon(pmol)
//			invalid = true;
//		} else{
//			if(isEtherLipid(pmol) || isGlyceroLipid(pmol) || 
//					isGlycerophosphoLipid(pmol)	 || isSphingoLipid(pmol)){
//				invalid = true;
//			}
//		}
		//Double.valueOf(weight.toString()) > 1500.0 || 
		return invalid;
	}	
	/**
	 * Functions are from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isGlycosylatedCompound(IAtomContainer molecule) throws SMARTSException {		
		SmartsPatternCDK glucuronidePattern = new SmartsPatternCDK("[#8;A;X2H1,X1-][#6](=O)-[#6;R1]-1-[#8;R1]-[#6;R1](-[#6,#7,#8,#16])-[#6;R1](-[#1,#7,#8,#16;R0])-[#6;R1](-[#1,#7,#8,#16;R0])-[#6;R1]-1-[#1,#7,#8,#16;R0]");
		SmartsPatternCDK glycosylMoietyPattern = new SmartsPatternCDK("["
				+ "$(CC1OC(O)C(O)C(O)C1O),"
				+ "$(CC1OC(O)C(O)CC1O),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[*,#1;OX2H1,$(NC(=O)C)])-[#6](!@-[#8])-[#6]-1!@-[#8]),"
				+ "$([#6]!@-[#6]-1-[#8]-[#6](!@-[#8])-[#6](!@-[OX2H1,$(NC(=O)C)])-[#6]-[#6]-1!@-[#8])"
				+ "]"
				);
		
		glucuronidePattern.match(molecule);
		glycosylMoietyPattern.match(molecule);		
		boolean isGlycosylated = glucuronidePattern.hasSMARTSPattern(molecule)>0 || glycosylMoietyPattern.hasSMARTSPattern(molecule)>0 ;
			
		return isGlycosylated;
	}	
	/**
	 * Functions are from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isGlutathioneConjugate(IAtomContainer molecule) throws SMARTSException{
		SmartsPatternCDK glutathioneConjugatePattern = new SmartsPatternCDK("[H][#7]([#6;A;H2X4][#6](-[#8])=O)-[#6](=O)[#6;A;H1X4]([#6;A;H2X4][#16][#6,#8,#16;A])[#7]([H])-[#6](=O)[#6;A;H2X4][#6;A;H2X4][#6;A;H1X4]([#7])[#6](-[#8])=O");	
		return glutathioneConjugatePattern.hasSMARTSPattern(molecule)>0;	
	}
	/**
	 * Functions are from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isSulfatedCompound(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK sulfatedRadicalPattern = new SmartsPatternCDK("[#6]-[#8;X2]S([#8;A;X2H1,X1-])(=[O;X1])=[O;X1]");	
		return sulfatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}
	/**
	 * Functions are from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isGlycinatedCompound(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK glycinatedRadicalPattern = new SmartsPatternCDK("[#6][#7;A;H1X3]!@-[#6;A;H2X4]!@-[#6;X3](!@-[#8;A;X2H1,X1-])=[O;X1]");	
		return glycinatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}	
	/**
	 * Functions are from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isTaurinatedCompound(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK taurinatedRadicalPattern = new SmartsPatternCDK("[#7;A;H1X3]!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-S(!@-[#8;A;X2H1,X1-])(!@=O)!@=O");	
		return taurinatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}
	/**
	 * Functions are from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isGlutamateConjugate(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK glutamatedRadicalPattern = new SmartsPatternCDK("[#7;A;H1X3][#6;A;H1X4](!@-[#6;A;H2X4]!@-[#6;A;H2X4]!@-[#6]([#8;A;X2H1,X1-])=O)!@-[#6]([#8;A;X2H1,X1-])=O");	
		return glutamatedRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}	
	/**
	 * Functions are from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isAcylcarnitineConjugate(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK acylcarnitineRadicalPattern = new SmartsPatternCDK("[#6]-[#6](=O)-[#8][#6;A;H1X4]([#6;A;H2X4][#6]([#8;A;X2H1,X1-])=O)[#6;A;H2X4][N;X4+]([#6;A;H3X4])([#6;A;H3X4])[#6;A;H3X4]");	
		return acylcarnitineRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}
	/**
	 * Functions are from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isCysteinylGlycineConjugate(IAtomContainer molecule) throws SMARTSException {
		SmartsPatternCDK cysteinylglycineRadicalPattern = new SmartsPatternCDK("[#6]-[#16;X2]-[#6]-[#6;X4](-[#7;X3])-[#6](=O)[#7;A;H1X3][#6;X4]-[#6]([#8;A;X2H1,X1-])=O");	
		return cysteinylglycineRadicalPattern.hasSMARTSPattern(molecule)>0;	
	}
	/**
	 * Functions are from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isAcylCoAConjugate(IAtomContainer molecule) throws SMARTSException{
		SmartsPatternCDK glutathioneConjugatePattern = new SmartsPatternCDK("[#6;R0]-[#6;R0](=O)-[#16]-[#6]-[#6]-[#7]-[#6](=O)-[#6]-[#6]-[#7]-[#6](=O)-[#6](-[#8])C([#6])([#6])[#6]-[#8]P([#8])(=O)[#8]P([#8])(=O)[#8]-[#6]-[#6]-1-[#8]-[#6](-[#6](-[#8])-[#6]-1-[#8]P([#8])([#8])=O)-[#7]-1-[#6]=[#7]-[#6]-2=[#6]-1-[#7]=[#6]-[#7]=[#6]-2-[#7]");	
		return glutathioneConjugatePattern.hasSMARTSPattern(molecule)>0;	
	}
	/**
	 * Functions are derived from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isTetrapyrrole(IAtomContainer molecule) throws SMARTSException{	
		SmartsPatternCDK tetrapyrrolePattern = new SmartsPatternCDK("["
				+ "$([#6]-1-[#6]-[#7]=[#6](-[#6]-1)\\[#6]=[#6]-1\\[#6]-[#6]-[#6](=[#7]-1)-[#6]-1-[#6]-[#6]-[#6](\\[#6]=[#6]-2\\[#6]-[#6]-[#6]=[#7]-2)=[#7]-1)"
				+ ",$([#6](~c1cccn1)~c1ccc(~[#6]~c2ccc(~[#6]~c3cccn3)n2)n1)"
				+ ",$([#6](~[#6]~1~[#6]~[#6]~[#6]~[#7]~1)~[#6]~1~[#6]~[#6]~[#6](~[#6]~[#6]~2~[#6]~[#6]~[#6](~[#6]~[#6]~3~[#6]~[#6]~[#6]~[#7]~3)~[#7]~2)~[#7]~1)"
				+ ",$([#6]-1-[#6]-[#7+]=[#6](-[#6]-1)-[#6]-1=[#6]-2-[#7]\\[#6](=[#6]/[#6]-3=[#7+]/[#6](/[#6]-[#6]3)=[#6]\\[#6]-3=[#6]-[#6]=[#6]-[#7]3)-[#6]=[#6]-2-[#6]-[#6]-1)"
				+ "]");	
		return tetrapyrrolePattern.hasSMARTSPattern(molecule)>0;
		
	}	
	/**
	 * Function is derived from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isSaccharide(IAtomContainer molecule) throws SMARTSException{	
		SmartsPatternCDK saccharidePattern = new SmartsPatternCDK("["
				+ "$([#6]-[#8;R0][#6;A;H1X4R1]1[#8][#6;A;H1X4R1]([#6;R0][#8,#7,#16;A])[#6;A;H1X4R1]([#8,#7,#16;A])[#6;A;H1X4R1]([#8,#7,#16;A])[#6;A;H1X4R1]1[#8,#7,#16;A])"
				+ ",$([#6]-[#8;R0][#6;A;H1X4R1]1[#8][#6;A;H1X4R1]([#6;R0][#8,#7,#16;A])[#6;A;H1X4R1]([#8,#7,#16;A])[#6;A;H1X4R1]1[#8,#7,#16;A])"
				+ "]");	
		return saccharidePattern.hasSMARTSPattern(molecule)>0;		
	}
	/**
	 * Function is derived from BioTransformer implemented by Yannick
	 * @param molecule
	 * @return
	 * @throws SMARTSException
	 */
	public static boolean isGlycerol_3_Phosphate_Inositol(IAtomContainer molecule) throws SMARTSException{	
		SmartsPatternCDK isGlycerol_3_Phosphate_Inositol = new SmartsPatternCDK("[#8][#6;A;H1X4R1]1[#6;A;H1X4R1]([#8])[#6;A;H1X4R1]([#8])[#6;A;H1X4R1]([#8]P([#8;X2H1,X1-])(=O)[#8]-[#6;H2X4]-[#6;H1X4](-[#6;H2X4]-[#8]-[#6]([#6,#1;A])=O)-[#8]-[#6]([#6,#1;A])=O)[#6;A;H1X4R1]([#8])[#6;A;H1X4R1]1[#8]");	
		return isGlycerol_3_Phosphate_Inositol.hasSMARTSPattern(molecule)>0;		
	}
}
