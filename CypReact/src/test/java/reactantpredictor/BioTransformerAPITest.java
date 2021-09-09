package reactantpredictor;

import biotransformerapi.BioTransformerAPI;
import org.junit.Test;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.silent.SilentChemObjectBuilder;
import org.openscience.cdk.smiles.SmilesParser;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.util.List;

public class BioTransformerAPITest {

    @Test
    public void testPredictEnzymeSubstrateSpecificity() throws Exception {
        IChemObjectBuilder builder = SilentChemObjectBuilder.getInstance();
        SmilesParser smiParser       = new SmilesParser(builder);
        IAtomContainer molecule = smiParser.parseSmiles("CC(=O)Nc1ccccc1");
        String[] enzymes = {"CYP1A2", "CYP2A6", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", "CYP2D6", "CYP2E1", "CYP3A4"};
        BioTransformerAPI bioTransformerAPI = new BioTransformerAPI();
        List<String> catalyzedEnzymes = bioTransformerAPI.predictEnzymeSubstrateSpecificity(molecule, enzymes);
        System.out.println("Catalyzed Enzymes - " + catalyzedEnzymes);
        //Catalyzed Enzymes - [CYP1A2, CYP2A6, CYP2D6]
        assertEquals("There should be 3 enzymes", 3, catalyzedEnzymes.size());
    }
}