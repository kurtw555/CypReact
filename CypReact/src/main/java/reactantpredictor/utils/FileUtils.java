package reactantpredictor.utils;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.io.iterator.IteratingSDFReader;
import org.openscience.cdk.silent.SilentChemObjectBuilder;

import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;

public class FileUtils {
    /**
     * Function that parses SMILES Strings in provided SDF file
     *
     * @param sdfFileName	Name of SDF file with SMILES Strings
     * @return	IAtomContainerSet of parsed SMILES Strings
     * @throws IOException    Invalid file type
     */
    public static IAtomContainerSet parseSdf(String sdfFileName) throws IOException {
        IAtomContainerSet containers = DefaultChemObjectBuilder.getInstance().newInstance(
                IAtomContainerSet.class);

        IChemObjectBuilder bldr = SilentChemObjectBuilder.getInstance();

        IteratingSDFReader sdfr = new IteratingSDFReader(new FileReader(sdfFileName), bldr);

        while (sdfr.hasNext()){
            IAtomContainer mol = sdfr.next();
            containers.addAtomContainer(mol);
        }

        sdfr.close();
        return containers;

    }
    /**
     * Function that saves an IAtomContainerSet of containers to an SDF file
     *
     * @param containers	IAtomContainerSet of containers to save
     * @param outputFileName	Name of output file to send containers to
     * @throws CDKException    CDK class exceptions
     * @throws IOException	Invalid file type
     */
    public static void saveAtomContainerSetToSDF(IAtomContainerSet containers, String outputFileName) throws CDKException, IOException{
        SDFWriter sdfWriter = new SDFWriter(new FileOutputStream(outputFileName));
        sdfWriter.write(containers);
        sdfWriter.close();
    }
}
