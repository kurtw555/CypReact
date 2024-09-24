
/**
 *
 * @author Leon Pu, Anisha Jauhari
 *
 */
package biotransformerapi;

import java.io.*;
import java.util.*;

import biotransformerapis.predictors.BioTransformerAPIs;
import org.apache.commons.lang3.StringUtils;
import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IAtomContainerSet;

import reactantpredictor.ReactantPred;
import reactantpredictor.SdfToSample;
import reactantpredictor.utils.FileUtils;
import weka.classifiers.Classifier;
import weka.core.Attribute;
import weka.core.DenseInstance;
import weka.core.Instance;
import weka.core.Instances;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import utils.Validator;


public class BioTransformerAPI implements BioTransformerAPIs {
	private static Logger logger = null;
	static {
			logger = LogManager.getLogger(BioTransformerAPI.class);
		}
		/**
		 * This function implements the BioTransformerAPIs interface. It returns an IAtomcontainerSet with enzymes tags that signify
		 * if the enzyme could catalyse the molecule or not.
		 * @param parameters - Contains molecule (IAtomContainerSet or Sdf file), enzymes. Make sure enzyme list is of form CYP1A2 etc.
		 * @return IAtomContainerSet or Sdf file.
		 * @throws Exception
		 */
	public Object predict(LinkedHashMap<String, Object> parameters) throws Exception {

		Validator validator = new Validator();
		Object results  = null;
		IAtomContainerSet inputMolecules = null;
		Validator.ioFormats inputFormat = Validator.ioFormats.valueOf( StringUtils.upperCase( (String) parameters.get("inputFormat")));
		Validator.ioFormats outputFormat = Validator.ioFormats.valueOf( StringUtils.upperCase( (String) parameters.get("outputFormat")));
		String[] enzymes = (String[]) parameters.get("properties");

		if(validator.validateParameters(parameters)) {
			if(inputFormat == Validator.ioFormats.IATOMCONTAINERSET) {
				inputMolecules = (IAtomContainerSet) parameters.get("input");
			}
			else if(inputFormat == Validator.ioFormats.SDFILE) {
				inputMolecules = FileUtils.parseSdf(String.valueOf(parameters.get("input")) );
			}

			IAtomContainer molsWithPredictions = null;
			try {
				logger.error("Values = "+ inputMolecules.getAtomContainer(0) + enzymes);
				molsWithPredictions = predictESSForMolecule(inputMolecules.getAtomContainer(0), enzymes);
				logger.error("MolWithPredictions - " + molsWithPredictions);
			} catch (Exception e) {
				e.printStackTrace();
			}

			IAtomContainerSet iAtomContainerSet = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
			iAtomContainerSet.addAtomContainer(molsWithPredictions);
			if(outputFormat == Validator.ioFormats.IATOMCONTAINERSET) {
				results = iAtomContainerSet;
			}
			else if(outputFormat == Validator.ioFormats.SDFILE) {
				String outputFileName = (String) parameters.get("output");

				if(molsWithPredictions != null) {
					FileUtils.saveAtomContainerSetToSDF(iAtomContainerSet, outputFileName);
					results = outputFileName;
				}
			}
		}

		return results;
	}

	/**
	 * This function predicts the Enzyme Substrate Specificity for the given molecule and returns the molecule with a true/false flag
	 * corresponding to the enzyme.
	 * @param molecule
	 * @param enzymes Make sure enzyme list is of form CYP1A2 etc.
	 * @return Molecule with enzyme tags.
	 * @throws Exception
	 */
	public IAtomContainer predictESSForMolecule (IAtomContainer molecule, String[] enzymes) throws Exception {
		List<String> catalyzedEnzymes =  predictEnzymeSubstrateSpecificity(molecule, enzymes);
		logger.info("CatalyzedEnzymes - " + catalyzedEnzymes);
		boolean isCatalyzed;
		Map<Object, Object> properties = new HashMap<>();
		for (String enzyme: enzymes) {
			isCatalyzed = false;
			if (catalyzedEnzymes.contains(enzyme)) {
				isCatalyzed = true;
			}
			properties.put(enzyme, isCatalyzed);
		}
		molecule.addProperties(properties);
		return molecule;
	}

	/**
	 * This function returns an arraylist of enzymes which catalysed the given molecule.
	 * @param molecule
	 * @param enzymes Make sure enzyme list is of form CYP1A2 etc.
	 * @return List of string of enzyme names which were able to catalyse the molecule.
	 * @throws Exception
	 */
	public List<String> predictEnzymeSubstrateSpecificity (IAtomContainer molecule, String[] enzymes) throws Exception {
		SdfToSample sdfTool = new SdfToSample();
		IAtomContainerSet inputMolecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		inputMolecules.addAtomContainer(molecule);
		Instances testSet = sdfTool.generateAllInstances(inputMolecules);
		boolean result = false;
		List<String> catalyzedEnzymeList = new ArrayList<>();
		for(String enzyme : enzymes) {
			String cyp = enzyme.substring(3);
			if(cyp.equals("2D6") || cyp.equals("2C9")){
				String supportfile = "/CYP" + cyp + "/" + "supportfile.csv";
				logger.error("SupportFile - "+ supportfile);
				result = makeEnsemblePrediction(cyp,testSet,supportfile);
			} else {
				String[] pathArray = ReactantPred.generatePath(cyp);
				logger.error("In Else - " + pathArray[0] + "  " + pathArray[1]);
				String model = pathArray[0];
				String supportfile = pathArray[1];
				result = makePrediction(cyp, model, testSet, supportfile);

			}
			if (result) {
				catalyzedEnzymeList.add(enzyme);
			}
		}
		return catalyzedEnzymeList;
	}

	public boolean predictReactant(IAtomContainer molecule, String[] enzymeList) throws Exception{
		SdfToSample sdfTool = new SdfToSample();
		IAtomContainerSet inputMolecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		inputMolecules.addAtomContainer(molecule);
		Instances testSet = sdfTool.generateAllInstances(inputMolecules);
		boolean result = false;
		for(int i = 0; i < enzymeList.length; i++){
			String cyp = enzymeList[i];
			if(cyp.equals("2D6") || cyp.equals("2C9")){
				String supportfile = "/CYP" + cyp + "/" + "supportfile.csv";
				result = makeEnsemblePrediction(cyp,testSet,supportfile);
			}
			else{
				String[] pathArray = ReactantPred.generatePath(cyp);
				String model = pathArray[0];
				String supportfile = pathArray[1];
				result = makePrediction(cyp, model, testSet, supportfile);
	
			}
			if(result == true){
				return result;
			}
		}
		return result;
	}
	
	public boolean predictReactant(IAtomContainer molecule, String cyp) throws Exception{
		//System.out.println("Running cypReact BioTransformerAPI");
		SdfToSample sdfTool = new SdfToSample();		
		IAtomContainerSet inputMolecules = DefaultChemObjectBuilder.getInstance().newInstance(IAtomContainerSet.class);
		inputMolecules.addAtomContainer(molecule);
		Instances testSet = sdfTool.generateAllInstances(inputMolecules);
		boolean result = false;
		if(cyp.equals("2D6") || cyp.equals("2C9")){
			String supportfile = "/CYP" + cyp + "/" + "supportfile.csv";
			System.out.println("Running Ensemble prediction");
			result = makeEnsemblePrediction(cyp,testSet,supportfile);
		}
		else{
			String[] pathArray = ReactantPred.generatePath(cyp);
			String model = pathArray[0];
			String supportfile = pathArray[1];
			System.out.println("Running single prediction");
			result = makePrediction(cyp, model, testSet, supportfile);

		}
		if(result == true){
			return result;
		}
		return result;
	}
	/**
	 * Predict whether the given molecules in the sdf/csv file are reactants or not for one CYP by using non-ensemble methods
	 * @param cyp
	 * @param model
	 * @param testSet
	 * @param supportfile
	 * @param predictedResult
	 * @return updated PredictedResults
	 * @throws Exception
	 */

	public boolean makePrediction(String cyp, String model, Instances testSet, String supportfile) throws Exception{
		InputStream model_Stream = getClass().getResourceAsStream(model);
		logger.error("Model Stream in Make PRedictions - " + model_Stream + "  " +cyp + model + supportfile);
		Classifier cls = (Classifier) weka.core.SerializationHelper.read(model_Stream);//Load the model
		InputStream supportPath_Stream = getClass().getResourceAsStream(supportfile);
		//File supfile = new File(supportPath_Stream);
		//FileReader spfr = new FileReader(supportPath_Stream);
		BufferedReader spbr = new BufferedReader(new InputStreamReader(supportPath_Stream));
		
		//Create attribute arraylist,mean arraylist, max arraylist, min arraylist
		ArrayList<String> attList = new ArrayList<String>();
		ArrayList<String> meanList = new ArrayList<String>();
		ArrayList<String> maxList = new ArrayList<String>();
		ArrayList<String> minList = new ArrayList<String>();
		//Write the output 
		int counter = 0;

		//attribute	mean	max	min. Skip the first line that contains all the titles
		String supLine = spbr.readLine();
		while((supLine = spbr.readLine())!=null){
			String[] elements = supLine.split(",");
			attList.add(elements[0]);
			meanList.add(elements[1]);
			maxList.add(elements[2]);
			minList.add(elements[3]);
		}
		int countR = 0;
		int countT = 0;
		Instances matched = matchAttributes(testSet,attList,meanList,maxList,minList);
		matched.setClassIndex(matched.numAttributes()-1);
		//System.out.println("--------------------------------------");
		for(int i = 0; i<matched.size();i++){
			counter++;
			Instance oneSample = matched.get(i);
			double result = cls.classifyInstance(oneSample);
			String pred = "";
		
			if(result==0.0){
				return false;
			}
			else if(result == 1.0) return true;

		}
		//System.out.println("N:" + countT + " R:" + countR );
		return false;
	}
	
	public boolean makeEnsemblePrediction(String cyp, Instances testSet, String supportfile) throws Exception{
		Double result = 0.0;
		//Load support file and creat the instances
		InputStream supportPath_Stream = getClass().getResourceAsStream(supportfile);
		BufferedReader spbr = new BufferedReader(new InputStreamReader(supportPath_Stream));		
		//Declare attributes arraylist,mean arraylist, max arraylist, min arraylist, will be used in normalization
		ArrayList<String> attList = new ArrayList<String>();
		ArrayList<String> meanList = new ArrayList<String>();
		ArrayList<String> maxList = new ArrayList<String>();
		ArrayList<String> minList = new ArrayList<String>();
		//update mean, max, min values of each Attribute.
		String supLine = spbr.readLine();//Skip the first line that contains all the titles
		while((supLine = spbr.readLine())!=null){
			String[] elements = supLine.split(",");
			attList.add(elements[0]);
			meanList.add(elements[1]);
			maxList.add(elements[2]);
			minList.add(elements[3]);
		}
		Instances matched = matchAttributes(testSet,attList,meanList,maxList,minList);
		matched.setClassIndex(matched.numAttributes()-1);
		Instance oneSample = matched.get(0);
		//Load the ensemble models
		if(cyp.equals("2D6")){
			for(int i = 0; i < 5; i++){
				//String model = supportFoldPath + "supportfiles/CYP" + cyp + "/model/" + cyp + "_NR_" + i +".model";
				String model = "/CYP" + cyp + "/" + cyp + "_NR_" + i + ".model";  
				InputStream model_Stream = getClass().getResourceAsStream(model);
				Classifier cls = (Classifier) weka.core.SerializationHelper.read(model_Stream);
				result = result + cls.classifyInstance(oneSample);
				//System.out.println("current result: " + result);
			}
			result = result / 5;
		}
		if(cyp.equals("2C9")){
			for(int i = 0; i < 6; i++){
				//String model = supportFoldPath + "supportfiles/CYP" + cyp + "/model/" + cyp + "_NR_" + i +".model";
				String model = "/CYP" + cyp + "/" + cyp + "_NR_" + i + ".model";  
				InputStream model_Stream = getClass().getResourceAsStream(model);
				Classifier cls = (Classifier) weka.core.SerializationHelper.read(model_Stream);
				result = result + cls.classifyInstance(oneSample);
				//System.out.println("current result: " + result);
			}
			result = result / 6;
		}
		if(result >= 0.5) return true;
		else return false;
		
	}	
	/**
	 * Reformat the instances of the input molecules so they can be input to the learned models by applying feature selection and normalization
	 * @param testSet: The instances of the user's input molecules
	 * @param attList: The attributes that are used in the learned model
	 * @param meanList: The mean values of each attribute
	 * @param maxList: The max values of each attribute
	 * @param minList: The min values of each attribute
	 * @return
	 * @throws Exception
	 */
	public Instances matchAttributes(Instances testSet, ArrayList<String> attList, ArrayList<String> meanList, ArrayList<String> maxList, ArrayList<String> minList) throws Exception{
		//Create the class atrribute
		List my_nominal_values = new ArrayList(3); 
		my_nominal_values.add("R"); 
		my_nominal_values.add("N"); 
		Attribute lastAttribute  = new Attribute(attList.get(attList.size()-1), my_nominal_values);
		
		ArrayList<Attribute> atts = new ArrayList<Attribute>();
		for(int i = 0; i<attList.size();i++){
			if(i<attList.size()-1){
				Attribute Attribute = new Attribute(attList.get(i));
				atts.add(Attribute);
			}
			else{
				
				atts.add(lastAttribute);
			}
		}
		
		int length = testSet.size();
		int numAttributes = attList.size();
		Instances matched = new Instances("Rel", atts, length);
		Instance sample = new DenseInstance(numAttributes);
		
		for(int j = 0; j< length; j++){
			ArrayList misFeatures = new ArrayList();
			Instance temp = testSet.get(j);
			for(int vidx = 0; vidx < numAttributes; vidx++){
				  Attribute att = atts.get(vidx);
				  if(vidx<numAttributes-1){
					  
					  String stAtt = att.toString();
					  String[] whatever = stAtt.split(" ");
					  
					  if(whatever.length!=3){
						  String combineSpaces = "";
						  for(int kk = 1; kk<whatever.length-1; kk++){
							  if(kk<whatever.length - 2){
								  combineSpaces = combineSpaces + whatever[kk]+ " ";
							  }
							  else combineSpaces = combineSpaces + whatever[kk];
						  }
						  whatever[1] = combineSpaces;
					  }
					 
					  if(whatever[1].contains("'")){
						  whatever[1] = whatever[1].replace("'", "");
					  }
					  Attribute convAtt = testSet.attribute(whatever[1]);
					  double vle = temp.value(convAtt);
					  if(!Double.isNaN(vle)){
						  double min = Double.parseDouble(minList.get(vidx));
						  double max = Double.parseDouble(maxList.get(vidx));
						  double normedVal = (vle-min)/(max-min);
						  sample.setValue(att, normedVal);
					  }
					  else{
						  //If the value of a attribute is not a number, use the mean value of that attribute to replace it.
						  misFeatures.add(whatever[1]);
						  double syn = Double.parseDouble(meanList.get(vidx));
						  sample.setValue(att, syn);

					  }
				  }
				  else{
					  //Set class with random value. Does not matter
					  sample.setValue(lastAttribute, "N");
					  
				  }
			}
			//Output those missing features
			if(misFeatures.size()>0){
				int numMisFeatures = misFeatures.size();
				String misFeatureString = "Molecule" + j + ": " +  numMisFeatures + 
						" features could not be computed from the structure representation, and were thus synthesized: " + misFeatures.get(0);
				for(int i = 1; i < numMisFeatures; i++){
					misFeatureString = misFeatureString + "," + misFeatures.get(i);
				}
				//System.out.println(misFeatureString);
			}
			
			
			matched.add(sample);
		}
		
		return matched;
		
	}
}
