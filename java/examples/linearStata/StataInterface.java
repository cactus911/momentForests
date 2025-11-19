package examples.linearStata;

import com.stata.sfi.Data;
import com.stata.sfi.SFIToolkit;
import Jama.Matrix;
import java.util.stream.Collectors;

import java.util.*;

/**
*
* @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
*/

public class StataInterface {
	public static int RunLinearModel(String[] args) {
		
		boolean debug = false;
		
		// Block of code to replace system.out.println to work with Stata
		System.setOut(new java.io.PrintStream(new java.io.OutputStream() {
		    private StringBuilder buffer = new StringBuilder();

		    @Override
		    public void write(int b) {
		        if (b == '\n') {
		            SFIToolkit.displayln(buffer.toString());
		            buffer.setLength(0);
		        } else {
		            buffer.append((char) b);
		        }
		    }
		}));
		
        int rc = 0;
        
        // CV grid
        List<Integer> cvGridMinLeaf = new ArrayList<>();
        List<Double> cvGridMinImprovement = new ArrayList<>();
        List<Integer> cvGridMaxDepth = new ArrayList<>();
        
        try {
            // Parse arguments into key-value pairs 
            Map<String, String> opts = new HashMap<>();
            String fullArgs = String.join(" ", args);
            String[] tokens = fullArgs.split(" (?=\\w+=)");

            for (String token : tokens) {
                String[] parts = token.split("=", 2);
                if (parts.length == 2) {
                    opts.put(parts[0].trim().toLowerCase(), parts[1].trim());
                }
            }

            // Check required arguments
            if (!opts.containsKey("y") || !opts.containsKey("x") || !opts.containsKey("z")) {
                SFIToolkit.errorln("Missing required argument(s): y= x= z=");
                return 9;
            }

            String yVar = opts.get("y");
            String[] xVars = opts.get("x").split(" ");
            String[] zVars = opts.get("z").split(" ");
            
            int varCount = Data.getVarCount();
            int N = (int) Data.getObsTotal();
            
            Matrix Y = new Matrix(N, 1);
            Matrix X = new Matrix(N, xVars.length);
            Matrix Z = new Matrix(N, varCount-1);
            Matrix balance = new Matrix(N, 1, 1.0); // default weight of 1.0
            
            // Variable names, skipping the y variable 
            List<String> zNames = new ArrayList<>();
            for (int i = 1; i <= varCount; i++) {
                String varName = Data.getVarName(i);
                if (!varName.equals(yVar)) {
                	zNames.add(varName);
                }
            }
            zNames.removeIf(name -> name == null || name.trim().isEmpty());
            String[] varNames = zNames.toArray(new String[0]);
                 
            if (debug) {
                SFIToolkit.displayln("varCount: " + varCount);
                SFIToolkit.displayln("varNames.length: " + varNames.length);
                SFIToolkit.displayln("Stored variable names:");
                SFIToolkit.displayln("--------------------");
                for(int i = 0; i < varNames.length; i++) {
                    SFIToolkit.displayln("varNames[" + i + "]: " + varNames[i]);
                }
            }

            // Load Y 
            int yIndex = Data.getVarIndex(yVar);
            for (int i = 0; i < N; i++) {
                double val = Data.getNum(yIndex, i + 1);
                if (Data.isValueMissing(val)) {
                    SFIToolkit.errorln("Missing Y at obs " + (i + 1));
                }
                Y.set(i, 0, val);
            }

            // Load X
            for (int j = 0; j < xVars.length; j++) {
                int xIndex = Data.getVarIndex(xVars[j]);
                for (int i = 0; i < N; i++) {
                    double val = Data.getNum(xIndex, i + 1);
                    if (Data.isValueMissing(val)) {
                        SFIToolkit.errorln("Missing X[" + j + "] at obs " + (i + 1));
                    }
                    X.set(i, j, val);
                }
            }

            // Load Z
            for (int j = 0; j < zNames.size(); j++) {
                int zIndex = Data.getVarIndex(zNames.get(j));
                for (int i = 0; i < N; i++) {
                    double val = Data.getNum(zIndex, i + 1);
                    if (Data.isValueMissing(val)) {
                        SFIToolkit.errorln("Missing Z[" + j + "] at obs " + (i + 1));
                    }
                    Z.set(i, j, val);
                }
            }
            int[] variableSearchIndex = new int[zVars.length];
            for (int j = 0; j < zVars.length; j++) {
            	variableSearchIndex[j] = zNames.indexOf(zVars[j]); 
            }                      
            
            if (debug) {
	            // This is supposed to skip the y-variable
	            SFIToolkit.displayln("Z matrix column indices:");
	            SFIToolkit.displayln("----------------");
	            for (int i = 0; i < variableSearchIndex.length; i++) {
	                SFIToolkit.displayln("Column " + (i + 1) + ": " + variableSearchIndex[i]);
	            }
            }
            
            // Optional: discrete variables
            String[] discreteVars = new String[0];
            if (opts.containsKey("discretevars")) {
                discreteVars = opts.get("discretevars").split(" ");
            }            
            
            // Construct LinearSpecification
            LinearStataSpecification spec = new LinearStataSpecification(X, Y, Z, balance);
            spec.setVarNames(varNames);
            spec.setVariableIndicesToSearchOver(variableSearchIndex);
            spec.setDiscreteVariables(zVars, Z, discreteVars, variableSearchIndex);
            
            if (debug) {
            	SFIToolkit.displayln("Discrete variables:");
            	for (int j = 0; j < spec.getDiscreteVector().length; j++) {
            	    SFIToolkit.displayln("Z[" + j + "] = " + varNames[j] + " -> " + spec.getDiscreteVector()[j]);
            	}
            }
            
			// July 15, 2025 attempt to automate value labels for discrete variables
			/*
			// Optional: value labels for discrete variables
			if (opts.containsKey("value_labels")) {
				Map<String, Map<Integer, String>> valueLabels = new HashMap<>();
			    String raw = opts.get("value_labels");
			    String[] varBlocks = raw.split(";");
			    for (String block : varBlocks) {
			        if (block.trim().isEmpty()) continue;
			        String[] varSplit = block.split("=", 2);
			        if (varSplit.length < 2) continue;

			        String varName = varSplit[0].trim();
			        String mappingStr = varSplit[1];
			        Map<Integer, String> labelMap = new HashMap<>();

			        for (String pair : mappingStr.split(",")) {
			            String[] kv = pair.split("=", 2);
			            if (kv.length < 2) continue;
			            try {
			                int val = Integer.parseInt(kv[0].trim());
			                String label = kv[1].trim();
			                labelMap.put(val, label);
			            } catch (NumberFormatException e) {
			                // skip malformed entries
			            }
			        }

			        valueLabels.put(varName, labelMap);
			    }
			    spec.setValueLabels(valueLabels);
			}
			*/
                   
            if (opts.containsKey("numtrees")) {
                int numTrees = Integer.parseInt(opts.get("numtrees"));
                spec.setNumTrees(numTrees);
            }
                      
            // Optional: set stratification column
            if (opts.containsKey("strata")) {
                String stratVars = opts.get("strata");
                if (stratVars != null && !stratVars.isEmpty()) {
                    String[] stratVarArray = stratVars.trim().split("\\s+"); // split by whitespace
                    List<Integer> stratIndices = new ArrayList<>();

                    for (String stratVar : stratVarArray) {
                        boolean found = false;
                        for (int j = 0; j < varNames.length; j++) {
                            if (varNames[j].equals(stratVar)) {
                                stratIndices.add(j);
                                found = true;
                                break;
                            }
                        }
                        if (!found) {
                            SFIToolkit.errorln("Strata variable '" + stratVar + "' not found in dataset.");
                        }
                    }

                    if (!stratIndices.isEmpty()) {
                        int[] indicesArray = stratIndices.stream().mapToInt(Integer::intValue).toArray();
                        spec.setStratificationIndex(indicesArray);
                    }
                }
            }
            
            // Optional: testing for homogeneity
            if (opts.containsKey("testhomogeneity")) {
                boolean doDetectHomogeneity = Boolean.parseBoolean(opts.get("testhomogeneity"));
                spec.setDetectHomogeneity(doDetectHomogeneity);
            }
            
            // Optional: cross-validation flag
            if (opts.containsKey("cv")) {
                boolean doCV = Boolean.parseBoolean(opts.get("cv"));
                spec.setCrossValidation(doCV);
            }
            
            // Optional: CV grid
            if (opts.containsKey("cvgrid_minleaf")) {
                cvGridMinLeaf = Arrays.stream(opts.get("cvgrid_minleaf").split(","))
                        .map(String::trim)
                        .map(Integer::parseInt)
                        .sorted()
                        .collect(Collectors.toList());
            } else {
                cvGridMinLeaf = Collections.singletonList(25);
            }
            spec.setCVGridMinLeaf(cvGridMinLeaf);
			
            if (opts.containsKey("cvgrid_minimp")) {
                cvGridMinImprovement = Arrays.stream(opts.get("cvgrid_minimp").split(","))
                        .map(String::trim)
                        .map(Double::parseDouble)
                        .sorted()
                        .collect(Collectors.toList());
            } else {
                cvGridMinImprovement = Collections.singletonList(0.1);
            }
            spec.setCVGridMinImprovement(cvGridMinImprovement);
			
            if (opts.containsKey("cvgrid_maxdepth")) {
                cvGridMaxDepth = Arrays.stream(opts.get("cvgrid_maxdepth").split(","))
                        .map(String::trim)
                        .map(Integer::parseInt)
                        .sorted(Comparator.reverseOrder())
                        .collect(Collectors.toList());
            } else {
                cvGridMaxDepth = Collections.singletonList(5);
            }
            spec.setCVGridMaxDepth(cvGridMaxDepth);
			
			// Optional: generating betas
			if (opts.containsKey("gen")) {
				String betaPrefix = opts.get("gen");
				spec.setBetaPrefixes(betaPrefix);
			}
			
			// Optional: proportion of observations to estimate tree structure 
			double propstructure = 0.35;
            if (opts.containsKey("propstructure")) {
            	propstructure = Double.parseDouble(opts.get("propstructure"));
            }  
            spec.setProportionObservationsToEstimateTreeStructure(propstructure);
            
			// Random seed
			Random rand = new Random();
			int seed = rand.nextInt();
            if (opts.containsKey("seed")) {
                seed = Integer.parseInt(opts.get("seed"));
            }
			
            // Run the model
            LinearStataMain.execute(spec, seed);

        } catch (Exception e) {
            java.io.StringWriter sw = new java.io.StringWriter();
            java.io.PrintWriter pw = new java.io.PrintWriter(sw);
            e.printStackTrace(pw);
            SFIToolkit.displayln(sw.toString());
            return 9;
        }

        return rc;
    }
}
