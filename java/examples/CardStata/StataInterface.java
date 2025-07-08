package examples.CardStata;

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
	public static int RunCardModel(String[] args) {
		
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
        
        // Default CV grid
        List<Integer> cvGridMinLeaf = Arrays.asList(25, 50, 100, 200);      
        List<Double> cvGridMinImprovement = Arrays.asList(0.1, 0.2, 0.4, 0.8, 1.6); 
        List<Integer> cvGridMaxDepth = Arrays.asList(7, 6, 5, 4, 3, 2, 1);   
        
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

            int N = (int) Data.getObsTotal();
            int Kx = xVars.length;
            int Kz = zVars.length;

            Matrix Y = new Matrix(N, 1);
            Matrix X = new Matrix(N, Kx);
            Matrix Z = new Matrix(N, Kz);
            Matrix balance = new Matrix(N, 1, 1.0); // default weight of 1.0

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
            for (int j = 0; j < Kx; j++) {
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
            for (int j = 0; j < Kz; j++) {
                int zIndex = Data.getVarIndex(zVars[j]);
                for (int i = 0; i < N; i++) {
                    double val = Data.getNum(zIndex, i + 1);
                    if (Data.isValueMissing(val)) {
                        SFIToolkit.errorln("Missing Z[" + j + "] at obs " + (i + 1));
                    }
                    Z.set(i, j, val);
                }
            }
            
            // Variable names 
            String[] varNames = new String[0];
            if (opts.containsKey("varnames")) {
                varNames = opts.get("varnames").split(" ");
            }
            
            // Discrete variables
            String[] discreteVars = new String[0];
            if (opts.containsKey("discretevars")) {
                discreteVars = opts.get("discretevars").split(" ");
            }
            
            // Construct CardSpecification
            CardSpecification spec = new CardSpecification(X, Y, Z, balance);
            spec.setVarNames(varNames);
            spec.setDiscreteVariables(zVars, Z, discreteVars);

            // Optional: set stratification column
            if (opts.containsKey("strata")) {
                String stratVar = opts.get("strata");
                int stratIndex = Data.getVarIndex(stratVar);
                if (stratIndex >= 0) {
                	spec.setStratificationIndex(stratIndex);
                } else {
                	SFIToolkit.errorln("Strata variable '" + stratVar + "' not found in dataset.");
                }
            }

            // Optional: number of trees
            if (opts.containsKey("numtrees")) {
                int numTrees = Integer.parseInt(opts.get("numtrees"));
                spec.setNumTrees(numTrees);
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
                        .collect(Collectors.toList());
			    spec.setCVGridMinLeaf(cvGridMinLeaf);
			}
			
			if (opts.containsKey("cvgrid_minimp")) {
			    cvGridMinImprovement = Arrays.stream(opts.get("cvgrid_minimp").split(","))
                        .map(String::trim)
                        .map(Double::parseDouble)
                        .collect(Collectors.toList());
			    spec.setCVGridMinImprovement(cvGridMinImprovement);
			}
			
			if (opts.containsKey("cvgrid_maxdepth")) {
			    cvGridMaxDepth = Arrays.stream(opts.get("cvgrid_maxdepth").split(","))
                        .map(String::trim)
                        .map(Integer::parseInt)
                        .collect(Collectors.toList());
				spec.setCVGridMaxDepth(cvGridMaxDepth);
			}
			
            // Run the model
            CardMain.execute(spec);

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
