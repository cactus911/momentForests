/*
 * The MIT License
 *
 * Copyright 2020 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 */
package examples.SimpleRCT;

import com.stata.sfi.*;
import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentPartitionObj;
import core.MomentSpecification;
import core.NaiveContainer;
import java.util.ArrayList;
import utility.pmUtility;
import java.util.TreeSet;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SimpleRCTMomentSpecification implements MomentSpecification {

    Jama.Matrix X;
    Jama.Matrix Y;
    int numObs;
    int numtrees;
    Jama.Matrix CVparameters;

    
    public void SimpleRCTMomentSpecification() {
        // this.numObs = numObs;
    }

    @Override
    public Double getPredictedY(Matrix xi, Jama.Matrix beta) {
        // pmUtility.prettyPrint(xi);
        if (beta != null) {
            double yhat = xi.get(0, 0) * beta.get(0, 0);
            return yhat;
        }
        return null;
    }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        int parsedVariables = Data.getParsedVarCount();
        int[] variableSearchIndex = new int[parsedVariables-2];
        for (int i=0; i<variableSearchIndex.length; i++)
            {
        	variableSearchIndex[i] = i + 1;
            }
        return variableSearchIndex;
    }

    @Override
    public MomentContinuousSplitObj getFminObjective(DataLens lens, int indexSplitVariable, double minProportionEachPartition, int minCountEachPartition) {
        return new MomentContinuousSplitObjRCT(indexSplitVariable, lens, minProportionEachPartition, minCountEachPartition);
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjRCT(partition, indexSplitVariable, lens);
    }

    @Override
    public ContainerMoment computeOptimalBeta(DataLens lens) {
        return new ContainerRCT(lens);
    }

    @Override
    public Matrix getY() {
        return Y;
    }

    @Override
    public Matrix getX() {
        return X;
    }

    @Override
    public Matrix getXoriginal() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
    @Override
    public Matrix cvparameters() {
        return CVparameters;
    }

    @Override
    public int numberoftrees() {
        return numtrees;
    }

    @Override
    public Boolean[] getDiscreteVector() {
        double parameter;
    	long obsEnd = Data.getObsParsedIn2();
        int  numvar = Data.getParsedVarCount();
    	Boolean[] whichVariablesDiscrete = new Boolean[numvar - 1];
        for(int i=0; i < whichVariablesDiscrete.length; i++){
        	parameter = Data.getNum(i+2,obsEnd - 1 - 6); // We need to skip those numtree and CV parameter rows
			if (parameter >= 1) {
				whichVariablesDiscrete[i] = true;
    		} else {
    			whichVariablesDiscrete[i] = false;
    		}
        }
        return whichVariablesDiscrete;
    }

    @Override
    public Matrix getBetaTruth(Matrix xi) {
        Jama.Matrix beta = new Jama.Matrix(1, 1);
        beta.set(0, 0, 0);
        beta.set(0, 0, xi.get(0, 1) + xi.get(0, 1) * (xi.get(0, 2) - 1));
        return beta;
    }

    @Override
    public Matrix getOutOfSampleX() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public NaiveContainer computeNaiveStatistics() {
        // this is to run separate OLS regressions for each treatment
        System.out.println("\nIndependent Outcome Estimates");
        for (int k = 0; k < 10; k++) {
            ArrayList<Jama.Matrix> typeX = new ArrayList<>();
            ArrayList<Double> typeY = new ArrayList<>();
            for (int i = 0; i < X.getRowDimension(); i++) {
                if (X.get(i, 1) == k) {
                    typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                    typeY.add(Y.get(i, 0));
                }
            }
            Jama.Matrix subX = new Jama.Matrix(typeX.size(), 1);
            Jama.Matrix subY = new Jama.Matrix(typeX.size(), 1);
            int countControl = 0;
                int countTreatment = 0;
            for (int i = 0; i < typeX.size(); i++) {
                for (int j = 0; j < 1; j++) {
                    subX.set(i, j, typeX.get(i).get(0, j));
                }
                subY.set(i, 0, typeY.get(i));
                if (subX.get(i, 0) == 0) {
                        countControl++;
                    } else {
                        countTreatment++;
                    }
            }
            // pmUtility.prettyPrint(pmUtility.concatMatrix(subY, subX));
            // Jama.Matrix[] bootOLS = pmUtility.bootstrapOLS(subX, subY, false, 500, 787);

            Jama.Matrix olsBeta = pmUtility.OLSsvd(subX, subY, false);
            Jama.Matrix olsVar = pmUtility.getOLSVariances(subY, subX, false);
            String sig = "";
            if (Math.abs(olsBeta.get(0, 0) / Math.sqrt(olsVar.get(0, 0))) > 1.98) {
                sig = "*";
            }
            // System.out.format("OLS Formula  Group %d: %g (%g) %s %n", k, olsBeta.get(0, 0), Math.sqrt(olsVar.get(0, 0)), sig);
            System.out.format("Bootstrapped Group %d: [%d, %d] %g (%g) %s %n", k, countControl, countTreatment, /* bootOLS[0].get(0, 0), bootOLS[1].get(0, 0),*/ sig);
        }

        boolean computeOracle = true;
        if (computeOracle) {
            // let's run the oracle estimator and see how that compares
            System.out.println("\nOracle Estimator");
            for (int k = 0; k < 4; k++) {
                ArrayList<Jama.Matrix> typeX = new ArrayList<>();
                ArrayList<Double> typeY = new ArrayList<>();
                for (int i = 0; i < X.getRowDimension(); i++) {
                    if (k == 0 && X.get(i, 1) < 4) {
                        typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                        typeY.add(Y.get(i, 0));
                    }
                    if (k == 1 && X.get(i, 1) >= 4 && X.get(i, 1) < 8) {
                        typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                        typeY.add(Y.get(i, 0));
                    }
                    if (k == 2 && X.get(i, 1) == 8) {
                        typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                        typeY.add(Y.get(i, 0));
                    }
                    if (k == 3 && X.get(i, 1) == 9) {
                        typeX.add(X.getMatrix(i, i, 0, X.getColumnDimension() - 1));
                        typeY.add(Y.get(i, 0));
                    }
                }
                Jama.Matrix subX = new Jama.Matrix(typeX.size(), 1);
                Jama.Matrix subY = new Jama.Matrix(typeX.size(), 1);
                int countControl = 0;
                int countTreatment = 0;
                for (int i = 0; i < typeX.size(); i++) {
                    for (int j = 0; j < 1; j++) {
                        subX.set(i, j, typeX.get(i).get(0, j));
                    }
                    subY.set(i, 0, typeY.get(i));
                    if (subX.get(i, 0) == 0) {
                        countControl++;
                    } else {
                        countTreatment++;
                    }
                }
                // pmUtility.prettyPrint(pmUtility.concatMatrix(subY, subX));
                Jama.Matrix olsBeta = pmUtility.OLS(subX, subY, false);
                Jama.Matrix olsVar = pmUtility.getOLSVariances(subY, subX, false);
                // Jama.Matrix[] bootOLS = pmUtility.bootstrapOLS(subX, subY, false, 5000, 787);

                String sig = "";
                if (Math.abs(olsBeta.get(0, 0) / Math.sqrt(olsVar.get(0, 0))) > 1.98) {
                    sig = "*";
                }
                // System.out.format("Group %d: [%d, %d] %g (%g) %s %n", k, countControl, countTreatment, olsBeta.get(0, 0), Math.sqrt(olsVar.get(0, 0)), sig);
                System.out.format("Bootstrapped Group %d: [%d, %d] %g (%g) %s %n", k, countControl, countTreatment, /* bootOLS[0].get(0, 0), bootOLS[1].get(0, 0),*/ sig);
            }
        }

        return null;
    }

    @Override
    public void loadData() {
        
	int varIndex_x, varIndex_y;
	double value_x, value_y;
	String msg_x, msg_y;
	TreeSet<Integer> exclusionTree = new TreeSet<>();    
        
        
	int nVariables = Data.getParsedVarCount(); // Get number of variables in varlist specified to javacall
	long firstObs = Data.getObsParsedIn1(); // Get first observation specified by an in restriction
	long lastObs = Data.getObsParsedIn2(); // Get last observation specified by an in restriction
	long nObs = Data.getObsTotal();
	int nObss = (int) nObs;  	   

	// find out missing values
	int counter = 0;
	for (long obs = firstObs; obs <= lastObs; obs++ ) {
	if (!Data.isParsedIfTrue(obs)) {
		exclusionTree.add(counter);
	   }
	   counter++;
	   }
        
        // should not count the last discrete/continuous parameter row
	// should not count the number of tree parameter row
	// should not count the last 6 CV parameter rows
	Y = new Jama.Matrix(counter - 1 - 1 - 6 - exclusionTree.size(), 1);
	X = new Jama.Matrix(counter - 1 - 1 - 6 - exclusionTree.size(), nVariables-1);
        
        numtrees = (int) Data.getNum(1,nObs - 1 - 5);
        
    	CVparameters = new Jama.Matrix(1,6);
	for (int cv = 5; cv >= 0; cv-- ) {
		CVparameters.set(0,5-cv, Data.getNum(1,nObs - cv));  
	}    
        
        // Loop over y to assign values 
	    for (long obs_y = 1; obs_y <= nObs - 1 - 1 - 6; obs_y++ ) {
                if (!exclusionTree.contains(obs_y-1)) {
	   	    int var_y = 1;
	            // get the real variable index for parsed variable -var-
	            varIndex_y = Data.mapParsedVarIndex(var_y);
	            value_y = Data.getNum(varIndex_y, obs_y);
			   
	 	    // Exit with error
	            if (Data.isValueMissing(value_y)) {
		        msg_y = "{err}missing values encountered" ;
			SFIToolkit.errorln(msg_y);
			}
			
	   	    int obss_y = (int) obs_y ;
	            Y.set(obss_y - 1 , 0 , value_y);
                }
	    } 		
        // Loop over x to assign values 		   
	    for(int var_x = 2; var_x <= nVariables; var_x++) {
	   	for (long obs_x = 1; obs_x<=nObs - 1 - 1 - 6; obs_x++ ) {
	            if (!exclusionTree.contains(obs_x-1)) {
		        // get the real variable index for parsed variable -var-
			varIndex_x = Data.mapParsedVarIndex(var_x);
			value_x = Data.getNum(varIndex_x, obs_x);
			   
		    // Exit with error
		    if (Data.isValueMissing(value_x)) {
			msg_x = "{err}missing values encountered" ;
		        SFIToolkit.errorln(msg_x);
			}

		    int obss_x = (int) obs_x ;
	            X.set(obss_x - 1 , var_x - 2 , value_x);
		   }
	   	}
	   }
        Data.addVarLong("beta_estimated");
        Data.addVarLong("beta_bootstrapped");
    }

    @Override
    public String getVariableName(int variableIndex) {
        if (variableIndex == 0) {
            return "Treatment Indicator";
        }
        if (variableIndex == 1) {
            return "Outcome";
        }
        return "Unknown";
    }

    @Override
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex) {
        return "Group " + fixedEffectIndex;
    }

    @Override
    public String formatTreeLeafOutput(Matrix beta, Matrix variance) {
        if (beta == null) {
            return "null (shouldn't be here!)";
        }
        double b = beta.get(0, 0);
        double se = Math.sqrt(variance.get(0, 0));
        String stars = "";
        NormalDistribution normal = new NormalDistribution(0, 1);
        if (Math.abs(b / se) > Math.abs(normal.inverse(0.90))) {
            stars = "*";
        }
        if (Math.abs(b / se) > Math.abs(normal.inverse(0.95))) {
            stars = "**";
        }
        if (Math.abs(b / se) > Math.abs(normal.inverse(0.99))) {
            stars = "***";
        }
        return String.format("%.2f (%.2f) %s", b, se, stars);
    }

}
