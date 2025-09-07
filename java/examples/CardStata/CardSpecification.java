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
package examples.CardStata;

// import JSci.maths.statistics.NormalDistribution;
import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.ContainerMoment;
import core.DataLens;
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentPartitionObj;
import core.MomentSpecification;
import java.util.*;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>, Nathan Jiang <jiang.n@wustl.edu>
 */
public class CardSpecification implements MomentSpecification {

    Jama.Matrix X;
    Jama.Matrix Y;
    Jama.Matrix Z;
    Jama.Matrix balancingVector; // is treatment status in the RCT setting
    int numObs;
    int numtrees;
    double proportionObservationsToEstimateTreeStructure;
    String[] varNames;
    int[] variableSearchIndex; 
    Boolean[] DiscreteVariables; 
    // July 15, 2015 attempt to automate value labels for discrete variables
    //Map<String, Map<Integer, String>> valueLabels = new HashMap<>();
    boolean failedEstimator = false;
    int numTrees = 100;
    boolean crossValidation = true;
    List<Integer> cvGridMinLeaf;
    List<Double> cvGridMinImprovement;
    List<Integer> cvGridMaxDepth;
    int stratificationIndex = -1;   
    String betaPrefixes = "";
    
    /**
     * We are going to control homogeneous parameters through these variables
     */
    final private boolean[] homogeneityIndex; // = new boolean[X.getColumnDimension()];
    final private Jama.Matrix homogeneousParameterVector; // = new Jama.Matrix(X.getColumnDimension(), 1); // this is a compact vector (only consists of the parameters we are imposing for homogeneity)
 
    public CardSpecification(Jama.Matrix X, Jama.Matrix Y, Jama.Matrix Z, Jama.Matrix balancingVector) {
        this.X = X;
        this.Y = Y;
        this.Z = Z;
        this.balancingVector = balancingVector;
        this.numObs = X.getRowDimension();

        homogeneityIndex = new boolean[X.getColumnDimension()];
        homogeneousParameterVector = new Jama.Matrix(X.getColumnDimension(), 1);
        resetHomogeneityIndex();
        failedEstimator = false;
        
        /*
        variableSearchIndex = new int[Z.getColumnDimension()];
        for (int i = 0; i < Z.getColumnDimension(); i++) {
            variableSearchIndex[i] = i;
        }
        */
    }

    @Override
    public boolean[] getHomogeneousIndex() {
        return homogeneityIndex;
    }

    @Override
    public double getHomogeneousParameter(int parameterIndex) {
        return homogeneousParameterVector.get(parameterIndex, 0);
    }

    @Override
    public int getNumMoments() {
        return X.getColumnDimension();
    }

    @Override
    public void setHomogeneousParameter(int parameterIndex, double value) {
        homogeneousParameterVector.set(parameterIndex, 0, value);
    }

    @Override
    public Matrix getBalancingVector() {
        return balancingVector;
    }

    @Override
    public double getGoodnessOfFit(double yi, Matrix xi, Matrix beta) {
        double fit = (xi.times(beta)).get(0, 0);
        double error = fit - yi;
        // System.out.println(fit+" "+yi);
        return error * error;
    }

    @Override
    public Double getPredictedY(Matrix xi, Jama.Matrix beta, Random rng) {
        /**
         * This may have to be adjusted when we impose homogeneity, depending on
         * what this is used for
         */
        // pmUtility.prettyPrint(xi);
        if (beta != null) {
            double yhat = (xi.times(beta)).get(0, 0);
            return yhat;
        }
        return null;
    }
    
    @Override
    public MomentContinuousSplitObj getFminObjective(DataLens lens, int indexSplitVariable, double minProportionEachPartition, int minCountEachPartition) {
        return new MomentContinuousSplitObjLinear(indexSplitVariable, lens, minProportionEachPartition, minCountEachPartition, this);
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjLinear(partition, indexSplitVariable, lens, this);
    }

    @Override
    public ContainerMoment computeOptimalBeta(DataLens lens, boolean allParametersHomogeneous) {
        ContainerCard l = new ContainerCard(lens, homogeneityIndex, homogeneousParameterVector, allParametersHomogeneous, this);
        l.computeBetaAndErrors();
        failedEstimator = l.didEstimatorFail();
        return l;
    }

    @Override
    public boolean didEstimatorFail() {
        return failedEstimator;
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
    public Matrix getZ() {
        return Z;
    }

    @Override
    public int numberoftrees() {
        return numtrees;
    }

    @Override
    public Matrix getBetaTruth(Matrix zi, Random rng) {
        // we don't know, this shouldn't be called in a real application
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }

    @Override
    public DataLens getOutOfSampleXYZ(int numObsOutOfSample, long rngSeed) {
        throw new UnsupportedOperationException("Not supported yet."); // Generated from nbfs://nbhost/SystemFileSystem/Templates/Classes/Code/GeneratedMethodBody
    }
    
    @Override
    public void loadData(long rngSeed) {
        // No-op: Data is already loaded in Stata and passed via constructor
    }
    
    @Override
    public String getVariableName(int variableIndex) {
        return varNames[variableIndex];
    }
    
    public void setVarNames(String[] varNames) {
        this.varNames = varNames;
    }
    
    @Override
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex) {
    	
    	/*
    	// Hard-coded version
        if (variableIndex == 7) {
            String[] regionNames = {"New England", "Mid Atlantic", "East North Central", "West North Central", "South Atlantic", "East South Central", "West South Central", "Mountain", "Pacific"};
            return regionNames[fixedEffectIndex - 1];
        }
        */
    	
    	/*
    	// July 15, 2025 attempt to automate
        if (variableSearchIndex == null || varNames == null || variableIndex >= variableSearchIndex.length)
            return String.valueOf(fixedEffectIndex);

        int originalIndex = variableSearchIndex[variableIndex];
        String varName = varNames[originalIndex];

        if (valueLabels.containsKey(varName)) {
            return valueLabels.get(varName).getOrDefault(fixedEffectIndex, String.valueOf(fixedEffectIndex));
        }

        return String.valueOf(fixedEffectIndex);
        */
    	
        return " " + fixedEffectIndex;
    }
    
    @Override
    public int[] getVariableIndicesToSearchOver() {
        return variableSearchIndex;
    }

    public void setVariableIndicesToSearchOver(int[] variableSearchIndex) {
    	this.variableSearchIndex = variableSearchIndex;
    }    
    
    @Override
    public Boolean[] getDiscreteVector() {
        return DiscreteVariables;
    }
    
    public void setDiscreteVariables(String[] zVars, Jama.Matrix Z, String[] discreteVarNames, int[] variableSearchIndex) {
        int numZVars = zVars.length;

        DiscreteVariables = new Boolean[Z.getColumnDimension()];  
        Arrays.fill(DiscreteVariables, false);  

        // Store variable names
        Set<String> isDiscrete = new HashSet<>();
        for (String name : discreteVarNames) {
        	isDiscrete.add(name);
        }
        
        // Check if variable name matches with list of discrete variables
        for (int j = 0; j < numZVars; j++) {
            String varName = zVars[j];
            int colIndex = variableSearchIndex[j];

            if (isDiscrete.contains(varName)) {
                DiscreteVariables[colIndex] = true;
                continue;
            }
        }
    }
    
    // July 15, 2015 attempt to automate value labels for discrete variables
    /*
    public void setValueLabels(Map<String, Map<Integer, String>> valueLabels) {
        this.valueLabels = valueLabels;
    }
    */
    
    @Override
    public String formatTreeLeafOutput(Matrix beta, Matrix variance) {
        if (beta == null) {
            return "null (shouldn't be here!)";
        }
        // double b = beta.get(0, 0);
//        double se = Math.sqrt(variance.get(0, 0));
//        String stars = "";
//        NormalDistribution normal = new NormalDistribution(0, 1);
//        if (Math.abs(b / se) > Math.abs(normal.inverse(0.90))) {
//            stars = "*";
//        }
//        if (Math.abs(b / se) > Math.abs(normal.inverse(0.95))) {
//            stars = "**";
//        }
//        if (Math.abs(b / se) > Math.abs(normal.inverse(0.99))) {
//            stars = "***";
//        }
//        return String.format("%.2f (%.2f) %s", b, se, stars);
        return pmUtility.stringPrettyPrintVector(beta);
    }

    /**
     * @return the homogeneityIndex
     */
    public boolean[] getHomogeneityIndex() {
        return homogeneityIndex;
    }

    @Override
    public void resetHomogeneityIndex() {
        for (int i = 0; i < homogeneityIndex.length; i++) {
            homogeneityIndex[i] = false;
        }
    }

    @Override
    public void setHomogeneousIndex(Integer i) {
        homogeneityIndex[i] = true;
    }

    /**
     * @return the homogeneousParameterVector
     */
    public Jama.Matrix getHomogeneousParameterVector() {
        return homogeneousParameterVector;
    }

    @Override
    public ContainerMoment getContainerMoment(DataLens lens) {
        return new ContainerCard(lens, homogeneityIndex, homogeneousParameterVector, false, this);
    }

    @Override
    public int getNumParams() {
        return X.getColumnDimension();
    }
    
    public void setNumTrees(int n) {
        this.numTrees = n;
    }

    public int getNumTrees() {
        return this.numTrees;
    }
    
    public void setProportionObservationsToEstimateTreeStructure(double b) {
        this.proportionObservationsToEstimateTreeStructure = b;
    }

    public double getProportionObservationsToEstimateTreeStructure() {
        return this.proportionObservationsToEstimateTreeStructure;
    }
    
    public void setCrossValidation(boolean b) {
        this.crossValidation = b;
    }

    public boolean doCrossValidation() {
        return this.crossValidation;
    }
    
    public List<Integer> getCVGridMinLeaf() {
        return this.cvGridMinLeaf;
    }
    
    public void setCVGridMinLeaf(List<Integer> grid) {
        this.cvGridMinLeaf = grid;
    }
    
    public List<Double> getCVGridMinImprovement() {
        return this.cvGridMinImprovement;
    }

    public void setCVGridMinImprovement(List<Double> grid) {
        this.cvGridMinImprovement = grid;
    }
    
    public List<Integer> getCVGridMaxDepth() {
        return this.cvGridMaxDepth;
    }
    
    public void setCVGridMaxDepth(List<Integer> grid) {
        this.cvGridMaxDepth = grid;
    }

    public void setStratificationIndex(int index) {
        this.stratificationIndex = index;
    }

    public int getStratificationIndex() {
        return this.stratificationIndex;
    }
    
    
    public void setBetaPrefixes(String varname) {
        this.betaPrefixes = varname;
    }
    
    @Override
    public String getBetaPrefixes() {
        return this.betaPrefixes;
    }
}
	