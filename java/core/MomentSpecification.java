/*
 * The MIT License
 *
 * Copyright 2018 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
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
package core;

import Jama.Matrix;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public interface MomentSpecification {

    public Double getPredictedY(Jama.Matrix xi, Jama.Matrix beta);

    public int[] getVariableIndicesToSearchOver();

    public MomentContinuousSplitObj getFminObjective(DataLens lens, int k, double minProportionEachPartition, int minCountEachPartition);

    public MomentPartitionObj getMomentPartitionObj(DataLens lens, int k, IntegerPartition get);

    public ContainerMoment computeOptimalBeta(DataLens lens);

    //  public void generateData(int numObs, Random rng, boolean addNoise);
    public Jama.Matrix getY(boolean residualizeY);

    public Jama.Matrix getX();

    public Jama.Matrix getZ();

    public Jama.Matrix getBalancingVector();

    public int numberoftrees();

    public Boolean[] getDiscreteVector();

    /**
     * Return the true \beta at a given vector z_i
     * 
     * @param zi Observable vector that determines the parameters.
     * @return 
     */
    public Matrix getBetaTruth(Matrix zi);

    public DataLens getOutOfSampleXYZ(int numObsOutOfSample, long rngSeed);

    public void loadData(long rngSeed);
    
    public String getVariableName(int variableIndex);

    public String getFixedEffectName(int variableIndex, int fixedEffectIndex);

    public String formatTreeLeafOutput(Jama.Matrix beta, Jama.Matrix variance);
    
    public void setHomogeneousParameter(int parameterIndex, double value);
    public double getHomogeneousParameter(int parameterIndex);

    public double getHomogeneousComponent(Jama.Matrix xi);

    public void resetHomogeneityIndex();

    public void setHomogeneousIndex(Integer i);
    public boolean[] getHomogeneousIndex();
    public Jama.Matrix getHomogeneousParameterVector();
    
    public Jama.Matrix residualizeX(Jama.Matrix Xp);

}
