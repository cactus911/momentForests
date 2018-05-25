/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import Jama.Matrix;
import java.util.Random;

/**
 *
 * @author Stephen P. Ryan
 */
public interface MomentSpecification {

    public Double getPredictedY(Matrix xi, Jama.Matrix beta);
    
    public int[] getVariableIndicesToSearchOver();

    public MomentContinuousSplitObj getFminObjective(Matrix nodeY, Matrix nodeX, int k, double minProportionEachPartition, int minCountEachPartition);

    public MomentPartitionObj getMomentPartitionObj(Matrix nodeX, Matrix nodeY, int k, IntegerPartition get);

    public ContainerMoment computeOptimalBeta(Matrix nodeY, Matrix nodeX);

    public void generateData(int numObs, Random rng, boolean addNoise);

    public Matrix getY();

    public Matrix getX();

    public Boolean[] getDiscreteVector();

    public Matrix getBetaTruth(Matrix xi);
    
    public Jama.Matrix getOutOfSampleX();
    
    public NaiveContainer computeNaiveStatistics();
    
    public void loadData();
    
    public String getVariableName(int variableIndex);
    
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex);

    public String formatTreeLeafOutput(Jama.Matrix beta, Jama.Matrix variance);

}
