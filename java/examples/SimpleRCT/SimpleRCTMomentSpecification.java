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

import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.ContainerMoment;
import core.IntegerPartition;
import core.MomentContinuousSplitObj;
import core.MomentPartitionObj;
import core.MomentSpecification;
import core.NaiveContainer;
import examples.RCT.ContainerRCT;
import examples.RCT.MomentContinuousSplitObjRCT;
import examples.RCT.MomentPartitionObjRCT;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SimpleRCTMomentSpecification implements MomentSpecification {
    
    Jama.Matrix X;
    Jama.Matrix Y;    
    
    int[] variableSearchIndex = {1};
    Boolean[] whichVariablesDiscrete = {true, true};

    public SimpleRCTMomentSpecification() {
    }

    @Override
    public Double getPredictedY(Matrix xi, Matrix beta) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int[] getVariableIndicesToSearchOver() {
        return variableSearchIndex;
    }

    @Override
    public MomentContinuousSplitObj getFminObjective(Matrix nodeY, Matrix nodeX, int indexSplitVariable, double minProportionEachPartition, int minCountEachPartition) {
        return new MomentContinuousSplitObjRCT(indexSplitVariable, nodeX, nodeY, minProportionEachPartition, minCountEachPartition);
    }

    @Override
    public MomentPartitionObj getMomentPartitionObj(Matrix nodeX, Matrix nodeY, int indexSplitVariable, IntegerPartition partition) {
        return new MomentPartitionObjRCT(partition, indexSplitVariable, nodeX, nodeY);
    }

    @Override
    public ContainerMoment computeOptimalBeta(Matrix nodeY, Matrix nodeX) {
        return new ContainerRCT(nodeX, nodeY);
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
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public int numberoftrees() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Boolean[] getDiscreteVector() {
        return whichVariablesDiscrete;
    }

    @Override
    public Matrix getBetaTruth(Matrix xi) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public Matrix getOutOfSampleX() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public NaiveContainer computeNaiveStatistics() {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public void loadData() {
        NormalDistribution normal = new NormalDistribution();
        int n = 500;
        Y = new Jama.Matrix(n,1);
        X = new Jama.Matrix(n,2);
        
        for(int i=0;i<n;i++) {
            X.set(i, 0, Math.floor(2*Math.random())); // treatment indicator
            X.set(i, 1, Math.floor(10*Math.random())); // group number
            if(X.get(i,0)==1) {
                Y.set(i, 0, 1.0 + normal.inverse(Math.random()));
            }
        }
    }

    @Override
    public String getVariableName(int variableIndex) {
        if(variableIndex==0) {
            return "Treatment Indicator";
        }
        if(variableIndex==1) {
            return "Group Indicator";
        }
        return "Unknown";
    }

    @Override
    public String getFixedEffectName(int variableIndex, int fixedEffectIndex) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }

    @Override
    public String formatTreeLeafOutput(Matrix beta, Matrix variance) {
        throw new UnsupportedOperationException("Not supported yet."); //To change body of generated methods, choose Tools | Templates.
    }
    
}
