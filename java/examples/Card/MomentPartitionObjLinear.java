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
package examples.Card;

import core.IntegerPartition;
import core.MomentPartitionObj;
import core.DataLens;
import core.MomentSpecification;
import core.SplitContainer;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentPartitionObjLinear extends MomentPartitionObj {

    SplitContainer container;
    MomentSpecification spec;
    
    public MomentPartitionObjLinear(IntegerPartition partition, int indexSplitVariable, DataLens lens, MomentSpecification spec) {
        this.partition = partition;
        this.indexSplitVariable = indexSplitVariable;
        this.lens = lens;
        this.spec = spec;
        
        initialize();
    }
    
    private void initialize() {
        container = getDataSplit();
        numObsLeft = container.getLeft().getNumObs();
        numObsRight = container.getRight().getNumObs();
    }
    
    @Override
    public int getEffectiveNumObsLeft() {
        return numObsLeft;
    }

    @Override
    public int getEffectiveNumObsRight() {
        return numObsRight;
    }
    
    @Override
    public double getSSE() {
        leftMSE = 0;
        rightMSE = 0;

        ContainerCard leftLinear = new ContainerCard(container.getLeft(), spec.getHomogeneousIndex(), spec.getHomogeneousParameterVector(), false); //This object will compute the beta and MSE for the left split
        ContainerCard rightLinear = new ContainerCard(container.getRight(), spec.getHomogeneousIndex(), spec.getHomogeneousParameterVector(), false);
        
        leftLinear.computeBetaAndErrors();
        rightLinear.computeBetaAndErrors();
        
        // pmUtility.prettyPrintVector(leftLinear.getBeta());
        // pmUtility.prettyPrintVector(rightLinear.getBeta());
        
        leftMSE = leftLinear.getGoodnessOfFit();
        rightMSE = rightLinear.getGoodnessOfFit();
                
        // System.out.println(numObsLeft+" "+numObsRight+" "+leftMSE+" "+rightMSE);
        // return (leftMSE + rightMSE) / X.getNumObs();
        return (leftMSE + rightMSE);
    }
    
    
    
}