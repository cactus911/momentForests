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
package examples.linear;

import core.IntegerPartition;
import core.MomentPartitionObj;
import core.DataLens;
import core.SplitContainer;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentPartitionObjLinear extends MomentPartitionObj {

    SplitContainer container;
    
    public MomentPartitionObjLinear(IntegerPartition partition, int indexSplitVariable, DataLens lens) {
        this.partition = partition;
        this.indexSplitVariable = indexSplitVariable;
        this.lens = lens;
        
        initialize();
    }
    
    private void initialize() {
        container = getDataSplit();
        numObsLeft = container.getLeft().getNumObs();
        numObsRight = container.getRight().getNumObs();
    }
    
    @Override
    public int getNumObsLeft() {
        // in the rct context, care about the minimum of count of 0's and 1's in each partition
        int count = 0;
        for (int i = 0; i < numObsLeft; i++) {
            if (container.getLeft().getX(i, 0) == 0) {
                count++;
            }
        }
        return Math.min(count, numObsLeft - count);
    }

    @Override
    public int getNumObsRight() {
        int count = 0;
        for (int i = 0; i < numObsRight; i++) {
            if (container.getRight().getX(i, 0) == 0) {
                count++;
            }
        }
        return Math.min(count, numObsRight - count);
    }
    
    @Override
    public double getMSE() {
        leftMSE = 0;
        rightMSE = 0;

        /**
         * RCT regresses outcome on indicator for treatment; there are no X's
         * There is a constant, and we measure the coefficient on the W
         * So in this implementation just use the first column to get OLS fits and errors, etc.
         */
        ContainerLinear leftRCT = new ContainerLinear(container.getLeft()); //This object will compute the beta and MSE for the left split
        ContainerLinear rightRCT = new ContainerLinear(container.getRight());
        
        leftMSE = leftRCT.getMSE();
        rightMSE = rightRCT.getMSE();
                
        // System.out.println(numObsLeft+" "+numObsRight+" "+leftMSE+" "+rightMSE);
        // return (leftMSE + rightMSE) / X.getNumObs();
        return (leftMSE + rightMSE);
    }
    
    
    
}