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

import core.MomentContinuousSplitObj;
import core.DataLens;
import core.SplitContainer;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentContinuousSplitObjLinear extends MomentContinuousSplitObj {

    SplitContainer container;
    DataLens lens;
    int minCount;
    double minProportion;
    boolean debugVerbose = false;

    public MomentContinuousSplitObjLinear(int indexSplitVariable, DataLens lens, double minProportion, int minCount) {
        this.indexSplitVariable = indexSplitVariable;
        this.lens = lens;
        this.minCount = minCount;
        this.minProportion = minProportion;

        // System.out.println("MomentContinuousSplitObjBart.java:26 -> min/max x_1 = "+pmUtility.min(X, 1)+" "+pmUtility.max(X,1));
    }

    @Override
    public int getEffectiveNumObsLeft() {
        // in the rct context, care about the minimum of count of 0's and 1's in each partition
        // here, total N is fine
        return numObsLeft;
    }

    @Override
    public int getEffectiveNumObsRight() {
        return numObsRight;
    }

    @Override
    public double getMSE() {

        leftMSE = 0;
        rightMSE = 0;

        ContainerLinear leftRCT = new ContainerLinear(container.getLeft()); //This object will compute the beta and MSE for the left split
        ContainerLinear rightRCT = new ContainerLinear(container.getRight());

        leftMSE = leftRCT.getSSE();
        rightMSE = rightRCT.getSSE();

        if(debugVerbose) {
            System.out.println("MSE = "+(leftMSE+rightMSE)+" n_Left: "+numObsLeft+" n_Right: "+numObsRight+" MSE_Left: "+leftMSE+" MSE_Right: "+rightMSE+" minCount: "+minCount);
        }
        
        if (getEffectiveNumObsLeft() < minCount) {
            if(debugVerbose) {
                System.out.println("Not enough n_left = "+getEffectiveNumObsLeft()+" ==> Triggered positive infinity");
            }
            return Double.POSITIVE_INFINITY;
        } else if (getEffectiveNumObsRight() < minCount) {
            if(debugVerbose) {
                System.out.println("Not enough n_right = "+getEffectiveNumObsRight()+" ==> Triggered positive infinity");
            }
            return Double.POSITIVE_INFINITY;
        }

        
        // return (leftMSE + rightMSE) / X.getNumObs();
        return (leftMSE + rightMSE);
    }

    
    @Override
    public double f_to_minimize(double splitPoint) {
        container = SplitContainer.getContinuousDataSplit(lens, splitPoint, indexSplitVariable); //This returns the data split into each leaf based on splitpoint and index of split variable
        numObsLeft = container.getLeft().getNumObs();
        numObsRight = container.getRight().getNumObs();
        return getMSE();
    }

}