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
package examples.logit;

import core.MomentContinuousSplitObj;
import core.DataLens;
import core.MomentSpecification;
import core.SplitContainer;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentContinuousSplitObjLogit extends MomentContinuousSplitObj {

    SplitContainer container; 
    DataLens lens;
    int minCount;
    double minProportion;
    boolean debugVerbose = false;
    MomentSpecification spec;

    public MomentContinuousSplitObjLogit(int indexSplitVariable, DataLens lens, double minProportion, int minCount, LogitMomentSpecification spec) {
        this.indexSplitVariable = indexSplitVariable;
        this.lens = lens;
        this.minCount = minCount;
        this.minProportion = minProportion;
        this.spec = spec;

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
    public double getSSE() {

        leftMSE = 0;
        rightMSE = 0;

        ContainerLogit leftLogit = new ContainerLogit(container.getLeft(), spec.getHomogeneousIndex(), spec.getHomogeneousParameterVector(), false); //This object will compute the beta and MSE for the left split
        ContainerLogit rightLogit = new ContainerLogit(container.getRight(), spec.getHomogeneousIndex(), spec.getHomogeneousParameterVector(), false);
        // System.out.println("Compute left");
        leftLogit.computeBetaAndErrors();
        // System.out.println("Compute right");
        rightLogit.computeBetaAndErrors();
        // System.out.println("Compute out");

        leftMSE = leftLogit.getGoodnessOfFit();
        rightMSE = rightLogit.getGoodnessOfFit();

        if (debugVerbose) {
            System.out.println("MSE = " + (leftMSE + rightMSE) + " n_Left: " + numObsLeft + " n_Right: " + numObsRight + " MSE_Left: " + leftMSE + " MSE_Right: " + rightMSE + " minCount: " + minCount);
        }

        if (getEffectiveNumObsLeft() < minCount) {
            if (debugVerbose) {
                System.out.println("Not enough n_left = " + getEffectiveNumObsLeft() + " ==> Triggered positive infinity");
            }
            leftMSE = Double.POSITIVE_INFINITY;
            // return Double.POSITIVE_INFINITY;
        }
        if (getEffectiveNumObsRight() < minCount) {
            if (debugVerbose) {
                System.out.println("Not enough n_right = " + getEffectiveNumObsRight() + " ==> Triggered positive infinity");
            }
            rightMSE = Double.POSITIVE_INFINITY;
            // return Double.POSITIVE_INFINITY;
        }

        // return (leftMSE + rightMSE) / X.getNumObs();
        return (leftMSE + rightMSE);
    }

    @Override
    public double f_to_minimize(double splitPoint) {
        container = SplitContainer.getContinuousDataSplit(lens, splitPoint, indexSplitVariable); //This returns the data split into each leaf based on splitpoint and index of split variable
        numObsLeft = container.getLeft().getNumObs();
        numObsRight = container.getRight().getNumObs();
        return getSSE();
    }

}
