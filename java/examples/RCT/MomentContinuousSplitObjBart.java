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
package examples.RCT;

import core.MomentContinuousSplitObj;
import core.SplitContainer;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentContinuousSplitObjBart extends MomentContinuousSplitObj {

    SplitContainer container;
    Jama.Matrix X;
    Jama.Matrix Y;
    int minCount;
    double minProportion;

    public MomentContinuousSplitObjBart(int indexSplitVariable, Jama.Matrix X, Jama.Matrix Y, double minProportion, int minCount) {
        this.indexSplitVariable = indexSplitVariable;
        this.X = X;
        this.Y = Y;
        this.minCount = minCount;
        this.minProportion = minProportion;

        // System.out.println("MomentContinuousSplitObjBart.java:26 -> min/max x_1 = "+pmUtility.min(X, 1)+" "+pmUtility.max(X,1));
    }

    @Override
    public int getNumObsLeft() {
        // in the rct context, care about the minimum of count of 0's and 1's in each partition
        int count = 0;
        for (int i = 0; i < container.getxLeft().getRowDimension(); i++) {
            if (container.getxLeft().get(i, 0) == 0) {
                count++;
            }
        }
        return Math.min(count, container.getxLeft().getRowDimension() - count);
    }

    @Override
    public int getNumObsRight() {
        int count = 0;
        for (int i = 0; i < container.getxRight().getRowDimension(); i++) {
            if (container.getxRight().get(i, 0) == 0) {
                count++;
            }
        }
        return Math.min(count, container.getxRight().getRowDimension() - count);
    }

    @Override
    public double getMSE() {

        leftMSE = 0;
        rightMSE = 0;

        /**
         * RCT regresses outcome on indicator for treatment; there are no X's
         * There is a constant, and we measure the coefficient on the W So in
         * this implementation just use the first column to get OLS fits and
         * errors, etc.
         */
        ContainerRCT leftRCT = new ContainerRCT(container.getxLeft(), container.getyLeft());
        ContainerRCT rightRCT = new ContainerRCT(container.getxRight(), container.getyRight());

        leftMSE = leftRCT.getMSE();
        rightMSE = rightRCT.getMSE();

        if (getNumObsLeft() < minCount || getNumObsRight() < minCount) {
            return Double.POSITIVE_INFINITY;
        }

        // System.out.println(numObsLeft+" "+numObsRight+" "+leftMSE+" "+rightMSE);
        // return (leftMSE + rightMSE) / X.getRowDimension();
        return (leftMSE + rightMSE);
    }

    @Override
    public double f_to_minimize(double splitPoint) {
        container = SplitContainer.getContinuousDataSplit(Y, X, splitPoint, indexSplitVariable);
        numObsLeft = container.getyLeft().getRowDimension();
        numObsRight = container.getyRight().getRowDimension();
        return getMSE();
    }

}
