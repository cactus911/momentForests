/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package examples.BartRCT;

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
