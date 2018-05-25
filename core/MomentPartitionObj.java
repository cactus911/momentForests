/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import Jama.Matrix;
import java.util.ArrayList;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public abstract class MomentPartitionObj {

    public IntegerPartition partition;
    public int indexSplitVariable;
    public Matrix X;
    public Matrix Y;

    public boolean verbose = false;
    public double leftMSE;
    public double rightMSE;
    public int numObsLeft;
    public int numObsRight;
    
    public abstract double getMSE();

    public double getRightMSE() {
        return rightMSE;
    }

    public double getLeftMSE() {
        return leftMSE;
    }

    public int getNumObsLeft() {
        return numObsLeft;
    }

    public int getNumObsRight() {
        return numObsRight;
    }

    public SplitContainer getDataSplit() {
        int countLeft = 0;
        int countRight = 0;
        int numX = X.getColumnDimension();

        ArrayList<Integer> leftList = partition.getLeft();

        for (int i = 0; i < X.getRowDimension(); i++) {
            if (leftList.contains((int) X.get(i, indexSplitVariable))) {
                countLeft++;
            } else {
                countRight++;
            }
        }

        Jama.Matrix xLeft = new Jama.Matrix(countLeft, numX);
        Jama.Matrix xRight = new Jama.Matrix(countRight, numX);
        Jama.Matrix yLeft = new Jama.Matrix(countLeft, 1);
        Jama.Matrix yRight = new Jama.Matrix(countRight, 1);
        countRight = 0;
        countLeft = 0;
        for (int i = 0; i < X.getRowDimension(); i++) {
            if (leftList.contains((int) X.get(i, indexSplitVariable))) {
                for (int j = 0; j < numX; j++) {
                    xLeft.set(countLeft, j, X.get(i, j));
                }
                yLeft.set(countLeft, 0, Y.get(i, 0));
                countLeft++;
            } else {
                for (int j = 0; j < numX; j++) {
                    xRight.set(countRight, j, X.get(i, j));
                }
                yRight.set(countRight, 0, Y.get(i, 0));
                countRight++;
            }
        }
        // System.out.println("Index split variable: "+indexSplitVariable);
        //pmUtility.prettyPrint(xLeft);
        //System.exit(0);
        return new SplitContainer(yLeft, xLeft, yRight, xRight);
    }

}
