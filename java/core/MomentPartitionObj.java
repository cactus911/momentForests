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
