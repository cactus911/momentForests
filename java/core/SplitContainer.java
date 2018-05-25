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
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SplitContainer {

    private final Jama.Matrix yLeft;
    private final Jama.Matrix yRight;
    private final Jama.Matrix xLeft;
    private final Jama.Matrix xRight;
    private final Jama.Matrix wLeft;
    private final Jama.Matrix wRight;

    public SplitContainer(Matrix yLeft, Matrix xLeft, Matrix yRight, Matrix xRight) {
        this.yLeft = yLeft;
        this.yRight = yRight; 
        this.xLeft = xLeft;
        this.xRight = xRight;
        wLeft = null;
        wRight = null;
    }

    public SplitContainer(Matrix yLeft, Matrix xLeft, Matrix yRight, Matrix xRight, Matrix wLeft, Matrix wRight) {
        this.yLeft = yLeft;
        this.yRight = yRight;
        this.xLeft = xLeft;
        this.xRight = xRight;
        this.wLeft = wLeft;
        this.wRight = wRight;
    }

    public Matrix getxLeft() {
        return xLeft;
    }

    public Matrix getyLeft() {
        return yLeft;
    }

    public Matrix getxRight() {
        return xRight;
    }

    public Matrix getyRight() {
        return yRight;
    }

    public double getMinimumProportionDataInPartition() {
        double total = (double)xLeft.getRowDimension()+(double)xRight.getRowDimension();
        double proportionLeft = (double)xLeft.getRowDimension() / total;
        // System.out.println(xLeft.getRowDimension()+" "+xRight.getRowDimension()+" "+total+" "+proportionLeft);
        return Math.min(proportionLeft, 1.0-proportionLeft);        
    }
    
    public int getMinimumCountInEachPartition() {
        return Math.min(xLeft.getRowDimension(), xRight.getRowDimension());
    }
    
    /**
     * @return the wLeft
     */
    public Jama.Matrix getwLeft() {
        return wLeft;
    }

    /**
     * @return the wRight
     */
    public Jama.Matrix getwRight() {
        return wRight;
    }
    
    public static SplitContainer getContinuousDataSplit(Jama.Matrix Y, Jama.Matrix X, double splitPoint, int indexSplitVariable) {
        // System.out.println("indexSplitVariable: " + indexSplitVariable + " splitPoint: " + splitPoint);
        int countLeft = 0;
        int countRight = 0;
        int numX = X.getColumnDimension();

        for (int i = 0; i < X.getRowDimension(); i++) {
            if (X.get(i, indexSplitVariable) < splitPoint) {
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
            if (X.get(i, indexSplitVariable) < splitPoint) {
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

//        System.out.println("left:");
//        if (yLeft.getRowDimension() > 10) {
//            pmUtility.prettyPrint(pmUtility.concatMatrix(yLeft, xLeft).getMatrix(0, 10, 0, xLeft.getColumnDimension()));
//        }
//        System.out.println("right:");
//        if (yRight.getRowDimension() > 10) {
//            pmUtility.prettyPrint(pmUtility.concatMatrix(yRight, xRight).getMatrix(0, 10, 0, xLeft.getColumnDimension()));
//        }

        return new SplitContainer(yLeft, xLeft, yRight, xRight);
    }

}
