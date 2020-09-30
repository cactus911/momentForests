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

import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SplitContainer {

    private final DataLens left;
    private final DataLens right;

    public SplitContainer(DataLens left, DataLens right) {
        this.left = left;
        this.right = right;
    }

    public DataLens getLeft() {
        return left;
    }

    public DataLens getRight() {
        return right;
    }

    public double getMinimumProportionDataInPartition() {
        double total = (double) left.getNumObs() + (double) right.getNumObs();
        double proportionLeft = (double) left.getNumObs() / total;
        // System.out.println(xLeft.getNumObs()+" "+xRight.getNumObs()+" "+total+" "+proportionLeft);
        return Math.min(proportionLeft, 1.0 - proportionLeft);
    }

    public int getMinimumCountInEachPartition() {
        return Math.min(left.getNumObs(), right.getNumObs());
    }

    public static SplitContainer getContinuousDataSplit(DataLens lens, double splitPoint, int indexSplitVariable) {
        // System.out.println("indexSplitVariable: " + indexSplitVariable + " splitPoint: " + splitPoint);
        int countLeft = 0;
        int countRight = 0;

        for (int i = 0; i < lens.getNumObs(); i++) {
            if (lens.getX(i, indexSplitVariable) < splitPoint) {
                countLeft++;
            } else {
                countRight++;
            }
        }

        int[] indicesObservationsForLeftSplit = new int[countLeft];
        int[] indicesObservationsForRightSplit = new int[countRight];

        countRight = 0;
        countLeft = 0;
        for (int i = 0; i < lens.getNumObs(); i++) {
            // System.out.print(lens.getX(i, indexSplitVariable)+" -> ");
            if (lens.getX(i, indexSplitVariable) < splitPoint) {
                indicesObservationsForLeftSplit[countLeft] = lens.dataIndex[i];
                countLeft++;
                // System.out.println("left");
            } else {
                indicesObservationsForRightSplit[countRight] = lens.dataIndex[i];
                countRight++;
                // System.out.println("right");
            }
        }
        
        DataLens left = new DataLens(lens, indicesObservationsForLeftSplit);
        DataLens right = new DataLens(lens, indicesObservationsForRightSplit);
        
//        System.out.println("Index to split "+indexSplitVariable+" split point: "+splitPoint);
//        System.out.println("Left:");
//        System.out.println(left);
//        System.out.println("Right:");
//        System.out.println(right);
//        System.exit(0);
        
        return new SplitContainer(left, right);
    }

}
