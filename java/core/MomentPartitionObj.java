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

import java.util.ArrayList;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public abstract class MomentPartitionObj {

    public IntegerPartition partition;
    public int indexSplitVariable;
    public DataLens lens;

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

        ArrayList<Integer> leftList = partition.getLeft();

        for (int i = 0; i < lens.getNumObs(); i++) {
            if (leftList.contains((int) lens.getX(i, indexSplitVariable))) { // If this obs value of the splitting variable is in the left partition
                countLeft++;
            } else {
                countRight++;
            }
        }

        int[] observationIndicesLeftSplit = new int[countLeft];
        int[] observationIndicesRightSplit = new int[countRight];

//        System.out.println("Entire data:");
//        System.out.println(lens);
        
        // Assign obs indices to either left or right list depending on split
        countRight = 0;
        countLeft = 0;
        for (int i = 0; i < lens.getNumObs(); i++) {
            if (leftList.contains((int) lens.getX(i, indexSplitVariable))) {
                observationIndicesLeftSplit[countLeft] = i;
                countLeft++;
            } else {
                observationIndicesRightSplit[countRight] = i;
                countRight++;
            }
        }
//        
//        for(int i : observationIndicesLeftSplit) {
//            System.out.print(" "+i);
//        }
//        System.out.println();
//        for(int i : observationIndicesRightSplit) {
//            System.out.print(" "+i);
//        }
//        System.out.println();
        
        //Create new data lens for left and right split
        DataLens left = lens.getDataLensSubset(observationIndicesLeftSplit);
        DataLens right = lens.getDataLensSubset(observationIndicesRightSplit);

        String leftVars = "";
        for (int i : partition.getLeft()) {
            leftVars = leftVars + i + " ";
        }
        String rightVars = "";
        for (int i : partition.getRight()) {
            rightVars = rightVars + i + " ";
        }

//        System.out.println("Left " + leftVars);
//        System.out.println(left);
//        System.out.println("Right " + rightVars);
//        System.out.println(right);
//        System.exit(0);

        return new SplitContainer(left, right);
    }

}
