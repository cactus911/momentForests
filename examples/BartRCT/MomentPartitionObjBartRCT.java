/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package examples.BartRCT;

import Jama.Matrix;
import core.IntegerPartition;
import core.MomentPartitionObj;
import core.SplitContainer;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class MomentPartitionObjBartRCT extends MomentPartitionObj {

    SplitContainer container;
    
    public MomentPartitionObjBartRCT(IntegerPartition partition, int indexSplitVariable, Matrix X, Matrix Y) {
        this.partition = partition;
        this.indexSplitVariable = indexSplitVariable;
        this.X = X;
        this.Y = Y;
        
        initialize();
    }
    
    private void initialize() {
        container = getDataSplit();
        numObsLeft = container.getyLeft().getRowDimension();
        numObsRight = container.getyRight().getRowDimension();
    }
    
    @Override
    public int getNumObsLeft() {
        // in the rct context, care about the minimum of count of 0's and 1's in each partition
        int count = 0;
        for (int i = 0; i < numObsLeft; i++) {
            if (container.getxLeft().get(i, 0) == 0) {
                count++;
            }
        }
        return Math.min(count, numObsLeft - count);
    }

    @Override
    public int getNumObsRight() {
        int count = 0;
        for (int i = 0; i < numObsRight; i++) {
            if (container.getxRight().get(i, 0) == 0) {
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
        ContainerRCT leftRCT = new ContainerRCT(container.getxLeft(), container.getyLeft());
        ContainerRCT rightRCT = new ContainerRCT(container.getxRight(), container.getyRight());
        
        leftMSE = leftRCT.getMSE();
        rightMSE = rightRCT.getMSE();
        
        // System.out.println(numObsLeft+" "+numObsRight+" "+leftMSE+" "+rightMSE);
        // return (leftMSE + rightMSE) / X.getRowDimension();
        return (leftMSE + rightMSE);
    }
    
    
    
    
    
}
