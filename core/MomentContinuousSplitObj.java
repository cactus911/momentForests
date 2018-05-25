/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import optimization.Fmin_methods;

/**
 *
 * @author Stephen P. Ryan
 */
public abstract class MomentContinuousSplitObj implements Fmin_methods {

    public int indexSplitVariable;
    
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

}
