/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import java.util.ArrayList;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class SplitRuleContainer implements Comparable<SplitRuleContainer> {

    ArrayList<Integer> group;
    int indexVariableSplit;
    double continuousVariableSplitPoint;
    boolean isDiscreteSplit;
    boolean left;

    public SplitRuleContainer(ArrayList<Integer> group, int indexVariableSplit) {
        this.group = group;
        this.isDiscreteSplit = true;
        this.indexVariableSplit = indexVariableSplit;
    }

    public SplitRuleContainer(int indexVariableSplit, double continuousVariableSplitPoint, boolean left) {
        this.indexVariableSplit = indexVariableSplit;
        this.continuousVariableSplitPoint = continuousVariableSplitPoint;
        this.isDiscreteSplit = false;
        this.left = left;
    }

    public boolean isLeft() {
        return left;
    }

    public double getContinuousVariableSplitPoint() {
        return continuousVariableSplitPoint;
    }

    public int getIndexVariableSplit() {
        return indexVariableSplit;
    }
    

    public ArrayList<Integer> getGroup() {
        return group;
    }

    public boolean isIsDiscreteSplit() {
        return isDiscreteSplit;
    }

    @Override
    public int compareTo(SplitRuleContainer otherRule) {
        if(otherRule.getIndexVariableSplit()<getIndexVariableSplit()) {
            return 1;
        }
        if(otherRule.getIndexVariableSplit()>getIndexVariableSplit()) {
            return -1;
        }
        if(otherRule.isDiscreteSplit && isDiscreteSplit) {
            if(otherRule.getGroup().get(0) < getGroup().get(0)) {
                return 1;
            }
            return -1;
        }
        if(otherRule.isLeft()) {
            return 1;
        }
        return -1;
    }

}
