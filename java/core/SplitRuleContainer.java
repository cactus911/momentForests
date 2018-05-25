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
