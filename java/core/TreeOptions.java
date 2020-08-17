/*
 * The MIT License
 *
 * Copyright 2020 Stephen P. Ryan <stephen.p.ryan@wustl.edu>.
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

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class TreeOptions {

    private double minProportion = 0.0001;
    private int minCount = 1;
    private double minMSEImprovement = 0.2;
    private int maxDepth = 100;

    public TreeOptions() {
    }

    public TreeOptions(double minProportion, int minCount, double minMSEImprovement, int maxDepth) {
        this.maxDepth = maxDepth;
        this.minMSEImprovement = minMSEImprovement;
        this.minCount = minCount;
        this.minProportion = minProportion;
    }

    /**
     * @return the minProportion
     */
    public double getMinProportion() {
        return minProportion;
    }

    /**
     * @param minProportion the minProportion to set
     */
    public void setMinProportion(double minProportion) {
        this.minProportion = minProportion;
    }

    /**
     * @return the minCount
     */
    public int getMinCount() {
        return minCount;
    }

    /**
     * @param minCount the minCount to set
     */
    public void setMinCount(int minCount) {
        this.minCount = minCount;
    }

    /**
     * @return the minMSEImprovement
     */
    public double getMinMSEImprovement() {
        return minMSEImprovement;
    }

    /**
     * @param minMSEImprovement the minMSEImprovement to set
     */
    public void setMinMSEImprovement(double minMSEImprovement) {
        this.minMSEImprovement = minMSEImprovement;
    }

    /**
     * @return the maxDepth
     */
    public int getMaxDepth() {
        return maxDepth;
    }

    /**
     * @param maxDepth the maxDepth to set
     */
    public void setMaxDepth(int maxDepth) {
        this.maxDepth = maxDepth;
    }

}
