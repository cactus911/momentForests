/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
class PValue implements Comparable<PValue> {
    
    int k;
    double p;

    /**
     * Container class for p-value and index of hypothesis
     * 
     * @param k Index of parameter we are testing
     * @param p Associated p-value for homogeneity test
     */
    public PValue(int k, double p) {
        this.k = k;
        this.p = p;
    }

    public int getK() {
        return k;
    }

    public double getP() {
        return p;
    }

    @Override
    public int compareTo(PValue pComp) {
        if(pComp.getP()<getP()) {
            return 1;
        }
        if(pComp.getP()==getP()) {
            return 0;
        }
        return -1;
    }

    @Override
    public String toString() {
        return "parameter index: "+k+" p-value: "+p;
    }
    
    
    
}
