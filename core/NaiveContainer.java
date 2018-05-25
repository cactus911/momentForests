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
public class NaiveContainer {

    double MSPE_OLS;
    int countSignificant;
    
    public NaiveContainer(int countSignificant, double MSPE_OLS) {
        this.countSignificant = countSignificant;
        this.MSPE_OLS = MSPE_OLS;
    }

    public int getCountSignificant() {
        return countSignificant;
    }

    public double getMSPE_OLS() {
        return MSPE_OLS;
    }
    
}

