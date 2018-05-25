/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package core;

import Jama.Matrix;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public abstract class ContainerMoment {

    public abstract Matrix getBeta();
    public abstract double getMSE();
    public abstract Jama.Matrix getVariance();
    
}
