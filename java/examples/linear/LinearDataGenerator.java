/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package examples.linear;

import JSci.maths.statistics.NormalDistribution;
import Jama.Matrix;
import core.MomentSpecification;
import java.util.Random;
import utility.pmUtility;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class LinearDataGenerator {

    private Jama.Matrix X;
    private Jama.Matrix Y; 
    private Jama.Matrix Z;

    public LinearDataGenerator(int numObs, MomentSpecification mySpecification, long randSeed, int dimensionX) {
        X = new Jama.Matrix(numObs, dimensionX);
        Z = new Jama.Matrix(numObs, mySpecification.getDiscreteVector().length);
        Y = new Jama.Matrix(numObs, 1);
        Random rng = new Random(randSeed);
        NormalDistribution normal = new NormalDistribution();
        for (int i = 0; i < numObs; i++) {
            // X.set(i, 0, normal.inverse(rng.nextDouble()));
            X.set(i, 0, 1.0);
            for(int k=1;k<dimensionX;k++) {
                X.set(i, k, Math.pow(normal.inverse(rng.nextDouble()), 2));
            }

            Z.set(i, 0, normal.inverse(rng.nextDouble()));
            // Z.set(i, 0, -6.0 + 10.0 * rng.nextDouble()); // uniform [-6,4]
            Z.set(i, 1, rng.nextDouble());

            double draw = rng.nextDouble();
            if (draw < 0.3) {
                Z.set(i, 2, 1);
            } else if (draw > 0.7) {
                Z.set(i, 2, 2);
            }
            // Z.set(i, 0, X.get(i, 0));
            // Z.set(i, 1, X.get(i, 1));
            Jama.Matrix beta = mySpecification.getBetaTruth(Z.getMatrix(i, i, 0, Z.getColumnDimension() - 1), rng); // Z1 and Z2 to compute beta
//                pmUtility.prettyPrintVector(beta);
            Jama.Matrix subX = X.getMatrix(i, i, 0, dimensionX-1);
            // pmUtility.prettyPrint(subX);
            Y.set(i, 0, (subX.times(beta)).get(0, 0) + 5.0*normal.inverse(rng.nextDouble()));
        }
        
//        pmUtility.prettyPrintVector(pmUtility.OLS(X,Y,false));
//        System.exit(0);
    }

    public Matrix getX() {
        return X;
    }

    public Matrix getY() {
        return Y;
    }

    public Matrix getZ() {
        return Z;
    }

}
