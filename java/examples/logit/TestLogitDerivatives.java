/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package examples.logit;

import utility.pmUtility;

/**
 *
 * @author stephen.p.ryan
 */
public class TestLogitDerivatives {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        TestLogitDerivatives go = new TestLogitDerivatives();
    }

    public TestLogitDerivatives() {
        Jama.Matrix X;
        Jama.Matrix theta;

        int numObs = 1000;
        int numX = 4;
        int numInsideOptions = 3;
        // compute matrix of derivatives numerically and analytically, check that my math is right
        X = Jama.Matrix.random(numInsideOptions, numX);
        theta = new Jama.Matrix(X.getColumnDimension(), 1);
        for (int i = 0; i < X.getColumnDimension(); i++) {
            theta.set(i, 0, 0.1 * i);
        }

        Jama.Matrix G = new Jama.Matrix(numInsideOptions, numX);
        for (int i = 0; i < numObs; i++) {
            G.plusEquals(G_i(theta, X));
        }
        G.timesEquals(1.0 / numObs);
        pmUtility.prettyPrint(G);
    }

    private double g_im(Jama.Matrix theta, Jama.Matrix X) {
        double gim = 0;
        return gim;
    }

    private double g_m(Jama.Matrix theta, Jama.Matrix X) {
        double gi = 0;
        return gi;
    }

    private Jama.Matrix g(Jama.Matrix theta, Jama.Matrix X) {
        Jama.Matrix g = new Jama.Matrix(X.getColumnDimension(), 1);
        return g;
    }

    private double p_ij(Jama.Matrix theta, Jama.Matrix X) {
        double pij = 0;
        return pij;
    }

    private Jama.Matrix p_i(Jama.Matrix theta, Jama.Matrix X) {
        Jama.Matrix pi = new Jama.Matrix(X.getRowDimension(), 1);
        for(int j=0;j<X.getRowDimension();j++) {
            pi.set(j, 0, p_ij(theta,X));
        }
        return pi;
    }

    private Jama.Matrix G_i(Jama.Matrix theta, Jama.Matrix X) {
        Jama.Matrix Gi = new Jama.Matrix(X.getRowDimension(), theta.getRowDimension());

        // pmUtility.prettyPrint(X);
        Jama.Matrix choiceProbabilities = p_i(theta, X);
        pmUtility.prettyPrintVector(choiceProbabilities);
        // moments are each choice * X_choice stacked on top of each other
        for (int choice = 0; choice < X.getRowDimension(); choice++) {
            // for each choice, go with the simple thing here of just 1'e
            for (int l = 0; l < theta.getRowDimension(); l++) {
                double sum = 0;
                for (int k = 0; k < X.getRowDimension(); k++) {
                    sum += choiceProbabilities.get(k, 0) * X.get(k, l);
                }
                Gi.set(choice, l, choiceProbabilities.get(choice, 0) * (X.get(choice, l) - sum));
            }
        }
        return Gi;
    }

}
