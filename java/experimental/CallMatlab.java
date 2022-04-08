/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package experimental;

import JSci.maths.statistics.NormalDistribution;
import com.mathworks.engine.MatlabEngine;
import utility.pmUtility;

/**
 *
 * @author stephen.p.ryan
 */
public class CallMatlab {

    public static void main(String[] args) {
        CallMatlab test = new CallMatlab();
    }

    public CallMatlab() {
        try {
            // MatlabEngine eng = MatlabEngine.startMatlab();

            String[] engines = MatlabEngine.findMatlab();
            for(String s : engines) {
                System.out.println(s);
            }
            MatlabEngine eng = MatlabEngine.connectMatlab(engines[0]);
            // Execute command on shared MATLAB session
            // eng.eval("plot(1:10); print('myPlot','-djpeg')");
            // eng.close();

            int numObs = 500;
            Jama.Matrix X = new Jama.Matrix(numObs, 5);
            for (int i = 0; i < X.getRowDimension(); i++) {
                X.set(i, 0, 1);
                X.set(i, 1, Math.random());
                X.set(i, 2, Math.pow(Math.random(), 2));
                X.set(i, 3, Math.exp(Math.random()));
                X.set(i, 4, 2.0 * Math.random() - 1.0);
            }
            Jama.Matrix beta = new Jama.Matrix(5, 1);
            beta.set(0, 0, 0.1);
            beta.set(1, 0, 0.1);
            beta.set(2, 0, 0.1);
            beta.set(3, 0, 0.2);
            beta.set(4, 0, 0.5);
            Jama.Matrix Y = X.times(beta);
            NormalDistribution normal = new NormalDistribution();
            for (int i = 0; i < Y.getRowDimension(); i++) {
                Y.set(i, 0, Y.get(i, 0) + 0.1 * normal.inverse(Math.random()));
            }

            eng.putVariable("X", X.getArray());

            // this is so cool!
            eng.putVariable("Y", Y.getArray());
            eng.eval("Aeq=ones(1,5)");
            eng.eval("beq=1"); // parameters sum to one
            eng.eval("lb=zeros(5,1)");
            eng.eval("ub=ones(5,1)");

            eng.eval("beta=lsqlin(X,Y,[],[],Aeq,beq,lb,ub)");

            double[] bhat = eng.getVariable("beta");
            pmUtility.prettyPrint(new Jama.Matrix(bhat, 1));

            eng.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
