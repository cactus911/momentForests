/*
 * gibbsLTE.java
 *
 * Created on May 10, 2004, 2:39 PM
 */

package mcmc;

import JSci.maths.statistics.NormalDistribution;
import JSci.maths.statistics.ProbabilityDistribution;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;

/**
 *
 * @author  Administrator
 */
public class gibbsLTEGeneralizedNoGUI {
    
    static boolean TUNE_VARIANCE = false;
    
    Jama.Matrix chain;
    double acceptance_ratio = 0.234;
    double[] variance;
    double[] savedVar;
    java.text.NumberFormat nf = java.text.NumberFormat.getInstance();
    boolean loadedVar = false;
    
    String[] varLabels;
    private int burnProgress = 0;
    boolean haveLabels = false;
    private double[] lowestPoint;
    private double[] point;
    double lowestValue;
    ArrayList stack = new ArrayList();
    
    // Parallel code parameters
    int numWorkers = 1;  // i should just always parallelize the evaluation of the objective function
    
    public gibbsLTEGeneralizedNoGUI(mcmcFunction function, int draws, int burnPeriod, double[] guess, String varFile, String[] varLabels) {
        this.varLabels = varLabels;
        haveLabels = true;
        try {
            ObjectInputStream in = new ObjectInputStream(new FileInputStream("data/"+varFile+".dat"));
            variance = new double[guess.length];
            savedVar = new double[guess.length];
            for(int i=0;i<guess.length;i++) {
                variance[i] = in.readDouble();
                savedVar[i] = variance[i];
            }
            in.close();
            // System.out.println("Read in variance vector: "+pmUtility.stringPrettyPrint(new Jama.Matrix(variance, 1)));
            loadedVar = true;
        } catch(Exception e) {
            // e.printStackTrace();
            System.out.println("Creating new variance storage file: "+varFile+".dat");
        }
        execute(function, draws, burnPeriod, guess);
        try {
            ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream("data/"+varFile+".dat"));
            for(int i=0;i<guess.length;i++) {
                out.writeDouble(variance[i]);
            }
            out.close();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    public gibbsLTEGeneralizedNoGUI(mcmcFunction function, int draws, int burnPeriod, double[] guess, String varFile) {
        try {
            ObjectInputStream in = new ObjectInputStream(new FileInputStream("data/"+varFile+".dat"));
            variance = new double[guess.length];
            savedVar = new double[guess.length];
            for(int i=0;i<guess.length;i++) {
                variance[i] = in.readDouble();
                savedVar[i] = variance[i];
            }
            in.close();
            // System.out.println("Read in variance vector: "+pmUtility.stringPrettyPrint(new Jama.Matrix(variance, 1)));
            loadedVar = true;
        } catch(Exception e) {
            // e.printStackTrace();
            System.out.println("Creating new variance storage file: "+varFile+".dat");
        }
        execute(function, draws, burnPeriod, guess);
        try {
            ObjectOutputStream out = new ObjectOutputStream(new FileOutputStream("data/"+varFile+".dat"));
            for(int i=0;i<guess.length;i++) {
                out.writeDouble(variance[i]);
            }
            out.close();
        } catch(Exception e) {
            e.printStackTrace();
        }
    }
    
    public gibbsLTEGeneralizedNoGUI(mcmcFunction function, int draws, int burnPeriod, double[] guess) {
        execute(function, draws, burnPeriod, guess);
    }
    
    /** Creates a new instance of gibbsLTE */
    public gibbsLTEGeneralizedNoGUI(mcmcFunction function, int draws, int burnPeriod, double[] guess, String[] varLabels) {
        this.varLabels = varLabels;
        haveLabels = true;
        execute(function, draws, burnPeriod, guess);
    }
    
    private void execute(final mcmcFunction function, int draws, int burnPeriodPassed, double[] guess) {
        final int burnPeriod = burnPeriodPassed;
        lowestPoint = guess;
        lowestValue = function.objectiveFunction(guess);
        
        final int numParams = guess.length;
        point = new double[numParams];
        if(!loadedVar) {
            variance = new double[numParams];
        }
        double count[] = new double[numParams];
        
        for(int i=0;i<numParams;i++) {
            point[i] = guess[i];
            if(!loadedVar) {
                variance[i] = Math.max(1.0, Math.abs(point[i]/10.0));
            }
            count[i] = 0;
        }
        
        double fold = function.objectiveFunction(point);
        
        // long t1 = System.currentTimeMillis();
        
        int ijk = 0;
        fold = function.objectiveFunction(point);
        // System.out.println("burnperiod: "+burnPeriod);
            
        ProbabilityDistribution normal = new NormalDistribution(0, 1);
        while(ijk<burnPeriod) {
            ijk++;
            double[] pointOld = new double[numParams];
            // pmUtility.prettyPrint(new Jama.Matrix(point, 1));
            for(int j=1;j<numParams;j++) {
                double pt = point[j]+variance[j]*normal.inverse(Math.random());
                pointOld[j] = point[j];
                point[j] = pt;
                
                double fnew = function.objectiveFunction(point);
                // System.out.print("old: "+(-fold)+" new: "+(-fnew)+" - "+pmUtility.stringPrettyPrint(new Jama.Matrix(point, 1))+" ");
                double diff = fnew - fold;
                double draw = Math.log(Math.random());
                if(draw<diff) {
                    // System.out.println("accept");
                    fold = fnew;
                } else {
                    // System.out.println("reject");
                    point[j] = pointOld[j];
                }
            }
        }
        
        // select variance using the replace one at a time method
        boolean go = true;
        
        if(!loadedVar) {
            for(int i=0;i<numParams;i++) {
                variance[i] = Math.max(1.0, Math.abs(point[i]/10.0));
            }
        } else {
            for(int i=0;i<numParams;i++) {
                variance[i] = savedVar[i];
            }
        }
        
        // pmUtility.prettyPrint(new Jama.Matrix(variance, 1));
        
        while(go && TUNE_VARIANCE) {
            for(int j=1;j<numParams;j++) {
                double acceptInside = 0;
                double rejectInside = 0;
                for(int i=0;i<100;i++) {
                    double pt = point[j]+variance[j]*normal.inverse(Math.random());
                    double po = point[j];
                    point[j] = pt;
                    double fnew = function.objectiveFunction(point);
                    double diff = fnew - fold;
                    double draw = Math.log(Math.random());
                    if(fnew>lowestValue) {
                        lowestValue = fnew;
                        lowestPoint = point;
                    }
                    if(draw<diff) {
                        acceptInside++;
                        fold = fnew;
                    } else {
                        rejectInside++;
                        point[j] = po;
                    }
                }
                double acceptRate = acceptInside/(acceptInside+rejectInside);
                // System.out.print("variance: "+variance[j]+" rate: "+acceptRate+" reject: "+rejectInside+" accept: "+acceptInside+"\n");
                // System.out.print("vector: "+pmUtility.stringPrettyPrint(new Jama.Matrix(point, 1)));
                double buffer = 0.20;
                double change = 0.15;
                if(acceptRate>(1.0+buffer)*acceptance_ratio) {
                    // if(acceptRate>0.45) {
                    variance[j] *= 1.15; // (1.0+change);
                    count[j] = Math.max(0, count[j]-1);
                }
                if(acceptRate<(1.0-buffer)*acceptance_ratio) {
                    // if(acceptRate<0.10) {
                    variance[j] *= 0.85; // (1.0-change);
                    count[j] = Math.max(0, count[j]-1);
                }
                if(acceptRate==0) {
                    variance[j] *= 0.5;
                }
                if(acceptRate==1) {
                    variance[j] *= 2;
                }
                if(acceptRate<=acceptance_ratio*(1.0+buffer) && acceptRate>=acceptance_ratio*(1.0-buffer)) {
                    // if(acceptRate<0.45 && acceptRate>0.10) {
                    count[j]++;
                }
            }
            go = false;
            for(int i=1;i<numParams;i++) {
                if(!go) {
                    if(count[i]>=3) {
                        go = false;
                    } else {
                        go = true;
                    }
                }
            }
        }
        // System.out.print("Post variance selection: "+pmUtility.stringPrettyPrint(new Jama.Matrix(point, 1))+"\n");
        // pmUtility.prettyPrint(new Jama.Matrix(variance, 1));
        
        // now run chain
        chain = new Jama.Matrix(draws, numParams);
        int accept = 0;
        int reject = 0;
        
        for(int i=0;i<draws;i++) {
            for(int j=1;j<numParams;j++) {
                // double pt = point[j]+cauchy.inverse(Math.random());
                double pt = point[j]+variance[j]*normal.inverse(Math.random());
                double po = point[j];
                point[j] = pt;
                double fnew = function.objectiveFunction(point);
                double diff = fnew - fold;
                double draw = Math.log(Math.random());
                if(fnew>lowestValue) {
                    lowestValue = fnew;
                    lowestPoint = point;
                }
                if(draw<diff) {
                    accept++;
                    fold = fnew;
                    // System.out.print("f: "+fnew+" vector: "+pmUtility.stringPrettyPrint(new Jama.Matrix(point, 1)));
                } else {
                    reject++;
                    point[j] = po;
                }
            }
            for(int k=0;k<numParams;k++) {
                chain.set(i, k, point[k]);
            }
        }
    }
    
    public Jama.Matrix getChain() {
        return chain;
    }
    
    public double[] getLowestPoint() {
        return lowestPoint;
    }

}
