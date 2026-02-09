/*
 * mcmcFunction.java
 *
 * Created on March 10, 2004, 6:03 PM
 */

package mcmc;

/** This interface defines the methods for the statistical function
 * needed in the LTE estimator method of Hong and Chernozhukov.
 * @author Stephen P. Ryan
 */
public interface mcmcFunction {
    /** This should return the value of the statistical function at the
     * vector of parameters <I>x</I>.
     * @param x Vector of parameters for the statistical function.
     * @return Returns the value of the statistical function at <I>x</I>.
     */    
    public double objectiveFunction(double[] x);
    /** This is a weight or prior probability density that is strictly
     * positive and continous over the parameter space.
     * @param x Vector of parameters to evaluate.
     * @return The probability associated with this vector of parameters.
     */    
    public double pi(double[] x); 
}
