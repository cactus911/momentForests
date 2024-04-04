/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package examples.linear;

/**
 *
 * @author stephen.p.ryan
 */
public class LinearMonteCarloTable {

    int n;
    double YMSE_unrestricted;
    double YMSE_SD_unrestricted;
    double YMSE_restricted;
    double YMSE_SD_restricted;

    double betaMSE_unrestricted;
    double betaMSE_restricted;
    double betaMSE_SD_unrestricted;
    double betaMSE_SD_restricted;

    double[] beta_mean;
    double[] beta_SD;
    double[] classificationRate;

    double YMSE_nonparametric;
    double YMSE_parametric;
    double[] mean_beta_parametric;
    double betaMSE_parametric;

    public LinearMonteCarloTable(int n, double YMSE_unrestricted, double YMSE_SD_unrestricted, double YMSE_restricted, double YMSE_SD_restricted,
            double betaMSE_unrestricted, double betaMSE_restricted, double betaMSE_SD_unrestricted, double betaMSE_SD_restricted, double[] beta_mean, double[] beta_SD, double[] classificationRate, 
            double YMSE_nonparametric, double YMSE_parametric, double[] mean_beta_parametric, double betaMSE_parametric) {
        this.n = n;
        this.YMSE_unrestricted = YMSE_unrestricted;
        this.YMSE_SD_unrestricted = YMSE_SD_unrestricted;
        this.YMSE_restricted = YMSE_restricted;
        this.YMSE_SD_restricted = YMSE_SD_restricted;
        this.betaMSE_unrestricted = betaMSE_unrestricted;
        this.betaMSE_restricted = betaMSE_restricted;
        this.betaMSE_SD_unrestricted = betaMSE_SD_unrestricted;
        this.betaMSE_SD_restricted = betaMSE_SD_restricted;
        this.beta_mean = beta_mean;
        this.beta_SD = beta_SD;
        this.classificationRate = classificationRate;
        this.YMSE_nonparametric = YMSE_nonparametric;
        this.YMSE_parametric = YMSE_parametric;
        this.mean_beta_parametric = mean_beta_parametric;
        this.betaMSE_parametric = betaMSE_parametric;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(String.format("%d & $MSE(Y)$ & %.3g & %.3g & %10.3g & (%.3g) & %.3g & (%.3g) \\\\ %n ", n, YMSE_nonparametric, YMSE_parametric, YMSE_unrestricted, YMSE_SD_unrestricted, YMSE_restricted, YMSE_SD_restricted));
        s.append(String.format("   & $MSE(\\beta)$ & & %.3g &  %.3g & (%.3g) & %.3g & (%.3g) \\\\ %n ", betaMSE_parametric, betaMSE_unrestricted, betaMSE_SD_unrestricted, betaMSE_restricted, betaMSE_SD_restricted));
        for (int k = 0; k < beta_SD.length; k++) {
            s.append(String.format("   & $\\beta_%d$ & & %.3g & & & %.3g & (%.3g) & %.3g\\%% \\\\ %n ", k, mean_beta_parametric[k], beta_mean[k], beta_SD[k], classificationRate[k] * 100.0));
        }
        return s.toString();
    }

}
