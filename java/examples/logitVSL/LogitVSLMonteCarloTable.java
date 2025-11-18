/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package examples.logitVSL;

/**
 *
 * @author stephen.p.ryan
 */
public class LogitVSLMonteCarloTable {

    int n;
    double YMSE_unrestricted;
    double YMSE_SD_unrestricted;
    double YMSE_restricted;
    double YMSE_SD_restricted;

    double betaMSE_unrestricted;
    double betaMSE_restricted;
    double betaMSE_SD_unrestricted;
    double betaMSE_SD_restricted;

    double beta1_mean;
    double beta1_SD;
    double beta2_mean;
    double beta2_SD;

    double classificationRate1;
    double classificationRate2;

    public LogitVSLMonteCarloTable(int n, double YMSE_unrestricted, double YMSE_SD_unrestricted, double YMSE_restricted, double YMSE_SD_restricted, double betaMSE_unrestricted, double betaMSE_restricted, double betaMSE_SD_unrestricted, double betaMSE_SD_restricted, double beta1_mean, double beta1_SD, double beta2_mean, double beta2_SD, double classificationRate1, double classificationRate2) {
        this.n = n;
        this.YMSE_unrestricted = YMSE_unrestricted;
        this.YMSE_SD_unrestricted = YMSE_SD_unrestricted;
        this.YMSE_restricted = YMSE_restricted;
        this.YMSE_SD_restricted = YMSE_SD_restricted;
        this.betaMSE_unrestricted = betaMSE_unrestricted;
        this.betaMSE_restricted = betaMSE_restricted;
        this.betaMSE_SD_unrestricted = betaMSE_SD_unrestricted;
        this.betaMSE_SD_restricted = betaMSE_SD_restricted;
        this.beta1_mean = beta1_mean;
        this.beta1_SD = beta1_SD;
        this.beta2_mean = beta2_mean;
        this.beta2_SD = beta2_SD;
        this.classificationRate1 = classificationRate1;
        this.classificationRate2 = classificationRate2;
    }

    @Override
    public String toString() {
        StringBuilder s = new StringBuilder();
        s.append(String.format("%d & $MSE(Y)$ & %10.3g & (%.3g) & %.3g & (%.3g) \\\\ %n ", n, YMSE_unrestricted, YMSE_SD_unrestricted, YMSE_restricted, YMSE_SD_restricted));
        s.append(String.format("   & $MSE(\\beta)$ & %.3g & (%.3g) & %.3g & (%.3g) \\\\ %n ", betaMSE_unrestricted, betaMSE_SD_unrestricted, betaMSE_restricted, betaMSE_SD_restricted));
        s.append(String.format("   & $\\beta_1$ &  &  & %.3g & (%.3g) & %.3g\\%% \\\\ %n ", beta1_mean, beta1_SD, classificationRate1 * 100.0));
        s.append(String.format("   & $\\beta_2$ &  &  & %.3g & (%.3g) & %.3g\\%% \\\\ %n ", beta2_mean, beta2_SD, classificationRate2 * 100.0));
        return s.toString();
    }

}
