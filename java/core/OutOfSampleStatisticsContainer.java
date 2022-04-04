/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package core;

/**
 *
 * @author stephen.p.ryan
 */
public class OutOfSampleStatisticsContainer {
    
    private double meanSquaredErrorParameters;
    private double outOfSampleMeasureY;

    public OutOfSampleStatisticsContainer(double meanSquaredErrorParameters, double outOfSampleMeasureY) {
        this.meanSquaredErrorParameters = meanSquaredErrorParameters;
        this.outOfSampleMeasureY = outOfSampleMeasureY;
    }

    public double getMeanSquaredErrorParameters() {
        return meanSquaredErrorParameters;
    }

    public double getOutOfSampleMeasureY() {
        return outOfSampleMeasureY;
    }
    
}
