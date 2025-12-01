/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package core;

import java.awt.Color;
import java.awt.Dimension;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;
import org.jfree.data.statistics.HistogramDataset;

/**
 *
 * @author stephen
 */
public class PlotPDF {
    
    static int numPlots = 0;
    
    static public void plotDistributionUnrestrictedTestStatistics(double[] values) {       
        HistogramDataset dataset = new HistogramDataset();
        
        // Convert ArrayList<Double> to double array
        
        
        // Add to dataset with automatic bin calculation
        int bins = (int) Math.ceil(Math.log(values.length) / Math.log(2) + 1);
        dataset.addSeries("PDF", values, bins);
        JFreeChart histogram = ChartFactory.createHistogram(
            "Probability Distribution Function",
            "Value",
            "Frequency",
            dataset,
            PlotOrientation.VERTICAL,
            true,
            true,
            false
        );
        
        // Customize the plot
        XYPlot plot = histogram.getXYPlot();
        plot.setBackgroundPaint(Color.WHITE);
        plot.setDomainGridlinePaint(Color.LIGHT_GRAY);
        plot.setRangeGridlinePaint(Color.LIGHT_GRAY);
        
        // Create panel and add to frame
        ChartPanel chartPanel = new ChartPanel(histogram);
        chartPanel.setPreferredSize(new Dimension(800, 600));
        JFrame f = new JFrame("PDF");
        f.setBounds(200+25*numPlots,200+25*numPlots,500,500);
        f.getContentPane().add(chartPanel);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.setVisible(true);
        
        numPlots++;
    }
    
}
