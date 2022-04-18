/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package core;

import java.awt.BorderLayout;
import javax.swing.JFrame;
import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.data.xy.XYSeriesCollection;

/**
 *
 * @author stephen.p.ryan
 */
public class ChartGenerator {

    public static void makeXYScatter(XYSeriesCollection xyc, String title, String xAxisLabel, String yAxisLabel) {
        ChartPanel cp = new ChartPanel(ChartFactory.createScatterPlot(title, xAxisLabel, yAxisLabel, xyc));
        
        NumberAxis axis = (NumberAxis)cp.getChart().getXYPlot().getRangeAxis();
        axis.setAutoRangeIncludesZero(false);
        
        JFrame f = new JFrame("Plot");
        f.setBounds(100,100,500,500);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().setLayout(new BorderLayout());
        f.getContentPane().add(cp, BorderLayout.CENTER);
        f.setVisible(true);
    }
    
    public static void makeXYLine(XYSeriesCollection xyc, String title, String xAxisLabel, String yAxisLabel) {
        ChartPanel cp = new ChartPanel(ChartFactory.createXYLineChart(title, xAxisLabel, yAxisLabel, xyc));
        
        NumberAxis axis = (NumberAxis)cp.getChart().getXYPlot().getRangeAxis();
        axis.setAutoRangeIncludesZero(false);
        
        JFrame f = new JFrame("Plot");
        f.setBounds(100,100,500,500);
        f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        f.getContentPane().setLayout(new BorderLayout());
        f.getContentPane().add(cp, BorderLayout.CENTER);
        f.setVisible(true);
    }
    
}
