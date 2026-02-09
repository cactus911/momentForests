/*
 * pmUtility.java
 *
 * Created on April 24, 2003, 4:58 PM
 */
package utility;

import Jama.Matrix;
import com.google.common.collect.HashMultiset;
import com.google.common.collect.Multiset;
import com.google.common.math.Quantiles;
import com.lowagie.text.pdf.DefaultFontMapper;
import com.lowagie.text.pdf.PdfContentByte;
import com.lowagie.text.pdf.PdfTemplate;
import com.lowagie.text.pdf.PdfWriter;
import java.awt.Graphics2D;
import java.awt.geom.Rectangle2D;
import java.io.FileOutputStream;
import java.util.*;
import java.text.*;
import org.jfree.chart.JFreeChart;

/**
 *
 * @author Stephen P. Ryan <stephen.p.ryan@wustl.edu>
 */
public class pmUtility {

    private static final NumberFormat nf = NumberFormat.getInstance();

    // this should be the newest version
    public static double getMean(ArrayList<Double> list) {
        double mean = 0;
        for (Double d : list) {
            mean += d;
        }
        if (list.size() > 0) {
            mean /= list.size();
        }
        return mean;
    }

    public static double mode(Matrix m, int colIndex) {
        Multiset<Double> freq = HashMultiset.create();
        
        double[][] data = m.getArray();
        for (double[] row : data) {
            freq.add(row[colIndex]);
        }
        
        return Collections.max(freq.entrySet(), 
            Comparator.comparingInt(Multiset.Entry::getCount)
        ).getElement();
    }

    public static Matrix getOLSVariances(Matrix yp, Matrix xp, boolean useIntercept) {
        Jama.Matrix y = yp.copy();
        Jama.Matrix x = xp.copy();

        if (useIntercept) {
            Jama.Matrix temp = new Jama.Matrix(x.getRowDimension(), x.getColumnDimension() + 1);
            for (int i = 0; i < x.getRowDimension(); i++) {
                for (int j = 0; j < x.getColumnDimension(); j++) {
                    temp.set(i, j + 1, x.get(i, j));
                    temp.set(i, 0, 1);
                }
            }
            x = temp.copy();
        }

//        System.out.println("Size y: "+y.getRowDimension()+" by "+y.getColumnDimension());
//        System.out.println("Size x: "+x.getRowDimension()+" by "+x.getColumnDimension());
        Jama.Matrix beta = OLSsvd(x.copy(), y.copy(), false);
        Jama.Matrix fits = getFitsAndErrors(y.copy(), x.copy(), beta.copy(), false);
        // fit in col 0 and error in col 1
        double sigma2 = 0;
        for (int i = 0; i < fits.getRowDimension(); i++) {
            sigma2 += Math.pow(fits.get(i, 2), 2);
        }
//        System.out.println("Raw SSE: "+sigma2);
//        System.out.println("After dividing by n-m: "+(sigma2/(x.getRowDimension()-x.getColumnDimension())));
//        System.out.println("After dividing by y rows: "+(sigma2/y.getRowDimension()));
        Jama.Matrix xprimex = (x.transpose()).times(x);
        Jama.SingularValueDecomposition svd = xprimex.svd();
        Jama.Matrix U = svd.getU();
        Jama.Matrix S = svd.getS();
        Jama.Matrix V = svd.getV();
        Jama.Matrix xpxinv = V.times(S.inverse());
        xpxinv = xpxinv.times(U.transpose());
        Jama.Matrix varianceMatrix = xpxinv.copy();

        sigma2 /= (x.getRowDimension() - x.getColumnDimension());
        varianceMatrix.timesEquals(sigma2);
        // System.out.println("Variance matrix:");
        // pmUtility.prettyPrint(varianceMatrix);
        Jama.Matrix variances = new Jama.Matrix(x.getColumnDimension(), 1);
        for (int i = 0; i < variances.getRowDimension(); i++) {
            variances.set(i, 0, varianceMatrix.get(i, i));
        }
        return variances;
    }

    public static double getVariance(ArrayList<Double> list) {
        double mean = getMean(list);
        double variance = 0;
        for (Double d : list) {
            variance += Math.pow(d - mean, 2);
        }
        if (list.size() > 1) {
            variance /= (list.size() - 1);
        }
        return variance;
    }

    public static Matrix getRowAverages(Matrix m) {
        Jama.Matrix avg = new Jama.Matrix(m.getRowDimension(), 1);
        for (int i = 0; i < m.getRowDimension(); i++) {
            avg.set(i, 0, getRowAverage(m, i));
        }
        return avg;
    }

    public static double getVariance(ArrayList<Double> list, double mean) {
        double variance = 0;
        for (Double d : list) {
            variance += Math.pow(d - mean, 2);
        }
        if (list.size() > 1) {
            variance /= (list.size() - 1);
        }
        return variance;
    }

    public static double meanAbsoluteDeviations(Matrix x, int column, double compareTo) {
        Jama.Matrix d = new Jama.Matrix(x.getRowDimension(), 1);
        for (int i = 0; i < x.getRowDimension(); i++) {
            d.set(i, 0, Math.abs(x.get(i, column) - compareTo));
        }
        return mean(d, 0);
    }

    public static double medianAbsoluteDeviations(Matrix x, int column, double compareTo) {
        Jama.Matrix d = new Jama.Matrix(x.getRowDimension(), 1);
        for (int i = 0; i < x.getRowDimension(); i++) {
            d.set(i, 0, Math.abs(x.get(i, column) - compareTo));
        }
        return median(d, 0);
    }

    private static double getRowAverage(Matrix m, int i) {
        double avg = 0;
        for (int j = 0; j < m.getColumnDimension(); j++) {
            avg += m.get(i, j);
        }
        if (m.getColumnDimension() > 0) {
            avg /= m.getColumnDimension();
        }
        return avg;
    }

    public static double sum(Matrix rdd, int i) {
        double sum = 0;
        for (int r = 0; r < rdd.getRowDimension(); r++) {
            sum += rdd.get(r, i);
        }
        return sum;
    }

    public static double sumSquaredElements(Matrix rdd, int i) {
        double sum = 0;
        for (int r = 0; r < rdd.getRowDimension(); r++) {
            sum += rdd.get(r, i) * rdd.get(r, i);
        }
        return sum;
    }

    public static double sumSquaredElements(Matrix rdd) {
        double sum = 0;
        for (int r = 0; r < rdd.getRowDimension(); r++) {
            sum += rdd.get(r, 0) * rdd.get(r, 0);
        }
        return sum;
    }

    public static Matrix[] bootstrapOLS(Matrix subX, Matrix subY, boolean useIntercept, int numBoots, long seed) {
        int addOne = 0;
        if (useIntercept) {
            addOne++;
        }
        Random rng = new Random(seed);
        Jama.Matrix pointEstimate = OLSsvd(subX, subY, useIntercept);
        Jama.Matrix bootEstimates = new Jama.Matrix(numBoots, subX.getColumnDimension() + addOne);
        for (int b = 0; b < numBoots; b++) {
            long bootSeed = rng.nextLong();
//            System.out.println("bootSeed in pmUtility is: "+bootSeed);
//            Jama.Matrix reX = resample(subX, bootSeed);
//            Jama.Matrix reY = resample(subY, bootSeed);
//            double meanTreatment = 0;
//            double meanControl = 0;
//            int countTreatment = 0;
//            int countControl = 0;
//            for(int i=0;i<reX.getRowDimension();i++) {
//                if(reX.get(i,0)==0) {
//                    meanControl+=reY.get(i,0);
//                    countControl++;
//                } else {
//                    meanTreatment += reY.get(i,0);
//                    countTreatment++;
//                }
//            }
//            meanControl /= countControl+0.0;
//            meanTreatment /= countTreatment+0.0;
            // System.out.println("Arithmetic difference: "+(meanTreatment-meanControl));
            // prettyPrint(concatMatrix(resample(subY, bootSeed), resample(subX, bootSeed)));
            Jama.Matrix ols = OLSsvd(resample(subX, bootSeed), resample(subY, bootSeed), useIntercept);
            for (int i = 0; i < ols.getRowDimension(); i++) {
                bootEstimates.set(b, i, ols.get(i, 0));
            }
        }

        // pmUtility.prettyPrintVector(bootEstimates);
        Jama.Matrix averagedEstimate = new Jama.Matrix(pointEstimate.getRowDimension(), 1);
        Jama.Matrix standardErrors = new Jama.Matrix(pointEstimate.getRowDimension(), 1);
//        System.out.print("[ ");
//        for(int b=0;b<numBoots;b++) {
//            System.out.format("%g ", bootEstimates.get(b,0));
//        }
//        System.out.print("]");
        for (int i = 0; i < standardErrors.getRowDimension(); i++) {
            standardErrors.set(i, 0, pmUtility.standardDeviation(bootEstimates, i));
            averagedEstimate.set(i, 0, pmUtility.mean(bootEstimates, i));
        }
        Jama.Matrix[] result = new Jama.Matrix[2];
        // result[0] = pointEstimate;
        result[0] = averagedEstimate;
        result[1] = standardErrors;
        return result;
    }

    public static Jama.Matrix resample(Matrix m, long seed) {
        Random rng = new Random(seed);
        Jama.Matrix result = m.copy();
        for (int i = 0; i < m.getRowDimension(); i++) {
            int index = rng.nextInt(m.getRowDimension());
            for (int j = 0; j < m.getColumnDimension(); j++) {
                result.set(i, j, m.get(index, j));
            }
        }
        return result;
    }

    public static double computePercentile(Matrix X, int k, int p) {
        ArrayList<Double> XList = new ArrayList<>();
        for (double du : pmUtility.getColumn(X, k).getRowPackedCopy()) {
            XList.add(du);
        }
        return Quantiles.percentiles().index(p).compute(XList);
    }

    /**
     * Creates a new instance of pmUtility
     */
    public pmUtility() {
    }

    public static void printChartToPDF(JFreeChart chart, int width, int height, String fileName) {
        try {
            com.lowagie.text.Document document = new com.lowagie.text.Document(new com.lowagie.text.Rectangle(width, height));
            PdfWriter writer = PdfWriter.getInstance(document, new FileOutputStream(fileName));
            document.addAuthor("Stephen P. Ryan");
            document.open();
            PdfContentByte cb = writer.getDirectContent();
            PdfTemplate tp = cb.createTemplate(width, height);
            Graphics2D g2 = tp.createGraphics(width, height, new DefaultFontMapper());
            Rectangle2D rectangle2D = new Rectangle2D.Double(0, 0, width, height);
            chart.draw(g2, rectangle2D);
            g2.dispose();
            cb.addTemplate(tp, 0, 0);
            document.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public static Jama.Matrix cutoutRow(Jama.Matrix start, int rowToRemove) {
        Jama.Matrix end = new Jama.Matrix(start.getRowDimension() - 1, start.getColumnDimension());
        for (int i = 0; i < rowToRemove; i++) {
            for (int j = 0; j < start.getColumnDimension(); j++) {
                end.set(i, j, start.get(i, j));
            }
        }
        for (int i = rowToRemove + 1; i < start.getRowDimension(); i++) {
            for (int j = 0; j < start.getColumnDimension(); j++) {
                end.set(i - 1, j, start.get(i, j));
            }
        }
        return end;
    }

    public static Jama.Matrix cutoutColumn(Jama.Matrix start, int columnToRemove) {
        Jama.Matrix end = new Jama.Matrix(start.getRowDimension(), start.getColumnDimension() - 1);
        for (int i = 0; i < start.getRowDimension(); i++) {
            for (int j = 0; j < columnToRemove; j++) {
                end.set(i, j, start.get(i, j));
            }
        }
        for (int i = 0; i < start.getRowDimension(); i++) {
            for (int j = columnToRemove + 1; j < start.getColumnDimension(); j++) {
                end.set(i, j - 1, start.get(i, j));
            }
        }
        return end;
    }

    public static double max(Jama.Matrix x, int col) {
        double value = x.get(0, col);
        for (int i = 1; i < x.getRowDimension(); i++) {
            if (x.get(i, col) > value) {
                value = x.get(i, col);
            }
        }
        return value;
    }

    public static double min(Jama.Matrix x, int col) {
        double value = x.get(0, col);
        for (int i = 1; i < x.getRowDimension(); i++) {
            if (x.get(i, col) < value) {
                value = x.get(i, col);
            }
        }
        return value;
    }

    public static Jama.Matrix getColumn(Jama.Matrix x, int index) {
        Jama.Matrix col = new Jama.Matrix(x.getRowDimension(), 1);
        for (int i = 0; i < x.getRowDimension(); i++) {
            col.set(i, 0, x.get(i, index));
        }
        return col;
    }

    public static Jama.Matrix getRow(Jama.Matrix x, int index) {
        Jama.Matrix row = new Jama.Matrix(1, x.getColumnDimension());
        for (int i = 0; i < x.getColumnDimension(); i++) {
            row.set(0, i, x.get(index, i));
        }
        return row;
    }

    public static Jama.Matrix concatMatrix(Jama.Matrix start, Jama.Matrix add) {
        Jama.Matrix result = new Jama.Matrix(start.getRowDimension(), start.getColumnDimension() + add.getColumnDimension());
        for (int i = 0; i < result.getRowDimension(); i++) {
            int counter = 0;
            for (int j = 0; j < start.getColumnDimension(); j++) {
                result.set(i, counter, start.get(i, j));
                counter++;
            }
            for (int j = 0; j < add.getColumnDimension(); j++) {
                result.set(i, counter, add.get(i, j));
                counter++;
            }
        }
        return result;
    }

    public static Jama.Matrix stackMatrix(Jama.Matrix top, Jama.Matrix bottom) {
        Jama.Matrix result = new Jama.Matrix(top.getRowDimension() + bottom.getRowDimension(), top.getColumnDimension());
        for (int j = 0; j < top.getColumnDimension(); j++) {
            for (int i = 0; i < top.getRowDimension(); i++) {
                result.set(i, j, top.get(i, j));
            }
            for (int i = 0; i < bottom.getRowDimension(); i++) {
                result.set(i + top.getRowDimension(), j, bottom.get(i, j));
            }
        }
        return result;
    }

    public static Jama.Matrix concatMatrixLOO(Jama.Matrix start, Jama.Matrix add) {
        Jama.Matrix result = new Jama.Matrix(start.getRowDimension(), start.getColumnDimension() + add.getColumnDimension() - 1);
        for (int i = 0; i < result.getRowDimension(); i++) {
            int counter = 0;
            for (int j = 0; j < start.getColumnDimension(); j++) {
                result.set(i, counter, start.get(i, j));
                counter++;
            }
            for (int j = 0; j < add.getColumnDimension() - 1; j++) {
                result.set(i, counter, add.get(i, j));
                counter++;
            }
        }
        return result;
    }

    public static Jama.Matrix getFitsAndErrors(Jama.Matrix y, Jama.Matrix x, Jama.Matrix beta, boolean useIntercept) {
        Jama.Matrix fits = new Jama.Matrix(y.getRowDimension(), 3);
        double error = 0;
        for (int i = 0; i < y.getRowDimension(); i++) {
            fits.set(i, 0, y.get(i, 0));
            Jama.Matrix xrow = x.getMatrix(i, i, 0, x.getColumnDimension() - 1);
            // pmUtility.prettyPrint(xrow);
            // pmUtility.prettyPrint(beta);
            fits.set(i, 1, xrow.times(beta).get(0, 0));
            fits.set(i, 2, fits.get(i, 1) - fits.get(i, 0));
            error += fits.get(i, 2) * fits.get(i, 2);
        }
        error /= y.getRowDimension();
        // System.out.println("regression sse: " + error);
        return fits;
    }

    public static Jama.Matrix addIntercept(Jama.Matrix x) {
        Jama.Matrix results = new Jama.Matrix(x.getRowDimension(), x.getColumnDimension() + 1);
        for (int i = 0; i < x.getRowDimension(); i++) {
            for (int j = 0; j < x.getColumnDimension(); j++) {
                results.set(i, j + 1, x.get(i, j));
                results.set(i, 0, 1);
            }
        }
        return results;
    }

    public static double mean(Jama.Matrix x, int columnToSort) {
        double mean = 0;
        for (int i = 0; i < x.getRowDimension(); i++) {
            mean += x.get(i, columnToSort);
        }
        mean /= (double) x.getRowDimension();
        return mean;
    }

    public static double IQR(Jama.Matrix x, int columnToSort) {
//        ArrayList<Double> sorter = new ArrayList<Double>();
//        double median = 0;
//        for (int i = 0; i < x.getRowDimension(); i++) {
//            sorter.add(new Double(x.get(i, columnToSort)));
//        }
//        Collections.sort(sorter);
//        double q1 = (sorter.get((int) Math.round(x.getRowDimension() / 4.0))).doubleValue();
//        double q3 = (sorter.get((int) Math.round(3.0 * x.getRowDimension() / 4.0))).doubleValue();
        double q3 = percentile(x, columnToSort, 0.75);
        double q1 = percentile(x, columnToSort, 0.25);
        return q3 - q1;
    }

    public static double median(Jama.Matrix x, int columnToSort) {
        ArrayList<Double> sorter = new ArrayList<Double>();
        double median = 0;
        for (int i = 0; i < x.getRowDimension(); i++) {
            sorter.add(x.get(i, columnToSort));
        }
        Collections.sort(sorter);
        median = (sorter.get((int) Math.round(x.getRowDimension() / 2.0) - 1)).doubleValue();
        return median;
    }

    public static double percentile(Jama.Matrix x, int columnToSort, double alpha) {
        ArrayList<Double> sorter = new ArrayList<Double>();
        double percentile = 0;
        for (int i = 0; i < x.getRowDimension(); i++) {
            sorter.add(x.get(i, columnToSort));
        }
        Collections.sort(sorter);
        int index = (int) Math.floor(x.getRowDimension() * alpha);
        // System.out.println("alpha: "+alpha+" index: "+index);
        percentile = sorter.get(Math.min(index, sorter.size() - 1));
        return percentile;
    }

    public static Jama.Matrix OLS(Jama.Matrix xp, Jama.Matrix yp, boolean useIntercept) {
        Jama.Matrix x = xp.copy();
        Jama.Matrix y = yp.copy();

        if (useIntercept) {
            Jama.Matrix temp = new Jama.Matrix(x.getRowDimension(), x.getColumnDimension() + 1);
            for (int i = 0; i < x.getRowDimension(); i++) {
                for (int j = 0; j < x.getColumnDimension(); j++) {
                    temp.set(i, j + 1, x.get(i, j));
                    temp.set(i, 0, 1);
                }
            }
            x = temp.copy();
        }

        Jama.Matrix beta = new Jama.Matrix(x.getColumnDimension(), 1);

        Jama.Matrix xprimex = (x.transpose()).times(x);
        Jama.Matrix xprimey = (x.transpose()).times(y);

        /*
         Jama.SingularValueDecomposition svd = xprimex.svd();
         Jama.Matrix U = svd.getU();
         Jama.Matrix S = svd.getS();
         Jama.Matrix V = svd.getV();
         Jama.Matrix xpxinv = V.times(S.inverse());
         xpxinv = xpxinv.times(U.transpose());
         */
        Jama.Matrix xpxinv = xprimex.inverse();

        beta = xpxinv.times(xprimey);

        return beta;
    }

    public static Jama.Matrix inverseSVD(Jama.Matrix x) {
        Jama.SingularValueDecomposition svd = x.svd();
        Jama.Matrix U = svd.getU();
        Jama.Matrix S = svd.getS();
        Jama.Matrix V = svd.getV();
        Jama.Matrix xinv = V.times(S.inverse());
        xinv = xinv.times(U.transpose());
        return xinv;
    }

    public static Jama.Matrix OLSsvd(Jama.Matrix xp, Jama.Matrix yp, boolean useIntercept) {
        Jama.Matrix x = xp.copy();
        Jama.Matrix y = yp.copy();

        if (useIntercept) {
            Jama.Matrix temp = new Jama.Matrix(x.getRowDimension(), x.getColumnDimension() + 1);
            for (int i = 0; i < x.getRowDimension(); i++) {
                for (int j = 0; j < x.getColumnDimension(); j++) {
                    temp.set(i, j + 1, x.get(i, j));
                    temp.set(i, 0, 1);
                }
            }
            x = temp.copy();
        }

        Jama.Matrix xprimex = (x.transpose()).times(x);
        Jama.Matrix xprimey = (x.transpose()).times(y);

        Jama.SingularValueDecomposition svd = xprimex.svd();
        Jama.Matrix U = svd.getU();
        Jama.Matrix S = svd.getS();
        Jama.Matrix V = svd.getV();
        Jama.Matrix xpxinv = V.times(S.inverse());
        xpxinv = xpxinv.times(U.transpose());

        Jama.Matrix beta = xpxinv.times(xprimey);

        return beta;
    }

    public static double standardDeviation(Jama.Matrix x) {
        double mean = 0;
        double result = 0;
        // find mean
        for (int i = 0; i < x.getRowDimension(); i++) {
            mean += x.get(i, 0);
        }
        mean /= x.getRowDimension();
        // find deviation
        for (int i = 0; i < x.getRowDimension(); i++) {
            result += (x.get(i, 0) - mean) * (x.get(i, 0) - mean);
        }
        result /= x.getRowDimension();
        return Math.sqrt(result);
    }

    public static double standardDeviation(Jama.Matrix x, int col) {
        double mean = 0;
        double result = 0;
        // find mean
        for (int i = 0; i < x.getRowDimension(); i++) {
            mean += x.get(i, col);
        }
        mean /= x.getRowDimension();
        // find deviation
        for (int i = 0; i < x.getRowDimension(); i++) {
            result += (x.get(i, col) - mean) * (x.get(i, col) - mean);
        }
        result /= x.getRowDimension();
        return Math.sqrt(result);
    }

    public static double variance(Jama.Matrix x, int col) {
        double mean = 0;
        double result = 0;
        // find mean
        for (int i = 0; i < x.getRowDimension(); i++) {
            mean += x.get(i, col);
        }
        mean /= x.getRowDimension();
        // find deviation
        for (int i = 0; i < x.getRowDimension(); i++) {
            result += (x.get(i, col) - mean) * (x.get(i, col) - mean);
        }
        result /= x.getRowDimension();
        return result;
    }

    // find a matrix x in the ArrayList list
    public static int findMatrixIndex(ArrayList list, Jama.Matrix x) {
        int value = 0;
        int size = list.size();
        Jama.Matrix y;
        for (int i = 0; i < size; i++) {
            y = (Jama.Matrix) list.get(i);
            boolean go = true;
            while (go) {
                double vy, vx;
                // pmArrayIndexer.prettyPrintVector(x);
                for (int j = 0; j < x.getRowDimension(); j++) {
                    vy = y.get(j, 0);  // we'll check element by element, starting on the left and going until the end
                    vx = x.get(j, 0);
                    if (vy != vx) {
                        // System.out.println("Computer does not match "+vx+" to "+vy);
                        go = false;
                    }
                }
                if (go) {
                    // System.out.print("x: ");
                    // prettyPrintVector(x);
                    // System.out.print("y: ");
                    // prettyPrintVector(y);
                    return i;
                }
            }
        }
        return -1;
    }

    public static String stringPrettyPrintVector(Jama.Matrix x) {
        nf.setMinimumFractionDigits(6);
        nf.setMaximumFractionDigits(6);
        String s = "";
        s = s.concat("{");

        if (x == null) {
            s = s.concat(" null ");
        } else {
            for (int i = 0; i < x.getRowDimension(); i++) {
                try {
                    s = s.concat(" " + nf.format(x.get(i, 0)));
                } catch (NullPointerException e) {
                    s = s.concat(" null");
                }
            }
        }
        s = s.concat(" }");
        return s;
    }

    public static String stringPrettyPrint(Jama.Matrix x) {
        String s = "";
        if (x != null) {
            nf.setMinimumFractionDigits(6);
            nf.setMaximumFractionDigits(6);

            for (int j = 0; j < x.getRowDimension(); j++) {
                s = s.concat("{");
                for (int i = 0; i < x.getColumnDimension(); i++) {
                    s = s.concat(" " + String.format("% .6f ", x.get(j, i)));
                }
                s = s.concat(" }");
                if (x.getRowDimension() > 1) {
                    s = s.concat("\n");
                }
            }
        } else {
            s = "{ null }";
        }
        return s;
    }

    public static String stringPrettyPrintVectorSci(Jama.Matrix x) {
        String s = "";
        if (x != null) {

            s = s.concat("{");
            for (int j = 0; j < x.getRowDimension(); j++) {
                s = s.concat(" " + String.format("% g ", x.get(j, 0)));
            }
            s = s.concat(" }");
        } else {
            s = "{ null }";
        }
        return s;
    }

    public static String stringPrettyPrintVector(Jama.Matrix x, int desiredFractionDigits) {
        nf.setMinimumFractionDigits(desiredFractionDigits);
        nf.setMaximumFractionDigits(desiredFractionDigits);
        String s = "";
        s = s.concat("{");

        for (int i = 0; i < x.getRowDimension(); i++) {
            s = s.concat(" " + nf.format(x.get(i, 0)));
        }
        s = s.concat(" }");
        return s;
    }

    public static String stringPrettyPrintVector(Jama.Matrix x, int desiredFractionDigits, int min) {
        nf.setMinimumFractionDigits(desiredFractionDigits);
        nf.setMaximumFractionDigits(desiredFractionDigits);
        nf.setMinimumIntegerDigits(min);
        String s = "";
        s = s.concat("{");

        for (int i = 0; i < x.getRowDimension(); i++) {
            s = s.concat(" " + nf.format(x.get(i, 0)));
        }
        s = s.concat(" }");
        return s;
    }

    public static Jama.Matrix sortOnFirstElement(Jama.Matrix x) {
        int size = x.getRowDimension();
        Jama.Matrix temp = new Jama.Matrix(size, 1);
        temp.set(0, 0, x.get(0, 0));
        ArrayList<Double> list = new ArrayList<Double>(size - 1);
        for (int i = 1; i < size; i++) {
            list.add(x.get(i, 0));
        }
        Collections.sort(list);
        for (int i = 0; i < size - 1; i++) {
            double value = (list.get(i)).doubleValue();
            temp.set(i + 1, 0, value);
        }
        return temp;
    }

    public static Jama.Matrix sortMatrixAscending(Jama.Matrix x) {
        int size = x.getRowDimension();
        Jama.Matrix temp = new Jama.Matrix(size, 1);
        ArrayList<Double> list = new ArrayList<Double>();
        for (int i = 0; i < size; i++) {
            list.add(x.get(i, 0));
        }
        Collections.sort(list);
        for (int i = 0; i < size; i++) {
            double value = (list.get(i)).doubleValue();
            temp.set(i, 0, value);
        }
        return temp;
    }

    public static Jama.Matrix reverseSortOnFirstElement(Jama.Matrix x) {
        Jama.Matrix temp = x.copy();
        int size = x.getRowDimension();
        ArrayList<Double> list = new ArrayList<Double>();
        for (int i = 1; i < size; i++) {
            list.add(x.get(i, 0));
        }
        Collections.sort(list);
        for (int i = 0; i < size - 1; i++) {
            double value = (list.get(i)).doubleValue();
            temp.set(size - i - 1, 0, value);
        }
        return temp;
    }

    public static void prettyPrint(Jama.Matrix x) {
        nf.setMinimumFractionDigits(6);
        nf.setMaximumFractionDigits(6);
        for (int i = 0; i < x.getRowDimension(); i++) {
            System.out.print("{");
            for (int j = 0; j < x.getColumnDimension(); j++) {
                // System.out.print(" " + nf.format(x.get(i, j)));
                System.out.format("% .6f ", x.get(i, j));
            }
            System.out.print(" }\n");
        }
    }

    public static void prettyPrint(Jama.Matrix x, int numberOfRowsFromBeginning) {
        nf.setMinimumFractionDigits(6);
        nf.setMaximumFractionDigits(6);
        for (int i = 0; i < numberOfRowsFromBeginning; i++) {
            System.out.print("{");
            for (int j = 0; j < x.getColumnDimension(); j++) {
                // System.out.print(" " + nf.format(x.get(i, j)));
                System.out.format("% .6f ", x.get(i, j));
            }
            System.out.print(" }\n");
        }

    }

    public static void prettyPrintVector(Jama.Matrix x) {
        if (x != null) {
            nf.setMinimumFractionDigits(6);
            nf.setMaximumFractionDigits(6);
            System.out.print("{");

            for (int i = 0; i < x.getRowDimension(); i++) {
                // System.out.print(" " + nf.format(x.get(i, 0)));
                System.out.format("% .6f ", x.get(i, 0));
            }
            System.out.println(" }");
        } else {
            System.out.println("{ null }");
        }
    }

    public static void prettyPrintMatrixDiag(Jama.Matrix x) {
        if (x != null) {
            nf.setMinimumFractionDigits(6);
            nf.setMaximumFractionDigits(6);
            System.out.print("{");

            for (int i = 0; i < x.getRowDimension(); i++) {
                // System.out.print(" " + nf.format(x.get(i, 0)));
                System.out.format("% .6f ", x.get(i, i));
            }
            System.out.println(" }");
        } else {
            System.out.println("{ null }");
        }
    }
}
