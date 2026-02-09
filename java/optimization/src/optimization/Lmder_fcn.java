package optimization;

public interface Lmder_fcn {

   void fcn(int m, int n, double x[], double fvec[],
            double fjac[][], int iflag[]);

}
