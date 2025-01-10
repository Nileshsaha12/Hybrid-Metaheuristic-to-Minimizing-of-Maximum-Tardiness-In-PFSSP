package NFlowShop;
import java.util.*;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;

class MatrixCreation{
    static int jobs;
    static int mcs;
    static int MAXCOLS = 510;

    static long timeseed;
    long[] seed=new long[2];
     int Matrix[][] ;
    MatrixCreation(){
            try{

                long sed = FileRead();
                long[] seed = {sed};
                int tmatrix[][]= new int[mcs][jobs];
                // Initialize and populate the matrix
                for (int i = 0; i < mcs; i++) {
                    for (int j = 0; j < jobs; j++) {
                        tmatrix[i][j] = MatrixGen.unif(seed);
                    }
                }

                // Print the original matrix
                //System.out.println("Original Matrix:");
               // printMatrix(matrix);

                // Transpose the matrix
               Matrix= transposeMatrix(tmatrix);

                // Print the transposed matrix
               // System.out.println("Transposed Matrix:");
              //  printMatrix(Matrix);


            }catch (Exception e){
                System.out.println(e);
            }


    }
    public   int[][] getMatrix(){
        return Matrix;
    }
    public static long FileRead() throws IOException {
        BufferedReader br = new BufferedReader(new FileReader("C:\\\\Users\\\\Nilesh Saha\\\\OneDrive\\\\Desktop\\\\4thyr\\\\NFlowShop\\\\tlrdata1.txt"));
        String line=br.readLine();

            StringTokenizer stk = new StringTokenizer(line);
            if (stk.hasMoreTokens()) {
                jobs = Integer.parseInt(stk.nextToken());
                System.out.println("Jobs: " + jobs);
                if (stk.hasMoreTokens()) {
                    mcs = Integer.parseInt(stk.nextToken());
                    System.out.println("Machines: " + mcs);
                    if (stk.hasMoreTokens()) {
                        timeseed = Long.parseLong(stk.nextToken());
                        System.out.println("Timeseed: " + timeseed);
                    } else {
                        System.err.println("Not enough tokens in line: " + line);
                    }
                } else {
                    System.err.println("Not enough tokens in line: " + line);
                }
            } else {
                System.err.println("Not enough tokens in line: " + line);
            }


        br.close(); // Close the BufferedReader
        return timeseed;
    }
    public static int[][] transposeMatrix(int[][] matrix) {
        int rowCount = matrix.length;
       // System.out.println(rowCount);
        int colCount = matrix[0].length;
       // System.out.println(colCount);
        int[][] transposedMatrix = new int[colCount][rowCount];

        for (int i = 0; i < colCount; i++) {

            for (int j = 0; j < rowCount; j++) {
                transposedMatrix[i][j] = matrix[j][i];

            }
        }
        return transposedMatrix;
    }

    public static void printMatrix(int[][] matrix) {
        for (int[] row : matrix) {
            for (int element : row) {
                System.out.print(element + "\t");
            }
            System.out.println();
        }
}

}
class MatrixGen{
    // Random number generator based on the Linear Congruential Generator (LCG)
    static int unif(long[] seed) {
        int low=1;
        int high=99;
        long m = 2147483647L, a = 16807L, b = 127773L, c = 2836L;
        int k;
        double val01;
        k = (int) (seed[0] / b);
        seed[0] = a * (seed[0] % b) - k * c;
        if (seed[0] < 0) {
            seed[0] = seed[0] + m;
        }
        val01 = (double) seed[0] / (double) m;
        return (low + (int) (val01 * (high - low + 1)));
    }
}

public class UniFuncMatrixGen {
   int[][] matrix;
    UniFuncMatrixGen(){
       MatrixCreation mc= new MatrixCreation();
       matrix= mc.getMatrix();
    }
    public int[][] returnMatrix(){
        return  matrix;
    }
}
