package NFlowShop;
import NFlowShop.UniFuncMatrixGen;
import NFlowShop.MakeDelayMatrix;
public class Main {

    public static void main(String[] args) {
        int[] sequence={4,3,0,1,2};
        UniFuncMatrixGen gn=new UniFuncMatrixGen();
        int[][] matrix = gn.returnMatrix();
        MakeDelayMatrix md=new MakeDelayMatrix(matrix);
        System.out.println("Matrix");
        for(int i=0;i<matrix.length;i++){
            for(int j=0;j<matrix[0].length;j++){
                System.out.print(matrix[i][j] + "\t");
            }
            System.out.println();
        }
        int Size=matrix.length;
        int mcs=matrix[0].length;
        int[][] delayedMatrix=md.returnDelayMatrix();
        System.out.println("Delay Matrix");
        for(int i=0;i<Size;i++){
            for(int j=0;j<Size;j++){
                System.out.print(delayedMatrix[i][j] + "\t");
            }
            System.out.println();
        }
        int[][] SpanMatrix=getMakeSpan(matrix,delayedMatrix,sequence);
        for(int i=0;i<SpanMatrix.length;i++){
            for(int j=0;j<SpanMatrix[0].length;j++){
                System.out.print(SpanMatrix[i][j] + "\t");
            }
            System.out.println();
        }
    }
    public static int[][] getMakeSpan(int[][] jbMatrix, int[][] delayMatrix, int[] jbSeq) {
        int jobs = jbSeq.length; // Use the length of the job sequence instead of jbMatrix
        int mcs = jbMatrix[0].length;
        int[][] SpanMatrix = new int[jobs][mcs]; // Create the SpanMatrix with the correct length

        for (int i = 0; i < jobs; i++) {
            int n = jbSeq[i];
            if (i == 0) {
                int a = 0;
                for (int j = 0; j < mcs; j++) {
                    SpanMatrix[i][j] = a + jbMatrix[n][j];
                    a = SpanMatrix[i][j];
                }
            } else {
                int prevJob = jbSeq[i - 1];
                int delay = delayMatrix[prevJob][n]; // Add delay from previous job to current job
                for (int j = 0; j < mcs; j++) {
                    if (j == 0) {
                        SpanMatrix[i][j] = SpanMatrix[i - 1][j] + jbMatrix[n][j] + delay;
                    } else {
                        SpanMatrix[i][j] = SpanMatrix[i][j-1] + jbMatrix[n][j];
                    }
                }
            }
        }
        return SpanMatrix;
    }
}
