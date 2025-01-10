package NFlowShop;
import NFlowShop.UniFuncMatrixGen;
import java.util.*;

class makeMatrix{
     int jobs;
     int mcs;
     int[][] jbmatrix;
     int[] jbSeq;
    makeMatrix(){
       try(Scanner sc=new Scanner(System.in)){
           System.out.println("Enter the no. of jobs...");
           int jobs = sc.nextInt();
           this.jobs=jobs;
           System.out.println("Enter the no. of machines...");
           int mcs = sc.nextInt();
           this.mcs=mcs;
           jbmatrix = new int[jobs][mcs];

           for (int i = 0; i < jobs; i++) {
               System.out.println("Enter the machines time for job no.:" + (i + 1));
               for (int j = 0; j < mcs; j++) {
                   this.jbmatrix[i][j] = sc.nextInt();
               }
           }
           System.out.println("The job matrix is :");
           for (int i = 0; i < jobs; i++) {
               for (int j = 0; j < mcs; j++) {
                   System.out.print(jbmatrix[i][j] + " ");
               }
               System.out.println();
           }
           jbSeq = new int[jobs];
           System.out.println("Enter the sequence of jobs according to the no. of jobs without repetition(Starting from 1) . ");

           for (int i = 0; i < jobs; i++) {
               this.jbSeq[i] = sc.nextInt();
           }
           System.out.println("The Sequence of Job is :");
           for(int i=0;i<jobs;i++){
               System.out.print(jbSeq[i]+" ");
           }
           System.out.println();
       }catch (Exception e){
           System.out.println(e);
       }
    }
    public int[][] getjbMatrix(){
        return jbmatrix;
    }
    public int[] getjbSeq(){
        return jbSeq;
    }
}
class Calculation{
    public  static int[][] getMakeSpan(int[][] jbMatrix,int[] jbSeq){
        int jobs=jbMatrix.length;
        int mcs=jbMatrix[0].length;
        int[][] SpanMatrix=new int[jobs][mcs];
        for(int i=0;i<jobs;i++){
            int n=jbSeq[i];
            if(i==0){
                int a=0;
                for(int j=0;j<mcs;j++){
                    SpanMatrix[i][j]=a+jbMatrix[n-1][j];
                    a=SpanMatrix[i][j];
                }
            }
            else{
               for(int j=0;j<mcs;j++){
                   if(j==0){
                       SpanMatrix[i][j]=SpanMatrix[i-1][j]+jbMatrix[n-1][j];
                   }
                   else{
                       if(SpanMatrix[i][j-1]>SpanMatrix[i-1][j]){
                           SpanMatrix[i][j]=SpanMatrix[i][j-1]+jbMatrix[n-1][j];
                       }
                       else{
                           SpanMatrix[i][j]=SpanMatrix[i-1][j]+jbMatrix[n-1][j];
                       }
                   }
               }
               }


        }
        System.out.println("The Span Matrix is");
        for(int i=0;i<jobs;i++){
            for(int j=0;j<mcs;j++){
                System.out.print(SpanMatrix[i][j] + " ");
            }
            System.out.println();
        }
        return  SpanMatrix;
    }
    public static int totalFlow(int[][] SpanMatrix){
        int tft=0;
        int mcs=SpanMatrix[0].length;
        for (int[] spanMatrix : SpanMatrix) {
            tft += spanMatrix[mcs - 1];
        }
        return tft;
    }
}

public class MatrixIn {
    public static  void main(String[] args){
        try  {
            makeMatrix mx=new makeMatrix();
            int[][] jbmatrix= mx.getjbMatrix();
            int[] jbSeq=mx.getjbSeq();
            int jobs=jbmatrix.length;
            int mcs=jbmatrix[0].length;
            int[][] SpanMatrix = Calculation.getMakeSpan(jbmatrix, jbSeq);
            int makeSpan= SpanMatrix[jobs-1][mcs-1];
            System.out.println("The makespan is: " + makeSpan);

            int TFT = Calculation.totalFlow(SpanMatrix);
            System.out.println("Total Flow Time is: " + TFT);

            float avgFT= (float) TFT /jobs;
            System.out.println("Average Flow time is: " +avgFT);

        } catch (Exception e) {
            System.out.println(e);
        }
    }
}
