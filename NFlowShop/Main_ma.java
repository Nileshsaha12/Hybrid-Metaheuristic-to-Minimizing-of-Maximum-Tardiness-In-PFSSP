package NFlowShop;
import NFlowShop.UniFuncMatrixGen;
import NFlowShop.MakeDelayMatrix;
import NFlowShop.Scheduling;

import java.util.List;

public class Main_ma {
    public static void main(String[] args) {
        double bestTardiness = Double.MAX_VALUE;
        List<Integer> bestOverallSchedule = null;
        int[][] matrix = {
                {54,79,16,66,58},
                {83,3,89,58,56},
                {15	,11	,49	,31	,20},
                {71	,99	,15	,68	,85}
        };
        MakeDelayMatrix md = new MakeDelayMatrix(matrix);
        System.out.println("Matrix");
        for (int i = 0; i < matrix.length; i++) {
            for (int j = 0; j < matrix[0].length; j++) {
                System.out.print(matrix[i][j] + "\t");
            }
            System.out.println();
        }
        int jobs = matrix.length;
        int mcs = matrix[0].length;
        int[][] delayedMatrix = md.returnDelayMatrix();
        System.out.println("Delay Matrix");
        for (int i = 0; i < jobs; i++) {
            for (int j = 0; j < jobs; j++) {
                System.out.print(delayedMatrix[i][j] + "\t");
            }
            System.out.println();
        }
        int[] dueDates = calculateDueDates(matrix);
        for(int i:dueDates){
            System.out.print(i+ " ");
        }
        System.out.println();

        Scheduling sc=new Scheduling(jobs,mcs,matrix,delayedMatrix,dueDates,100,200,50,2000,0.9,0.15,2000,0.98 );
        sc.optimize();
        System.out.println(sc.getBestSchedule());
        System.out.println(sc.getBestTardiness());

    }
    /**
     * Calculates the due dates for all jobs based on their processing times.
     */
    private static int[] calculateDueDates(int[][] processingTimes) {
        int numJobs = processingTimes.length;
        int numMachines = processingTimes[0].length;
        int[] dueDates = new int[numJobs];

        for (int job = 0; job < numJobs; job++) {
            for (int machine = 0; machine < numMachines; machine++) {
                dueDates[job] += processingTimes[job][machine];
            }
        }

        return dueDates;
    }

}
