package NFlowShop;

import java.util.*;
public class Main {
    public static void main(String[] args) {
        long startTime = System.currentTimeMillis();

        // Problem setup (same as before)
        UniFuncMatrixGen gn = new UniFuncMatrixGen();
        int[][] matrix = gn.returnMatrix();
        MakeDelayMatrix md = new MakeDelayMatrix(matrix);
        int[][] delayedMatrix = md.returnDelayMatrix();
        int[] dueDates = calculateDueDates(matrix);
        int jobs = matrix.length;
        int mcs = matrix[0].length;

        // Algorithm parameters (same as before)
        int InitialpopulationSize = 500;
        int GApopulation=200;
        int maxGenerationsGA = 250;
        int maxIterationsJAYA = 1500;
        int maxIterationsSA = 50000;
        double crossoverRate = 0.8;
        double mutationRate = 0.05;
        double initialTemperature = 1000.0;
        double coolingRate = 0.90;

        // Results tracking
        int numIterations = 10;
        double[] iterationResults = new double[numIterations];
        double bestTardiness = Double.MAX_VALUE;
        List<Integer> bestOverallSchedule = null;

        for (int i = 0; i < numIterations; i++) {
            System.out.println("\n--- Iteration " + (i+1) + " ---");
            Scheduling scheduler = new Scheduling(jobs, mcs, matrix, delayedMatrix, dueDates,
                    InitialpopulationSize,GApopulation, maxGenerationsGA, maxIterationsJAYA, maxIterationsSA,
                    crossoverRate, mutationRate, initialTemperature, coolingRate);
            // Run optimization
            scheduler.optimize();





            // Get iteration results
            double currentTardiness = scheduler.getBestTardiness();
            iterationResults[i] = currentTardiness;

            // Update best overall
            if (currentTardiness < bestTardiness) {
                bestTardiness = currentTardiness;
                bestOverallSchedule = new ArrayList<>(scheduler.getBestSchedule());
            }

            System.out.printf("Iteration %d Tardiness: %.2f%n", (i+1), currentTardiness);
        }

        // Calculate average
        double averageTardiness = Arrays.stream(iterationResults).average().orElse(0);

        // Final results
        System.out.println("\n=== Final Results ===");
        System.out.printf("Average Maximum Tardiness: %.2f%n", averageTardiness);
        System.out.printf("Best Maximum Tardiness: %.2f%n", bestTardiness);
        System.out.println("Best Schedule: " + bestOverallSchedule);
        System.out.println("Execution Time: " + (System.currentTimeMillis() - startTime) + " ms");
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
