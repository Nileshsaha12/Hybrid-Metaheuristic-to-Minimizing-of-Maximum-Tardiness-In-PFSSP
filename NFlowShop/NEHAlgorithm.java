import java.util.Arrays;
import java.util.Comparator;

public class NEHAlgorithm {

    // Function to calculate makespan for a given sequence of jobs
    public static int calculateMakespan(int[][] processingTimes, int[] jobSequence) {
        int numMachines = processingTimes[0].length;
        int numJobs = jobSequence.length;

        int[][] completionTimes = new int[numJobs + 1][numMachines + 1];

        for (int i = 1; i <= numJobs; i++) {
            int job = jobSequence[i - 1];
            for (int j = 1; j <= numMachines; j++) {
                completionTimes[i][j] = Math.max(completionTimes[i - 1][j], completionTimes[i][j - 1]) + processingTimes[job][j - 1];
            }
        }

        return completionTimes[numJobs][numMachines];
    }

    // Function to apply NEH algorithm and return the best job sequence
    public static int[] nehAlgorithm(int[][] processingTimes) {
        int numJobs = processingTimes.length;
        Integer[] jobIndices = new Integer[numJobs];
        int[] jobSequence = new int[numJobs];

        for (int i = 0; i < numJobs; i++) {
            jobIndices[i] = i;
        }

        Arrays.sort(jobIndices, Comparator.comparingInt(i -> -Arrays.stream(processingTimes[i]).sum()));

        jobSequence[0] = jobIndices[0];
        for (int i = 1; i < numJobs; i++) {
            int bestPosition = -1;
            int bestMakespan = Integer.MAX_VALUE;

            for (int j = 0; j <= i; j++) {
                int[] tempSequence = new int[i + 1];
                for (int k = 0, l = 0; k <= i; k++) {
                    if (k == j) {
                        tempSequence[k] = jobIndices[i];
                    } else {
                        tempSequence[k] = jobSequence[l++];
                    }
                }

                int makespan = calculateMakespan(processingTimes, tempSequence);
                if (makespan < bestMakespan) {
                    bestMakespan = makespan;
                    bestPosition = j;
                }
            }

            for (int j = i; j > bestPosition; j--) {
                jobSequence[j] = jobSequence[j - 1];
            }
            jobSequence[bestPosition] = jobIndices[i];
        }

        return jobSequence;
    }

    public static void main(String[] args) {
        // Example processing times matrix
        int[][] processingTimes = {
                {2,4,2,3},
                {1,3,2,4},
                {5,1,1,6},
                {3,2,3,5},
                {2,1,3,3}
        };

        int[] bestSequence = nehAlgorithm(processingTimes);

        System.out.println("Best job sequence: " + Arrays.toString(bestSequence));
        System.out.println("Minimum makespan: " + calculateMakespan(processingTimes, bestSequence));
    }
}
