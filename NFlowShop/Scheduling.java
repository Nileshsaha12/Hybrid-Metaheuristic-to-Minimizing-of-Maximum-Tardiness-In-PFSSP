package NFlowShop;
import java.util.*;
import java.util.stream.Collectors;
public class Scheduling {
    // Define problem parameters
    int numJobs;  // Example number of jobs
    int machine;
    int populationSize = 50;
    int maxGenerations = 100;
    double crossoverRate = 0.8;
    double mutationRate = 0.1;
    int[][] processingTimes; // Assume 3 machines for simplicity
    int[] dueDates;
    int[][] delayMatrix;

    Scheduling(int jobs, int mcs, int[][] jbmatrix, int[][] delayMtrix, int[] dueTime) {
        numJobs = jobs;
        machine = mcs;
        processingTimes = jbmatrix;
        dueDates = dueTime;
        delayMatrix = delayMtrix;
        List<List<Integer>> population = generateInitialPopulation(populationSize, numJobs);
        //print
        for (List<Integer> e :
                population) {
            System.out.println(e);
        }
        List<Integer> globalBest = population.get(0);
        for (int i = 0; i < globalBest.size(); i++) {
            System.out.print(globalBest.get(i));
        }
        double globalBestTardiness = calculateMaxTardiness(globalBest, processingTimes, dueDates, delayMatrix);
        for (int generation = 0; generation < maxGenerations; generation++) {
            List<List<Integer>> newPopulation = new ArrayList<>();

            // Select parents and apply crossover
            for (int i = 0; i < populationSize; i += 2) {
                List<Integer> parent1 = population.get(i);
                List<Integer> parent2 = population.get(i + 1);
                List<Integer> child1 = crossover(parent1, parent2);
                List<Integer> child2 = crossover(parent2, parent1);

                // Apply mutation
                if (Math.random() < mutationRate) {
                    mutate(child1);
                }
                if (Math.random() < mutationRate) {
                    mutate(child2);
                }

                newPopulation.add(child1);
                newPopulation.add(child2);
            }
            for (List<Integer> t :
                    newPopulation) {
                System.out.println(t);
            }
            // Apply Jaya refinement
//            for (int i = 0; i < newPopulation.size(); i++) {
//                newPopulation.set(i, jayaAlgorithm(newPopulation.get(i), globalBest, processingTimes, dueDates, delayMatrix));
//            }
//            System.out.println("after jaya");
//            for (List<Integer> t :
//                    newPopulation) {
//                System.out.println(t);
//            }

            // Update the population with the best solution found so far
            population = newPopulation;
            for (List<Integer> individual : population) {
                double tardiness = calculateMaxTardiness(individual, processingTimes, dueDates, delayMatrix);
                System.out.println(tardiness);
                if (tardiness < globalBestTardiness) {
                    globalBest = individual;
                    globalBestTardiness = tardiness;
                }
            }
        }
        System.out.println("Best Schedule: " + globalBest);
        System.out.println("Best Maximum Tardiness: " + globalBestTardiness);
    }

    public static List<List<Integer>> generateInitialPopulation(int populationSize, int numJobs) {
        List<List<Integer>> population = new ArrayList<>();

        for (int i = 0; i < populationSize; i++) {
            List<Integer> individual = new ArrayList<>();
            for (int j = 0; j < numJobs; j++) {
                individual.add(j);  // Add job index to the individual list
            }
            Collections.shuffle(individual);  // Randomly shuffle the jobs
            population.add(individual);
        }

        return population;
    }

    public static double calculateMaxTardiness(List<Integer> schedule, int[][] processingTimes, int[] dueDates, int[][] delayMatrix) {
        int[][] completionTimes = calculateCompletionTimes(schedule, processingTimes, delayMatrix);
        int numJobs = schedule.size();
        double maxTardiness = 0;

        for (int j = 0; j < numJobs; j++) {
            int job = schedule.get(j);
            double tardiness = Math.max(0, completionTimes[completionTimes.length - 1][j] - dueDates[job]);
            if (tardiness > maxTardiness) {
                maxTardiness = tardiness;
            }
        }

        return maxTardiness;
    }

    public static int[][] calculateCompletionTimes(List<Integer> schedule, int[][] processingTimes, int[][] delayMatrix) {
        int jobs = schedule.size(); // Use the length of the job sequence instead of jbMatrix
        int mcs = processingTimes[0].length;
        int[][] SpanMatrix = new int[jobs][mcs]; // Create the SpanMatrix with the correct length

        for (int i = 0; i < jobs; i++) {
            int n = schedule.get(i);
            if (i == 0) {
                int a = 0;
                for (int j = 0; j < mcs; j++) {
                    SpanMatrix[i][j] = a + processingTimes[n][j];
                    a = SpanMatrix[i][j];
                }
            } else {
                int prevJob = schedule.get(i - 1);
                int delay = delayMatrix[prevJob][n]; // Add delay from previous job to current job
                for (int j = 0; j < mcs; j++) {
                    if (j == 0) {
                        SpanMatrix[i][j] = SpanMatrix[i - 1][j] + processingTimes[n][j] + delay;
                    } else {
                        SpanMatrix[i][j] = SpanMatrix[i][j - 1] + processingTimes[n][j];
                    }
                }
            }
        }
        return SpanMatrix;

    }

    public static List<Integer> crossover(List<Integer> parent1, List<Integer> parent2) {
        Random rd = new Random();
        System.out.println(parent1);
        System.out.println(parent2);
        List<Integer> child = new ArrayList<Integer>(Collections.nCopies(parent1.size(), -1));

        ArrayList<Integer> notPresent = new ArrayList<>();


        int length = parent1.size();
        int lowerbound = rd.nextInt(length);
        System.out.println(lowerbound);
        int upperbound = rd.nextInt(length - lowerbound) + lowerbound;
        System.out.println(upperbound);

        for (int i = lowerbound; i < upperbound; i++) {
            child.set(i, parent1.get(i));
        }

        // Collect elements from parent2 that are not in the copied segment
        for (int i = 0; i < parent2.size(); i++) {
            if (!child.subList(lowerbound, upperbound).contains(parent2.get(i))) {
                notPresent.add(parent2.get(i));
            }
        }
        System.out.println(notPresent.toString());

        int count = 0;
        for (int i = 0; i < child.size(); i++) {
            if (i < lowerbound || i >= upperbound) {
                child.set(i, notPresent.get(count));
                count++;
                System.out.println("c"+count);

            }
        }
        System.out.println("Crossover");
        System.out.println((child));
        return  child;
    }

    public static void mutate(List<Integer> individual) {
                int size = individual.size();
                int index1 = (int) (Math.random() * size);
                int index2 = (int) (Math.random() * size);

                // Swap two jobs
                int temp = individual.get(index1);
                individual.set(index1, individual.get(index2));
                individual.set(index2, temp);
        System.out.println("Mutation");
        System.out.println(individual);
            }
//    public static List<Integer> jayaAlgorithm(List<Integer> individual, List<Integer> globalBest, int[][] processingTimes,int[] dueDates,int[][] delayMatrix) {
//        int size = individual.size();
//
//        System.out.println();
//        for (int i = 0; i < size; i++) {
//            // Move individual towards the global best
//            if (Math.random() < 0.5) {
//                individual.set(i, globalBest.get(i));
//            }
//        }
//        //print
//        System.out.println(individual);
//        // Optionally, recompute fitness for the individual after Jaya update
//        double individualTardiness = calculateMaxTardiness(individual, processingTimes,dueDates,delayMatrix);
//        System.out.println(individualTardiness);
//
//        return individual;
//    }
        }






