package NFlowShop;
import java.util.*;

public class Scheduling {
    // Problem parameters
     int numJobs;
     int numMachines;
     int[][] processingTimes;
     int[] dueDates;
     int[][] delayMatrix;

    // Algorithm parameters
     int populationSize;
     int maxGenerationsGA;
     int maxIterationsJAYA;
     int maxIterationsSA;
     double crossoverRate;
     double mutationRate;
     double initialTemperature;
     double coolingRate;
     List<Integer> bestOverallSchedule;
     double bestOverallTardiness = Double.MAX_VALUE;
    public Scheduling(int jobs, int mcs, int[][] jbmatrix, int[][] delayMatrix, int[] dueDates,
                      int populationSize, int maxGenerationsGA, int maxIterationsJAYA, int maxIterationsSA,
                      double crossoverRate, double mutationRate, double initialTemperature, double coolingRate) {
        this.numJobs = jobs;
        this.numMachines = mcs;
        this.processingTimes = jbmatrix;
        this.delayMatrix = delayMatrix;
        this.dueDates = dueDates;
        this.populationSize = populationSize;
        this.maxGenerationsGA = maxGenerationsGA;
        this.maxIterationsJAYA = maxIterationsJAYA;
        this.maxIterationsSA = maxIterationsSA;
        this.crossoverRate = crossoverRate;
        this.mutationRate = mutationRate;
        this.initialTemperature = initialTemperature;
        this.coolingRate = coolingRate;
    }
    public List<Integer> getBestSchedule() {
        return new ArrayList<>(bestOverallSchedule);
    }

    public double getBestTardiness() {
        return bestOverallTardiness;
    }
    private void updateBest(List<Integer> schedule) {
        double tardiness = calculateMaxTardiness(schedule, processingTimes, dueDates, delayMatrix);
        if (tardiness < bestOverallTardiness) {
            bestOverallSchedule = new ArrayList<>(schedule);
            bestOverallTardiness = tardiness;
        }
    }

    public void optimize() {
        // Phase 1: Genetic Algorithm
        bestOverallTardiness = Double.MAX_VALUE;

        List<Integer> gaBest = runGA();
        updateBest(gaBest);

        List<Integer> jayaBest = runJAYA(gaBest);
        updateBest(jayaBest);

        List<Integer> saBest = runSA(jayaBest);
        updateBest(saBest);
    }



    private List<Integer> runGA() {
        List<List<Integer>> population = generateInitialPopulation(populationSize, numJobs);
        List<Integer> globalBest = new ArrayList<>(population.get(0));
        List<Integer> globalWorst = new ArrayList<>(population.get(0));

        double bestTardiness = calculateMaxTardiness(globalBest, processingTimes, dueDates, delayMatrix);
        double worstTardiness = bestTardiness;

        for (int generation = 0; generation < maxGenerationsGA; generation++) {
            List<List<Integer>> newPopulation = new ArrayList<>();
            for (int i = 0; i < populationSize; i += 2) {
                // Selection (Tournament Selection)
                List<Integer> parent1 = tournamentSelection(population);
                List<Integer> parent2 = tournamentSelection(population);

                // Crossover
                List<Integer> child1 = crossover(parent1, parent2);
                List<Integer> child2 = crossover(parent2, parent1);

                // Mutation
                if (Math.random() < mutationRate) mutate(child1);
                if (Math.random() < mutationRate) mutate(child2);

                newPopulation.add(child1);
                newPopulation.add(child2);
            }
            population = newPopulation;

            // Update global best
            for (List<Integer> individual : population) {
                double tardiness = calculateMaxTardiness(individual, processingTimes, dueDates, delayMatrix);
                if (tardiness < bestTardiness) {
                    globalBest = new ArrayList<>(individual);
                    bestTardiness = tardiness;
                }
                if (tardiness > worstTardiness) {
                    globalWorst = new ArrayList<>(individual);
                    worstTardiness = tardiness;
                }
            }
        }
        return globalBest;
    }

    private List<Integer> tournamentSelection(List<List<Integer>> population) {
        Random rand = new Random();
        List<Integer> best = null;
        double bestTardiness = Double.MAX_VALUE;
        for (int i = 0; i < 5; i++) { // Tournament size = 5
            int index = rand.nextInt(populationSize);
            List<Integer> candidate = population.get(index);
            double tardiness = calculateMaxTardiness(candidate, processingTimes, dueDates, delayMatrix);
            if (tardiness < bestTardiness) {
                best = candidate;
                bestTardiness = tardiness;
            }
        }
        return new ArrayList<>(best);
    }
    private List<Integer> runJAYA(List<Integer> initialSolution) {
        // Create a small population for JAYA (e.g., 50 individuals)
        List<List<Integer>> population = generateJAYAPopulation(initialSolution, 50);
        List<Integer> globalBest = new ArrayList<>(initialSolution);
        List<Integer> globalWorst = new ArrayList<>(initialSolution);
        double bestTardiness = calculateMaxTardiness(globalBest, processingTimes, dueDates, delayMatrix);
        double worstTardiness = bestTardiness;

        // Find initial best and worst in the population
        for (List<Integer> individual : population) {
            double tardiness = calculateMaxTardiness(individual, processingTimes, dueDates, delayMatrix);
            if (tardiness < bestTardiness) {
                globalBest = new ArrayList<>(individual);
                bestTardiness = tardiness;
            }
            if (tardiness > worstTardiness) {
                globalWorst = new ArrayList<>(individual);
                worstTardiness = tardiness;
            }
        }

        // JAYA iterations
        for (int iteration = 0; iteration < maxIterationsJAYA; iteration++) {
            for (List<Integer> individual : population) {
                List<Integer> newSolution = new ArrayList<>(individual);
                Random rand = new Random();

                // Update each job position using JAYA formula
                for (int i = 0; i < numJobs; i++) {
                    double r1 = rand.nextDouble();
                    double r2 = rand.nextDouble();
                    int newValue = (int) (individual.get(i)
                            + r1 * (globalBest.get(i) - individual.get(i))
                            - r2 * (globalWorst.get(i) - individual.get(i)));
                    newSolution.set(i, Math.abs(newValue) % numJobs); // Ensure valid job index
                }

                // Repair duplicates
                newSolution = repairSolution(newSolution);

                // Calculate tardiness
                double newTardiness = calculateMaxTardiness(newSolution, processingTimes, dueDates, delayMatrix);
                if (newTardiness < calculateMaxTardiness(individual, processingTimes, dueDates, delayMatrix)) {
                    individual.clear();
                    individual.addAll(newSolution);
                }
            }

            // Update global best and worst after each iteration
            for (List<Integer> individual : population) {
                double tardiness = calculateMaxTardiness(individual, processingTimes, dueDates, delayMatrix);
                if (tardiness < bestTardiness) {
                    globalBest = new ArrayList<>(individual);
                    bestTardiness = tardiness;
                }
                if (tardiness > worstTardiness) {
                    globalWorst = new ArrayList<>(individual);
                    worstTardiness = tardiness;
                }
            }
        }

        return globalBest;
    }
    private List<List<Integer>> generateJAYAPopulation(List<Integer> initialSolution, int populationSize) {
        List<List<Integer>> population = new ArrayList<>();
        population.add(new ArrayList<>(initialSolution)); // Include the initial solution

        for (int i = 1; i < populationSize; i++) {
            List<Integer> perturbed = new ArrayList<>(initialSolution);
            // Swap two random jobs to create a neighbor
            int idx1 = (int) (Math.random() * numJobs);
            int idx2 = (int) (Math.random() * numJobs);
            Collections.swap(perturbed, idx1, idx2);
            population.add(perturbed);
        }
        return population;
    }
    private List<Integer> repairSolution(List<Integer> solution) {
        // Fix duplicates (e.g., using set to track missing jobs)
        Set<Integer> missing = new HashSet<>();
        for (int i = 0; i < numJobs; i++) missing.add(i);
        for (int job : solution) missing.remove(job);

        List<Integer> repaired = new ArrayList<>(solution);
        for (int i = 0; i < numJobs; i++) {
            if (Collections.frequency(repaired, repaired.get(i)) > 1) {
                repaired.set(i, missing.iterator().next());
                missing.remove(repaired.get(i));
            }
        }
        return repaired;
    }

    private List<Integer> runSA(List<Integer> initialSolution) {
        List<Integer> currentSolution = new ArrayList<>(initialSolution);
        List<Integer> bestSolution = new ArrayList<>(currentSolution);
        double currentTardiness = calculateMaxTardiness(currentSolution, processingTimes, dueDates, delayMatrix);
        double bestTardiness = currentTardiness;
        double temperature = initialTemperature;
        Random rand = new Random();

        for (int iteration = 0; iteration < maxIterationsSA; iteration++) {
            // Generate neighbors without causing infinite loops
            List<Integer> neighborForward = new ArrayList<>(currentSolution);
            List<Integer> neighborBackward = new ArrayList<>(currentSolution);

            // Forward shift: move job to a later position
            if (numJobs > 1) {
                int idx1 = rand.nextInt(numJobs - 1); // Ensures idx1 can have a valid idx2 > idx1
                int idx2 = idx1 + 1 + rand.nextInt(numJobs - idx1 - 1);
                Integer job = neighborForward.remove(idx1);
                neighborForward.add(idx2, job);
            }

            // Backward shift: move job to an earlier position
            if (numJobs > 1) {
                int idx3 = 1 + rand.nextInt(numJobs - 1); // Ensures idx3 > 0
                int idx4 = rand.nextInt(idx3);
                Integer job2 = neighborBackward.remove(idx3);
                neighborBackward.add(idx4, job2);
            }

            // Evaluate neighbors
            double tardinessForward = calculateMaxTardiness(neighborForward, processingTimes, dueDates, delayMatrix);
            double tardinessBackward = calculateMaxTardiness(neighborBackward, processingTimes, dueDates, delayMatrix);

            // Determine the best neighbor
            List<Integer> bestNeighbor;
            double bestNeighborTardiness;
            if (tardinessForward < tardinessBackward) {
                bestNeighbor = neighborForward;
                bestNeighborTardiness = tardinessForward;
            } else {
                bestNeighbor = neighborBackward;
                bestNeighborTardiness = tardinessBackward;
            }

            // Metropolis acceptance criterion
            double delta = bestNeighborTardiness - currentTardiness;
            if (delta < 0 || Math.exp(-delta / temperature) > rand.nextDouble()) {
                currentSolution = bestNeighbor;
                currentTardiness = bestNeighborTardiness;

                if (currentTardiness < bestTardiness) {
                    bestSolution = new ArrayList<>(currentSolution);
                    bestTardiness = currentTardiness;
                }
            }

            // Cool down temperature
            temperature *= coolingRate;
        }

        return bestSolution;
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
        // Calculate completion times
        int[][] completionTimes = calculateCompletionTimes(schedule, processingTimes, delayMatrix);
        int numJobs = schedule.size();
        double maxTardiness = 0;

        for (int j = 0; j < numJobs; j++) {
            int job = schedule.get(j); // Get the job ID from the schedule
            // Access the completion time of the last machine for this job
            int lastMachineIndex = processingTimes[0].length - 1;
            double tardiness = Math.max(0, completionTimes[j][lastMachineIndex] - dueDates[job]);
            if (tardiness > maxTardiness) {
                maxTardiness = tardiness;
            }
        }

        return maxTardiness;
    }


    public static int[][] calculateCompletionTimes(List<Integer> schedule, int[][] processingTimes, int[][] delayMatrix) {
        // Validate the schedule
        for (int job : schedule) {
            if (job < 0 || job >= processingTimes.length) {
                throw new IllegalArgumentException("Invalid job index in schedule: " + job);
            }
        }

        int jobs = schedule.size();
        int mcs = processingTimes[0].length;
        int[][] SpanMatrix = new int[jobs][mcs];

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
//        for(int[] x:SpanMatrix){
//            for(int y:x){
//                System.out.print(y+" ");
//            }
//            System.out.println();
//        }
        return SpanMatrix;
    }
    public static List<Integer> crossover(List<Integer> parent1, List<Integer> parent2) {
        Random rd = new Random();
        // System.out.println("Parent 1: " + parent1);
       // System.out.println("Parent 2: " + parent2);

        int length = parent1.size();
        List<Integer> child = new ArrayList<>(Collections.nCopies(length, -1));

        // Select a random segment from parent1
        int lowerbound = rd.nextInt(length);
        int upperbound = rd.nextInt(length - lowerbound) + lowerbound;

        // Copy the segment from parent1 to child
        for (int i = lowerbound; i < upperbound; i++) {
            child.set(i, parent1.get(i));
        }

        // Collect elements from parent2 that are not in the copied segment
        List<Integer> notPresent = new ArrayList<>();
        for (int i = 0; i < length; i++) {
            int element = parent2.get(i);
            if (!child.contains(element)) {
                notPresent.add(element);
            }
        }

        // Fill the remaining positions in child with elements from notPresent
        int count = 0;
        for (int i = 0; i < length; i++) {
            if (child.get(i) == -1) { // Check if the position is empty
                if (count < notPresent.size()) { // Ensure we don't exceed notPresent size
                    child.set(i, notPresent.get(count));
                    count++;
                } else {
                    // If notPresent is exhausted, fill with remaining elements from parent2
                    for (int j = 0; j < length; j++) {
                        if (!child.contains(parent2.get(j))) {
                            child.set(i, parent2.get(j));
                            break;
                        }
                    }
                }
            }
        }

       // System.out.println("Crossover Child: " + child);
        return child;
    }

    public static void mutate(List<Integer> individual) {
        int size = individual.size();

        // Randomly select two indices
        int index1 = (int) (Math.random() * size);
        int index2 = (int) (Math.random() * size);

        // Ensure index1 is less than index2
        if (index1 > index2) {
            int temp = index1;
            index1 = index2;
            index2 = temp;
        }

        // Reverse the subset of the list between index1 and index2
        Collections.reverse(individual.subList(index1, index2 + 1));

        // System.out.println("Mutation: " + individual);
    }
}