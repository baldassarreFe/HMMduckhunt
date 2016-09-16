import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.Map;
import java.util.Random;
import java.util.stream.Collectors;

class Player {

	/** states */
	private static final int N = 2;
	/** emissions */
	private static final int M = Constants.COUNT_MOVE;
	/** Max sequence length */
	private static final int MAX_T = 100;

	/** do not try to learn if the seq is shorter */
	private static final int WAIT_UNTIL_UPDATING = 60;
	/** update guesses every X steps */
	private static final int UPDATE_GUESSES_EVERY = 5;
	/** do not attempt to shoot if the seq is shorter */
	private static final int WAIT_UNTIL_SHOOTING = 60;

	/** Constant action for not shooting */
	public static final Action ACTION_DONT_SHOOT = new Action(-1, -1);
	private static final double MIN_ABSOLUTE_LOG_PROB_TO_GUESS = -260;
	private static final double MIN_RELATIVE_PROB_TO_GUESS = 0.7;
	private static final double MIN_ABSOLUTE_LOG_PROB_TO_GUESS_BLACK_STORK = -800;
	private static final int MIN_ROUND_TO_SHOOT = 1;

	// TODO better use a fixed value or one that varies linearly with T?
	private static final double MIN_PROB_TO_SHOOT = 0.97;
	private static double minProbToShoot(int T) {
		// double startProb = .4;
		// double endProb = 0;
		// double result = ((double) T - WAIT_UNTIL_SHOOTING) * (endProb -
		// startProb) / (MAX_T - WAIT_UNTIL_SHOOTING)
		// + startProb;
		// return result;
		return MIN_PROB_TO_SHOOT;
	}

	// debug info
	private static final boolean DEBUG_LINES = true;
	private static final boolean DEBUG_BIRD_UPDATES = false;
	private static final boolean DEBUG_SHOOT = false;
	private static final boolean DEBUG_GUESSING_UPDATES = false;
	private static final boolean DEBUG_ROUND_GUESSING = true;
	private static final boolean DEBUG_ROUND_REVEAL = true;
	private static final boolean DEBUG_TOTAL_SCORES = true;
	
//	private static final boolean DEBUG_LINES = false;
//	private static final boolean DEBUG_BIRD_UPDATES = false;
//	private static final boolean DEBUG_SHOOT = false;
//	private static final boolean DEBUG_GUESSING_UPDATES = false;
//	private static final boolean DEBUG_ROUND_GUESSING = false;
//	private static final boolean DEBUG_ROUND_REVEAL = false;
//	private static final boolean DEBUG_TOTAL_SCORES = true;

	// Matrices for the HMM of the birds in the current round
	// A = transition
	// B = emissions
	// Pi = initial state
	// O = observations so far
	// currentObsLength = length of the sequence until death
	// [bird number] [i] [j]
	double[][][] Abirds;
	double[][][] Bbirds;
	double[][] Pibirds;
	int[][] Obirds;
	int[] Tbirds;
	boolean[] riskyShots;

	// Matrices for the best HMM that are able to identify a species across the
	// rounds, keep more than one per each
	// [species id] [history index] [i] [j]
	double[][][] Aspecies = new double[Constants.COUNT_SPECIES][][];
	double[][][] Bspecies = new double[Constants.COUNT_SPECIES][][];
	double[][] Pispecies = new double[Constants.COUNT_SPECIES][];

	// info variables
	int currentRound = -1;
	int currentObsLength = 0;
	int numBirds;
	int lastGuessUpdate;
	int[] guessesMade;

	// overall statistics
	int totalBirds = 0;
	int totalShots = 0;
	int totalHits = 0;
	int totalGuessAttempts = 0;
	int totalGuesses = 0;
	Map<Integer, ArrayList<int[]>> sequencesPerSpecies;

	public Player() {
		sequencesPerSpecies = new HashMap<Integer, ArrayList<int[]>>();
		for (int s = 0; s < Constants.COUNT_SPECIES; s++) {
			sequencesPerSpecies.put(s, new ArrayList<int[]>());
		}
		for (int s = 0; s < Constants.COUNT_SPECIES; s++) {
			Aspecies[s] = null;
			Bspecies[s] = null;
			Pispecies[s] = null;
		}
	}

	/**
	 * Shoot!
	 *
	 * This is the function where you start your work.
	 *
	 * You will receive a variable pState, which contains information about all
	 * birds, both dead and alive. Each bird contains all past moves.
	 *
	 * The state also contains the scores for all players and the number of time
	 * steps elapsed since the last time this function was called.
	 *
	 * @param pState
	 *            the GameState object with observations etc
	 * @param pDue
	 *            time before which we must have returned
	 * @return the prediction of a bird we want to shoot at, or
	 *         ACTION_DONT_SHOOT to pass
	 */
	public Action shoot(GameState pState, Deadline pDue) {
		// is it a new round?
		if (pState.getRound() != currentRound) {
			initRound(pState);
		}

		// update currentObsLength
		int newTurns = pState.getNumNewTurns();
		currentObsLength += newTurns;

		// get new observations for alive birds
		for (int b = 0; b < numBirds; b++) {
			Bird bird = pState.getBird(b);
			for (int t = Tbirds[b]; t < bird.getSeqLength(); t++) {
				if (bird.wasAlive(t)) {
					int temp = 0;
					temp = bird.getObservation(t);
					Obirds[b][t] = temp;
					Tbirds[b]++;
				}
			}
		}

		// update matrices of this round
		if (currentObsLength >= WAIT_UNTIL_UPDATING) {
			updateBirdModels(pState);

			// choose whether to shoot or not
			if (currentRound >= MIN_ROUND_TO_SHOOT && currentObsLength >= WAIT_UNTIL_SHOOTING) {
				if (currentObsLength - lastGuessUpdate > UPDATE_GUESSES_EVERY) {
					updateGuesses();
				}
				return hitMeWithYourBestShot(pState);
			}
		}

		return ACTION_DONT_SHOOT;
	}

	// TODO is this a good method to do this thing??
	private double[] getLastStateVector(double[][] emissionMatrix, int lastEmission) {
		double[] returnVector = new double[emissionMatrix.length];
		double totalValue = 0;
		for (int i = 0; i < emissionMatrix.length; i++) {
			double value = emissionMatrix[i][lastEmission];
			returnVector[i] = value;
			totalValue += value;
		}
		for (int i = 0; i < returnVector.length; i++) {
			returnVector[i] /= totalValue;
		}
		return returnVector;
	}

	// Consider every bird, except for the ones that might be a STORK or whose
	// species is unknown, calculate the next emission vector, take the most
	// probable action and compare it with the most probable action of
	// other birds, eventually shoot the best action
	private Action hitMeWithYourBestShot(GameState pState) {
		double bestProb = 0.0;
		int bestBird = -1;
		int bestEmission = -1;
		for (int b = 0; b < numBirds; b++) {
			if (riskyShots[b])
				continue;
			if (pState.getBird(b).isAlive()) {
				// TODO pick the best between the 2 emissions
				double[] emissions = HiddenMarkovModel.nextEmissionDistributionKnowingSequence(
						Arrays.copyOf(Obirds[b], Tbirds[b]), Abirds[b], Bbirds[b], Pibirds[b]);

				// double[] lastState = getLastStateVector(Bbirds[b],
				// Obirds[b][Tbirds[b] - 1]);
				// double[] emissions2 =
				// HiddenMarkovModel.nextEmissions(Abirds[b], Bbirds[b],
				// lastState);

				int move = MatrixHelper.argMax(emissions);
				if (emissions[move] > bestProb) {
					bestProb = emissions[move];
					bestBird = b;
					bestEmission = move;
				}
			}
		}

		Action action;
		String debug = String.format("Round %d\tT %d\tbird %d\tmove %d\tconfidence %f", currentRound, currentObsLength,
				bestBird, bestEmission, bestProb);

		// if (bestProb > (currentObsLength/MAX_T) * MIN_PROB_TO_SHOOT) {
		if (bestBird != -1 && bestProb >= minProbToShoot(currentObsLength)) {
			totalShots++;
			action = new Action(bestBird, bestEmission);
			debug += " <-- SHOOT";
		} else {
			action = ACTION_DONT_SHOOT;
		}
		if (DEBUG_SHOOT)
			System.err.println(debug);
		return action;
	}

	// Use the new emissions to update the matrices (Baum Welch)
	// 1. Start with a new randomized model, use baum welch to estimate it
	// 2. Test the old model and the new against the sequence
	// 3. Keep the best model
	private void updateBirdModels(GameState pState) {
		// TODO parallelize over the birds?
		for (int b = 0; b < numBirds; b++) {
			Bird bird = pState.getBird(b);
			if (bird.isAlive()) {
				int T = Tbirds[b];
				int[] O = Arrays.copyOf(Obirds[b], T);

				double[][] Anew = MatrixHelper.newRowStochasticMatrix(N, N);
				double[][] Bnew = MatrixHelper.newRowStochasticMatrix(N, M);
				double[] Pinew = MatrixHelper.newRowStochasticArray(N);

				// TODO try different ratios
				// double[] statsNew =
				// HiddenMarkovModel.baumWelchAvoidOverfitting(O, Anew, Bnew,
				// Pinew, 0.7);
				double[] statsNew = HiddenMarkovModel.baumWelch(O, Anew, Bnew, Pinew);
				int itersNew = (int) statsNew[0];
				double logProbNew = statsNew[1];

				double[][] Aold = Abirds[b];
				double[][] Bold = Bbirds[b];
				double[] Piold = Pibirds[b];

				double[] statsOld = HiddenMarkovModel.baumWelch(O, Aold, Bold, Piold);
				// double[] statsOld =
				// HiddenMarkovModel.baumWelchAvoidOverfitting(O, Aold, Bold,
				// Piold,0.8);
				int itersOld = (int) statsOld[0];
				double logProbOld = statsOld[1];

				String debug = String.format("Round %d\tbird %d\tT %d\tlogProbNew %f (%d)\tlogProbOld %f (%d)\t",
						currentRound, b, T, logProbNew, itersNew, logProbOld, itersOld);

				// which one is the best model?
				if (Double.isFinite(logProbOld) && Double.isFinite(logProbNew)) {
					if (logProbNew > logProbOld) {
						Abirds[b] = Anew;
						Bbirds[b] = Bnew;
						Pibirds[b] = Pinew;
						debug += "kept new";
					} else {
						debug += "kept old";
					}
				} else if (Double.isFinite(logProbOld) && !Double.isFinite(logProbNew)) {
					debug += "kept old, new is NaN";
				} else if (!Double.isFinite(logProbOld) && Double.isFinite(logProbNew)) {
					Abirds[b] = Anew;
					Bbirds[b] = Bnew;
					Pibirds[b] = Pinew;
					debug += "kept new, old was NaN";
				} else {
					MatrixHelper.rowStochasticMatrix(Anew);
					MatrixHelper.rowStochasticMatrix(Bnew);
					MatrixHelper.rowStochasticArray(Pinew);
					debug += "randomized new";
				}
				if (DEBUG_BIRD_UPDATES)
					System.err.println(debug);
			}
		}
	}

	// At the beginning of every round initialize the matrices A,B,Pi
	private void initRound(GameState pState) {
		currentObsLength = 0;
		numBirds = pState.getNumBirds();
		Abirds = new double[numBirds][][];
		Bbirds = new double[numBirds][][];
		Pibirds = new double[numBirds][];
		Obirds = new int[numBirds][MAX_T];
		Tbirds = new int[numBirds];
		guessesMade = new int[numBirds];
		riskyShots = new boolean[numBirds];
		lastGuessUpdate = 0;
		currentRound++;
		totalBirds += numBirds;

		for (int b = 0; b < numBirds; b++) {
			Abirds[b] = MatrixHelper.newAlmostDiagonalMatrix(N);
			Bbirds[b] = MatrixHelper.newRowStochasticMatrix(N, M);
			Pibirds[b] = MatrixHelper.newRowStochasticArray(N);
		}

		if (DEBUG_LINES)
			System.err.println("\n------------------------------------------------------------------- ROUND "
					+ currentRound + " ---------");
	}

	/**
	 * Guess the species! This function will be called at the end of each round,
	 * to give you a chance to identify the species of the birds for extra
	 * points.
	 *
	 * Fill the vector with guessesMade for the all birds. Use SPECIES_UNKNOWN
	 * to avoid guessing.
	 *
	 * @param pState
	 *            the GameState object with observations etc
	 * @param pDue
	 *            time before which we must have returned
	 * @return a vector with guessesMade for all the birds
	 */
	public int[] guess(GameState pState, Deadline pDue) {
		if (DEBUG_LINES)
			System.err.println("\n---------------------------------------------------------------- GUESSING "
					+ currentRound + " ---------");

		// first round, no info available, random choice
		if (currentRound == 0) {
			// TODO put the most likely one here
			Random r = new Random();
			for (int b = 0; b < numBirds; b++) {
				guessesMade[b] = r.nextInt(Constants.SPECIES_BLACK_STORK);
			}
			return guessesMade;
		}

		updateGuesses();

		// use guesses unknown to get info on unseen species,
		// but do not fill the guesses vector with random choices,
		// fill no more than two slots per every unknown species,
		// otherwise we'll keep filling every hole with BLACK STORK
		// until we see one
		if (currentRound != 9) {
			for (int s = 0; s < 2 * Constants.COUNT_SPECIES; s++) {
				if (Aspecies[s / 2] == null) {
					for (int b = 0; b < numBirds; b++) {
						if (guessesMade[b] == Constants.SPECIES_UNKNOWN) {
							guessesMade[b] = s / 2;
							break;
						}
					}
				}
			}
		}

		return guessesMade;
	}

	private void updateGuesses() {
		lastGuessUpdate = currentObsLength;
		Arrays.fill(riskyShots, false);
		guessesMade = new int[numBirds];
		Arrays.fill(guessesMade, Constants.SPECIES_UNKNOWN);

		// for every bird
		for (int b = 0; b < numBirds; b++) {
			int[] O = Arrays.copyOf(Obirds[b], Tbirds[b]);
			double[] absoluteProbs = new double[Constants.COUNT_SPECIES];

			// for every species, keep the highest prob from the species models
			for (int s = 0; s < Constants.COUNT_SPECIES; s++) {
				// for every saved model in history, keep the best
				absoluteProbs[s] = Double.NEGATIVE_INFINITY;
				if (Aspecies[s] == null) {
					continue;
				}
				double[][] A = Aspecies[s];
				double[][] B = Bspecies[s];
				// try with every Pi that has one 1 and the rest 0,
				// keep only the best result
				for (int p = 0; p < N; p++) {
					double[] Pi = new double[N];
					Pi[p] = 1;
					double logProbTemp = HiddenMarkovModel.alphaPass(O, A, B, Pi);
					if (logProbTemp > absoluteProbs[s]) {
						absoluteProbs[s] = logProbTemp;
					}
				}
				if (DEBUG_GUESSING_UPDATES || (DEBUG_ROUND_GUESSING
						&& currentObsLength >= MAX_T - UPDATE_GUESSES_EVERY))
					System.err.printf("Bird %d\tT %d\tspecies %d\tconfidence %f\n", b, Tbirds[b], s, absoluteProbs[s]);
			}

			double[] relativeProbs = MatrixHelper.normalize(Arrays.stream(absoluteProbs).map(Math::exp).toArray());
			int bestGuess = MatrixHelper.argMax(relativeProbs);
			if (absoluteProbs[bestGuess] > MIN_ABSOLUTE_LOG_PROB_TO_GUESS
					&& relativeProbs[bestGuess] > MIN_RELATIVE_PROB_TO_GUESS) {
				guessesMade[b] = bestGuess;
			}

			if (guessesMade[b] == Constants.SPECIES_UNKNOWN
					|| absoluteProbs[Constants.SPECIES_BLACK_STORK] > MIN_ABSOLUTE_LOG_PROB_TO_GUESS_BLACK_STORK) {
				riskyShots[b] = true;
			}

			if (DEBUG_GUESSING_UPDATES || (DEBUG_ROUND_GUESSING
					&& currentObsLength >= MAX_T - UPDATE_GUESSES_EVERY)) {
				System.err.println("Relative prob: " + Arrays.stream(relativeProbs)
						.mapToObj(p -> String.format("%.10f", p)).collect(Collectors.joining(", ", "[", "]")));
				System.err.printf(" --->\tspecies %d\n", guessesMade[b]);
			}
		}
		
		if (DEBUG_GUESSING_UPDATES) {
			System.err.printf("Risky shots:");
			for (int b = 0; b < riskyShots.length; b++) {
				if (riskyShots[b]) {
					System.err.printf(" " + b);
				}
			}
			System.err.println();
		}
	}

	/**
	 * If you hit the bird you were trying to shoot, you will be notified
	 * through this function.
	 *
	 * @param pState
	 *            the GameState object with observations etc
	 * @param pBird
	 *            the bird you hit
	 * @param pDue
	 *            time before which we must have returned
	 */
	public void hit(GameState pState, int pBird, Deadline pDue) {
		totalHits++;
		if (DEBUG_SHOOT)
			System.err.printf("Round %d\tT %d\tbird %d\tHIT BIRD!!! (hits/shots = %d/%d)\n", currentRound,
					Tbirds[pBird], pBird, totalHits, totalShots);
	}

	/**
	 * If you made any guessesMade, you will find out the true species of those
	 * birds through this function.
	 *
	 * @param pState
	 *            the GameState object with observations etc
	 * @param pSpecies
	 *            the vector with species
	 * @param pDue
	 *            time before which we must have returned
	 */
	public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {
		if (DEBUG_LINES)
			System.err.println("\n--------------------------------------------------------------- REVEALING "
					+ currentRound + " ---------");
		String corrects = "            ";
		int guesses = 0;
		int guessAttempts = 0;
		for (int i = 0; i < guessesMade.length; i++) {
			if (guessesMade[i] == Constants.SPECIES_UNKNOWN)
				corrects += "    ";
			else {
				totalGuessAttempts++;
				guessAttempts++;
				if (guessesMade[i] == pSpecies[i]) {
					guesses++;
					totalGuesses++;
					corrects += "+  ";
				} else {
					corrects += "-  ";
				}
			}
		}
		if (DEBUG_ROUND_REVEAL) {
			System.err.println("My guess:  " + Arrays.toString(guessesMade));
			System.err.println("The truth: " + Arrays.toString(pSpecies));
			System.err.println(corrects);
		}

		for (int b = 0; b < numBirds; b++) {
			if (pSpecies[b] != -1) {
				sequencesPerSpecies.get(pSpecies[b]).add(Arrays.copyOf(Obirds[b], Tbirds[b]));
			}
		}

		for (int s = 0; s < Constants.COUNT_SPECIES; s++) {
			if (sequencesPerSpecies.get(s).isEmpty())
				continue;

			// train a new model
			double[][] Anew = MatrixHelper.newRowStochasticMatrix(N, N);
			double[][] Bnew = MatrixHelper.newRowStochasticMatrix(N, M);
			double[] Pinew = MatrixHelper.newRowStochasticArray(N);

			double[] statsNew = HiddenMarkovModel
					.baumWelchMultiSequences(sequencesPerSpecies.get(s).toArray(new int[][] {}), Anew, Bnew, Pinew);
			int itersNew = (int) statsNew[0];
			double logProbNew = statsNew[1];

			double logProbOld = Double.NEGATIVE_INFINITY;
			if (Aspecies[s] != null) {
				// train the old model
				double[][] Aold = Aspecies[s];
				double[][] Bold = Bspecies[s];
				double[] Piold = Pispecies[s];

				double[] statsOld = HiddenMarkovModel
						.baumWelchMultiSequences(sequencesPerSpecies.get(s).toArray(new int[][] {}), Aold, Bold, Piold);
				int itersOld = (int) statsOld[0];
				logProbOld = statsOld[1];
			}

			String debug = String.format("Species %d retrained on %d sequences\told confidence %f\tnew confidence %f\t",
					s, sequencesPerSpecies.get(s).size(), logProbOld, logProbNew);

			// which one is the best model?
			if (Double.isFinite(logProbOld) && Double.isFinite(logProbNew)) {
				if (logProbNew > logProbOld) {
					Aspecies[s] = Anew;
					Bspecies[s] = Bnew;
					Pispecies[s] = Pinew;
					debug += "kept new";
				} else {
					debug += "kept old";
				}
			} else if (Double.isFinite(logProbOld) && !Double.isFinite(logProbNew)) {
				debug += "kept old";
			} else if (!Double.isFinite(logProbOld) && Double.isFinite(logProbNew)) {
				Aspecies[s] = Anew;
				Bspecies[s] = Bnew;
				Pispecies[s] = Pinew;
				debug += "kept new";
			} else {
				Aspecies[s] = null;
				Bspecies[s] = null;
				Pispecies[s] = null;
				debug += "threw away both";
			}
			if (DEBUG_ROUND_REVEAL) {
				System.err.println(debug);
			}
		}

		if (DEBUG_ROUND_REVEAL) {
			int scoreByGuess = guesses - (guessAttempts - guesses);
			System.err.printf("Guesses:\t%d/%d\t%.2f%%\t\n", guesses, guessAttempts,
					(double) guesses / guessAttempts * 100);
			System.err.println("Score: " + scoreByGuess);
		}

		if (DEBUG_LINES) {
			System.err.println(
					"\n\n#########################################################################################");
			System.err.println(
					"######################################## END ROUND ######################################");
			System.err.println(
					"#########################################################################################\n");
		}

		if (currentRound == 9 && DEBUG_TOTAL_SCORES) {
			printTotalScores();
		}

	}

	private void printTotalScores() {
		int scoreByHits = totalHits - (totalShots - totalHits);
		int scoreByGuess = totalGuesses - (totalGuessAttempts - totalGuesses);
		System.err.printf("Shots:  \t%d/%d\t%.2f%%\t%d\n", totalHits, totalShots, (double) totalHits / totalShots * 100,
				scoreByHits);
		System.err.printf("Guesses:\t%d/%d\t%.2f%%\t%d\n", totalGuesses, totalGuessAttempts,
				(double) totalGuesses / totalGuessAttempts * 100, scoreByGuess);
		System.err.println("Score: " + (scoreByHits + scoreByGuess));
	}
}
