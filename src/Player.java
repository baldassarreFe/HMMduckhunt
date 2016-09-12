import java.util.Arrays;
import java.util.Random;

class Player {

	/** states */
	private static final int N = 5;
	/** emissions */
	private static final int M = Constants.COUNT_MOVE;
	/** Max sequence length */
	private static final int MAX_T = 100;

	/** do not try to learn if the seq is shorter */
	private static final int WAIT_UNTIL_UPDATING = 50;
	/** do not attempt to shoot if the seq is shorter */
	private static final int WAIT_UNTIL_SHOOTING = 70;

	/** Constant action for not shooting */
	public static final Action ACTION_DONT_SHOOT = new Action(-1, -1);
	private static final double MIN_PROB_TO_SHOOT = Math.log(0.6);
	private static final double MIN_PROB_TO_GUESS = Math.log(0.6);

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

	// Matrices for the best HMM that are able to identify a species across the
	// rounds
	// [species id] [i] [j]
	double[][][] Aspecies;
	double[][][] Bspecies;
	double[][] Pispecies;

	// info variables
	int currentRound = -1;
	int currentObsLength = 0;
	private int numBirds;

	// overall statistics
	int totalBirds = 0;
	int totalShots = 0;
	int totalHits = 0;

	public Player() {
		// Initialize the matrices for the species
		Aspecies = new double[Constants.COUNT_SPECIES][N][N];
		Bspecies = new double[Constants.COUNT_SPECIES][N][M];
		Pispecies = new double[Constants.COUNT_SPECIES][N];

		// TODO probably we can do something smarter for B
		for (int s = 0; s < Constants.COUNT_SPECIES; s++) {
			MatrixHelper.rowStochasticMatrix(Aspecies[s]);
			MatrixHelper.rowStochasticMatrix(Bspecies[s]);
			MatrixHelper.rowStochasticArray(Pispecies[s]);
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
			if (bird.isAlive()) {
				for (int t = Tbirds[b]; t < Tbirds[b] + newTurns; t++) {
					Obirds[b][t] = bird.getObservation(t);
				}
				Tbirds[b] += newTurns;
			}
		}

		// update matrices of this round
		if (currentObsLength >= WAIT_UNTIL_UPDATING) {
			reEstimateBirdModels(pState, newTurns);

			// choose whether to shoot or not
			if (currentObsLength >= WAIT_UNTIL_SHOOTING)
				return makeShot(pState);
		}

		return ACTION_DONT_SHOOT;
	}

	private Action makeShot(GameState pState) {
		double bestLogProb = Double.NEGATIVE_INFINITY;
		int bestBird = -1;
		int bestEmission = -1;
		for (int b = 0; b < numBirds; b++) {
			if (pState.getBird(b).isAlive()) {
				for (int o = 0; o < M; o++) {
					int[] O = Arrays.copyOf(Obirds[b], Tbirds[b] + 1);
					O[Tbirds[b]] = o;
					double logProb = HiddenMarkovModel.probabilityLogOfEmissionSequence(O, Abirds[b], Bbirds[b],
							Pibirds[b]);
					//System.err.printf("Evaluating:\tbird %d\tmove %d\tconfidence %f\n", b, o, logProb);
					if (logProb > bestLogProb) {
						bestLogProb = logProb;
						bestEmission = o;
						bestBird = b;
					}
				}
			}
		}

		//System.err.printf("Best guess:\tbird %d\tmove %d\tconfidence %f\n", bestBird, bestEmission, bestLogProb);

		if (bestLogProb > MIN_PROB_TO_SHOOT) {
			return new Action(bestBird, bestEmission);
		} else {
			return ACTION_DONT_SHOOT;
		}
	}

	// Use the new emissions to update the matrices (Baum Welch)
	// 1. Start with a new randomized model, use baum welch to estimate it
	// 2. Test the old model and the new against the sequence
	// 3. Keep the best model
	private void reEstimateBirdModels(GameState pState, int newTurns) {

		// TODO parallelize over the birds?
		for (int b = 0; b < numBirds; b++) {
			Bird bird = pState.getBird(b);
			if (bird.isAlive()) {
				double[][] Anew = new double[N][N];
				double[][] Bnew = new double[N][M];
				double[] Pinew = new double[N];
//				MatrixHelper.rowStochasticMatrix(Anew);
//				MatrixHelper.rowStochasticMatrix(Bnew);
				MatrixHelper.almostDiagonalMatrix(Anew);
				MatrixHelper.almostDiagonalMatrix(Bnew);
				MatrixHelper.rowStochasticArray(Pinew);

				int T = Tbirds[b];
				int[] O = Arrays.copyOf(Obirds[b], T);

				double[] stats = HiddenMarkovModel.baumWelch(O, Anew, Bnew, Pinew);
				int iters = (int) stats[0];
				double logProbNew = stats[1];

				double[][] Aold = Abirds[b];
				double[][] Bold = Bbirds[b];
				double[] Piold = Pibirds[b];

				double logProbOld = HiddenMarkovModel.probabilityLogOfEmissionSequence(O, Aold, Bold, Piold);

				System.err.printf("Round %d\tbird %d\tT %d\titers %d\tlogProbNew %f\tlogProbOld %f\t", currentRound, b, T, iters, logProbNew, logProbOld);

				// best model?
				if (Double.isNaN(logProbOld) || logProbNew > logProbOld) {
					Abirds[b] = Anew;
					Bbirds[b] = Bnew;
					Pibirds[b] = Pinew;
					System.err.println("kept new");
				} else {
					System.err.println("kept old");
				}
			}
		}
	}

	// At the beginning of every round initialize the matrices A,B,Pi
	private void initRound(GameState pState) {
		currentObsLength = 0;
		numBirds = 1;//pState.getNumBirds();
		Abirds = new double[numBirds][N][N];
		Bbirds = new double[numBirds][N][M];
		Pibirds = new double[numBirds][N];
		Obirds = new int[numBirds][MAX_T];
		Tbirds = new int[numBirds];
		currentRound++;
		totalBirds += numBirds;

		for (int b = 0; b < numBirds; b++) {
			MatrixHelper.rowStochasticMatrix(Abirds[b]);
			MatrixHelper.rowStochasticMatrix(Bbirds[b]);
			MatrixHelper.rowStochasticArray(Pibirds[b]);
		}
	}

	/**
	 * Guess the species! This function will be called at the end of each round,
	 * to give you a chance to identify the species of the birds for extra
	 * points.
	 *
	 * Fill the vector with guesses for the all birds. Use SPECIES_UNKNOWN to
	 * avoid guessing.
	 *
	 * @param pState
	 *            the GameState object with observations etc
	 * @param pDue
	 *            time before which we must have returned
	 * @return a vector with guesses for all the birds
	 */
	public int[] guess(GameState pState, Deadline pDue) {
		// for every bird evaluate the likelihood of having generated the
		// sequence,
		// considering A, B from the species model and Pi from the specific bird
		int[] guesses = new int[numBirds];
		
		if (currentRound==0) {			
			Arrays.fill(guesses, Constants.SPECIES_UNKNOWN);
		}
		
		for (int b = 0; b < numBirds; b++) {
			double bestLogProb = MIN_PROB_TO_GUESS;
			for (int s = 0; s < Constants.COUNT_SPECIES; s++) {
				int O[] = Arrays.copyOf(Obirds[b], Tbirds[b]);
				double logProb = HiddenMarkovModel.probabilityLogOfEmissionSequence(O, Aspecies[s], Bspecies[s],
						Pibirds[b]);
				//System.err.printf("Evaluating:\tbird %d\tspecies %d\tconfidence %f\n", b, s, logProb);
				if (!Double.isNaN(logProb) && logProb > bestLogProb) {
					bestLogProb = logProb;
					guesses[b] = s;
				}					
			}
		}
		return guesses;
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
		System.err.println("HIT BIRD!!!");
		totalHits++;
	}

	/**
	 * If you made any guesses, you will find out the true species of those
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
		System.err.println(Arrays.toString(pSpecies));
		for (int b = 0; b < numBirds; b++) {
			
		}
	}
}
