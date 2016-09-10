
class Player {
	
	// states
	private static final int N = 5;
	// emissions
	private static final int M = Constants.COUNT_MOVE;
	// do not try to learn if the seq is shorter
	private static final int MIN_T = 30;
	
	// Matrices for the HMM of the birds in the current round
	// A = transition
	// B = emissions
	// Pi = initial state
	// [bird number] [i] [j]
	double Abirds[][][];
	double Bbirds[][][];
	double Pibirds[][];
	
	// Matrices for the best HMM that are able to identify a species across the rounds
	// [species id] [i] [j]
	double Aspecies[][][];
	double Bspecies[][][];
	double Pispecies[][];
	
	// info variables
	int currentRound = 0;
	int T = 0;

    public Player() {
    	Aspecies = new double[Constants.COUNT_SPECIES][N][N];
    	Bspecies = new double[Constants.COUNT_SPECIES][N][M];
    	Pispecies = new double[Constants.COUNT_SPECIES][N];
    }

    /**
     * Shoot!
     *
     * This is the function where you start your work.
     *
     * You will receive a variable pState, which contains information about all
     * birds, both dead and alive. Each bird contains all past moves.
     *
     * The state also contains the scores for all players and the number of
     * time steps elapsed since the last time this function was called.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return the prediction of a bird we want to shoot at, or cDontShoot to pass
     */
    public Action shoot(GameState pState, Deadline pDue) {
    	
    	// is it a new round?
    	if(pState.getRound()!=currentRound){
    		initRound();
    	}
    	
    	T += pState.getNumNewTurns();
    	if (T < MIN_T)
    		return cDontShoot;

        // This line would predict that bird 0 will move right and shoot at it.
        // return Action(0, MOVE_RIGHT);
    }

    /**
     * Guess the species!
     * This function will be called at the end of each round, to give you
     * a chance to identify the species of the birds for extra points.
     *
     * Fill the vector with guesses for the all birds.
     * Use SPECIES_UNKNOWN to avoid guessing.
     *
     * @param pState the GameState object with observations etc
     * @param pDue time before which we must have returned
     * @return a vector with guesses for all the birds
     */
    public int[] guess(GameState pState, Deadline pDue) {
        /*
         * Here you should write your clever algorithms to guess the species of
         * each bird. This skeleton makes no guesses, better safe than sorry!
         */

        int[] lGuess = new int[pState.getNumBirds()];
        for (int i = 0; i < pState.getNumBirds(); ++i)
            lGuess[i] = Constants.SPECIES_UNKNOWN;
        return lGuess;
    }

    /**
     * If you hit the bird you were trying to shoot, you will be notified
     * through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pBird the bird you hit
     * @param pDue time before which we must have returned
     */
    public void hit(GameState pState, int pBird, Deadline pDue) {
        System.err.println("HIT BIRD!!!");
    }

    /**
     * If you made any guesses, you will find out the true species of those
     * birds through this function.
     *
     * @param pState the GameState object with observations etc
     * @param pSpecies the vector with species
     * @param pDue time before which we must have returned
     */
    public void reveal(GameState pState, int[] pSpecies, Deadline pDue) {
    }

    public static final Action cDontShoot = new Action(-1, -1);
}
