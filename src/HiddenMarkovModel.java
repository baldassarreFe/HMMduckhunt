import java.util.Arrays;

public class HiddenMarkovModel {

	private static final double NO_DIV_BY_ZERO = Double.MIN_NORMAL;
	/** baum welch iterations */
	private static final long MAX_ITERATIONS = 30;

	public static double[] nextEmissions(double[][] A, double[][] B, double[] Pi) {
		return currentEmissions(B, nextStates(A, Pi));
	}

	public static double[] nextStates(double[][] A, double[] Pi) {
		return MatrixHelper.multiply(new double[][] { Pi }, A)[0];
	}

	public static double[] currentEmissions(double[][] B, double[] Pi) {
		int M = B[0].length;
		double[] result = new double[M];
		for (int currentEmission = 0; currentEmission < M; currentEmission++) {
			result[currentEmission] = MatrixHelper.multiply(Pi, MatrixHelper.getColumn(B, currentEmission));
		}
		return result;
	}

	public static double nextMostLikelyEmission(double[][] A, double[][] B, double[] Pi) {
		return MatrixHelper.indexOfHigherElement(currentEmissions(B, nextStates(A, Pi)));
	}

	public static double nextMostLikelyState(double[][] A, double[] Pi) {
		return MatrixHelper.indexOfHigherElement(MatrixHelper.multiply(new double[][] { Pi }, A)[0]);
	}

	public static double currentMostLikelyEmission(double[][] B, double[] Pi) {
		return MatrixHelper.indexOfHigherElement(currentEmissions(B, Pi));
	}

	public static double[] nextStateDistributionKnowingSequence(int[] O, double[][] A, double[][] B, double[] Pi) {
		int T = O.length;
		int N = A.length;
		@SuppressWarnings("unused")
		int M = B[0].length;

		double[][] alpha = new double[T][N];
		double[] c = new double[T];
		double[][] beta = new double[T][N];
		double[][] gamma = new double[T][N];
		double[][][] digamma = new double[T][N][N];
		alphaPass(O, A, B, Pi, alpha, c);
		betaPass(O, A, B, c, beta);
		digammaGamma(O, A, B, alpha, beta, digamma, gamma);

		// gamma is defined as:
		// gamma[t][i] = P(X[t] = i | O[0:T-1], model)
		return MatrixHelper.getRow(gamma, T - 1);
	}

	public static double[] nextEmissionDistributionKnowingSequence(int[] O, double[][] A, double[][] B, double[] Pi) {
		double[] nextState = nextStateDistributionKnowingSequence(O, A, B, Pi);
		return currentEmissions(B, nextState);
	}

	public static int[] estimateStateSequence(int[] O, double[][] A, double[][] B, double[] Pi) {
		int T = O.length;
		int result[] = new int[T];
		@SuppressWarnings("unused")
		double probability = viterbi(O, A, B, Pi, result);
		return result;
	}

	public static double probabilityOfEmissionSequence(int[] O, double[][] A, double[][] B, double[] Pi) {
		int N = A.length;
		int T = O.length;
		double[][] alpha = new double[T][N];
		double[] c = new double[T];
		double result;

		alphaPass(O, A, B, Pi, alpha, c);
		result = 1;
		for (int t = 0; t < T; t++) {
			result *= c[t];
		}
		result = 1.0 / result;

		return result;
	}

	public static double probabilityLogOfEmissionSequence(int[] O, double[][] A, double[][] B, double[] Pi) {
		return alphaPass(O, A, B, Pi);
	}

	public static double viterbi(int[] O, double[][] A, double[][] B, double[] Pi, int[] X) {
		int N = A.length;
		int T = O.length;

		double delta[][] = new double[T][N];
		int indexes[][] = new int[T][N];

		for (int i = 0; i < N; i++) {
			delta[0][i] = Pi[i] * B[i][O[0]];
			// useless because indexes[0] is never used
			indexes[0][i] = 0;
		}

		for (int t = 1; t < T; t++) {
			for (int i = 0; i < N; i++) {
				double max = 0;
				for (int j = 0; j < N; j++) {
					double temp = delta[t - 1][j] * A[j][i];
					if (temp > max) {
						max = temp;
						indexes[t][i] = j;
					}
				}
				delta[t][i] = max * B[i][O[t]];
			}
		}

		double probability = 0.0;
		int maxIndex = 0;
		for (int i = 0; i < N; i++) {
			if (probability < delta[T - 1][i]) {
				probability = delta[T - 1][i];
				maxIndex = i;
			}
		}

		for (int t = T - 1; t >= 0; t--) {
			X[t] = maxIndex;
			maxIndex = indexes[t][maxIndex];
		}

		return probability;
	}

	public static double viterbiNoLogs(int[] O, double[][] A, double[][] B, double[] Pi, int[] X) {
		int N = A.length;
		int T = O.length;

		double delta[][] = new double[T][N];
		int indexes[][] = new int[T][N];

		for (int i = 0; i < N; i++) {
			delta[0][i] = Math.log(Pi[i] * B[i][O[0]]);
			// useless because indexes[0] is never used
			indexes[0][i] = 0;
		}

		for (int t = 1; t < T; t++) {
			for (int i = 0; i < N; i++) {
				double max = Double.NEGATIVE_INFINITY;
				for (int j = 0; j < N; j++) {
					double temp = delta[t - 1][j] + Math.log(A[j][i]);
					if (temp > max) {
						max = temp;
						indexes[t][i] = j;
					}
				}
				delta[t][i] = max + Math.log(B[i][O[t]]);
			}
		}

		double probability = Double.NEGATIVE_INFINITY;
		int maxIndex = 0;
		for (int i = 0; i < N; i++) {
			if (probability < delta[T - 1][i]) {
				probability = delta[T - 1][i];
				maxIndex = i;
			}
		}

		for (int t = T - 1; t >= 0; t--) {
			X[t] = maxIndex;
			maxIndex = indexes[t][maxIndex];
		}

		return probability;
	}

	// public static double giveMeAnyModel(int[] O, double[][] A, double[][] B,
	// double[] Pi) {
	// // initialization
	// int N = A.length;
	// int M = B[0].length;
	// int T = O.length;
	//
	// double logProb;
	// double bestLogProb = Double.NEGATIVE_INFINITY;
	//
	// // no normalization new
	// double[][] Anew = MatrixHelper.newRowStochasticMatrix(N, N);
	// double[][] Bnew = MatrixHelper.newRowStochasticMatrix(N, M);
	// double[] Pinew = MatrixHelper.newRowStochasticArray(N);
	// logProb = Math.log(baumWelchNoNormalization(O, Anew, Bnew, Pinew)[1]);
	// if (logProb > bestLogProb) {
	// bestLogProb = logProb;
	// A = Anew
	// }
	//
	// double[][] Aold = MatrixHelper.copy(A);
	// double[][] Bold = MatrixHelper.copy(B);
	// double[] Piold = MatrixHelper.copy(Pi);
	// }

	/**
	 * @param O
	 * @param A
	 * @param B
	 * @param Pi
	 * @return
	 */
	public static double[] baumWelch(int[] O, double[][] A, double[][] B, double[] Pi) {
		// initialization
		int N = A.length;
		@SuppressWarnings("unused")
		int M = B[0].length;
		int T = O.length;

		double oldProbLog;
		double logProbAfter, probAfter;
		double logProbBefore = Double.NEGATIVE_INFINITY;
		long iters = 0;
		do {
			oldProbLog = logProbBefore;

			double[][] alpha = new double[T][N];
			double[] c = new double[T];
			// alpha pass also gives log( P(O|lambda) ) of the model that we're
			// going to update
			logProbBefore = alphaPass(O, A, B, Pi, alpha, c);

			double[][] beta = new double[T][N];
			betaPass(O, A, B, c, beta);

			double[][][] digamma = new double[T][N][N];
			double[][] gamma = new double[T][N];
			digammaGamma(O, A, B, alpha, beta, digamma, gamma);

			reEstimate(O, A, B, Pi, digamma, gamma);

			// compute log( P(O|lambda) ) for the new model
			// computing it again with a new alpha pass gives a different
			// result because A,B,Pi has changed, the latter should be more
			// correct right?
			// logProbAfter = alphaPass(O, A, B, Pi);
			// probAfter = alphaPassNoNormalization(O, A, B, Pi);
			logProbAfter = logProbBefore;
		} while (++iters < MAX_ITERATIONS && logProbAfter > oldProbLog);

		// System.err.println("\t\t\t" + logProbBefore + "\t" + logProbAfter +
		// "\t" + probAfter);
		return new double[] { iters, logProbAfter };
	}

	public static double[] baumWelchClustered(int[] O, double[][] A, double[][] B, double[] Pi, int cluserSize) {
		int[][] Osplit = new int[O.length / cluserSize][];
		for (int i = 0; i < Osplit.length; i++) {
			Osplit[i] = Arrays.copyOfRange(O, cluserSize * i, cluserSize * (1 + i) - 1);
		}
		return baumWelchMultiSequences(Osplit, A, B, Pi);
	}

	/**
	 * @param Otrain
	 * @param A
	 * @param B
	 * @param Pi
	 * @return
	 */
	public static double[] baumWelchAvoidOverfitting(int[] O, double[][] A, double[][] B, double[] Pi,
			double partitionRatio) {
		// initialization
		int N = A.length;
		@SuppressWarnings("unused")
		int M = B[0].length;

		int[] Otrain = Arrays.copyOf(O, (int) (O.length * partitionRatio));
		int[] Ovalidate = Arrays.copyOfRange(O, Otrain.length, O.length);
		int T = Otrain.length;

		double oldProbLog;
		double logProb2;
		double logProb = Double.NEGATIVE_INFINITY;
		long iters = 0;
		do {
			oldProbLog = logProb;

			double[][] alpha = new double[T][N];
			double[] c = new double[T];
			// alpha pass also gives log( P(O|lambda) ) of the model that we're
			// going to update
			logProb = alphaPass(Otrain, A, B, Pi, alpha, c);

			double[][] beta = new double[T][N];
			betaPass(Otrain, A, B, c, beta);

			double[][][] digamma = new double[T][N][N];
			double[][] gamma = new double[T][N];
			digammaGamma(Otrain, A, B, alpha, beta, digamma, gamma);

			reEstimate(Otrain, A, B, Pi, digamma, gamma);

			// compute log( P(O|lambda) ) for the new model
			// computing it again with a new alpha pass gives a different
			// result because A,B,Pi has changed, the latter should be more
			// correct right?
			logProb2 = alphaPass(O, A, B, Pi);
		} while (++iters < MAX_ITERATIONS && logProb2 > oldProbLog);

		// System.err.println(logProb + "\t" + logProb2);
		return new double[] { iters, logProb2 };
	}

	public static double alphaPass(int[] O, double[][] A, double[][] B, double[] Pi) {
		int N = A.length;
		@SuppressWarnings("unused")
		int M = B[0].length;
		int T = O.length;
		double[][] alpha = new double[T][N];
		double[] c = new double[T];
		return alphaPass(O, A, B, Pi, alpha, c);
	}

	public static double alphaPass(int[] O, double[][] A, double[][] B, double[] Pi, double[][] alpha, double[] c) {
		int N = A.length;
		int T = O.length;

		c[0] = 0;
		for (int i = 0; i < N; i++) {
			alpha[0][i] = Pi[i] * B[i][O[0]];
			c[0] += alpha[0][i];
		}
		// HACK
		if (Double.isNaN(c[0]) || c[0] == 0) {
			c[0] = NO_DIV_BY_ZERO;
		}
		c[0] = 1.0 / c[0];
		for (int i = 0; i < N; i++) {
			alpha[0][i] *= c[0];
		}

		for (int t = 1; t < T; t++) {
			c[t] = 0;
			for (int i = 0; i < N; i++) {
				alpha[t][i] = 0;
				for (int j = 0; j < N; j++) {
					alpha[t][i] += alpha[t - 1][j] * A[j][i];
				}
				alpha[t][i] *= B[i][O[t]];
				c[t] += alpha[t][i];
			}
			// HACK
			if (Double.isNaN(c[t]) || c[t] == 0) {
				c[t] = NO_DIV_BY_ZERO;
			}
			c[t] = 1.0 / c[t];
			for (int i = 0; i < N; i++) {
				alpha[t][i] *= c[t];
			}
		}

		double logProb = 0;
		for (int t = 0; t < T; t++) {
			logProb += Math.log(c[t]);
		}
		return -logProb;
	}

	private static void betaPass(int[] O, double[][] A, double[][] B, double[] c, double[][] beta) {
		int N = A.length;
		int T = O.length;

		for (int i = 0; i < N; i++) {
			beta[T - 1][i] = c[T - 1];
		}

		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < N; i++) {
				beta[t][i] = 0;
				for (int j = 0; j < N; j++) {
					beta[t][i] += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
				}
				beta[t][i] *= c[t];
			}
		}
	}

	private static void digammaGamma(int[] O, double[][] A, double[][] B, double[][] alpha, double[][] beta,
			double[][][] digamma, double[][] gamma) {
		int N = A.length;
		int T = O.length;

		for (int t = 0; t < T - 1; t++) {
			double denom = 0;
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					denom += alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
				}
			}
			for (int i = 0; i < N; i++) {
				gamma[t][i] = 0;
				for (int j = 0; j < N; j++) {
					digamma[t][i][j] = alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j] / denom;
					gamma[t][i] += digamma[t][i][j];
				}
			}
		}

		// special case T-1
		double denom = 0;
		for (int i = 0; i < N; i++) {
			denom += alpha[T - 1][i];
		}
		for (int i = 0; i < N; i++) {
			gamma[T - 1][i] = alpha[T - 1][i] / denom;
		}
	}

	private static void reEstimate(int[] O, double[][] A, double[][] B, double[] Pi, double[][][] digamma,
			double[][] gamma) {

		int N = A.length;
		int M = B[0].length;
		int T = O.length;

		for (int i = 0; i < N; i++) {
			// Pi
			Pi[i] = gamma[0][i];
		}
		for (int i = 0; i < N; i++) {
			// A
			for (int j = 0; j < N; j++) {
				double numerA = 0;
				double denomA = 0;
				for (int t = 0; t < T - 1; t++) {
					numerA += digamma[t][i][j];
					denomA += gamma[t][i];
				}
				A[i][j] = numerA / denomA;
			}
		}

		for (int i = 0; i < N; i++) {
			// B
			for (int j = 0; j < M; j++) {
				double numerB = 0;
				double denomB = 0;
				for (int t = 0; t < T; t++) {
					if (O[t] == j) {
						numerB += gamma[t][i];
					}
					denomB += gamma[t][i];
				}
				B[i][j] = numerB / denomB;
			}
		}
	}

	/**
	 * @param O
	 * @param A
	 * @param B
	 * @param Pi
	 * @return
	 */
	public static double[] baumWelchNoNormalization(int[] O, double[][] A, double[][] B, double[] Pi) {
		// initialization
		int N = A.length;
		@SuppressWarnings("unused")
		int M = B[0].length;
		int T = O.length;

		double oldProbLog;
		double logProb2;
		double logProb = 0.0;
		long iters = 0;
		do {
			oldProbLog = logProb;

			double[][] alpha = new double[T][N];
			// alpha pass also gives P(O|lambda) of the model that we're
			// going to update
			logProb = alphaPassNoNormalization(O, A, B, Pi, alpha);

			double[][] beta = new double[T][N];
			betaPassNoNormalization(O, A, B, beta);

			double[][][] digamma = new double[T][N][N];
			double[][] gamma = new double[T][N];
			digammaGammaNoNormalization(O, A, B, alpha, beta, digamma, gamma);

			reEstimateNoNormalization(O, A, B, Pi, digamma, gamma);

			// compute P(O|lambda) for the new model
			// computing it again with a new alpha pass gives a different
			// result because A,B,Pi has changed, the latter should be more
			// correct right?
			logProb2 = alphaPassNoNormalization(O, A, B, Pi);
		} while (++iters < MAX_ITERATIONS && logProb2 > oldProbLog);

		// System.err.println(logProb + "\t" + logProb2);
		return new double[] { iters, logProb2 };
	}

	private static double alphaPassNoNormalization(int[] O, double[][] A, double[][] B, double[] Pi) {
		int N = A.length;
		int T = O.length;
		return alphaPassNoNormalization(O, A, B, Pi, new double[T][N]);
	}

	private static double alphaPassNoNormalization(int[] O, double[][] A, double[][] B, double[] Pi, double[][] alpha) {
		int N = A.length;
		int T = O.length;

		for (int i = 0; i < N; i++) {
			alpha[0][i] = Pi[i] * B[i][O[0]];
		}

		for (int t = 1; t < T; t++) {
			for (int i = 0; i < N; i++) {
				alpha[t][i] = 0;
				for (int j = 0; j < N; j++) {
					alpha[t][i] += alpha[t - 1][j] * A[j][i];
				}
				alpha[t][i] *= B[i][O[t]];
			}
		}
		double result = 0;
		for (int i = 0; i < N; i++) {
			result += alpha[T - 1][i];
		}
		return result;
	}

	private static void betaPassNoNormalization(int[] O, double[][] A, double[][] B, double[][] beta) {
		int N = A.length;
		int T = O.length;

		for (int i = 0; i < N; i++) {
			beta[T - 1][i] = 1;
		}

		for (int t = T - 2; t >= 0; t--) {
			for (int i = 0; i < N; i++) {
				beta[t][i] = 0;
				for (int j = 0; j < N; j++) {
					beta[t][i] += A[i][j] * B[j][O[t + 1]] * beta[t + 1][j];
				}
			}
		}
	}

	private static void digammaGammaNoNormalization(int[] O, double[][] A, double[][] B, double[][] alpha,
			double[][] beta, double[][][] digamma, double[][] gamma) {

		int N = A.length;
		int T = O.length;

		double denom = 0;
		for (int i = 0; i < N; i++) {
			denom += alpha[T - 1][i];
		}
		for (int t = 0; t < T - 1; t++) {
			for (int i = 0; i < N; i++) {
				for (int j = 0; j < N; j++) {
					digamma[t][i][j] = alpha[t][i] * A[i][j] * B[j][O[t + 1]] * beta[t + 1][j] / denom;
				}
			}
		}
		for (int t = 0; t < T - 1; t++) {
			for (int i = 0; i < N; i++) {
				gamma[t][i] = 0;
				for (int j = 0; j < N; j++) {
					gamma[t][i] += digamma[t][i][j];
				}
			}
		}
	}

	private static void reEstimateNoNormalization(int[] O, double[][] A, double[][] B, double[] Pi,
			double[][][] digamma, double[][] gamma) {
		int N = A.length;
		int M = B[0].length;
		int T = O.length;

		for (int i = 0; i < N; i++) {
			// Pi
			Pi[i] = gamma[0][i];

			// A and B denom
			double denom = 0;
			for (int t = 0; t < T - 1; t++) {
				denom += gamma[t][i];
			}

			// A
			for (int j = 0; j < N; j++) {
				double numerA = 0;
				for (int t = 0; t < T - 1; t++) {
					numerA += digamma[t][i][j];
				}
				A[i][j] = numerA / denom;
			}

			// B
			for (int j = 0; j < M; j++) {
				double numerB = 0;
				for (int t = 0; t < T; t++) {
					if (O[t] == j) {
						numerB += gamma[t][i];
					}
				}
				B[i][j] = numerB / denom;
			}
		}
	}

	/**
	 * @param O
	 * @param A
	 * @param B
	 * @param Pi
	 * @return
	 */
	public static double[] baumWelchMultiSequences(int[][] O, double[][] A, double[][] B, double[] Pi) {
		// initialization
		int N = A.length;
		@SuppressWarnings("unused")
		int M = B[0].length;
		int NUM_SEQ = O.length;

		double oldProbLog;
		double logProb = Double.NEGATIVE_INFINITY;
		long iters = 0;
		do {
			oldProbLog = logProb;

			// only certain sequences are relevant
			int validSequences = 0;
			double[][][][] digammas = new double[NUM_SEQ][][][];
			double[][][] gammas = new double[NUM_SEQ][][];
			int[][] validOs = new int[NUM_SEQ][];

			for (int seqIndex = 0; seqIndex < NUM_SEQ; seqIndex++) {
				int T = O[seqIndex].length;
				double[][] alpha = new double[T][N];
				double[] c = new double[T];

				logProb = alphaPass(O[seqIndex], A, B, Pi, alpha, c);
				if (!Double.isFinite(logProb)) {
					continue;
				}

				double[][] beta = new double[T][N];
				betaPass(O[seqIndex], A, B, c, beta);

				validOs[validSequences] = O[seqIndex];
				digammas[validSequences] = new double[T][N][N];
				gammas[validSequences] = new double[T][N];
				digammaGamma(O[seqIndex], A, B, alpha, beta, digammas[validSequences], gammas[validSequences]);
				validSequences++;
			}

			reEstimateMultiSequences(Arrays.copyOf(validOs, validSequences), A, B, Pi,
					Arrays.copyOf(digammas, validSequences), Arrays.copyOf(gammas, validSequences));
		} while (++iters < MAX_ITERATIONS && logProb > oldProbLog);

		// System.err.println(logProb + "\t" + logProb2);
		return new double[] { iters, logProb };
	}

	private static void reEstimateMultiSequences(int[][] O, double[][] A, double[][] B, double[] Pi,
			double[][][][] digamma, double[][][] gamma) {

		int N = A.length;
		int M = B[0].length;
		int NUM_SEQ = O.length;

		// Pi
		for (int i = 0; i < N; i++) {
			for (int seqIndex = 0; seqIndex < NUM_SEQ; seqIndex++) {
				Pi[i] += gamma[seqIndex][0][i];
			}
			Pi[i] /= NUM_SEQ;
		}

		// A
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				double numerA = 0;
				double denomA = 0;
				for (int seqIndex = 0; seqIndex < NUM_SEQ; seqIndex++) {
					int T = O[seqIndex].length;
					for (int t = 0; t < T - 1; t++) {
						numerA += digamma[seqIndex][t][i][j];
						denomA += gamma[seqIndex][t][i];
					}
				}
				A[i][j] = numerA / denomA;
			}
		}

		// B
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < M; j++) {
				double numerB = 0;
				double denomB = 0;
				for (int seqIndex = 0; seqIndex < NUM_SEQ; seqIndex++) {
					int T = O[seqIndex].length;
					for (int t = 0; t < T; t++) {
						if (O[seqIndex][t] == j) {
							numerB += gamma[seqIndex][t][i];
						}
						denomB += gamma[seqIndex][t][i];
					}
				}
				B[i][j] = numerB / denomB;
			}
		}
	}

}
